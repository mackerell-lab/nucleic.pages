import test from 'node:test';
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { BUNDLED_SURVEY_ENCODING, expandBundledSurveyColumns, verifySurveyBundle } from '../core/bundled-survey-codec.js';
import { decodeSurveyRows } from '../core/survey-codec.js';
import { PACKED_SURVEY_ENCODING, encodePackedSurvey, expandPackedSurvey } from '../core/packed-survey-codec.js';

const hash = value => createHash('sha256').update(JSON.stringify(value)).digest('hex');
const noBundle = () => assert.fail('Unexpected bundle request');
const table = values => ({ encoding: BUNDLED_SURVEY_ENCODING, build_id: 'test', row_count: values.length,
  columns: { id: values.map((_, i) => `residue:${i}|survey|o2`), observation_id: values.map((_, i) => `residue:${i}`),
    term_id: values.map(() => 'o2'), value: values, status: values.map(value => value === null ? 'missing_atom' : 'ok') } });
// Independent Node IEEE writer (also constructs malicious nonfinite transports).
function numeric(values) {
  const raw = Buffer.alloc(values.length * 8), lanes = Buffer.alloc(raw.length);
  values.forEach((value, index) => raw.writeDoubleLE(value, index * 8));
  for (let index = 0; index < values.length; index++) for (let byte = 0; byte < 8; byte++) lanes[byte * values.length + index] = raw[index * 8 + byte];
  return { encoding: 'float64-le-shuffled-base64-1', count: values.length, data: lanes.toString('base64') };
}

test('Survey Float64, nulls and factored identities survive JSON with exact rows', async () => {
  const values = [0, -0, Number.MIN_VALUE, -Number.MIN_VALUE, Number.MAX_VALUE, -Number.MAX_VALUE, 1.2345678901234567, null];
  const input = table(values), encoded = await encodePackedSurvey(input, noBundle);
  assert.equal(encoded.encoding, PACKED_SURVEY_ENCODING);
  assert.deepEqual(encoded.columns.id, { encoding: 'rna-survey-id-suffix-1', column: 'observation_id', suffix: '|survey|o2' });
  assert.deepEqual(encoded.columns.term_id, { encoding: 'rna-survey-constant-1', value: 'o2' });
  assert.deepEqual(encoded.columns.value.values, numeric(values.map(value => value === null ? 0 : value)));
  const expanded = await expandPackedSurvey(JSON.parse(JSON.stringify(encoded)), noBundle);
  assert.deepEqual(expanded, await expandBundledSurveyColumns(input, noBundle));
  const rows = decodeSurveyRows(expanded);
  values.forEach((value, index) => assert(Object.is(rows[index].value, value)));
  assert.deepEqual(rows.map(row => row.id), input.columns.id);
  assert.deepEqual(decodeSurveyRows(expanded, ['value']), values.map(value => ({ value })));
});

test('Survey factoring uses authenticated observation strings and retains column order', async () => {
  const input = table([1.5, null]), observations = input.columns.observation_id;
  const column = hash(observations), bundle = { encoding: 'rna-survey-column-bundle-1', columns: { [column]: observations } }, reference = hash(bundle);
  input.columns.observation_id = { bundle: reference, column };
  const calls = [], reader = async key => { calls.push(key); return (await verifySurveyBundle(key, bundle)).bundle; };
  const encoded = await encodePackedSurvey(input, reader);
  const expanded = await expandPackedSurvey(encoded, reader);
  assert.deepEqual(calls, [reference, reference]);
  assert.deepEqual(Object.keys(expanded.columns), Object.keys(input.columns));
  assert.deepEqual(expanded.columns.id, ['residue:0|survey|o2', 'residue:1|survey|o2']);
  await assert.rejects(expandPackedSurvey(encoded, key => verifySurveyBundle(key, { ...bundle, columns: { [column]: ['evil', 'data'] } })), /hash mismatch/);
});

test('Ineligible identities and varying terms remain exact fallback columns', async () => {
  for (const mutate of [
    data => { data.columns.id[0] = 'external-id'; },
    data => { data.columns.term_id[1] = 'another-term'; },
    data => { data.columns.observation_id[0] = null; },
    data => { delete data.columns.observation_id; },
    data => { data.columns.observation_id[0] = ''; },
    data => { data.columns.term_id = ['', '']; },
  ]) {
    const input = table([1, null]); mutate(input);
    const encoded = await encodePackedSurvey(input, noBundle);
    assert.equal(encoded.columns.id, input.columns.id);
    assert.deepEqual(await expandPackedSurvey(encoded, noBundle), await expandBundledSurveyColumns(input, noBundle));
  }
  for (const values of [[], [null, null]]) {
    const input = table(values);
    assert.deepEqual(await expandPackedSurvey(await encodePackedSurvey(input, noBundle), noBundle), await expandBundledSurveyColumns(input, noBundle));
  }
});

test('Packed Survey rejects unknown descriptors, cycles and malformed columns before I/O', async () => {
  const base = await encodePackedSurvey(table([1, null]), noBundle);
  const mutations = [
    d => { d.encoding = 'unknown'; }, d => { d.row_count = -1; }, d => { d.row_count = 2.5; },
    d => { d.columns = []; }, d => { d.columns.id.extra = true; }, d => { d.columns.id.suffix = 12; },
    d => { d.columns.id.column = 'id'; }, d => { d.columns.id.column = 'term_id'; },
    d => { d.columns.id.encoding = 'future'; }, d => { delete d.columns.observation_id; },
    d => { delete d.columns.term_id; }, d => { d.columns.term_id.value = ''; },
    d => { d.columns.observation_id = d.columns.id; }, d => { d.columns.term_id.value = 12; },
    d => { d.columns.status = d.columns.term_id; }, d => { d.columns.status = ['ok']; },
    d => { d.columns.status = Array(2); }, d => { d.columns.status = [undefined, null]; },
    d => { d.columns.observation_id = { bundle: '../bad', column: 'a'.repeat(64) }; },
    d => { d.columns.value.encoding = 'future'; }, d => { d.columns.value.extra = true; },
  ];
  for (const mutate of mutations) {
    const input = structuredClone(base); mutate(input);
    await assert.rejects(expandPackedSurvey(input, noBundle), /Survey/);
  }
  const badSource = structuredClone(base); badSource.columns.observation_id[0] = null;
  await assert.rejects(expandPackedSurvey(badSource, noBundle), /source must contain nonempty strings/);
});

test('Survey masks and Float64 descriptors fail closed for corruption', async () => {
  const base = await encodePackedSurvey(table([1, null]), noBundle);
  const mutations = [
    d => { d.nulls = [1, 1]; }, d => { d.nulls = [1, 0]; }, d => { d.nulls = [-1]; }, d => { d.nulls = [2]; },
    d => { d.nulls = [0.5]; }, d => { d.nulls = Array(1); }, d => { d.nulls = null; },
    d => { d.nulls = [0]; }, d => { d.values = numeric([1, -0]); }, d => { d.values = numeric([1, Infinity]); },
    d => { d.values = numeric([NaN, 0]); }, d => { d.values.count = 1; }, d => { d.values.data = 'invalid'; },
    d => { d.values.extra = true; },
  ];
  for (const mutate of mutations) {
    const input = structuredClone(base); mutate(input.columns.value);
    await assert.rejects(expandPackedSurvey(input, noBundle), /mask|Masked|coordinate|finite|values descriptor/i);
  }
});

test('Survey missing fields and invalid numeric fallback inputs are explicitly unsupported', async () => {
  for (const value of [undefined, NaN, Infinity, '1', {}, []]) {
    await assert.rejects(encodePackedSurvey(table([value]), noBundle), /Survey|finite/);
    await assert.rejects(expandPackedSurvey({ ...table([value]), encoding: PACKED_SURVEY_ENCODING }, noBundle), /Survey|finite/);
  }
  for (const missing of [{ value: [0] }, { absent: [] }, [], null, { value: '0' }]) {
    const input = { ...table([null]), missing };
    await assert.rejects(encodePackedSurvey(input, noBundle), /missing/);
    await assert.rejects(expandPackedSurvey({ ...input, encoding: PACKED_SURVEY_ENCODING }, noBundle), /missing/);
  }
  const empty = { ...table([null]), missing: { value: [] } };
  assert.deepEqual(await expandPackedSurvey(await encodePackedSurvey(empty, noBundle), noBundle), await expandBundledSurveyColumns(empty, noBundle));
});


test('Derived Survey IDs enforce the complete source construction after dependency resolution', async () => {
  const base = await encodePackedSurvey(table([1, null]), noBundle);
  for (const mutate of [
    d => { d.columns.id.suffix = '|wrong|o2'; },
    d => { d.columns.id.suffix = '|survey|another'; },
    d => { d.columns.term_id = ['o2', 'different']; },
    d => { d.columns.observation_id[0] = ''; },
  ]) {
    const input = structuredClone(base); mutate(input);
    await assert.rejects(expandPackedSurvey(input, noBundle), /suffix|nonempty/);
  }
});
