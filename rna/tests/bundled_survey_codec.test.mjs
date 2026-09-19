import test from 'node:test';
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { decodeSurveyRows, SURVEY_COLUMNAR_ENCODING } from '../core/survey-codec.js';
import { BUNDLED_SURVEY_ENCODING, SURVEY_BUNDLE_ENCODING, verifySurveyBundle,
  expandBundledSurveyColumns } from '../core/bundled-survey-codec.js';

const hash = value => createHash('sha256').update(JSON.stringify(value)).digest('hex');
function bundleFor(...arrays) {
  const bundle = { encoding: SURVEY_BUNDLE_ENCODING,
    columns: Object.fromEntries(arrays.map(values => [hash(values), values]).sort(([a], [b]) => a.localeCompare(b))) };
  return { bundle, reference: hash(bundle), ref: values => ({ bundle: hash(bundle), column: hash(values) }) };
}
const table = columns => ({ encoding: BUNDLED_SURVEY_ENCODING, row_count: 2, build_id: 'fixture', columns });

test('Bundle verification authenticates values and names and reports UTF8 bytes', async () => {
  const { bundle, reference } = bundleFor(['α', 'RNA'], [null, { nested: [3, null] }]);
  assert.deepEqual(await verifySurveyBundle(reference, bundle), {
    bundle, byteLength: Buffer.byteLength(JSON.stringify(bundle))
  });
  const corrupt = structuredClone(bundle);
  Object.values(corrupt.columns)[0][0] = 'changed';
  await assert.rejects(verifySurveyBundle(reference, corrupt), /hash mismatch/);
  const renamed = { ...bundle, columns: { ['f'.repeat(64)]: Object.values(bundle.columns)[0] } };
  await assert.rejects(verifySurveyBundle(reference, renamed), /hash mismatch/);
  for (const malformed of [null, [], {}, { encoding: SURVEY_BUNDLE_ENCODING, columns: [] },
    { encoding: SURVEY_BUNDLE_ENCODING, columns: { wrong: [] } },
    { encoding: SURVEY_BUNDLE_ENCODING, columns: { ['a'.repeat(64)]: null } }]) {
    await assert.rejects(verifySurveyBundle(hash(malformed), malformed), /Invalid Survey/);
  }
  await assert.rejects(verifySurveyBundle('../escape', bundle), /content hash/);
});

test('Mixed bundled and inline columns preserve precise numbers, nested values and metadata', async () => {
  const ids = ['a', 'b'], nested = [{ atoms: ['N1', null] }, null];
  const fixture = bundleFor(ids, nested);
  const input = table({ id: fixture.ref(ids), residue_id: fixture.ref(ids), nested: fixture.ref(nested),
    value: [1.2345678901234567, -0.0000000000000023] });
  input.provenance = { build: 'test' };
  let reads = 0;
  const expanded = await expandBundledSurveyColumns(JSON.parse(JSON.stringify(input)), async reference => {
    reads++; return (await verifySurveyBundle(reference, fixture.bundle)).bundle;
  });
  assert.equal(reads, 1);
  assert.equal(expanded.encoding, SURVEY_COLUMNAR_ENCODING);
  assert.equal(expanded.build_id, input.build_id);
  assert.deepEqual(expanded.provenance, input.provenance);
  assert.deepEqual(decodeSurveyRows(expanded), [
    { id: 'a', residue_id: 'a', nested: nested[0], value: input.columns.value[0] },
    { id: 'b', residue_id: 'b', nested: null, value: input.columns.value[1] }
  ]);
  assert.deepEqual(decodeSurveyRows(expanded, ['id']), [{ id: 'a' }, { id: 'b' }]);
});

test('All descriptors and inline lengths fail before network reads', async () => {
  const fixture = bundleFor(['a', 'b']);
  let reads = 0;
  for (const invalid of [null, 'text', 0, {}, { bundle: fixture.reference },
    { ...fixture.ref(['a', 'b']), extra: true }, { bundle: '../escape', column: hash([]) },
    { bundle: fixture.reference, column: '../escape' }, [1]]) {
    await assert.rejects(expandBundledSurveyColumns(table({ id: fixture.ref(['a', 'b']), hidden: invalid }),
      async () => { reads++; return fixture.bundle; }), /reference|length mismatch/);
  }
  for (const row_count of [-1, 1.5, Number.MAX_SAFE_INTEGER + 1, null]) {
    await assert.rejects(expandBundledSurveyColumns({ ...table({}), row_count }, async () => {}), /row count/);
  }
  for (const columns of [null, [], 'text']) {
    await assert.rejects(expandBundledSurveyColumns(table(columns), async () => {}), /columns/);
  }
  assert.equal(reads, 0);
});

test('Referenced lengths and missing columns fail even when projected away', async () => {
  const fixture = bundleFor(['a', 'b'], [1]);
  await assert.rejects(expandBundledSurveyColumns(table({ id: fixture.ref(['a', 'b']), value: fixture.ref([1]) }),
    async () => fixture.bundle), /length mismatch/);
  await assert.rejects(expandBundledSurveyColumns(table({ id: fixture.ref(['absent', 'column']) }),
    async () => fixture.bundle), /Missing Survey bundle column/);
});

test('Distinct bundle reads overlap at most four and are never persistently cached', async () => {
  const fixtures = Array.from({ length: 9 }, (_, i) => bundleFor([i, null]));
  const assets = new Map(fixtures.map(item => [item.reference, item.bundle]));
  const input = table(Object.fromEntries(fixtures.map((item, i) => [`field${i}`, item.ref([i, null])])));
  let active = 0, maximum = 0, reads = 0;
  const read = async reference => {
    reads++; active++; maximum = Math.max(maximum, active);
    await new Promise(resolve => setTimeout(resolve, 1));
    active--; return assets.get(reference);
  };
  const expanded = await expandBundledSurveyColumns(input, read);
  assert.equal(maximum, 4);
  assert.equal(active, 0);
  assert.equal(reads, 9);
  assert.deepEqual(expanded.columns.field8, [8, null]);
  await expandBundledSurveyColumns(input, read);
  assert.equal(reads, 18);
  await assert.rejects(expandBundledSurveyColumns(input, async () => { throw new Error('HTTP failure'); }), /HTTP failure/);
  assert.deepEqual(await expandBundledSurveyColumns(input, read), expanded);
});

test('Prototype-like field names stay ordinary own columns', async () => {
  const fixture = bundleFor(['a', 'b']);
  const input = table(JSON.parse(JSON.stringify(Object.fromEntries([
    ['__proto__', fixture.ref(['a', 'b'])], ['constructor', [null, { nested: true }]]
  ]))));
  const expanded = await expandBundledSurveyColumns(input, async () => fixture.bundle);
  assert.equal(Object.getPrototypeOf(expanded.columns), Object.prototype);
  assert(Object.hasOwn(expanded.columns, '__proto__'));
  assert.deepEqual(expanded.columns.__proto__, ['a', 'b']);
  assert.deepEqual(decodeSurveyRows(expanded)[0], JSON.parse('{"__proto__":"a","constructor":null}'));
});
