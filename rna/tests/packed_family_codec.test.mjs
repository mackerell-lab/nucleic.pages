import test from 'node:test';
import assert from 'node:assert/strict';
import { BUNDLED_FAMILY_ENCODING, expandBundledFamilyColumns } from '../core/bundled-family-codec.js';
import { decodeFamilyRows } from '../core/survey-codec.js';
import { encodeFloat64Column } from '../core/packed-coordinate-codec.js';
import { PACKED_FAMILY_ENCODING, PACKED_FAMILY_VALUES_ENCODING, encodePackedFamily, expandPackedFamily } from '../core/packed-family-codec.js';

const table = values => ({ encoding: BUNDLED_FAMILY_ENCODING, build_id: 'test', row_count: values.length,
  columns: { id: values.map((_, index) => `residue-${index}`), values, statuses: values.map(() => ({ chi: 'available' })) } });
const roundtrip = data => expandPackedFamily(JSON.parse(JSON.stringify(encodePackedFamily(data))));

test('Packed family values preserve Float64 identity, absent keys and explicit nulls', async () => {
  const values = [0, -0, Number.MIN_VALUE, -Number.MIN_VALUE, Number.MAX_VALUE, -Number.MAX_VALUE, 1.2345678901234567];
  const rows = values.map((value, index) => ({ chi: value, ...(index % 2 ? {} : { alpha: index ? null : 0.25 }) }));
  const input = { ...table(rows), missing: { statuses: [2] }, provenance: { source: 'fixture' } };
  const expanded = roundtrip(input);
  assert.deepEqual(expanded, input);
  for (let index = 0; index < rows.length; index++) assert(Object.is(expanded.columns.values[index].chi, values[index]));
  assert.equal(expanded.columns.values[2].alpha, null);
  assert(!Object.hasOwn(expanded.columns.values[1], 'alpha'));
  assert.deepEqual(decodeFamilyRows(await expandBundledFamilyColumns(expanded, () => assert.fail('No bundle expected'))),
    decodeFamilyRows(await expandBundledFamilyColumns(input, () => assert.fail('No bundle expected'))));
});

test('Independent IEEE bytes and Node Float64 writer verify parameter lane order', () => {
  const encoded = encodePackedFamily(table([{ chi: 1 }, { chi: -2.5 }]));
  assert.equal(encoded.columns.values.parameters.chi.values.data, Buffer.from('000000000000000000000000f0043fc0', 'hex').toString('base64'));
  const numbers = [Number.MIN_VALUE, -0, 1e250, -1.25, 1.2345678901234567];
  const raw = Buffer.alloc(numbers.length * 8), shuffled = Buffer.alloc(raw.length);
  numbers.forEach((value, index) => raw.writeDoubleLE(value, index * 8));
  for (let lane = 0; lane < 8; lane++) for (let index = 0; index < numbers.length; index++) shuffled[lane * numbers.length + index] = raw[index * 8 + lane];
  const data = { ...table(numbers.map(value => ({ chi: value }))), encoding: PACKED_FAMILY_ENCODING };
  data.columns.values = { encoding: PACKED_FAMILY_VALUES_ENCODING, parameters: {
    chi: { values: { encoding: 'float64-le-shuffled-base64-1', count: numbers.length, data: shuffled.toString('base64') }, nulls: [], missing: [] }
  } };
  expandPackedFamily(data).columns.values.forEach((row, index) => assert(Object.is(row.chi, numbers[index])));
});

test('Empty objects, safe prototype-like parameters and top-level missing values survive', () => {
  const rows = [JSON.parse('{"__proto__":1.25,"constructor":null}'), {}, JSON.parse('{"__proto__":null}')];
  const input = { ...table(rows), missing: { values: [1] } }, output = roundtrip(input);
  assert.deepEqual(output, input);
  assert.equal(Object.getPrototypeOf(output.columns.values[0]), Object.prototype);
  assert(Object.hasOwn(output.columns.values[0], '__proto__'));
  assert.deepEqual(roundtrip(table([{}, {}])), table([{}, {}]));
  assert.deepEqual(roundtrip(table([])), table([]));
});

test('Validated plain arrays, dictionaries and bundle refs remain supported fallback forms', () => {
  const dict = { dictionary: [{ chi: null }, { chi: 3 }], indices: [0, 1] };
  const ref = { bundle: 'a'.repeat(64), column: 'b'.repeat(64) };
  for (const column of [dict, ref]) {
    const input = table([{}, {}]); input.columns.values = column;
    const packed = encodePackedFamily(input);
    assert.equal(packed.columns.values, column);
    assert.deepEqual(expandPackedFamily(packed), input);
  }
  const plain = table([{ chi: 1 }, { chi: null }]);
  assert.deepEqual(expandPackedFamily({ ...plain, encoding: PACKED_FAMILY_ENCODING }), plain);
  const input = table([{}, {}]); input.columns.id = ref;
  assert.equal(encodePackedFamily(input).columns.id, ref);
});

test('Family parameter order follows the source and dictionary values retain numeric validation', () => {
  const input = table([{ gamma: 1, alpha: 2, beta: 3 }, { gamma: null, alpha: 5, beta: 6 }]);
  assert.equal(JSON.stringify(roundtrip(input)), JSON.stringify(input));
  for (const value of [null, [], 1, { chi: undefined }, { chi: NaN }, { chi: Infinity }, { chi: '1' }]) {
    const data = table([{}]); data.columns.values = { dictionary: [value], indices: [0] };
    assert.throws(() => encodePackedFamily(data), /Family values/);
    assert.throws(() => expandPackedFamily({ ...data, encoding: PACKED_FAMILY_ENCODING }), /Family values/);
  }
});

test('Eligible values reject sparse rows, undefined, nonfinite and nonnumeric parameters', () => {
  for (const rows of [[null], [undefined], [1], [[]], [new Date()], Array(1), [{ chi: undefined }],
    [{ chi: NaN }], [{ chi: Infinity }], [{ chi: -Infinity }], [{ chi: '1' }], [{ chi: {} }]]) {
    assert.throws(() => encodePackedFamily(table(rows)), /Family values|Sparse packed family/);
    assert.throws(() => expandPackedFamily({ ...table(rows), encoding: PACKED_FAMILY_ENCODING }), /Family values|Sparse packed family/);
  }
});

test('Packed masks require ordered disjoint in-range indices and zero numeric placeholders', () => {
  const base = encodePackedFamily(table([{ chi: 1 }, { chi: null }, {}]));
  const mutations = [
    c => { c.nulls = [1, 1]; }, c => { c.missing = [2, 1]; }, c => { c.nulls = [-1]; },
    c => { c.nulls = [3]; }, c => { c.nulls = [0.5]; }, c => { c.nulls = Array(1); },
    c => { c.missing = [1]; }, c => { c.nulls = null; }, c => { c.nulls = [0]; },
    c => { c.values = encodeFloat64Column([1, -0, 0]); }, c => { c.extra = true; },
    c => { c.values.count = 2; }, c => { c.values.data = 'malformed'; },
  ];
  for (const mutate of mutations) {
    const input = structuredClone(base); mutate(input.columns.values.parameters.chi);
    assert.throws(() => expandPackedFamily(input), /mask|Masked|parameter descriptor|coordinate/i);
  }
});

test('Whole table validation rejects malformed unrelated fields and misplaced packed descriptors', () => {
  const base = encodePackedFamily(table([{ chi: 1 }]));
  const mutations = [
    d => { d.columns.id = []; }, d => { d.columns.id = Array(1); },
    d => { d.columns.id = d.columns.values; }, d => { d.columns.id = { bundle: '../escape', column: 'a'.repeat(64) }; },
    d => { d.columns.id = { dictionary: ['a'], indices: [1] }; },
    d => { d.columns.id = { dictionary: ['a'], indices: Array(1) }; },
    d => { d.columns.values.encoding = 'other'; }, d => { d.columns.values.extra = true; },
    d => { d.columns.values.parameters = []; }, d => { d.row_count = -1; },
    d => { d.missing = { id: [1] }; }, d => { d.missing = { absent: [0] }; }, d => { d.missing = { id: Array(1) }; },
  ];
  for (const mutate of mutations) {
    const input = structuredClone(base); mutate(input);
    assert.throws(() => expandPackedFamily(input), /family|Family/);
  }
  assert.throws(() => encodePackedFamily(base), /encoding/);
  assert.throws(() => expandPackedFamily(table([{}])), /encoding/);
});
