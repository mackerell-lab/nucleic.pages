import test from 'node:test';
import assert from 'node:assert/strict';
import { COORDINATE_COLUMNAR_ENCODING, decodeCoordinateRows } from '../core/survey-codec.js';
import { PACKED_COORDINATE_ENCODING, FLOAT64_COLUMN_ENCODING,
  encodePackedCoordinates, expandPackedCoordinates } from '../core/packed-coordinate-codec.js';

const table = (columns, count = columns.x?.length ?? 0) => ({ encoding: COORDINATE_COLUMNAR_ENCODING,
  row_count: count, build_id: 'test', columns });
const packed = (column, count = column.count) => ({ encoding: PACKED_COORDINATE_ENCODING, row_count: count, columns: { x: column } });
const descriptor = (base64, count = 1) => ({ encoding: FLOAT64_COLUMN_ENCODING, count, data: base64 });

test('Known independent IEEE-754 bytes establish little endian and byte-lane order', () => {
  // 1.0 = 00 00 00 00 00 00 f0 3f; -2.5 = 00 00 00 00 00 00 04 c0.
  const shuffled = Buffer.from('000000000000000000000000f0043fc0', 'hex').toString('base64');
  const encoded = encodePackedCoordinates(table({ x: [1, -2.5] }));
  assert.deepEqual(encoded.columns.x, descriptor(shuffled, 2));
  assert.deepEqual(expandPackedCoordinates(packed(descriptor(shuffled, 2))).columns.x, [1, -2.5]);
});

test('JSON roundtrip preserves every Float64 bit including signed zero and finite extremes', () => {
  const values = [0, -0, Number.MIN_VALUE, -Number.MIN_VALUE, Number.MAX_VALUE, -Number.MAX_VALUE,
    1.2345678901234567, -1e-250, 1e250, Number.MIN_SAFE_INTEGER, Number.MAX_SAFE_INTEGER];
  const input = { ...table({ x: values, y: [...values].reverse(), z: values.map(value => -value),
    id: values.map((_, index) => `coordinate-${index}`), status: values.map(() => 'available') }),
    missing: { status: [2] }, provenance: { axis_unit: 'Å' } };
  const restored = expandPackedCoordinates(JSON.parse(JSON.stringify(encodePackedCoordinates(input))));
  for (const axis of ['x', 'y', 'z']) for (let index = 0; index < values.length; index++) {
    assert(Object.is(restored.columns[axis][index], input.columns[axis][index]), `${axis}/${index} lost Float64 identity`);
  }
  assert.deepEqual(restored, input);
  assert.deepEqual(decodeCoordinateRows(restored), decodeCoordinateRows(input));
});

test('Null and absent coordinates use unchanged array fallbacks with missing maps intact', () => {
  const input = { ...table({ x: [1, null, 3], y: [4, 5, 6], z: [7, 8, 9], id: ['a', 'b', 'c'] }), missing: { y: [1] } };
  const encoded = encodePackedCoordinates(input);
  assert.equal(encoded.columns.x, input.columns.x);
  assert.equal(encoded.columns.y, input.columns.y);
  assert.equal(encoded.columns.z.encoding, FLOAT64_COLUMN_ENCODING);
  const restored = expandPackedCoordinates(JSON.parse(JSON.stringify(encoded)));
  assert.deepEqual(restored, input);
  const rows = decodeCoordinateRows(restored);
  assert.equal(rows[1].x, null);
  assert(!Object.hasOwn(rows[1], 'y'));
});

test('Sparse columns cannot silently become null during JSON transport', () => {
  for (const columns of [{ x: [1, , 3] }, { x: [1, 2, 3], id: ['a', , 'c'] }]) {
    const source = table(columns, 3);
    assert.throws(() => encodePackedCoordinates(source), /Sparse coordinate transport column/);
    assert.throws(() => expandPackedCoordinates({ ...source, encoding: PACKED_COORDINATE_ENCODING }), /Sparse coordinate transport column/);
  }
});

test('Chunked base64 supports a complete 10000-row partition and empty columns', () => {
  const x = Array.from({ length: 10000 }, (_, index) => (index - 5000) * Math.PI);
  const source = table({ x });
  const result = expandPackedCoordinates(encodePackedCoordinates(source));
  assert.deepEqual(result, source);
  const independent = Buffer.alloc(8 * x.length);
  for (let index = 0; index < x.length; index++) independent.writeDoubleLE(x[index], index * 8);
  const bytes = Buffer.from(encodePackedCoordinates(source).columns.x.data, 'base64');
  for (let lane = 0; lane < 8; lane++) for (let index = 0; index < x.length; index++) assert.equal(bytes[lane * x.length + index], independent[index * 8 + lane]);
  assert.deepEqual(expandPackedCoordinates(encodePackedCoordinates(table({ x: [], y: [], z: [] }))), table({ x: [], y: [], z: [] }));
});

test('Malformed descriptors, base64, row counts and unsupported tags fail closed', () => {
  const valid = descriptor('AAAAAAAA8D8=');
  for (const invalid of [null, {}, { ...valid, count: -1 }, { ...valid, count: 0.5 }, { ...valid, count: 2 },
    { ...valid, count: 2 ** 32 }, { ...valid, encoding: 'float32' }, { ...valid, extra: true },
    { ...valid, data: '' }, { ...valid, data: 'AAAAAAAA8D8' }, { ...valid, data: 'AAAAAAAA8D8=\n' },
    { ...valid, data: 'AAAAAAAA8D8-' }, { ...valid, data: 'AAAAAAAA8D9=' }, { ...valid, data: 123 }]) {
    assert.throws(() => expandPackedCoordinates(packed(invalid, 1)), /coordinate|Coordinate/);
  }
  for (const row_count of [-1, 0.5, null, Number.MAX_SAFE_INTEGER, 2 ** 32]) {
    assert.throws(() => expandPackedCoordinates({ ...packed(valid), row_count }), /row count/);
  }
  assert.throws(() => expandPackedCoordinates({ ...packed(valid), columns: { status: valid } }), /Only x\/y\/z/);
  assert.throws(() => expandPackedCoordinates({ ...packed(valid), encoding: COORDINATE_COLUMNAR_ENCODING }), /encoding/);
  assert.throws(() => encodePackedCoordinates(packed(valid)), /encoding/);
});

test('Nonfinite binary values and nonnumeric coordinate fallbacks are rejected', () => {
  // Independent LE bytes for +Infinity, -Infinity, quiet NaN.
  for (const bytes of ['000000000000f07f', '000000000000f0ff', '000000000000f87f']) {
    assert.throws(() => expandPackedCoordinates(packed(descriptor(Buffer.from(bytes, 'hex').toString('base64')))), /finite numbers/);
  }
  for (const value of [NaN, Infinity, -Infinity, '1', false, undefined, {}]) {
    assert.throws(() => encodePackedCoordinates(table({ x: [value] })), /finite numbers/);
    assert.throws(() => expandPackedCoordinates({ ...packed(descriptor('', 1)), columns: { x: [value] } }), /finite numbers/);
  }
});

test('All columns and missing maps are validated even if coordinate data itself is valid', () => {
  const input = table({ x: [1], id: ['a'] });
  const encoded = encodePackedCoordinates(input);
  assert.throws(() => encodePackedCoordinates({ ...input, columns: { ...input.columns, id: [] } }), /length mismatch/);
  assert.throws(() => expandPackedCoordinates({ ...encoded, columns: { ...encoded.columns, id: [] } }), /length mismatch/);
  for (const missing of [null, [], { absent: [0] }, { x: [1] }, { x: [-1] }, { x: [0.5] }, { x: null }, { x: Array(1) }]) {
    assert.throws(() => encodePackedCoordinates({ ...input, missing }), /missing-field/);
    assert.throws(() => expandPackedCoordinates({ ...encoded, missing }), /missing-field/);
  }
  for (const columns of [null, [], 'text']) {
    assert.throws(() => encodePackedCoordinates({ ...input, columns }), /columns/);
    assert.throws(() => expandPackedCoordinates({ ...encoded, columns }), /columns/);
  }
});
