import { COORDINATE_COLUMNAR_ENCODING } from './survey-codec.js';

export const PACKED_COORDINATE_ENCODING = 'rna-coordinate-float64-1';
export const FLOAT64_COLUMN_ENCODING = 'float64-le-shuffled-base64-1';

const axes = new Set(['x', 'y', 'z']);
const record = value => value !== null && typeof value === 'object' && !Array.isArray(value);
const countIsValid = value => Number.isSafeInteger(value) && value >= 0 && value <= 0xffffffff;
const define = (object, key, value) => Object.defineProperty(object, key, { value, enumerable: true });

function validateTable(data, encoding) {
  if (data?.encoding !== encoding) throw new Error('Unsupported packed coordinate encoding');
  if (!countIsValid(data.row_count)) throw new Error('Invalid packed coordinate row count');
  if (!record(data.columns)) throw new Error('Invalid packed coordinate columns');
  if (Object.hasOwn(data, 'missing')) {
    if (!record(data.missing)) throw new Error('Invalid coordinate missing-field map');
    for (const [key, indices] of Object.entries(data.missing)) {
      if (!Object.hasOwn(data.columns, key) || !Array.isArray(indices)) {
        throw new Error('Invalid coordinate missing-field index');
      }
      for (const index of indices) if (!Number.isSafeInteger(index) || index < 0 || index >= data.row_count) {
        throw new Error('Invalid coordinate missing-field index');
      }
    }
  }
}

/** Explicit nulls use array fallback; missing fields use the separate map. */
function validateArray(values, count, coordinate) {
  if (!Array.isArray(values) || values.length !== count) throw new Error('Coordinate column length mismatch');
  let eligible = true;
  for (let index = 0; index < count; index++) {
    if (!Object.hasOwn(values, index)) throw new Error('Sparse coordinate transport column cannot preserve JSON identity');
    if (coordinate) {
      if (values[index] === null) { eligible = false; continue; }
      if (typeof values[index] !== 'number' || !Number.isFinite(values[index])) throw new Error('Coordinate values must be finite numbers or null');
    }
  }
  return eligible;
}

function pack(values) {
  const count = values.length, bytes = new Uint8Array(count * 8), view = new DataView(bytes.buffer);
  for (let index = 0; index < count; index++) view.setFloat64(index * 8, values[index], true);
  const shuffled = new Uint8Array(bytes.length);
  for (let lane = 0; lane < 8; lane++) for (let index = 0; index < count; index++) shuffled[lane * count + index] = bytes[index * 8 + lane];
  const chunks = [];
  for (let offset = 0; offset < shuffled.length; offset += 0x8000) {
    chunks.push(String.fromCharCode(...shuffled.subarray(offset, offset + 0x8000)));
  }
  return { encoding: FLOAT64_COLUMN_ENCODING, count, data: btoa(chunks.join('')) };
}

function unpack(descriptor, count) {
  if (!record(descriptor) || Object.keys(descriptor).length !== 3
      || !['encoding', 'count', 'data'].every(key => Object.hasOwn(descriptor, key))
      || descriptor.encoding !== FLOAT64_COLUMN_ENCODING || !countIsValid(descriptor.count) || descriptor.count !== count) {
    throw new Error('Invalid packed coordinate descriptor');
  }
  const byteLength = count * 8;
  if (!Number.isSafeInteger(byteLength) || typeof descriptor.data !== 'string'
      || descriptor.data.length !== 4 * Math.ceil(byteLength / 3)
      || !/^[A-Za-z0-9+/]*={0,2}$/.test(descriptor.data)) throw new Error('Invalid packed coordinate base64 length or characters');
  let binary;
  try { binary = atob(descriptor.data); } catch { throw new Error('Invalid packed coordinate base64'); }
  if (binary.length !== byteLength || btoa(binary) !== descriptor.data) throw new Error('Noncanonical packed coordinate base64 or byte length');
  const bytes = new Uint8Array(byteLength), view = new DataView(bytes.buffer), values = new Array(count);
  for (let lane = 0; lane < 8; lane++) for (let index = 0; index < count; index++) bytes[index * 8 + lane] = binary.charCodeAt(lane * count + index);
  for (let index = 0; index < count; index++) {
    const value = view.getFloat64(index * 8, true);
    if (!Number.isFinite(value)) throw new Error('Packed coordinates must decode to finite numbers');
    values[index] = value;
  }
  return values;
}

export function encodePackedCoordinates(data) {
  validateTable(data, COORDINATE_COLUMNAR_ENCODING);
  // Finish all validation before allocating binary columns.
  const eligible = new Set();
  for (const [key, values] of Object.entries(data.columns)) {
    if (validateArray(values, data.row_count, axes.has(key)) && axes.has(key) && !(data.missing?.[key]?.length)) eligible.add(key);
  }
  const columns = {};
  for (const [key, values] of Object.entries(data.columns)) define(columns, key, eligible.has(key) ? pack(values) : values);
  return { ...data, encoding: PACKED_COORDINATE_ENCODING, columns };
}

export function expandPackedCoordinates(data) {
  validateTable(data, PACKED_COORDINATE_ENCODING);
  const columns = {};
  for (const [key, descriptor] of Object.entries(data.columns)) {
    let values = descriptor;
    if (!Array.isArray(descriptor)) {
      if (!axes.has(key)) throw new Error('Only x/y/z coordinate columns may be packed');
      values = unpack(descriptor, data.row_count);
    }
    validateArray(values, data.row_count, axes.has(key));
    define(columns, key, values);
  }
  return { ...data, encoding: COORDINATE_COLUMNAR_ENCODING, columns };
}
