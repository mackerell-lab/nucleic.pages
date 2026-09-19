import { BUNDLED_FAMILY_ENCODING } from './bundled-family-codec.js';
import { encodeFloat64Column, decodeFloat64Column } from './packed-coordinate-codec.js';

export const PACKED_FAMILY_ENCODING = 'rna-family-float64-1';
export const PACKED_FAMILY_VALUES_ENCODING = 'rna-family-values-float64-1';

const record = value => value !== null && typeof value === 'object' && !Array.isArray(value)
  && [Object.prototype, null].includes(Object.getPrototypeOf(value));
const exact = (value, keys) => record(value) && Object.keys(value).length === keys.length && keys.every(key => Object.hasOwn(value, key));
const hash = value => typeof value === 'string' && /^[a-f0-9]{64}$/.test(value);
const define = (object, key, value) => Object.defineProperty(object, key, { value, enumerable: true });

function dense(values, count) {
  if (!Array.isArray(values) || values.length !== count) throw new Error('Packed family column length mismatch');
  for (let index = 0; index < count; index++) if (!Object.hasOwn(values, index)) throw new Error('Sparse packed family column');
}

function validateValueRows(rows) {
  for (const row of rows) {
    if (!record(row)) throw new Error('Family values require record objects');
    for (const value of Object.values(row)) if (value !== null && (typeof value !== 'number' || !Number.isFinite(value))) {
      throw new Error('Family values require finite numbers or null');
    }
  }
}

function validateFallback(column, count, isValues) {
  if (Array.isArray(column)) {
    dense(column, count);
    if (isValues) validateValueRows(column);
    return;
  }
  if (exact(column, ['bundle', 'column']) && hash(column.bundle) && hash(column.column)) return;
  if (exact(column, ['dictionary', 'indices']) && Array.isArray(column.dictionary)) {
    dense(column.dictionary, column.dictionary.length); dense(column.indices, count);
    if (isValues) validateValueRows(column.dictionary);
    for (const index of column.indices) if (!Number.isSafeInteger(index) || index < 0 || index >= column.dictionary.length) {
      throw new Error('Invalid packed family dictionary index');
    }
    return;
  }
  throw new Error('Invalid packed family column descriptor');
}

function validateTable(data, encoding) {
  if (data?.encoding !== encoding) throw new Error('Unsupported packed family encoding');
  if (!Number.isSafeInteger(data.row_count) || data.row_count < 0 || data.row_count > 0xffffffff) throw new Error('Invalid packed family row count');
  if (!record(data.columns)) throw new Error('Invalid packed family columns');
  if (Object.hasOwn(data, 'missing')) {
    if (!record(data.missing)) throw new Error('Invalid family missing-field map');
    for (const [key, indices] of Object.entries(data.missing)) {
      if (!Object.hasOwn(data.columns, key) || !Array.isArray(indices)) throw new Error('Invalid family missing-field index');
      for (const index of indices) if (!Number.isSafeInteger(index) || index < 0 || index >= data.row_count) throw new Error('Invalid family missing-field index');
    }
  }
}

function packValues(rows, count) {
  const parameters = {}, names = [...new Set(rows.flatMap(row => Object.keys(row)))];
  for (const name of names) {
    const values = new Array(count).fill(0), nulls = [], missing = [];
    for (let index = 0; index < count; index++) {
      if (!Object.hasOwn(rows[index], name)) missing.push(index);
      else if (rows[index][name] === null) nulls.push(index);
      else values[index] = rows[index][name];
    }
    define(parameters, name, { values: encodeFloat64Column(values), nulls, missing });
  }
  return { encoding: PACKED_FAMILY_VALUES_ENCODING, parameters };
}

function validateMask(mask, count) {
  if (!Array.isArray(mask)) throw new Error('Invalid family values mask');
  let previous = -1;
  for (const index of mask) {
    if (!Number.isSafeInteger(index) || index <= previous || index >= count) throw new Error('Invalid family values mask index');
    previous = index;
  }
}

function expandValues(descriptor, count) {
  if (!exact(descriptor, ['encoding', 'parameters']) || descriptor.encoding !== PACKED_FAMILY_VALUES_ENCODING || !record(descriptor.parameters)) {
    throw new Error('Invalid packed family values descriptor');
  }
  const parameters = [];
  for (const [name, column] of Object.entries(descriptor.parameters)) {
    if (!exact(column, ['values', 'nulls', 'missing'])) throw new Error('Invalid packed family parameter descriptor');
    validateMask(column.nulls, count); validateMask(column.missing, count);
    const nulls = new Set(column.nulls), missing = new Set(column.missing);
    if (column.missing.some(index => nulls.has(index))) throw new Error('Overlapping family values masks');
    const values = decodeFloat64Column(column.values, count);
    for (const index of [...column.nulls, ...column.missing]) {
      if (!Object.is(values[index], 0)) throw new Error('Masked family numeric slots must be zero');
    }
    parameters.push({ name, values, nulls, missing });
  }
  const rows = Array.from({ length: count }, () => ({}));
  for (const { name, values, nulls, missing } of parameters) for (let index = 0; index < count; index++) {
    if (!missing.has(index)) define(rows[index], name, nulls.has(index) ? null : values[index]);
  }
  return rows;
}

export function encodePackedFamily(data) {
  validateTable(data, BUNDLED_FAMILY_ENCODING);
  for (const [key, column] of Object.entries(data.columns)) validateFallback(column, data.row_count, key === 'values');
  const columns = {};
  for (const [key, column] of Object.entries(data.columns)) define(columns, key,
    key === 'values' && Array.isArray(column) ? packValues(column, data.row_count) : column);
  return { ...data, encoding: PACKED_FAMILY_ENCODING, columns };
}

export function expandPackedFamily(data) {
  validateTable(data, PACKED_FAMILY_ENCODING);
  // Check every non-packed column before decoding any numeric arrays.
  for (const [key, column] of Object.entries(data.columns)) {
    if (key === 'values' && record(column) && Object.hasOwn(column, 'encoding')) continue;
    validateFallback(column, data.row_count, key === 'values');
  }
  const columns = {};
  for (const [key, column] of Object.entries(data.columns)) define(columns, key,
    key === 'values' && record(column) && Object.hasOwn(column, 'encoding') ? expandValues(column, data.row_count) : column);
  return { ...data, encoding: BUNDLED_FAMILY_ENCODING, columns };
}
