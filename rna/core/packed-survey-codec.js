import { BUNDLED_SURVEY_ENCODING, expandBundledSurveyColumns } from './bundled-survey-codec.js';
import { encodeFloat64Column, decodeFloat64Column } from './packed-coordinate-codec.js';

export const PACKED_SURVEY_ENCODING = 'rna-survey-float64-1';
export const PACKED_SURVEY_VALUES_ENCODING = 'rna-survey-values-float64-1';
export const PACKED_SURVEY_ID_ENCODING = 'rna-survey-id-suffix-1';
export const PACKED_SURVEY_CONSTANT_ENCODING = 'rna-survey-constant-1';

const record = value => value !== null && typeof value === 'object' && !Array.isArray(value)
  && [Object.prototype, null].includes(Object.getPrototypeOf(value));
const exact = (value, keys) => record(value) && Object.keys(value).length === keys.length && keys.every(key => Object.hasOwn(value, key));
const hash = value => typeof value === 'string' && /^[a-f0-9]{64}$/.test(value);
const define = (object, key, value) => Object.defineProperty(object, key, { value, enumerable: true });

function dense(values, count) {
  if (!Array.isArray(values) || values.length !== count) throw new Error('Packed Survey column length mismatch');
  for (let index = 0; index < count; index++) if (!Object.hasOwn(values, index) || values[index] === undefined) {
    throw new Error('Sparse or undefined packed Survey column');
  }
}

function validateTable(data, encoding) {
  if (data?.encoding !== encoding) throw new Error('Unsupported packed Survey encoding');
  if (!Number.isSafeInteger(data.row_count) || data.row_count < 0 || data.row_count > 0xffffffff) throw new Error('Invalid packed Survey row count');
  if (!record(data.columns)) throw new Error('Invalid packed Survey columns');
  if (Object.hasOwn(data, 'missing')) {
    if (!record(data.missing)) throw new Error('Invalid packed Survey missing-field map');
    for (const [key, indices] of Object.entries(data.missing)) {
      if (!Object.hasOwn(data.columns, key) || !Array.isArray(indices) || indices.length) {
        // The existing scalar row decoder has no absent-field representation.
        throw new Error('Packed Survey does not support missing fields');
      }
    }
  }
}

function validateValues(values) {
  for (const value of values) if (value !== null && (typeof value !== 'number' || !Number.isFinite(value))) {
    throw new Error('Survey values require finite numbers or null');
  }
}

function validateFallback(column, count, key) {
  if (Array.isArray(column)) {
    dense(column, count);
    if (key === 'value') validateValues(column);
  } else if (!(exact(column, ['bundle', 'column']) && hash(column.bundle) && hash(column.column))) {
    throw new Error('Invalid packed Survey column descriptor');
  }
}

function expandValues(column, count) {
  if (!exact(column, ['encoding', 'values', 'nulls']) || column.encoding !== PACKED_SURVEY_VALUES_ENCODING
      || !Array.isArray(column.nulls)) throw new Error('Invalid packed Survey values descriptor');
  let previous = -1;
  for (const index of column.nulls) {
    if (!Number.isSafeInteger(index) || index <= previous || index >= count) throw new Error('Invalid packed Survey null mask index');
    previous = index;
  }
  const values = decodeFloat64Column(column.values, count);
  for (const index of column.nulls) {
    if (!Object.is(values[index], 0)) throw new Error('Masked Survey numeric slots must be positive zero');
    values[index] = null;
  }
  return values;
}

function validateResolved(data) {
  for (const [key, values] of Object.entries(data.columns)) {
    dense(values, data.row_count);
    if (key === 'value') validateValues(values);
  }
  return data;
}

/** readBundle must authenticate the existing content-addressed Survey bundles. */
export async function encodePackedSurvey(data, readBundle) {
  validateTable(data, BUNDLED_SURVEY_ENCODING);
  for (const [key, column] of Object.entries(data.columns)) validateFallback(column, data.row_count, key);
  const resolved = validateResolved(await expandBundledSurveyColumns(data, readBundle));
  const columns = {}, term = resolved.columns.term_id?.[0];
  const constantTerm = data.row_count > 0 && typeof term === 'string' && term.length > 0
    && resolved.columns.term_id.every(value => value === term);
  const suffix = constantTerm ? `|survey|${term}` : null;
  const suffixIds = constantTerm && Array.isArray(resolved.columns.id) && Array.isArray(resolved.columns.observation_id)
    && resolved.columns.observation_id.every((value, index) => typeof value === 'string' && value.length > 0 && resolved.columns.id[index] === value + suffix);
  for (const [key, column] of Object.entries(data.columns)) {
    let packed = column;
    if (key === 'value') {
      const nulls = [], values = resolved.columns.value.map((value, index) => {
        if (value === null) { nulls.push(index); return 0; }
        return value;
      });
      packed = { encoding: PACKED_SURVEY_VALUES_ENCODING, values: encodeFloat64Column(values), nulls };
    } else if (key === 'id' && suffixIds) packed = { encoding: PACKED_SURVEY_ID_ENCODING, column: 'observation_id', suffix };
    else if (key === 'term_id' && constantTerm) packed = { encoding: PACKED_SURVEY_CONSTANT_ENCODING, value: term };
    define(columns, key, packed);
  }
  return { ...data, encoding: PACKED_SURVEY_ENCODING, columns };
}

/** Resolve all authenticated columns before deriving IDs or allowing projection. */
export async function expandPackedSurvey(data, readBundle) {
  validateTable(data, PACKED_SURVEY_ENCODING);
  const columns = {}, deferred = new Map();
  // Validate all local descriptors before any bundle request. Derived IDs have
  // exactly one permitted source; recursive expressions and cycles are invalid.
  for (const [key, column] of Object.entries(data.columns)) {
    if (record(column) && Object.hasOwn(column, 'encoding')) {
      if (key === 'value') define(columns, key, expandValues(column, data.row_count));
      else if (key === 'term_id' && exact(column, ['encoding', 'value'])
          && column.encoding === PACKED_SURVEY_CONSTANT_ENCODING && typeof column.value === 'string' && column.value.length > 0) {
        define(columns, key, Array(data.row_count).fill(column.value));
      } else if (key === 'id' && exact(column, ['encoding', 'column', 'suffix'])
          && column.encoding === PACKED_SURVEY_ID_ENCODING && column.column === 'observation_id'
          && typeof column.suffix === 'string' && Object.hasOwn(data.columns, 'observation_id')
          && Object.hasOwn(data.columns, 'term_id')) deferred.set(key, column);
      else throw new Error('Invalid packed Survey derived column descriptor');
    } else {
      validateFallback(column, data.row_count, key);
      define(columns, key, column);
    }
  }
  const resolved = validateResolved(await expandBundledSurveyColumns({ ...data, encoding: BUNDLED_SURVEY_ENCODING, columns }, readBundle));
  const output = {};
  for (const key of Object.keys(data.columns)) {
    let values = resolved.columns[key];
    if (deferred.has(key)) {
      const descriptor = deferred.get(key), source = resolved.columns[descriptor.column];
      if (!source.every(value => typeof value === 'string' && value.length > 0)) throw new Error('Packed Survey ID source must contain nonempty strings');
      if (!resolved.columns.term_id.every(value => typeof value === 'string' && value.length > 0 && descriptor.suffix === `|survey|${value}`)) {
        throw new Error('Packed Survey ID suffix must match every term');
      }
      values = source.map(value => value + descriptor.suffix);
    }
    define(output, key, values);
  }
  return { ...resolved, columns: output };
}
