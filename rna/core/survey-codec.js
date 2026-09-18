/** Lossless columnar transport for large Survey scalar partitions. */
export const SURVEY_COLUMNAR_ENCODING = 'rna-survey-columnar-1';
export const COORDINATE_COLUMNAR_ENCODING = 'rna-coordinate-columnar-1';
export const FAMILY_COLUMNAR_ENCODING = 'rna-family-columnar-1';
export const INTERACTION_COLUMNAR_ENCODING = 'rna-interaction-columnar-1';

function encodeDictionaryColumns(rows, encoding, buildId = null) {
  if (!Array.isArray(rows)) throw new TypeError('RNA rows must be an array');
  const keys = [...new Set(rows.flatMap(row => Object.keys(row)))].sort();
  const columns = {}, missing = {};
  for (const key of keys) {
    const values = rows.map(row => row[key] ?? null);
    const absent = rows.flatMap((row, index) => Object.hasOwn(row, key) ? [] : [index]);
    if (absent.length) missing[key] = absent;
    const dictionary = [], indices = [], lookup = new Map();
    for (const value of values) {
      const identity = JSON.stringify(value);
      let index = lookup.get(identity);
      if (index === undefined) { index = dictionary.length; lookup.set(identity, index); dictionary.push(value); }
      indices.push(index);
    }
    const encoded = { dictionary, indices };
    columns[key] = JSON.stringify(encoded).length < JSON.stringify(values).length ? encoded : values;
  }
  return { encoding, ...(buildId ? { build_id: buildId } : {}), row_count: rows.length, columns,
    ...(Object.keys(missing).length ? { missing } : {}) };
}

function decodeDictionaryColumns(data, encoding, label) {
  if (!data || data.encoding !== encoding || !data.columns) throw new Error(`Unsupported RNA ${label} encoding`);
  const keys = Object.keys(data.columns), count = data.row_count;
  if (!Number.isSafeInteger(count) || count < 0) throw new Error(`Invalid RNA ${label} row count`);
  const rows = Array.from({ length: count }, () => ({}));
  for (const key of keys) {
    const column = data.columns[key];
    if (Array.isArray(column)) {
      if (column.length !== count) throw new Error(`${label} column length mismatch`);
      for (let index = 0; index < count; index++) rows[index][key] = column[index];
    } else if (column && Array.isArray(column.dictionary) && Array.isArray(column.indices)) {
      if (column.indices.length !== count || column.indices.some(index => !Number.isSafeInteger(index) || index < 0 || index >= column.dictionary.length)) throw new Error(`${label} dictionary column mismatch`);
      for (let index = 0; index < count; index++) rows[index][key] = column.dictionary[column.indices[index]];
    } else throw new Error(`${label} column is invalid`);
  }
  for (const [key, indices] of Object.entries(data.missing ?? {})) {
    if (!Object.hasOwn(data.columns, key) || !Array.isArray(indices) || indices.some(index => !Number.isSafeInteger(index) || index < 0 || index >= count)) throw new Error(`Invalid ${label} missing-field index`);
    for (const index of indices) delete rows[index][key];
  }
  return rows;
}

export function encodeSurveyRows(rows, buildId = null) {
  if (!Array.isArray(rows)) throw new TypeError('Survey rows must be an array');
  const keys = [...new Set(rows.flatMap(row => Object.keys(row)))].sort();
  const columns = Object.fromEntries(keys.map(key => [key, rows.map(row => row[key] ?? null)]));
  return { encoding: SURVEY_COLUMNAR_ENCODING, ...(buildId ? { build_id: buildId } : {}), row_count: rows.length, columns };
}

export function decodeSurveyRows(data, fields = null) {
  if (Array.isArray(data)) {
    if (!fields) return data;
    const wanted = new Set(fields);
    return data.map(row => Object.fromEntries(Object.entries(row).filter(([key]) => wanted.has(key))));
  }
  if (!data || data.encoding !== SURVEY_COLUMNAR_ENCODING || !data.columns) throw new Error('Unsupported RNA Survey encoding');
  const available = Object.keys(data.columns);
  const wanted = fields ? new Set(fields) : null;
  const keys = wanted ? available.filter(key => wanted.has(key)) : available;
  const count = data.row_count ?? data.columns[available[0]]?.length ?? 0;
  if (!keys.every(key => Array.isArray(data.columns[key]) && data.columns[key].length === count)) throw new Error('Survey column length mismatch');
  return Array.from({ length: count }, (_, index) => Object.fromEntries(keys.map(key => [key, data.columns[key][index]])));
}

export function encodeCoordinateRows(rows, buildId = null) {
  if (!Array.isArray(rows)) throw new TypeError('Coordinate rows must be an array');
  const keys = [...new Set(rows.flatMap(row => Object.keys(row)))].sort();
  const columns = Object.fromEntries(keys.map(key => [key, rows.map(row => row[key] ?? null)]));
  // Absent JSON fields differ from explicit null, especially for identity metadata.
  const missing = Object.fromEntries(keys.map(key => [key, rows.flatMap((row, index) => Object.hasOwn(row, key) ? [] : [index])]).filter(([, indices]) => indices.length));
  return { encoding: COORDINATE_COLUMNAR_ENCODING, ...(buildId ? { build_id: buildId } : {}), row_count: rows.length, columns,
    ...(Object.keys(missing).length ? { missing } : {}) };
}

export function decodeCoordinateRows(data) {
  if (Array.isArray(data)) return data;
  if (!data || data.encoding !== COORDINATE_COLUMNAR_ENCODING || !data.columns) throw new Error('Unsupported RNA coordinate encoding');
  const keys = Object.keys(data.columns), count = data.row_count;
  if (!Number.isSafeInteger(count) || count < 0) throw new Error('Invalid coordinate row count');
  if (!keys.every(key => Array.isArray(data.columns[key]) && data.columns[key].length === count)) throw new Error('Coordinate column length mismatch');
  const rows = Array.from({ length: count }, (_, index) => Object.fromEntries(keys.map(key => [key, data.columns[key][index]])));
  for (const [key, indices] of Object.entries(data.missing ?? {})) {
    if (!Object.hasOwn(data.columns, key) || !Array.isArray(indices) || indices.some(index => !Number.isSafeInteger(index) || index < 0 || index >= count)) throw new Error('Invalid coordinate missing-field index');
    for (const index of indices) delete rows[index][key];
  }
  return rows;
}

export function encodeFamilyRows(rows, buildId = null) { return encodeDictionaryColumns(rows, FAMILY_COLUMNAR_ENCODING, buildId); }
export function decodeFamilyRows(data) { return decodeDictionaryColumns(data, FAMILY_COLUMNAR_ENCODING, 'family'); }
export function encodeInteractionRows(rows, buildId = null) { return encodeDictionaryColumns(rows, INTERACTION_COLUMNAR_ENCODING, buildId); }
export function decodeInteractionRows(data) { return decodeDictionaryColumns(data, INTERACTION_COLUMNAR_ENCODING, 'interaction'); }
