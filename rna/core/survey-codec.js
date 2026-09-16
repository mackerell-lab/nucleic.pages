/** Lossless columnar transport for large Survey scalar partitions. */
export const SURVEY_COLUMNAR_ENCODING = 'rna-survey-columnar-1';

export function encodeSurveyRows(rows, buildId = null) {
  if (!Array.isArray(rows)) throw new TypeError('Survey rows must be an array');
  const keys = [...new Set(rows.flatMap(row => Object.keys(row)))].sort();
  const columns = Object.fromEntries(keys.map(key => [key, rows.map(row => row[key] ?? null)]));
  return { encoding: SURVEY_COLUMNAR_ENCODING, ...(buildId ? { build_id: buildId } : {}), row_count: rows.length, columns };
}

export function decodeSurveyRows(data) {
  if (Array.isArray(data)) return data;
  if (!data || data.encoding !== SURVEY_COLUMNAR_ENCODING || !data.columns) throw new Error('Unsupported RNA Survey encoding');
  const keys = Object.keys(data.columns), count = data.row_count ?? data.columns[keys[0]]?.length ?? 0;
  if (!keys.every(key => Array.isArray(data.columns[key]) && data.columns[key].length === count)) throw new Error('Survey column length mismatch');
  return Array.from({ length: count }, (_, index) => Object.fromEntries(keys.map(key => [key, data.columns[key][index]])));
}
