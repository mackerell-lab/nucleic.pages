import { SURVEY_COLUMNAR_ENCODING } from './survey-codec.js';

export const SHARED_SURVEY_ENCODING = 'rna-survey-shared-columns-1';

function validateColumns(data) {
  if (!Number.isSafeInteger(data?.row_count) || data.row_count < 0) throw new Error('Invalid shared Survey row count');
  if (!data.columns || typeof data.columns !== 'object' || Array.isArray(data.columns)) throw new Error('Invalid shared Survey columns');
}

/** The writer returns an opaque reference; the release builder owns paths/hashes. */
export async function shareSurveyColumns(data, writeColumn) {
  if (data?.encoding !== SURVEY_COLUMNAR_ENCODING) throw new Error('Expected columnar Survey input');
  validateColumns(data);
  const columns = {};
  for (const [key, values] of Object.entries(data.columns)) {
    if (!Array.isArray(values) || values.length !== data.row_count) throw new Error('Survey column length mismatch');
    const reference = await writeColumn(values);
    if (typeof reference !== 'string' || !reference) throw new Error('Invalid shared Survey reference');
    Object.defineProperty(columns, key, { value: { reference }, enumerable: true });
  }
  return { ...data, encoding: SHARED_SURVEY_ENCODING, columns };
}

/** Resolve and validate every column before any field projection by the caller.
 * Requests are deduplicated only within this operation; no permanent raw cache.
 */
export async function expandSharedSurveyColumns(data, readColumn) {
  if (data?.encoding !== SHARED_SURVEY_ENCODING) throw new Error('Unsupported shared Survey encoding');
  validateColumns(data);
  const columns = {}, resolved = new Map();
  for (const [key, descriptor] of Object.entries(data.columns)) {
    if (!descriptor || typeof descriptor.reference !== 'string' || !descriptor.reference
        || Object.keys(descriptor).length !== 1) throw new Error('Invalid shared Survey reference');
    const reference = descriptor.reference;
    if (!resolved.has(reference)) resolved.set(reference, await readColumn(reference));
    const values = resolved.get(reference);
    if (!Array.isArray(values) || values.length !== data.row_count) throw new Error('Shared Survey column length mismatch');
    Object.defineProperty(columns, key, { value: values, enumerable: true });
  }
  return { ...data, encoding: SURVEY_COLUMNAR_ENCODING, columns };
}
