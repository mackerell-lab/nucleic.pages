import { SURVEY_COLUMNAR_ENCODING } from './survey-codec.js';

export const SHARED_SURVEY_ENCODING = 'rna-survey-shared-columns-1';

export async function verifySharedColumn(reference, values) {
  if (!/^[a-f0-9]{64}$/.test(reference)) throw new Error('Invalid shared Survey content hash');
  if (!Array.isArray(values)) throw new Error('Shared Survey column requires an array');
  const digest = await globalThis.crypto.subtle.digest('SHA-256', new TextEncoder().encode(JSON.stringify(values)));
  const actual = Array.from(new Uint8Array(digest), byte => byte.toString(16).padStart(2, '0')).join('');
  if (actual !== reference) throw new Error('Shared Survey column hash mismatch');
  return values;
}

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
  for (const descriptor of Object.values(data.columns)) {
    if (!descriptor || typeof descriptor.reference !== 'string' || !descriptor.reference
        || Object.keys(descriptor).length !== 1) throw new Error('Invalid shared Survey reference');
    resolved.set(descriptor.reference, null);
  }
  const references = [...resolved.keys()];
  let cursor = 0, failure = null;
  const worker = async () => {
    while (!failure && cursor < references.length) {
      const reference = references[cursor++];
      try {
        const values = await readColumn(reference);
        if (!Array.isArray(values) || values.length !== data.row_count) throw new Error('Shared Survey column length mismatch');
        resolved.set(reference, values);
      } catch (error) { failure = error; throw error; }
    }
  };
  // Bound outstanding transfers without imposing one network round trip per column.
  await Promise.all(Array.from({ length: Math.min(4, references.length) }, worker));
  for (const [key, descriptor] of Object.entries(data.columns)) {
    const values = resolved.get(descriptor.reference);
    Object.defineProperty(columns, key, { value: values, enumerable: true });
  }
  return { ...data, encoding: SURVEY_COLUMNAR_ENCODING, columns };
}
