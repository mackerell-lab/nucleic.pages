import { SURVEY_COLUMNAR_ENCODING } from './survey-codec.js';

export const BUNDLED_SURVEY_ENCODING = 'rna-survey-bundled-columns-1';
export const SURVEY_BUNDLE_ENCODING = 'rna-survey-column-bundle-1';

const isHash = value => typeof value === 'string' && /^[a-f0-9]{64}$/.test(value);
const isRecord = value => value !== null && typeof value === 'object' && !Array.isArray(value);

function validateBundle(bundle) {
  if (!isRecord(bundle) || bundle.encoding !== SURVEY_BUNDLE_ENCODING || !isRecord(bundle.columns)) {
    throw new Error('Invalid Survey column bundle');
  }
  for (const [hash, values] of Object.entries(bundle.columns)) {
    if (!isHash(hash) || !Array.isArray(values)) throw new Error('Invalid Survey bundle column');
  }
}

/** Authenticate the entire bundle, including its column names and values. */
export async function verifySurveyBundle(reference, payload) {
  if (!isHash(reference)) throw new Error('Invalid Survey bundle content hash');
  validateBundle(payload);
  const bytes = new TextEncoder().encode(JSON.stringify(payload));
  const digest = await globalThis.crypto.subtle.digest('SHA-256', bytes);
  const actual = Array.from(new Uint8Array(digest), byte => byte.toString(16).padStart(2, '0')).join('');
  if (actual !== reference) throw new Error('Survey bundle hash mismatch');
  return { bundle: payload, byteLength: bytes.byteLength };
}

/** Resolve all columns before downstream projection. The reader authenticates
 * bundles; this operation retains them only until the table is expanded.
 */
export async function expandBundledSurveyColumns(data, readBundle) {
  if (data?.encoding !== BUNDLED_SURVEY_ENCODING) throw new Error('Unsupported bundled Survey encoding');
  if (!Number.isSafeInteger(data.row_count) || data.row_count < 0) throw new Error('Invalid bundled Survey row count');
  if (!isRecord(data.columns)) throw new Error('Invalid bundled Survey columns');

  const resolved = new Map();
  // Validate every descriptor and inline column before issuing any request.
  for (const descriptor of Object.values(data.columns)) {
    if (Array.isArray(descriptor)) {
      if (descriptor.length !== data.row_count) throw new Error('Bundled Survey column length mismatch');
    } else {
      if (!isRecord(descriptor) || Object.keys(descriptor).length !== 2
          || !Object.hasOwn(descriptor, 'bundle') || !Object.hasOwn(descriptor, 'column')
          || !isHash(descriptor.bundle) || !isHash(descriptor.column)) {
        throw new Error('Invalid bundled Survey reference');
      }
      resolved.set(descriptor.bundle, null);
    }
  }

  const references = [...resolved.keys()];
  let cursor = 0, failed = false;
  const worker = async () => {
    while (!failed && cursor < references.length) {
      const reference = references[cursor++];
      try {
        const bundle = await readBundle(reference);
        validateBundle(bundle);
        resolved.set(reference, bundle);
      } catch (error) { failed = true; throw error; }
    }
  };
  await Promise.all(Array.from({ length: Math.min(4, references.length) }, worker));

  const columns = {};
  for (const [key, descriptor] of Object.entries(data.columns)) {
    let values = descriptor;
    if (!Array.isArray(descriptor)) {
      const bundle = resolved.get(descriptor.bundle);
      if (!Object.hasOwn(bundle.columns, descriptor.column)) throw new Error('Missing Survey bundle column');
      values = bundle.columns[descriptor.column];
    }
    if (!Array.isArray(values) || values.length !== data.row_count) throw new Error('Bundled Survey column length mismatch');
    Object.defineProperty(columns, key, { value: values, enumerable: true });
  }
  return { ...data, encoding: SURVEY_COLUMNAR_ENCODING, columns };
}
