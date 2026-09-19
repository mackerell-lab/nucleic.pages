import { FAMILY_COLUMNAR_ENCODING } from './survey-codec.js';

export const BUNDLED_FAMILY_ENCODING = 'rna-family-bundled-columns-1';
export const FAMILY_BUNDLE_ENCODING = 'rna-family-column-bundle-1';

const isHash = value => typeof value === 'string' && /^[a-f0-9]{64}$/.test(value);
const isRecord = value => value !== null && typeof value === 'object' && !Array.isArray(value);
const exactKeys = (value, names) => isRecord(value) && Object.keys(value).length === names.length && names.every(name => Object.hasOwn(value, name));

function validateColumn(column, count = null) {
  if (Array.isArray(column)) {
    if (count !== null && column.length !== count) throw new Error('Bundled family column length mismatch');
    return;
  }
  if (!exactKeys(column, ['dictionary', 'indices']) || !Array.isArray(column.dictionary) || !Array.isArray(column.indices)) {
    throw new Error('Invalid bundled family dictionary column');
  }
  if (count !== null && column.indices.length !== count) throw new Error('Bundled family dictionary length mismatch');
  if (column.indices.some(index => !Number.isSafeInteger(index) || index < 0 || index >= column.dictionary.length)) {
    throw new Error('Invalid bundled family dictionary index');
  }
}

function validateBundle(bundle) {
  if (!isRecord(bundle) || bundle.encoding !== FAMILY_BUNDLE_ENCODING || !isRecord(bundle.columns)) {
    throw new Error('Invalid family column bundle');
  }
  for (const [hash, column] of Object.entries(bundle.columns)) {
    if (!isHash(hash)) throw new Error('Invalid family bundle column hash');
    validateColumn(column);
  }
}

/** Authenticate complete bundle names and values, including dictionary indices. */
export async function verifyFamilyBundle(reference, payload) {
  if (!isHash(reference)) throw new Error('Invalid family bundle content hash');
  validateBundle(payload);
  const bytes = new TextEncoder().encode(JSON.stringify(payload));
  const digest = await globalThis.crypto.subtle.digest('SHA-256', bytes);
  const actual = Array.from(new Uint8Array(digest), byte => byte.toString(16).padStart(2, '0')).join('');
  if (actual !== reference) throw new Error('Family bundle hash mismatch');
  return { bundle: payload, byteLength: bytes.byteLength };
}

/** Expand transport references before decoding any family rows. The caller
 * authenticates returned bundles; no permanent bundle cache is retained here.
 */
export async function expandBundledFamilyColumns(data, readBundle) {
  if (data?.encoding !== BUNDLED_FAMILY_ENCODING) throw new Error('Unsupported bundled family encoding');
  if (!Number.isSafeInteger(data.row_count) || data.row_count < 0) throw new Error('Invalid bundled family row count');
  if (!isRecord(data.columns)) throw new Error('Invalid bundled family columns');
  const resolved = new Map();
  for (const column of Object.values(data.columns)) {
    if (isRecord(column) && (Object.hasOwn(column, 'bundle') || Object.hasOwn(column, 'column'))) {
      if (!exactKeys(column, ['bundle', 'column']) || !isHash(column.bundle) || !isHash(column.column)) {
        throw new Error('Invalid bundled family reference');
      }
      resolved.set(column.bundle, null);
    } else validateColumn(column, data.row_count);
  }
  if (Object.hasOwn(data, 'missing')) {
    if (!isRecord(data.missing)) throw new Error('Invalid bundled family missing-field map');
    for (const [key, indices] of Object.entries(data.missing)) {
      if (!Object.hasOwn(data.columns, key) || !Array.isArray(indices)
          || indices.some(index => !Number.isSafeInteger(index) || index < 0 || index >= data.row_count)) {
        throw new Error('Invalid bundled family missing-field index');
      }
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
    let column = descriptor;
    if (isRecord(descriptor) && Object.hasOwn(descriptor, 'bundle')) {
      const bundle = resolved.get(descriptor.bundle);
      if (!Object.hasOwn(bundle.columns, descriptor.column)) throw new Error('Missing family bundle column');
      column = bundle.columns[descriptor.column];
    }
    validateColumn(column, data.row_count);
    Object.defineProperty(columns, key, { value: column, enumerable: true });
  }
  return { ...data, encoding: FAMILY_COLUMNAR_ENCODING, columns };
}
