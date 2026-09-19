/** Pack and fully verify family numeric transport without activating a release. */
import assert from 'node:assert/strict';
import { mkdir, readFile, readdir, realpath, writeFile } from 'node:fs/promises';
import { gzipSync, gunzipSync } from 'node:zlib';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { performance } from 'node:perf_hooks';
import { createHash } from 'node:crypto';
import { OutputScope, sha256 } from './output_scope.mjs';
import { decodeFamilyRows } from '../core/survey-codec.js';
import { BUNDLED_FAMILY_ENCODING, verifyFamilyBundle, expandBundledFamilyColumns } from '../core/bundled-family-codec.js';
import { PACKED_FAMILY_ENCODING, encodePackedFamily, expandPackedFamily } from '../core/packed-family-codec.js';

const [sourceArgument, outputArgument, ...extra] = process.argv.slice(2);
if (!sourceArgument || !outputArgument || extra.length) {
  throw new Error('Usage: node build_packed_family_candidate.mjs RELEASE_MANIFEST NEW_OUTPUT_DIRECTORY');
}
const startedAt = new Date().toISOString(), started = performance.now();
const sourceFile = await realpath(sourceArgument), sourceRoot = path.dirname(sourceFile);
const sourceBytes = await readFile(sourceFile), source = JSON.parse(sourceBytes);
assert.equal(source.schema_version, 'rna-explorer-1');
assert.equal(source.molecule_type, 'RNA');
assert.ok(source.build_id && Array.isArray(source.families) && source.families.length, 'Source build and families');
assert.equal(new Set(source.families.map(family => family.id)).size, source.families.length, 'Unique family IDs');
assert.ok(source.family_bundles && typeof source.family_bundles === 'object' && !Array.isArray(source.family_bundles), 'Family bundle registry');
const inside = (root, file) => file === root || file.startsWith(`${root}${path.sep}`);
const output = path.resolve(outputArgument), scope = new OutputScope([output]);
const assets = await realpath(fileURLToPath(new URL('../../assets', import.meta.url)));
if ([assets, sourceRoot].some(root => inside(root, output) || inside(output, root))) {
  throw new Error('Candidate output must not overlap published assets or the source release');
}

function safeRelative(relative) {
  if (typeof relative !== 'string' || !relative || path.isAbsolute(relative) || relative.includes('\\')
      || relative.split('/').some(part => !part || part === '.' || part === '..')) throw new Error('Unsafe family resource path');
  if (!relative.startsWith('families/') || !relative.endsWith('.json.gz')) throw new Error('Expected family gzip resource path');
  return relative;
}

const paths = new Set();
for (const descriptor of [...source.families, ...Object.values(source.family_bundles)]) {
  safeRelative(descriptor.path);
  assert.ok(!paths.has(descriptor.path), 'Unique resource path');
  paths.add(descriptor.path);
}
for (const descriptor of source.families) {
  assert.match(descriptor.id, /^[A-Za-z0-9_-]+$/);
  assert.equal(descriptor.encoding, BUNDLED_FAMILY_ENCODING, 'Source family encoding');
  assert.ok(Number.isSafeInteger(descriptor.row_count) && descriptor.row_count >= 0, 'Safe source row count');
}
await scope.resolve(output);
await realpath(path.dirname(output));
await mkdir(output); // Exclusive acquisition; preserve existing and failed candidates.

async function readResource(root, descriptor) {
  const file = path.resolve(root, safeRelative(descriptor.path));
  const actual = await realpath(file);
  if (actual !== file || !inside(root, actual)) throw new Error('Resource path or symlink escaped its release');
  const compressed = await readFile(file);
  assert.equal(compressed.length, descriptor.bytes, `Compressed bytes: ${descriptor.path}`);
  assert.equal(sha256(compressed), descriptor.sha256, `Compressed checksum: ${descriptor.path}`);
  const raw = gunzipSync(compressed);
  assert.equal(raw.length, descriptor.uncompressed_bytes, `Raw bytes: ${descriptor.path}`);
  return { compressed, data: JSON.parse(raw) };
}

async function writeResource(relative, compressed) {
  const file = await scope.resolve(path.join(output, safeRelative(relative)));
  await mkdir(path.dirname(file), { recursive: true });
  await scope.resolve(file);
  await writeFile(file, compressed, { flag: 'wx' });
}

// Explicitly check IEEE numeric identity, in addition to full structural equality.
function checkNumbers(actual, expected) {
  if (typeof expected === 'number') {
    assert.ok(Object.is(actual, expected), 'Exact numeric identity');
    return 1;
  }
  if (expected === null || typeof expected !== 'object') return 0;
  let count = 0;
  for (const key of Object.keys(expected)) count += checkNumbers(actual[key], expected[key]);
  return count;
}

// Object insertion order is not a scientific difference. Hash canonical rows
// individually to avoid materializing a second complete decoded JSON document.
function semanticHash(rows) {
  const digest = createHash('sha256');
  const canonicalKeys = (_key, value) => value !== null && typeof value === 'object' && !Array.isArray(value)
    ? Object.fromEntries(Object.keys(value).sort().map(key => [key, value[key]])) : value;
  digest.update('[');
  for (let index = 0; index < rows.length; index++) {
    if (index) digest.update(',');
    digest.update(JSON.stringify(rows[index], canonicalKeys));
  }
  digest.update(']');
  return digest.digest('hex');
}

const bundles = structuredClone(source.family_bundles), bundleValidations = [];
for (const [reference, descriptor] of Object.entries(bundles)) {
  const original = await readResource(sourceRoot, descriptor);
  await verifyFamilyBundle(reference, original.data);
  await writeResource(descriptor.path, original.compressed);
  const reread = await readResource(output, descriptor);
  assert.ok(reread.compressed.equals(original.compressed), 'Copied bundle compressed bytes are identical');
  await verifyFamilyBundle(reference, reread.data);
  assert.deepEqual(reread.data, original.data, 'Copied bundle payload is identical');
  bundleValidations.push({ reference, path: descriptor.path, bytes: descriptor.bytes, sha256: descriptor.sha256, exact_copy: true });
}

async function expandFamily(data, root) {
  const cache = new Map();
  return expandBundledFamilyColumns(data, async reference => {
    assert.ok(Object.hasOwn(bundles, reference), 'Registered family bundle');
    if (!cache.has(reference)) {
      const resource = await readResource(root, bundles[reference]);
      cache.set(reference, (await verifyFamilyBundle(reference, resource.data)).bundle);
    }
    return cache.get(reference);
  });
}

const families = [], validations = [];
const totals = { checked_transport_numbers: 0, checked_decoded_row_numbers: 0,
  packed_value_families: 0, unchanged_value_families: 0,
  encode_ms: 0, serialize_compress_write_ms: 0, reread_expand_compare_ms: 0 };
for (const descriptor of source.families) {
  const familyStarted = performance.now();
  const { data: original } = await readResource(sourceRoot, descriptor);
  assert.equal(original.encoding, BUNDLED_FAMILY_ENCODING);
  assert.equal(original.build_id, source.build_id, 'Source family build identity');
  assert.equal(original.row_count, descriptor.row_count, 'Source family row count');
  const encodeStarted = performance.now(), packed = encodePackedFamily(original);
  const encodeMs = performance.now() - encodeStarted;
  assert.equal(packed.encoding, PACKED_FAMILY_ENCODING);
  assert.equal(packed.build_id, source.build_id, 'Candidate retains source build identity');
  assert.equal(packed.row_count, original.row_count);
  const valuesPacked = Array.isArray(original.columns.values) && !Array.isArray(packed.columns.values);
  const writeStarted = performance.now();
  const raw = Buffer.from(JSON.stringify(packed)), compressed = gzipSync(raw, { level: 9 });
  await writeResource(descriptor.path, compressed);
  const destination = { ...descriptor, encoding: PACKED_FAMILY_ENCODING,
    bytes: compressed.length, uncompressed_bytes: raw.length, sha256: sha256(compressed) };
  const writeMs = performance.now() - writeStarted;
  const verifyStarted = performance.now();
  const { data: reread } = await readResource(output, destination);
  assert.equal(reread.encoding, PACKED_FAMILY_ENCODING);
  const expanded = expandPackedFamily(reread);
  assert.deepEqual(expanded, original, 'Complete bundled family transport metadata and columns');
  const transportNumbers = checkNumbers(expanded, original);
  const expectedTransport = await expandFamily(original, sourceRoot);
  const actualTransport = await expandFamily(expanded, output);
  assert.deepEqual(actualTransport, expectedTransport, 'Complete expanded family transport');
  const expectedRows = decodeFamilyRows(expectedTransport), actualRows = decodeFamilyRows(actualTransport);
  assert.deepEqual(actualRows, expectedRows, `All decoded family rows and fields: ${descriptor.id}`);
  assert.equal(actualRows.length, descriptor.row_count);
  const rowNumbers = checkNumbers(actualRows, expectedRows);
  const semanticSha256 = semanticHash(expectedRows);
  assert.equal(semanticHash(actualRows), semanticSha256, 'Decoded family semantic identity');
  const verifyMs = performance.now() - verifyStarted;
  const references = new Set(Object.values(original.columns).filter(column => column?.bundle).map(column => column.bundle));
  families.push(destination);
  validations.push({ family: descriptor.id, path: descriptor.path, row_count: original.row_count,
    source_bytes: descriptor.bytes, candidate_bytes: destination.bytes, saved_bytes: descriptor.bytes - destination.bytes,
    source_sha256: descriptor.sha256, candidate_sha256: destination.sha256,
    values_packed: valuesPacked, checked_transport_numbers: transportNumbers, checked_decoded_row_numbers: rowNumbers,
    full_transport_equal: true, decoded_rows_equal: true, semantic_sha256: semanticSha256,
    cold_resource_count: 1 + references.size,
    timing_ms: { encode: encodeMs, serialize_compress_write: writeMs, reread_expand_compare: verifyMs, family_total: performance.now() - familyStarted } });
  totals.checked_transport_numbers += transportNumbers; totals.checked_decoded_row_numbers += rowNumbers;
  totals[valuesPacked ? 'packed_value_families' : 'unchanged_value_families']++;
  totals.encode_ms += encodeMs; totals.serialize_compress_write_ms += writeMs; totals.reread_expand_compare_ms += verifyMs;
  console.error(`Verified ${validations.length}/${source.families.length} packed families: ${descriptor.id}`);
}

const candidate = {
  schema_version: 'rna-family-packed-candidate-1', family_only: true, build_id: source.build_id,
  source_manifest: { path: sourceFile, sha256: sha256(sourceBytes) }, generated_at: new Date().toISOString(),
  limitation: 'Family-only candidate retaining source identity. Not a complete release or activation manifest. Candidate and validation files are local reports.',
  families, family_bundles: bundles,
};
const candidateJson = JSON.stringify(candidate, null, 2) + '\n';
await writeFile(path.join(output, 'candidate.json'), candidateJson, { flag: 'wx' });
assert.equal(await readFile(path.join(output, 'candidate.json'), 'utf8'), candidateJson, 'Exact candidate index bytes');
const files = [];
async function list(directory) {
  for (const entry of await readdir(directory, { withFileTypes: true })) {
    const file = path.join(directory, entry.name);
    if (entry.isDirectory()) await list(file);
    else if (entry.isFile()) files.push(path.relative(output, file));
    else throw new Error('Unexpected non-regular candidate resource');
  }
}
await list(output);
assert.deepEqual(files.sort(), [...paths, 'candidate.json'].sort(), 'Exact candidate inventory before report');
const size = descriptors => descriptors.reduce((sum, item) => sum + item.bytes, 0);
const sourceFamilyBytes = size(source.families), candidateFamilyBytes = size(families), bundleBytes = size(Object.values(bundles));
const report = {
  schema_version: 'rna-family-packed-candidate-validation-1', family_only: true,
  started_at: startedAt, completed_at: new Date().toISOString(), elapsed_ms: performance.now() - started,
  source_manifest: candidate.source_manifest, source_build_id: source.build_id, output_directory: output,
  family_count: families.length, bundle_count: Object.keys(bundles).length,
  row_count: validations.reduce((sum, item) => sum + item.row_count, 0),
  all_families_verified: true, all_transport_equal: true, all_decoded_rows_equal: true, all_bundles_exact_copies: true,
  source_family_bytes: sourceFamilyBytes, candidate_payload_bytes: candidateFamilyBytes,
  source_bundle_bytes: bundleBytes, candidate_bundle_bytes: bundleBytes,
  source_family_resource_bytes: sourceFamilyBytes + bundleBytes, candidate_family_resource_bytes: candidateFamilyBytes + bundleBytes,
  saved_resource_bytes: sourceFamilyBytes - candidateFamilyBytes,
  saved_resource_percent: 100 * (sourceFamilyBytes - candidateFamilyBytes) / (sourceFamilyBytes + bundleBytes),
  candidate_index_bytes: Buffer.byteLength(candidateJson), candidate_bytes_including_index: candidateFamilyBytes + bundleBytes + Buffer.byteLength(candidateJson),
  source_resource_count: paths.size, candidate_resource_count: paths.size,
  maximum_bundle_uncompressed_bytes: Math.max(0, ...Object.values(bundles).map(bundle => bundle.uncompressed_bytes)),
  minimum_family_cold_resource_count: Math.min(...validations.map(item => item.cold_resource_count)),
  maximum_family_cold_resource_count: Math.max(...validations.map(item => item.cold_resource_count)),
  totals,
  validation: 'Authenticated source and destination gzip/raw sizes and hashes, exact bundle copies, Object.is every numeric transport and decoded-row value, full transport metadata and every decoded row/field after gzip reread',
  limitation: 'Resource/index sizes exclude this validation report. Timings include validation on the shared host and are not browser performance evidence. Full release and browser acceptance remain separate gates. No active assets were modified.',
  families: validations, bundles: bundleValidations,
};
await writeFile(path.join(output, 'validation.json'), JSON.stringify(report, null, 2) + '\n', { flag: 'wx' });
const { families: verifiedFamilies, bundles: verifiedBundles, ...summary } = report;
console.log(JSON.stringify(summary, null, 2));
