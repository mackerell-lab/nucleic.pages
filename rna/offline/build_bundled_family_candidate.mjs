/** Build and fully verify family transport only; never activate a release. */
import assert from 'node:assert/strict';
import { mkdir, readFile, readdir, realpath, writeFile } from 'node:fs/promises';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { gzipSync, gunzipSync } from 'node:zlib';
import { OutputScope, sha256 } from './output_scope.mjs';
import { FAMILY_COLUMNAR_ENCODING, decodeFamilyRows } from '../core/survey-codec.js';
import {
  BUNDLED_FAMILY_ENCODING, FAMILY_BUNDLE_ENCODING, verifyFamilyBundle, expandBundledFamilyColumns,
} from '../core/bundled-family-codec.js';

const [sourceArgument, outputArgument, ...extra] = process.argv.slice(2);
if (!sourceArgument || !outputArgument || extra.length) {
  throw new Error('Usage: node build_bundled_family_candidate.mjs RELEASE_MANIFEST NEW_OUTPUT_DIRECTORY');
}
const MAX_BUNDLE_BYTES = 30 * 1024 * 1024;
const sourceFile = await realpath(sourceArgument), sourceRoot = path.dirname(sourceFile);
const sourceBytes = await readFile(sourceFile), source = JSON.parse(sourceBytes);
assert.equal(source.schema_version, 'rna-explorer-1');
assert.equal(source.molecule_type, 'RNA');
assert.ok(source.build_id && Array.isArray(source.families) && source.families.length, 'Source build and families');
assert.equal(new Set(source.families.map(family => family.id)).size, source.families.length, 'Unique family IDs');
const inside = (root, target) => target === root || target.startsWith(`${root}${path.sep}`);
const output = path.resolve(outputArgument);
const assets = await realpath(fileURLToPath(new URL('../../assets', import.meta.url)));
if ([sourceRoot, assets].some(root => inside(root, output) || inside(output, root))) {
  throw new Error('Candidate output must not overlap published assets or the source release');
}
const scope = new OutputScope([output]);
await scope.resolve(output);
// The caller supplies an existing parent; only acquire a new isolated directory.
await realpath(path.dirname(output));
await mkdir(output);
await mkdir(path.join(output, 'families', 'bundles'), { recursive: true });
const startedAt = new Date().toISOString();

function safeRelative(relative) {
  if (typeof relative !== 'string' || !relative || path.isAbsolute(relative) || relative.includes('\\')
      || relative.split('/').some(part => !part || part === '.' || part === '..')) throw new Error('Unsafe resource path');
  return relative;
}

async function readResource(root, descriptor) {
  const file = path.resolve(root, safeRelative(descriptor.path));
  const actual = await realpath(file);
  if (actual !== file || !inside(root, actual)) throw new Error('Resource path or symlink escaped its release');
  const compressed = await readFile(file);
  assert.equal(compressed.length, descriptor.bytes, `Compressed bytes: ${descriptor.path}`);
  assert.equal(sha256(compressed), descriptor.sha256, `Compressed checksum: ${descriptor.path}`);
  const raw = gunzipSync(compressed);
  assert.equal(raw.length, descriptor.uncompressed_bytes, `Raw bytes: ${descriptor.path}`);
  return JSON.parse(raw);
}

async function readFamily(descriptor) {
  assert.match(descriptor.id, /^[A-Za-z0-9_-]+$/);
  assert.equal(descriptor.encoding, FAMILY_COLUMNAR_ENCODING, 'Source family encoding');
  const data = await readResource(sourceRoot, descriptor);
  assert.equal(data.encoding, FAMILY_COLUMNAR_ENCODING);
  assert.equal(data.build_id, source.build_id, `Source build: ${descriptor.id}`);
  assert.equal(data.row_count, descriptor.row_count, `Source row count: ${descriptor.id}`);
  assert.ok(Number.isSafeInteger(data.row_count) && data.row_count >= 0, 'Safe source row count');
  assert.ok(data.columns && typeof data.columns === 'object' && !Array.isArray(data.columns), 'Source family columns');
  return data;
}

async function writeResource(relative, data) {
  const raw = Buffer.from(JSON.stringify(data)), compressed = gzipSync(raw, { level: 9 });
  const file = await scope.resolve(path.join(output, safeRelative(relative)));
  await writeFile(file, compressed, { flag: 'wx' });
  return { path: relative, bytes: compressed.length, uncompressed_bytes: raw.length, sha256: sha256(compressed) };
}

// Retain identities and one source location per column, not every family's rows.
const contents = new Map(), familyReferences = new Map();
for (const descriptor of source.families) {
  const data = await readFamily(descriptor), references = new Map();
  for (const [field, column] of Object.entries(data.columns)) {
    const raw = JSON.stringify(column), reference = sha256(raw);
    if (!contents.has(reference)) contents.set(reference, {
      bytes: Buffer.byteLength(raw), source: descriptor, field, families: new Set(),
    });
    contents.get(reference).families.add(descriptor.id);
    references.set(field, reference);
  }
  familyReferences.set(descriptor.id, references);
  console.error(`Indexed family ${descriptor.id}`);
}

const groups = new Map();
for (const [reference, record] of contents) {
  if (record.families.size < 2) continue;
  const membership = JSON.stringify([...record.families].sort());
  if (!groups.has(membership)) groups.set(membership, []);
  groups.get(membership).push(reference);
}
const emptyBundleBytes = Buffer.byteLength(JSON.stringify({ encoding: FAMILY_BUNDLE_ENCODING, columns: {} }));
const bundlePlans = [];
for (const [membership, references] of [...groups].sort(([a], [b]) => a.localeCompare(b))) {
  let current = [], currentBytes = emptyBundleBytes;
  for (const reference of references.sort()) {
    const entryBytes = Buffer.byteLength(JSON.stringify(reference)) + 1 + contents.get(reference).bytes;
    assert.ok(emptyBundleBytes + entryBytes <= MAX_BUNDLE_BYTES, 'Shared column exceeds the bundle byte limit');
    if (current.length && currentBytes + 1 + entryBytes > MAX_BUNDLE_BYTES) {
      bundlePlans.push({ references: current, bytes: currentBytes, families: JSON.parse(membership) });
      current = []; currentBytes = emptyBundleBytes;
    }
    currentBytes += entryBytes + (current.length ? 1 : 0);
    current.push(reference);
  }
  if (current.length) bundlePlans.push({ references: current, bytes: currentBytes, families: JSON.parse(membership) });
}

const bundles = Object.create(null), bundleForColumn = new Map();
for (const plan of bundlePlans) {
  // A membership group's columns normally share one source; release it each bundle.
  const bySource = new Map(), columns = Object.create(null);
  for (const reference of plan.references) {
    const record = contents.get(reference);
    if (!bySource.has(record.source.id)) bySource.set(record.source.id, []);
    bySource.get(record.source.id).push(reference);
  }
  for (const references of bySource.values()) {
    const data = await readFamily(contents.get(references[0]).source);
    for (const reference of references) {
      const column = data.columns[contents.get(reference).field];
      assert.equal(sha256(JSON.stringify(column)), reference, 'Source column identity');
      columns[reference] = column;
    }
  }
  const ordered = Object.fromEntries(Object.entries(columns).sort(([a], [b]) => a.localeCompare(b)));
  const payload = { encoding: FAMILY_BUNDLE_ENCODING, columns: ordered }, raw = JSON.stringify(payload);
  assert.equal(Buffer.byteLength(raw), plan.bytes, 'Exact bundle envelope accounting');
  assert.ok(Buffer.byteLength(raw) <= MAX_BUNDLE_BYTES, 'Bundle byte limit');
  const reference = sha256(raw);
  const descriptor = {
    ...await writeResource(`families/bundles/${reference}.json.gz`, payload),
    content_sha256: reference, column_count: plan.references.length, families: plan.families,
  };
  await verifyFamilyBundle(reference, await readResource(output, descriptor));
  bundles[reference] = descriptor;
  for (const column of plan.references) bundleForColumn.set(column, reference);
}

const families = [], validations = [];
for (const descriptor of source.families) {
  const original = await readFamily(descriptor), columns = Object.create(null);
  for (const [field, column] of Object.entries(original.columns)) {
    const reference = familyReferences.get(descriptor.id).get(field);
    assert.equal(sha256(JSON.stringify(column)), reference, 'Indexed source column identity');
    columns[field] = bundleForColumn.has(reference)
      ? { bundle: bundleForColumn.get(reference), column: reference } : column;
  }
  const bundleReferences = new Set(Object.values(columns).filter(column => column?.bundle).map(column => column.bundle));
  assert.ok(bundleReferences.size <= 2, `At most three cold resources per family: ${descriptor.id}`);
  const payload = { ...original, encoding: BUNDLED_FAMILY_ENCODING, columns };
  const destination = {
    ...descriptor, ...await writeResource(`families/${descriptor.id}.json.gz`, payload), encoding: BUNDLED_FAMILY_ENCODING,
  };
  const reread = await readResource(output, destination);
  const expanded = await expandBundledFamilyColumns(reread, async reference => {
    assert.ok(Object.hasOwn(bundles, reference), 'Registered family bundle');
    return (await verifyFamilyBundle(reference, await readResource(output, bundles[reference]))).bundle;
  });
  assert.deepEqual(expanded, original, `Complete transport metadata and columns: ${descriptor.id}`);
  const expectedRows = decodeFamilyRows(original), actualRows = decodeFamilyRows(expanded);
  assert.deepEqual(actualRows, expectedRows, `All decoded rows and fields: ${descriptor.id}`);
  assert.equal(actualRows.length, descriptor.row_count, 'Decoded family row count');
  const semanticSha256 = sha256(JSON.stringify(expectedRows));
  assert.equal(sha256(JSON.stringify(actualRows)), semanticSha256, 'Decoded family semantic identity');
  families.push(destination);
  validations.push({ family: descriptor.id, row_count: actualRows.length, decoded_equal: true,
    semantic_sha256: semanticSha256, cold_resource_count: 1 + bundleReferences.size });
  console.error(`Verified ${validations.length}/${source.families.length} families: ${descriptor.id}`);
}

const candidate = {
  schema_version: 'rna-family-bundled-candidate-1', family_only: true, build_id: source.build_id,
  source_manifest: { path: sourceFile, sha256: sha256(sourceBytes) }, generated_at: new Date().toISOString(),
  limitation: 'Family-only candidate with original build identity. Not a complete release or activation manifest. Candidate and validation files are local reports.',
  families, family_bundles: bundles,
};
const candidateJson = JSON.stringify(candidate, null, 2) + '\n';
await writeFile(path.join(output, 'candidate.json'), candidateJson, { flag: 'wx' });
assert.equal(await readFile(path.join(output, 'candidate.json'), 'utf8'), candidateJson, 'Exact candidate index bytes on disk');
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
assert.deepEqual(files.sort(), [...families, ...Object.values(bundles)].map(item => item.path).concat('candidate.json').sort(), 'Exact candidate inventory before report');
const size = descriptors => descriptors.reduce((sum, item) => sum + item.bytes, 0);
const sourceFamilyBytes = size(source.families), payloadBytes = size(families), bundleBytes = size(Object.values(bundles));
const report = {
  schema_version: 'rna-family-bundled-candidate-validation-1', started_at: startedAt, completed_at: new Date().toISOString(),
  family_only: true, source_build_id: source.build_id, source_manifest: candidate.source_manifest,
  output_directory: output, family_count: families.length, bundle_count: Object.keys(bundles).length,
  row_count: validations.reduce((sum, item) => sum + item.row_count, 0), all_families_verified: true,
  source_family_bytes: sourceFamilyBytes, candidate_payload_bytes: payloadBytes, candidate_bundle_bytes: bundleBytes,
  candidate_family_resource_bytes: payloadBytes + bundleBytes, saved_resource_bytes: sourceFamilyBytes - payloadBytes - bundleBytes,
  candidate_index_bytes: Buffer.byteLength(candidateJson), candidate_bytes_including_index: payloadBytes + bundleBytes + Buffer.byteLength(candidateJson),
  source_resource_count: source.families.length, candidate_resource_count: families.length + Object.keys(bundles).length,
  maximum_bundle_uncompressed_bytes: Math.max(0, ...Object.values(bundles).map(bundle => bundle.uncompressed_bytes)),
  bundle_uncompressed_byte_limit: MAX_BUNDLE_BYTES,
  minimum_family_cold_resource_count: Math.min(...validations.map(item => item.cold_resource_count)),
  maximum_family_cold_resource_count: Math.max(...validations.map(item => item.cold_resource_count)),
  validation: 'Source and destination compressed bytes/checksums, content-addressed bundles, full transport metadata, all decoded rows and fields, semantic SHA256 per family',
  limitation: 'Candidate resource and index bytes exclude this report. Browser cache/request and complete-release validation are separate gates. No active assets were modified.',
  families: validations,
};
await writeFile(path.join(output, 'validation.json'), JSON.stringify(report, null, 2) + '\n', { flag: 'wx' });
const { families: verifiedFamilies, ...summary } = report;
console.log(JSON.stringify(summary, null, 2));
