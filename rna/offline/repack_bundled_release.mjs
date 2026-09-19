/** Stage a complete lossless bundled release without changing published assets. */
import assert from 'node:assert/strict';
import { mkdir, readFile, readdir, realpath, writeFile, access } from 'node:fs/promises';
import { gzipSync, gunzipSync } from 'node:zlib';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { OutputScope, sha256 } from './output_scope.mjs';
import {
  decodeSurveyRows, decodeCoordinateRows, decodeFamilyRows, decodeInteractionRows,
  SURVEY_COLUMNAR_ENCODING, COORDINATE_COLUMNAR_ENCODING, FAMILY_COLUMNAR_ENCODING, INTERACTION_COLUMNAR_ENCODING,
} from '../core/survey-codec.js';
import { BUNDLED_SURVEY_ENCODING, expandBundledSurveyColumns, verifySurveyBundle } from '../core/bundled-survey-codec.js';

const [sourceArgument, candidateArgument, outputArgument, buildId, ...extra] = process.argv.slice(2);
if (!sourceArgument || !candidateArgument || !outputArgument || !buildId || extra.length) {
  throw new Error('Usage: node repack_bundled_release.mjs SOURCE_MANIFEST SCALAR_CANDIDATE NEW_OUTPUT_DIRECTORY NEW_BUILD_ID');
}
if (!/^[A-Za-z0-9][A-Za-z0-9_.-]*$/.test(buildId)) throw new Error('Unsafe new build ID');
const inside = (root, file) => file === root || file.startsWith(`${root}${path.sep}`);
const sourceFile = await realpath(sourceArgument), sourceRoot = path.dirname(sourceFile);
const candidateFile = await realpath(candidateArgument), candidateRoot = path.dirname(candidateFile);
const sourceBytes = await readFile(sourceFile), candidateBytes = await readFile(candidateFile);
const source = JSON.parse(sourceBytes), candidate = JSON.parse(candidateBytes);
assert.equal(source.schema_version, 'rna-explorer-1');
assert.equal(source.molecule_type, 'RNA');
assert.equal(source.partial, false, 'Repack requires a complete source release');
assert.ok(source.build_id && source.build_id !== buildId, 'New immutable build identity');
assert.equal(candidate.scalar_only, true);
assert.equal(candidate.build_id, source.build_id, 'Candidate source identity');
assert.equal(candidate.source_manifest.sha256, sha256(sourceBytes), 'Candidate source manifest identity');
assert.deepEqual(candidate.survey.terms, source.survey.terms, 'Unchanged term definitions');
assert.deepEqual(candidate.survey.opening_bins, source.survey.opening_bins, 'Unchanged opening bins');
assert.deepEqual(Object.keys(candidate.survey.scalars.terms).sort(), Object.keys(source.survey.scalars.terms).sort(), 'Complete scalar registry');

const output = path.resolve(outputArgument), parent = path.dirname(output);
const assets = await realpath(fileURLToPath(new URL('../../assets', import.meta.url)));
if ([assets, sourceRoot, candidateRoot].some(root => inside(root, output) || inside(output, root))) {
  throw new Error('Output must not overlap published assets or either source');
}
const scope = new OutputScope([output]), reportScope = new OutputScope([parent]);
await scope.resolve(output); // Checks existing ancestors before any directory creation.
await mkdir(parent, { recursive: true });
await scope.resolve(output);
const reportFile = `${output}-repack.json`, inventoryFile = `${output}-inventory.json`;
for (const file of [reportFile, inventoryFile]) {
  await reportScope.resolve(file);
  try { await access(file); throw new Error(`Report already exists: ${file}`); }
  catch (error) { if (error.code !== 'ENOENT') throw error; }
}
await mkdir(output); // Exclusive output ownership; no overwrite or resume.
const startedAt = new Date().toISOString();
const resourceRecords = [], seenPaths = new Set();

function safeRelative(relative) {
  if (typeof relative !== 'string' || !relative || path.isAbsolute(relative)
      || relative.includes('\\') || relative.split('/').some(part => !part || part === '.' || part === '..')) {
    throw new Error('Unsafe release resource path');
  }
  return relative;
}

async function readResource(root, descriptor) {
  const file = path.resolve(root, safeRelative(descriptor.path));
  const resolved = await realpath(file);
  if (!inside(root, resolved) || resolved !== file) throw new Error('Resource path or symlink escaped release scope');
  const compressed = await readFile(file);
  assert.equal(sha256(compressed), descriptor.sha256, `Resource checksum: ${descriptor.path}`);
  if (descriptor.bytes != null) assert.equal(compressed.length, descriptor.bytes, `Resource bytes: ${descriptor.path}`);
  const raw = descriptor.path.endsWith('.gz') ? gunzipSync(compressed) : compressed;
  if (descriptor.uncompressed_bytes != null) assert.equal(raw.length, descriptor.uncompressed_bytes, `Resource raw bytes: ${descriptor.path}`);
  const data = JSON.parse(raw);
  if (descriptor.encoding != null) assert.equal(data.encoding, descriptor.encoding, `Resource encoding: ${descriptor.path}`);
  return { compressed, raw, data };
}

async function writeResource(descriptor, compressed, rawLength) {
  const relative = safeRelative(descriptor.path);
  assert.ok(!seenPaths.has(relative), `Unique destination resource: ${relative}`);
  seenPaths.add(relative);
  await scope.write(path.join(output, relative), compressed);
  return { ...descriptor, bytes: compressed.length, uncompressed_bytes: rawLength, sha256: sha256(compressed) };
}

async function decode(data, root, bundles = {}) {
  if (data?.encoding === BUNDLED_SURVEY_ENCODING) {
    data = await expandBundledSurveyColumns(data, async reference => {
      assert.ok(Object.hasOwn(bundles, reference), 'Registered bundle reference');
      return (await verifySurveyBundle(reference, (await readResource(root, bundles[reference])).data)).bundle;
    });
  }
  if (data?.encoding === SURVEY_COLUMNAR_ENCODING) return decodeSurveyRows(data);
  if (data?.encoding === COORDINATE_COLUMNAR_ENCODING) return decodeCoordinateRows(data);
  if (data?.encoding === FAMILY_COLUMNAR_ENCODING) return decodeFamilyRows(data);
  if (data?.encoding === INTERACTION_COLUMNAR_ENCODING) return decodeInteractionRows(data);
  if (data?.encoding) throw new Error(`Unsupported release transport: ${data.encoding}`);
  return data;
}

const manifest = structuredClone(source);
manifest.build_id = buildId;
manifest.generated_at = new Date().toISOString();
manifest.survey.bundles = {};
// Bundle identities are independent of release ID; validate and copy exact bytes.
for (const [reference, descriptor] of Object.entries(candidate.survey.bundles)) {
  assert.equal(descriptor.path, `survey/bundles/${reference}.json.gz`, 'Canonical bundle path');
  assert.equal(descriptor.content_sha256, reference, 'Bundle content descriptor');
  const { compressed, raw, data } = await readResource(candidateRoot, descriptor);
  await verifySurveyBundle(reference, data);
  assert.equal(Object.keys(data.columns).length, descriptor.column_count, 'Bundle column count');
  const destination = await writeResource(descriptor, compressed, raw.length);
  const reread = await readResource(output, destination);
  await verifySurveyBundle(reference, reread.data);
  assert.deepEqual(reread.data, data, 'Complete bundle equality');
  manifest.survey.bundles[reference] = destination;
  resourceRecords.push({ path: destination.path, sha256: destination.sha256, bytes: destination.bytes,
    uncompressed_bytes: destination.uncompressed_bytes, kind: 'bundle', unchanged_bytes: true });
}

let rewritten = 0, copied = 0, validatedRows = 0;
async function repack(descriptor, replacement = null) {
  const original = await readResource(sourceRoot, descriptor);
  if (original.data?.encoding) assert.equal(original.data.build_id, source.build_id, `Source build identity: ${descriptor.path}`);
  const input = replacement ? await readResource(candidateRoot, replacement) : original;
  if (replacement) {
    assert.equal(input.data.encoding, BUNDLED_SURVEY_ENCODING);
    assert.equal(input.data.build_id, source.build_id, 'Scalar candidate build identity');
    assert.equal(replacement.row_count, descriptor.row_count, 'Scalar candidate descriptor row count');
    const expanded = await expandBundledSurveyColumns(input.data, async reference => {
      assert.ok(Object.hasOwn(manifest.survey.bundles, reference), 'Registered scalar bundle');
      return (await verifySurveyBundle(reference, (await readResource(output, manifest.survey.bundles[reference])).data)).bundle;
    });
    assert.deepEqual(expanded, original.data, `All original scalar transport metadata and columns: ${descriptor.path}`);
  }
  let raw = input.raw, compressed = input.compressed;
  if (input.data && Object.hasOwn(input.data, 'build_id')) {
    assert.equal(input.data.build_id, source.build_id, 'Input build identity');
    raw = Buffer.from(JSON.stringify({ ...input.data, build_id: buildId }));
    compressed = gzipSync(raw, { level: 9 });
    rewritten++;
  } else {
    assert.ok(!input.data?.encoding, 'Encoded assets require explicit build identity');
    copied++;
  }
  const destination = await writeResource({ ...descriptor, ...(replacement ? { encoding: BUNDLED_SURVEY_ENCODING } : {}) }, compressed, raw.length);
  const reread = await readResource(output, destination);
  if (reread.data?.encoding) assert.equal(reread.data.build_id, buildId, 'Destination build identity');
  if (input.data?.encoding) {
    const { build_id: oldBuild, ...before } = input.data;
    const { build_id: newBuild, ...after } = reread.data;
    assert.deepEqual(after, before, `Transport metadata and columns: ${descriptor.path}`);
  } else assert.deepEqual(reread.data, original.data, `Unencoded asset equality: ${descriptor.path}`);
  const expected = await decode(original.data, sourceRoot, source.survey.bundles);
  const actual = await decode(reread.data, output, manifest.survey.bundles);
  assert.deepEqual(actual, expected, `Every decoded row and field: ${descriptor.path}`);
  const rows = Array.isArray(actual) ? actual.length : null;
  if (descriptor.row_count != null) assert.equal(rows, descriptor.row_count, `Decoded row count: ${descriptor.path}`);
  if (rows != null) validatedRows += rows;
  resourceRecords.push({ path: destination.path, sha256: destination.sha256, bytes: destination.bytes,
    uncompressed_bytes: destination.uncompressed_bytes, row_count: rows, decoded_equal: true,
    semantic_sha256: sha256(JSON.stringify(actual)), rewritten_build_id: Boolean(input.data?.encoding) });
  if ((copied + rewritten) % 25 === 0) console.error(`Verified ${copied + rewritten} release resources`);
  return destination;
}

manifest.metadata = await repack(source.metadata);
manifest.provenance.decisions = await repack(source.provenance.decisions);
for (let index = 0; index < source.families.length; index++) manifest.families[index] = await repack(source.families[index]);
for (const [kind, descriptor] of Object.entries(source.relations)) manifest.relations[kind] = await repack(descriptor);
for (const [term, descriptor] of Object.entries(source.survey.scalars.terms)) {
  manifest.survey.scalars.terms[term] = await repack(descriptor, candidate.survey.scalars.terms[term]);
}
for (const [group, value] of Object.entries(source.survey.coordinates.groups)) {
  if (value.partitions) {
    manifest.survey.coordinates.groups[group].partitions = [];
    for (const descriptor of value.partitions) manifest.survey.coordinates.groups[group].partitions.push(await repack(descriptor));
  } else manifest.survey.coordinates.groups[group] = await repack(value);
}

manifest.provenance.source_release = {
  build_id: source.build_id, manifest_sha256: sha256(sourceBytes),
  build_stages: source.provenance.build_stages ?? {},
  ...(source.provenance.source_release ? { source_release: source.provenance.source_release } : {}),
};
delete manifest.provenance.build_stages;
manifest.provenance.repack = {
  operation: 'lossless_bundled_survey_transport', source_build_id: source.build_id,
  source_manifest_sha256: sha256(sourceBytes), scalar_candidate_sha256: sha256(candidateBytes),
  scalar_encoding: BUNDLED_SURVEY_ENCODING,
  code_sha256: Object.fromEntries(await Promise.all([
    './repack_bundled_release.mjs', './output_scope.mjs', '../core/survey-codec.js', '../core/bundled-survey-codec.js',
  ].map(async file => [file, sha256(await readFile(new URL(file, import.meta.url)))]))),
  note: 'Storage transformation only. Scientific rows and source selection are unchanged; source stage history is retained separately.',
};
await scope.json(path.join(output, 'manifest.json'), manifest);

function descriptorsIn(value, result = []) {
  if (!value || typeof value !== 'object') return result;
  if (typeof value.path === 'string') result.push(value);
  for (const item of Object.values(value)) if (item && typeof item === 'object') descriptorsIn(item, result);
  return result;
}
const referenced = descriptorsIn(manifest);
assert.equal(new Set(referenced.map(item => item.path)).size, referenced.length, 'Unique manifest resource references');
assert.deepEqual(referenced.map(item => item.path).sort(), [...seenPaths].sort(), 'Complete manifest resource references');
const actualFiles = [];
async function walk(directory) {
  for (const entry of await readdir(directory, { withFileTypes: true })) {
    const file = path.join(directory, entry.name);
    if (entry.isSymbolicLink()) throw new Error('Symlink in staged release');
    if (entry.isDirectory()) await walk(file);
    else if (entry.isFile()) actualFiles.push(path.relative(output, file));
    else throw new Error('Non-regular file in staged release');
  }
}
await walk(output);
assert.deepEqual(actualFiles.sort(), [...seenPaths, 'manifest.json'].sort(), 'Exact staged file inventory');
const manifestRaw = await readFile(path.join(output, 'manifest.json'));
assert.deepEqual(JSON.parse(manifestRaw), manifest, 'Manifest disk roundtrip');
const inventory = {
  schema_version: 'rna-release-repack-inventory-1', build_id: buildId,
  manifest_sha256: sha256(manifestRaw), resource_count: referenced.length, file_count: actualFiles.length,
  bytes: resourceRecords.reduce((sum, item) => sum + item.bytes, manifestRaw.length),
  files: [...resourceRecords.map(({ path, sha256, bytes, uncompressed_bytes }) => ({ path, sha256, bytes, uncompressed_bytes })),
    { path: 'manifest.json', sha256: sha256(manifestRaw), bytes: manifestRaw.length }].sort((a, b) => a.path.localeCompare(b.path)),
};
const report = {
  schema_version: 'rna-release-repack-validation-1', started_at: startedAt, completed_at: new Date().toISOString(),
  source_manifest: { path: sourceFile, sha256: sha256(sourceBytes), build_id: source.build_id },
  scalar_candidate: { path: candidateFile, sha256: sha256(candidateBytes) },
  output_directory: output, build_id: buildId, partial: false,
  resource_count: referenced.length, file_count: actualFiles.length, bytes: inventory.bytes,
  rewritten_build_assets: rewritten, copied_unencoded_assets: copied, bundle_count: Object.keys(manifest.survey.bundles).length,
  validated_rows: validatedRows, all_decoded_equal: true, exact_inventory: true,
  manifest_sha256: inventory.manifest_sha256,
  limitation: 'Storage and full decoded equality validated. Run validateRelease on this staged manifest and browser acceptance separately. No published assets or activation pointer were modified.',
  resources: resourceRecords,
};
await writeFile(inventoryFile, JSON.stringify(inventory, null, 2) + '\n', { flag: 'wx' });
await writeFile(reportFile, JSON.stringify(report, null, 2) + '\n', { flag: 'wx' });
const { resources, ...summary } = report;
console.log(JSON.stringify({ ...summary, inventory_file: inventoryFile, report_file: reportFile }, null, 2));
