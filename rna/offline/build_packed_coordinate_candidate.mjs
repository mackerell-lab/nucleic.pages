/** Build and verify every coordinate partition without activating a release. */
import assert from 'node:assert/strict';
import { mkdir, readFile, readdir, realpath, writeFile } from 'node:fs/promises';
import { gzipSync, gunzipSync } from 'node:zlib';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { performance } from 'node:perf_hooks';
import { OutputScope, sha256 } from './output_scope.mjs';
import { COORDINATE_COLUMNAR_ENCODING, decodeCoordinateRows } from '../core/survey-codec.js';
import { PACKED_COORDINATE_ENCODING, encodePackedCoordinates, expandPackedCoordinates } from '../core/packed-coordinate-codec.js';

const [sourceArgument, outputArgument, ...extra] = process.argv.slice(2);
if (!sourceArgument || !outputArgument || extra.length) {
  throw new Error('Usage: node build_packed_coordinate_candidate.mjs RELEASE_MANIFEST NEW_OUTPUT_DIRECTORY');
}
const startedAt = new Date().toISOString(), started = performance.now();
const sourceFile = await realpath(sourceArgument), sourceRoot = path.dirname(sourceFile);
const sourceBytes = await readFile(sourceFile), source = JSON.parse(sourceBytes);
assert.equal(source.schema_version, 'rna-explorer-1');
assert.equal(source.molecule_type, 'RNA');
assert.ok(source.build_id && source.survey?.coordinates?.groups, 'Source build and coordinate groups');
const inside = (root, file) => file === root || file.startsWith(`${root}${path.sep}`);
const output = path.resolve(outputArgument), scope = new OutputScope([output]);
const assets = await realpath(fileURLToPath(new URL('../../assets', import.meta.url)));
if ([assets, sourceRoot].some(root => inside(root, output) || inside(output, root))) {
  throw new Error('Candidate output must not overlap published assets or the source release');
}

function safeRelative(relative) {
  if (typeof relative !== 'string' || !relative || path.isAbsolute(relative) || relative.includes('\\')
      || relative.split('/').some(part => !part || part === '.' || part === '..')) throw new Error('Unsafe coordinate resource path');
  if (!relative.startsWith('survey/coordinates/') || !relative.endsWith('.json.gz')) throw new Error('Expected coordinate gzip resource path');
  return relative;
}

const jobs = [], paths = new Set();
for (const [group, descriptor] of Object.entries(source.survey.coordinates.groups)) {
  assert.match(group, /^[A-Za-z0-9_-]+$/);
  const partitions = descriptor.partitions ?? [descriptor];
  assert.ok(Array.isArray(partitions) && partitions.length, `Coordinate partitions: ${group}`);
  let groupRows = 0;
  for (const [index, partition] of partitions.entries()) {
    safeRelative(partition.path);
    assert.ok(!paths.has(partition.path), 'Unique source coordinate resource');
    paths.add(partition.path);
    assert.equal(partition.encoding, COORDINATE_COLUMNAR_ENCODING, 'Source coordinate encoding');
    assert.ok(Number.isSafeInteger(partition.row_count) && partition.row_count >= 0 && partition.row_count <= 10000, 'Coordinate partition row limit');
    groupRows += partition.row_count;
    jobs.push({ group, index, descriptor: partition, partitioned: Boolean(descriptor.partitions) });
  }
  assert.equal(groupRows, descriptor.row_count, `Coordinate group row count: ${group}`);
}
assert.ok(jobs.length, 'At least one coordinate partition');
await scope.resolve(output);
await realpath(path.dirname(output));
await mkdir(output); // Exclusive acquisition; failed/existing candidates are retained.

async function readResource(root, descriptor) {
  const file = path.resolve(root, safeRelative(descriptor.path));
  const actual = await realpath(file);
  if (file !== actual || !inside(root, actual)) throw new Error('Coordinate resource path or symlink escaped its release');
  const compressed = await readFile(file);
  assert.equal(compressed.length, descriptor.bytes, `Compressed bytes: ${descriptor.path}`);
  assert.equal(sha256(compressed), descriptor.sha256, `Compressed checksum: ${descriptor.path}`);
  const raw = gunzipSync(compressed);
  assert.equal(raw.length, descriptor.uncompressed_bytes, `Uncompressed bytes: ${descriptor.path}`);
  const data = JSON.parse(raw);
  assert.equal(data.encoding, descriptor.encoding, `Payload encoding: ${descriptor.path}`);
  return data;
}

const coordinates = structuredClone(source.survey.coordinates);
const validations = [];
const totals = { packed_axis_columns: 0, fallback_axis_columns: 0, checked_axis_elements: 0, checked_finite_axis_numbers: 0,
  fallback_null_values: 0, fallback_absent_values: 0, fallback_nonnumeric_values: 0, fallback_nonfinite_values: 0,
  fallback_missing_field_indices: 0, encode_ms: 0, serialize_compress_write_ms: 0, reread_expand_compare_ms: 0 };

for (const job of jobs) {
  const { group, index, descriptor } = job, partitionStarted = performance.now();
  const original = await readResource(sourceRoot, descriptor);
  assert.equal(original.build_id, source.build_id, 'Source partition build identity');
  assert.equal(original.row_count, descriptor.row_count, 'Source partition row count');
  const encodeStarted = performance.now();
  const packed = await encodePackedCoordinates(original);
  const encodeMs = performance.now() - encodeStarted;
  assert.equal(packed.encoding, PACKED_COORDINATE_ENCODING);
  assert.equal(packed.build_id, source.build_id, 'Candidate retains source build identity');
  assert.equal(packed.row_count, original.row_count);
  const packedAxes = [], fallbackAxes = [];
  const fallback = { null_values: 0, absent_values: 0, nonnumeric_values: 0, nonfinite_values: 0, missing_field_indices: 0 };
  for (const axis of ['x', 'y', 'z']) {
    assert.ok(Array.isArray(original.columns[axis]), `Original ${axis} array`);
    if (Array.isArray(packed.columns[axis])) {
      fallbackAxes.push(axis);
      assert.deepEqual(packed.columns[axis], original.columns[axis], 'Fallback arrays remain exact');
      fallback.missing_field_indices += (original.missing?.[axis] ?? []).length;
      for (let row = 0; row < original.row_count; row++) {
        const value = original.columns[axis][row];
        if (!Object.hasOwn(original.columns[axis], row)) fallback.absent_values++;
        else if (value === null) fallback.null_values++;
        else if (typeof value !== 'number') fallback.nonnumeric_values++;
        else if (!Number.isFinite(value)) fallback.nonfinite_values++;
      }
    } else {
      assert.equal(packed.columns[axis].encoding, 'float64-le-shuffled-base64-1');
      assert.equal(packed.columns[axis].count, original.row_count);
      packedAxes.push(axis);
    }
  }
  const writeStarted = performance.now();
  const raw = Buffer.from(JSON.stringify(packed)), compressed = gzipSync(raw, { level: 9 });
  const file = await scope.resolve(path.join(output, descriptor.path));
  await mkdir(path.dirname(file), { recursive: true });
  await scope.resolve(file);
  await writeFile(file, compressed, { flag: 'wx' });
  const destination = { ...descriptor, encoding: PACKED_COORDINATE_ENCODING,
    bytes: compressed.length, uncompressed_bytes: raw.length, sha256: sha256(compressed) };
  const writeMs = performance.now() - writeStarted;
  const verifyStarted = performance.now();
  const reread = await readResource(output, destination);
  const expanded = await expandPackedCoordinates(reread);
  assert.deepEqual(expanded, original, 'Complete coordinate transport metadata and columns');
  let finiteNumbers = 0;
  for (const axis of ['x', 'y', 'z']) {
    for (let row = 0; row < original.row_count; row++) {
      assert.equal(Object.hasOwn(expanded.columns[axis], row), Object.hasOwn(original.columns[axis], row), 'Coordinate array presence');
      assert.ok(Object.is(expanded.columns[axis][row], original.columns[axis][row]), `Exact ${axis} at ${descriptor.path}/${row}`);
      if (typeof original.columns[axis][row] === 'number' && Number.isFinite(original.columns[axis][row])) finiteNumbers++;
    }
  }
  const expectedRows = decodeCoordinateRows(original), actualRows = decodeCoordinateRows(expanded);
  assert.deepEqual(actualRows, expectedRows, `All decoded coordinate rows and fields: ${descriptor.path}`);
  assert.equal(actualRows.length, descriptor.row_count);
  const semanticSha256 = sha256(JSON.stringify(expectedRows));
  assert.equal(sha256(JSON.stringify(actualRows)), semanticSha256, 'Decoded coordinate semantic identity');
  const verifyMs = performance.now() - verifyStarted;
  if (job.partitioned) coordinates.groups[group].partitions[index] = destination;
  else coordinates.groups[group] = destination;
  validations.push({ group, path: descriptor.path, row_count: original.row_count,
    source_bytes: descriptor.bytes, candidate_bytes: destination.bytes, saved_bytes: descriptor.bytes - destination.bytes,
    source_sha256: descriptor.sha256, candidate_sha256: destination.sha256,
    packed_axes: packedAxes, fallback_axes: fallbackAxes, fallback,
    object_is_checked_elements: 3 * original.row_count, object_is_checked_finite_numbers: finiteNumbers,
    full_transport_equal: true, decoded_rows_equal: true, semantic_sha256: semanticSha256,
    timing_ms: { encode: encodeMs, serialize_compress_write: writeMs, reread_expand_compare: verifyMs, partition_total: performance.now() - partitionStarted } });
  totals.packed_axis_columns += packedAxes.length; totals.fallback_axis_columns += fallbackAxes.length;
  totals.checked_axis_elements += 3 * original.row_count; totals.checked_finite_axis_numbers += finiteNumbers;
  for (const [key, value] of Object.entries(fallback)) totals[`fallback_${key}`] += value;
  totals.encode_ms += encodeMs; totals.serialize_compress_write_ms += writeMs; totals.reread_expand_compare_ms += verifyMs;
  if (validations.length % 25 === 0 || validations.length === jobs.length) console.error(`Verified ${validations.length}/${jobs.length} coordinate partitions`);
}

const candidate = {
  schema_version: 'rna-coordinate-packed-candidate-1', coordinate_only: true, build_id: source.build_id,
  source_manifest: { path: sourceFile, sha256: sha256(sourceBytes) }, generated_at: new Date().toISOString(),
  limitation: 'Coordinate-only candidate retaining source identity. Not a complete release or activation manifest. Candidate and validation files are local reports.',
  survey: { coordinates, opening_bins: source.survey.opening_bins },
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
const sourceResourceBytes = validations.reduce((sum, item) => sum + item.source_bytes, 0);
const candidateResourceBytes = validations.reduce((sum, item) => sum + item.candidate_bytes, 0);
const report = {
  schema_version: 'rna-coordinate-packed-candidate-validation-1', coordinate_only: true,
  started_at: startedAt, completed_at: new Date().toISOString(), elapsed_ms: performance.now() - started,
  source_manifest: candidate.source_manifest, source_build_id: source.build_id, output_directory: output,
  group_count: Object.keys(coordinates.groups).length, partition_count: jobs.length,
  row_count: validations.reduce((sum, item) => sum + item.row_count, 0),
  all_partitions_verified: true, all_transport_equal: true, all_decoded_rows_equal: true,
  source_coordinate_bytes: sourceResourceBytes, candidate_coordinate_bytes: candidateResourceBytes,
  saved_resource_bytes: sourceResourceBytes - candidateResourceBytes,
  saved_resource_percent: 100 * (sourceResourceBytes - candidateResourceBytes) / sourceResourceBytes,
  candidate_index_bytes: Buffer.byteLength(candidateJson), candidate_bytes_including_index: candidateResourceBytes + Buffer.byteLength(candidateJson),
  source_resource_count: jobs.length, candidate_resource_count: jobs.length,
  maximum_partition_row_count: Math.max(...validations.map(item => item.row_count)),
  totals,
  validation: 'Authenticated source and destination gzip/raw sizes and hashes, Object.is every x/y/z element, full transport metadata and every decoded row/field after gzip reread',
  limitation: 'Resource/index sizes exclude this validation report. Timings include validation on the shared host and are not browser performance evidence. Full release and browser acceptance remain separate gates. No active assets were modified.',
  partitions: validations,
};
await writeFile(path.join(output, 'validation.json'), JSON.stringify(report, null, 2) + '\n', { flag: 'wx' });
const { partitions: verifiedPartitions, ...summary } = report;
console.log(JSON.stringify(summary, null, 2));
