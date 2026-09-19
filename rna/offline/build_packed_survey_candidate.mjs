/** Build and verify exact packed Survey scalars without activating a release. */
import assert from 'node:assert/strict';
import { mkdir, readFile, readdir, realpath, writeFile } from 'node:fs/promises';
import { gzipSync, gunzipSync } from 'node:zlib';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { performance } from 'node:perf_hooks';
import { createHash } from 'node:crypto';
import { OutputScope, sha256 } from './output_scope.mjs';
import { decodeSurveyRows } from '../core/survey-codec.js';
import { BUNDLED_SURVEY_ENCODING, verifySurveyBundle, expandBundledSurveyColumns } from '../core/bundled-survey-codec.js';
import { PACKED_SURVEY_ENCODING, encodePackedSurvey, expandPackedSurvey } from '../core/packed-survey-codec.js';

const [sourceArgument, outputArgument, ...extra] = process.argv.slice(2);
if (!sourceArgument || !outputArgument || extra.length) {
  throw new Error('Usage: node build_packed_survey_candidate.mjs RELEASE_MANIFEST NEW_OUTPUT_DIRECTORY');
}
const startedAt = new Date().toISOString(), started = performance.now();
const sourceFile = await realpath(sourceArgument), sourceRoot = path.dirname(sourceFile);
const sourceBytes = await readFile(sourceFile), source = JSON.parse(sourceBytes);
assert.equal(source.schema_version, 'rna-explorer-1');
assert.equal(source.molecule_type, 'RNA');
const isRecord = value => value !== null && typeof value === 'object' && !Array.isArray(value);
assert.ok(typeof source.build_id === 'string' && source.build_id, 'Source build identity');
assert.ok(isRecord(source.survey?.scalars?.terms), 'Survey scalar registry');
assert.ok(isRecord(source.survey.bundles), 'Survey bundle registry');
assert.ok(Array.isArray(source.survey.terms), 'Survey term definitions');
const definitions = new Map(source.survey.terms.map(term => [term.term_id, term]));
assert.equal(definitions.size, source.survey.terms.length, 'Unique term definitions');
const terms = Object.entries(source.survey.scalars.terms);
assert.ok(terms.length, 'Nonempty Survey scalar registry');
const inside = (root, file) => file === root || file.startsWith(`${root}${path.sep}`);
const output = path.resolve(outputArgument), scope = new OutputScope([output]);
const assets = await realpath(fileURLToPath(new URL('../../assets', import.meta.url)));
if ([assets, sourceRoot].some(root => inside(root, output) || inside(output, root))) {
  throw new Error('Candidate output must not overlap published assets or the source release');
}
function safeRelative(relative) {
  if (typeof relative !== 'string' || !relative || path.isAbsolute(relative) || relative.includes('\\')
      || relative.split('/').some(part => !part || part === '.' || part === '..')) throw new Error('Unsafe Survey resource path');
  if (!/^survey\/(scalars|bundles)\/.+\.json\.gz$/.test(relative)) throw new Error('Expected Survey scalar or bundle gzip path');
  return relative;
}
const paths = new Set();
for (const descriptor of [...terms.map(([, descriptor]) => descriptor), ...Object.values(source.survey.bundles)]) {
  safeRelative(descriptor.path);
  assert.ok(!paths.has(descriptor.path), 'Unique resource path');
  paths.add(descriptor.path);
  assert.ok(Number.isSafeInteger(descriptor.bytes) && descriptor.bytes > 0, 'Safe compressed byte count');
  assert.ok(Number.isSafeInteger(descriptor.uncompressed_bytes) && descriptor.uncompressed_bytes > 0, 'Safe raw byte count');
  assert.match(descriptor.sha256, /^[a-f0-9]{64}$/);
}
for (const [termId, descriptor] of terms) {
  assert.match(termId, /^[A-Za-z0-9_-]+$/);
  assert.ok(definitions.has(termId), 'Registered Survey term definition');
  assert.equal(descriptor.encoding, BUNDLED_SURVEY_ENCODING, 'Source Survey encoding');
  assert.ok(Number.isSafeInteger(descriptor.row_count) && descriptor.row_count >= 0, 'Safe source row count');
}
await scope.resolve(output);
await realpath(path.dirname(output));
await mkdir(output); // Exclusive acquisition; never overwrite a previous candidate.
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
function semanticHash(rows) {
  const digest = createHash('sha256');
  const canonicalKeys = (_key, value) => isRecord(value)
    ? Object.fromEntries(Object.keys(value).sort().map(key => [key, value[key]])) : value;
  digest.update('[');
  for (let index = 0; index < rows.length; index++) {
    if (index) digest.update(',');
    digest.update(JSON.stringify(rows[index], canonicalKeys));
  }
  digest.update(']');
  return digest.digest('hex');
}
const bundles = structuredClone(source.survey.bundles), bundleValidations = [];
for (const [reference, descriptor] of Object.entries(bundles)) {
  assert.equal(descriptor.content_sha256, reference, 'Bundle descriptor content identity');
  const original = await readResource(sourceRoot, descriptor);
  await verifySurveyBundle(reference, original.data);
  assert.equal(Object.keys(original.data.columns).length, descriptor.column_count, 'Bundle column count');
  await writeResource(descriptor.path, original.compressed);
  const reread = await readResource(output, descriptor);
  assert.ok(reread.compressed.equals(original.compressed), 'Copied bundle compressed bytes are identical');
  await verifySurveyBundle(reference, reread.data);
  assert.deepEqual(reread.data, original.data, 'Copied bundle payload is identical');
  bundleValidations.push({ reference, path: descriptor.path, bytes: descriptor.bytes, sha256: descriptor.sha256, exact_copy: true });
}
function bundleReader(root) {
  const cache = new Map();
  return async reference => {
    assert.ok(Object.hasOwn(bundles, reference), 'Registered Survey bundle');
    if (!cache.has(reference)) cache.set(reference, (async () => {
      const resource = await readResource(root, bundles[reference]);
      return (await verifySurveyBundle(reference, resource.data)).bundle;
    })());
    return cache.get(reference);
  };
}
const scalars = {}, validations = [];
const totals = { checked_expanded_transport_numbers: 0, checked_decoded_row_numbers: 0,
  packed_value_terms: 0, factored_id_terms: 0, constant_term_terms: 0,
  encode_ms: 0, serialize_compress_write_ms: 0, reread_expand_compare_ms: 0 };
for (const [termId, descriptor] of terms) {
  const termStarted = performance.now();
  const { data: original } = await readResource(sourceRoot, descriptor);
  assert.equal(original.encoding, BUNDLED_SURVEY_ENCODING);
  assert.equal(original.build_id, source.build_id, 'Source Survey build identity');
  assert.equal(original.row_count, descriptor.row_count, 'Source Survey row count');
  const readSourceBundle = bundleReader(sourceRoot);
  const expectedTransport = await expandBundledSurveyColumns(original, readSourceBundle);
  assert.ok(expectedTransport.columns.term_id.every(value => value === termId), 'Every source row belongs to its declared term');
  const encodeStarted = performance.now(), packed = await encodePackedSurvey(original, readSourceBundle);
  const encodeMs = performance.now() - encodeStarted;
  assert.equal(packed.encoding, PACKED_SURVEY_ENCODING);
  assert.equal(packed.build_id, source.build_id, 'Candidate retains source build identity');
  assert.equal(packed.row_count, original.row_count);
  const retained = table => Object.fromEntries(Object.entries(table).filter(([key]) => key !== 'encoding' && key !== 'columns'));
  assert.deepEqual(retained(packed), retained(original), 'All outer metadata preserved');
  for (const [key, column] of Object.entries(original.columns)) {
    if (!['value', 'id', 'term_id'].includes(key)) assert.deepEqual(packed.columns[key], column, `Unchanged column/reference: ${key}`);
  }
  const writeStarted = performance.now();
  const raw = Buffer.from(JSON.stringify(packed)), compressed = gzipSync(raw, { level: 9 });
  await writeResource(descriptor.path, compressed);
  const destination = { ...descriptor, encoding: PACKED_SURVEY_ENCODING,
    bytes: compressed.length, uncompressed_bytes: raw.length, sha256: sha256(compressed) };
  const writeMs = performance.now() - writeStarted, verifyStarted = performance.now();
  const { data: reread } = await readResource(output, destination);
  const actualTransport = await expandPackedSurvey(reread, bundleReader(output));
  assert.deepEqual(actualTransport, expectedTransport, 'Complete expanded Survey transport metadata and columns');
  const transportNumbers = checkNumbers(actualTransport, expectedTransport);
  const expectedRows = decodeSurveyRows(expectedTransport), actualRows = decodeSurveyRows(actualTransport);
  assert.deepEqual(actualRows, expectedRows, `Every decoded Survey row and field: ${termId}`);
  assert.equal(actualRows.length, descriptor.row_count);
  const rowNumbers = checkNumbers(actualRows, expectedRows), semanticSha256 = semanticHash(expectedRows);
  assert.equal(semanticHash(actualRows), semanticSha256, 'Decoded Survey semantic identity');
  const references = new Set(Object.values(original.columns).filter(column => column?.bundle).map(column => column.bundle));
  const actualReferences = new Set(Object.values(packed.columns).filter(column => column?.bundle).map(column => column.bundle));
  assert.deepEqual(actualReferences, references, 'Per-term bundle dependencies preserved');
  const valuesPacked = !Array.isArray(packed.columns.value) && packed.columns.value?.encoding === 'rna-survey-values-float64-1';
  const idFactored = packed.columns.id?.encoding === 'rna-survey-id-suffix-1';
  const termConstant = packed.columns.term_id?.encoding === 'rna-survey-constant-1';
  scalars[termId] = destination;
  validations.push({ term_id: termId, path: descriptor.path, row_count: descriptor.row_count,
    null_count: expectedRows.filter(row => row.value === null).length,
    source_bytes: descriptor.bytes, candidate_bytes: destination.bytes, saved_bytes: descriptor.bytes - destination.bytes,
    source_sha256: descriptor.sha256, candidate_sha256: destination.sha256,
    values_packed: valuesPacked, id_factored: idFactored, term_constant: termConstant,
    checked_expanded_transport_numbers: transportNumbers, checked_decoded_row_numbers: rowNumbers,
    expanded_transport_equal: true, decoded_rows_equal: true, semantic_sha256: semanticSha256,
    cold_resource_count: 1 + references.size,
    timing_ms: { encode: encodeMs, serialize_compress_write: writeMs, reread_expand_compare: performance.now() - verifyStarted,
      term_total: performance.now() - termStarted } });
  totals.checked_expanded_transport_numbers += transportNumbers; totals.checked_decoded_row_numbers += rowNumbers;
  totals.packed_value_terms += Number(valuesPacked); totals.factored_id_terms += Number(idFactored); totals.constant_term_terms += Number(termConstant);
  totals.encode_ms += encodeMs; totals.serialize_compress_write_ms += writeMs;
  totals.reread_expand_compare_ms += validations.at(-1).timing_ms.reread_expand_compare;
  console.error(`Verified ${validations.length}/${terms.length} packed Survey terms: ${termId}`);
}
const candidate = {
  schema_version: 'rna-survey-packed-candidate-1', survey_only: true, build_id: source.build_id,
  source_manifest: { path: sourceFile, sha256: sha256(sourceBytes) }, generated_at: new Date().toISOString(),
  limitation: 'Scalar-only candidate retaining source identity. Not a complete release or activation manifest. Coordinate descriptors are preserved metadata; coordinate resources are not copied.',
  survey: { ...source.survey, scalars: { ...source.survey.scalars, terms: scalars }, bundles },
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
assert.ok((await readFile(sourceFile)).equals(sourceBytes), 'Source manifest remained byte-identical throughout the build');
const size = descriptors => descriptors.reduce((sum, item) => sum + item.bytes, 0);
const sourceScalarBytes = size(terms.map(([, descriptor]) => descriptor)), candidateScalarBytes = size(Object.values(scalars));
const bundleBytes = size(Object.values(bundles));
const report = {
  schema_version: 'rna-survey-packed-candidate-validation-1', survey_only: true,
  started_at: startedAt, completed_at: new Date().toISOString(), elapsed_ms: performance.now() - started,
  source_manifest: candidate.source_manifest, source_build_id: source.build_id, source_manifest_stable: true, output_directory: output,
  term_count: terms.length, bundle_count: Object.keys(bundles).length,
  row_count: validations.reduce((sum, item) => sum + item.row_count, 0), null_count: validations.reduce((sum, item) => sum + item.null_count, 0),
  all_terms_verified: true, all_expanded_transport_equal: true, all_decoded_rows_equal: true, all_bundles_exact_copies: true,
  source_scalar_bytes: sourceScalarBytes, candidate_payload_bytes: candidateScalarBytes,
  source_bundle_bytes: bundleBytes, candidate_bundle_bytes: bundleBytes,
  source_survey_resource_bytes: sourceScalarBytes + bundleBytes, candidate_survey_resource_bytes: candidateScalarBytes + bundleBytes,
  saved_resource_bytes: sourceScalarBytes - candidateScalarBytes,
  saved_resource_percent: 100 * (sourceScalarBytes - candidateScalarBytes) / (sourceScalarBytes + bundleBytes),
  candidate_index_bytes: Buffer.byteLength(candidateJson), candidate_bytes_including_index: candidateScalarBytes + bundleBytes + Buffer.byteLength(candidateJson),
  source_resource_count: paths.size, candidate_resource_count: paths.size,
  minimum_term_cold_resource_count: Math.min(...validations.map(item => item.cold_resource_count)),
  maximum_term_cold_resource_count: Math.max(...validations.map(item => item.cold_resource_count)), totals,
  validation: 'Authenticated source and destination gzip/raw sizes and hashes, exact bundle copies, every decoded field/status/id, Object.is every expanded-transport and decoded-row number, source manifest stability',
  limitation: 'Counts include repeated observations across terms. Numeric scopes overlap and must not be summed. Sizes exclude this validation report. Shared-host build timings are not browser evidence. Full release, scientific validation and browser acceptance remain separate gates. No active assets were modified.',
  terms: validations, bundles: bundleValidations,
};
await writeFile(path.join(output, 'validation.json'), JSON.stringify(report, null, 2) + '\n', { flag: 'wx' });
const { terms: verifiedTerms, bundles: verifiedBundles, ...summary } = report;
console.log(JSON.stringify(summary, null, 2));
