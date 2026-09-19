/** Build and verify a scalar-only storage experiment; never activate a release. */
import assert from 'node:assert/strict';
import { mkdir, readFile, realpath, writeFile } from 'node:fs/promises';
import { createHash } from 'node:crypto';
import { gzipSync, gunzipSync } from 'node:zlib';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { decodeSurveyRows, SURVEY_COLUMNAR_ENCODING } from '../core/survey-codec.js';
import { shareSurveyColumns, expandSharedSurveyColumns, verifySharedColumn } from '../core/shared-survey-codec.js';

const [manifestArgument, outputArgument, ...extra] = process.argv.slice(2);
if (!manifestArgument || !outputArgument || extra.length) {
  throw new Error('Usage: node build_shared_survey_candidate.mjs RELEASE_MANIFEST NEW_OUTPUT_DIRECTORY');
}
const hash = value => createHash('sha256').update(value).digest('hex');
const inside = (root, file) => file === root || file.startsWith(`${root}${path.sep}`);
const manifestFile = await realpath(manifestArgument);
const sourceRoot = path.dirname(manifestFile);
const manifestBytes = await readFile(manifestFile);
const manifest = JSON.parse(manifestBytes);
const terms = Object.entries(manifest.survey?.scalars?.terms ?? {});
if (!manifest.build_id || !terms.length) throw new Error('Expected a release with scalar terms and build_id');
const declaredTerms = manifest.survey.terms.map(term => term.id ?? term.term_id);
assert.deepEqual([...terms.map(([term]) => term)].sort(), [...declaredTerms].sort(), 'All declared terms must have scalar data');

// Resolve the existing parent before ownership acquisition to reject symlink aliases.
const outputRequested = path.resolve(outputArgument);
const outputParent = await realpath(path.dirname(outputRequested));
const output = path.join(outputParent, path.basename(outputRequested));
const publishedAssets = await realpath(fileURLToPath(new URL('../../assets', import.meta.url)));
if (inside(publishedAssets, output) || inside(sourceRoot, output)) {
  throw new Error('Candidate output must be outside published assets and the source release');
}
await mkdir(output); // Exclusive: an existing output is never overwritten or resumed.
await mkdir(path.join(output, 'survey', 'columns'), { recursive: true });
await mkdir(path.join(output, 'survey', 'scalars'));

async function readCompressed(root, descriptor) {
  if (!descriptor || typeof descriptor.path !== 'string' || path.isAbsolute(descriptor.path)) {
    throw new Error('Expected a relative resource path');
  }
  const file = await realpath(path.resolve(root, descriptor.path));
  if (!inside(root, file)) throw new Error('Resource escaped its release directory');
  const compressed = await readFile(file);
  assert.equal(compressed.length, descriptor.bytes, `Compressed size: ${descriptor.path}`);
  assert.equal(hash(compressed), descriptor.sha256, `Compressed checksum: ${descriptor.path}`);
  const raw = gunzipSync(compressed);
  assert.equal(raw.length, descriptor.uncompressed_bytes, `Uncompressed size: ${descriptor.path}`);
  return JSON.parse(raw);
}

async function writeCompressed(relativePath, value) {
  const raw = Buffer.from(JSON.stringify(value));
  const compressed = gzipSync(raw, { level: 9 });
  await writeFile(path.join(output, relativePath), compressed, { flag: 'wx' });
  return { path: relativePath, bytes: compressed.length, uncompressed_bytes: raw.length, sha256: hash(compressed) };
}

const sharedColumns = Object.create(null);
const scalarTerms = Object.create(null);
const validations = [];
let sourceBytes = 0, sourceRows = 0, references = 0;
const startedAt = new Date().toISOString();
for (const [term, descriptor] of terms) {
  if (!/^[a-zA-Z0-9_-]+$/.test(term)) throw new Error(`Unsafe term id: ${term}`);
  const original = await readCompressed(sourceRoot, descriptor);
  assert.equal(original.encoding, SURVEY_COLUMNAR_ENCODING);
  assert.equal(original.build_id, manifest.build_id, `Source build id: ${term}`);
  assert.equal(original.row_count, descriptor.row_count, `Source row count: ${term}`);
  const shared = await shareSurveyColumns(original, async values => {
    // The content identity is exactly the JSON array consumed by the browser.
    const reference = hash(JSON.stringify(values));
    references++;
    if (!Object.hasOwn(sharedColumns, reference)) {
      sharedColumns[reference] = {
        ...await writeCompressed(`survey/columns/${reference}.json.gz`, values),
        row_count: values.length, content_sha256: reference,
      };
    }
    return reference;
  });
  const destination = {
    ...descriptor,
    ...await writeCompressed(`survey/scalars/${term}.json.gz`, shared),
    encoding: shared.encoding,
  };
  scalarTerms[term] = destination;
  // Verify bytes read back from disk, not only the in-memory transform.
  const reread = await readCompressed(output, destination);
  const expanded = await expandSharedSurveyColumns(reread, async reference => {
    if (!Object.hasOwn(sharedColumns, reference)) throw new Error('Unknown shared-column reference');
    return verifySharedColumn(reference, await readCompressed(output, sharedColumns[reference]));
  });
  assert.deepEqual(expanded, original, `All transport metadata and columns: ${term}`);
  const expectedRows = decodeSurveyRows(original);
  const candidateRows = decodeSurveyRows(expanded);
  assert.deepEqual(candidateRows, expectedRows, `All decoded rows and fields: ${term}`);
  const semanticSha256 = hash(JSON.stringify(expectedRows));
  assert.equal(hash(JSON.stringify(candidateRows)), semanticSha256);
  sourceBytes += descriptor.bytes;
  sourceRows += expectedRows.length;
  validations.push({ term, row_count: expectedRows.length, decoded_equal: true, semantic_sha256: semanticSha256 });
  if (validations.length % 10 === 0 || validations.length === terms.length) {
    console.error(`Verified ${validations.length}/${terms.length} scalar terms`);
  }
}

const descriptorBytes = Object.values(scalarTerms).reduce((sum, item) => sum + item.bytes, 0);
const columnBytes = Object.values(sharedColumns).reduce((sum, item) => sum + item.bytes, 0);
const candidate = {
  schema_version: 'rna-survey-scalar-candidate-1', scalar_only: true,
  build_id: manifest.build_id,
  source_manifest: { path: manifestFile, sha256: hash(manifestBytes) },
  generated_at: new Date().toISOString(),
  limitation: 'Scalar-only experiment. Preserves source build identity; not a complete release or an activation manifest.',
  survey: {
    terms: manifest.survey.terms, opening_bins: manifest.survey.opening_bins,
    scalars: { ...manifest.survey.scalars, terms: scalarTerms }, shared_columns: sharedColumns,
  },
};
const candidateJson = JSON.stringify(candidate, null, 2) + '\n';
await writeFile(path.join(output, 'candidate.json'), candidateJson, { flag: 'wx' });
const report = {
  schema_version: 'rna-survey-scalar-candidate-validation-1',
  started_at: startedAt, completed_at: new Date().toISOString(),
  source_build_id: manifest.build_id, source_manifest: candidate.source_manifest,
  output_directory: output, scalar_only: true, all_terms_verified: true,
  term_count: terms.length, row_count: sourceRows,
  original_scalar_bytes: sourceBytes,
  candidate_scalar_descriptor_bytes: descriptorBytes,
  candidate_shared_column_bytes: columnBytes,
  candidate_scalar_resource_bytes: descriptorBytes + columnBytes,
  candidate_index_bytes: Buffer.byteLength(candidateJson),
  candidate_bytes_including_index: descriptorBytes + columnBytes + Buffer.byteLength(candidateJson),
  saved_scalar_resource_bytes: sourceBytes - descriptorBytes - columnBytes,
  shared_column_count: Object.keys(sharedColumns).length,
  column_reference_count: references,
  original_resource_count: terms.length,
  candidate_resource_count: terms.length + Object.keys(sharedColumns).length,
  validation: 'Compressed checksums, content hashes, transport metadata, all decoded rows and fields, per-term semantic SHA256',
  limitation: 'Size accounting includes scalar resources and candidate index, excludes this validation report. Request overhead and browser behavior require separate validation. No active assets were modified.',
  terms: validations,
};
await writeFile(path.join(output, 'validation.json'), JSON.stringify(report, null, 2) + '\n', { flag: 'wx' });
const { terms: verifiedTerms, ...summary } = report;
console.log(JSON.stringify(summary, null, 2));
