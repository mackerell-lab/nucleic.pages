/** Consolidate shared scalar columns into term-membership bundles; never activate. */
import assert from 'node:assert/strict';
import { mkdir, readFile, realpath, writeFile } from 'node:fs/promises';
import { createHash } from 'node:crypto';
import { gzipSync, gunzipSync } from 'node:zlib';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { decodeSurveyRows } from '../core/survey-codec.js';
import { SHARED_SURVEY_ENCODING, expandSharedSurveyColumns, verifySharedColumn } from '../core/shared-survey-codec.js';
import { BUNDLED_SURVEY_ENCODING, expandBundledSurveyColumns, verifySurveyBundle } from '../core/bundled-survey-codec.js';

const [sourceArgument, outputArgument, ...extra] = process.argv.slice(2);
if (!sourceArgument || !outputArgument || extra.length) {
  throw new Error('Usage: node build_bundled_survey_candidate.mjs SHARED_CANDIDATE NEW_OUTPUT_DIRECTORY');
}
const hash = value => createHash('sha256').update(value).digest('hex');
const inside = (root, file) => file === root || file.startsWith(`${root}${path.sep}`);
const sourceFile = await realpath(sourceArgument);
const sourceRoot = path.dirname(sourceFile);
const sourceBytes = await readFile(sourceFile);
const source = JSON.parse(sourceBytes);
const terms = Object.entries(source.survey?.scalars?.terms ?? {});
if (!source.scalar_only || !source.build_id || !terms.length || !source.survey.shared_columns) {
  throw new Error('Expected a scalar-only shared-column candidate with build_id');
}
assert.deepEqual(terms.map(([term]) => term).sort(), source.survey.terms.map(term => term.id ?? term.term_id).sort(), 'All declared terms must have scalar data');

const outputRequested = path.resolve(outputArgument);
const outputParent = await realpath(path.dirname(outputRequested));
const output = path.join(outputParent, path.basename(outputRequested));
const publishedAssets = await realpath(fileURLToPath(new URL('../../assets', import.meta.url)));
if (inside(publishedAssets, output) || inside(sourceRoot, output)) {
  throw new Error('Candidate output must be outside published assets and the source candidate');
}
await mkdir(output); // Acquire exclusive ownership; never overwrite or resume output.
await mkdir(path.join(output, 'survey', 'bundles'), { recursive: true });
await mkdir(path.join(output, 'survey', 'scalars'));
const startedAt = new Date().toISOString();

async function readCompressed(root, descriptor) {
  if (!descriptor || typeof descriptor.path !== 'string' || path.isAbsolute(descriptor.path)) throw new Error('Expected a relative resource path');
  const file = await realpath(path.resolve(root, descriptor.path));
  if (!inside(root, file)) throw new Error('Resource escaped its candidate directory');
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

async function readColumn(reference) {
  assert.match(reference, /^[a-f0-9]{64}$/);
  assert.ok(Object.hasOwn(source.survey.shared_columns, reference), 'Registered shared column');
  const descriptor = source.survey.shared_columns[reference];
  assert.equal(descriptor.content_sha256, reference, 'Source content identity');
  const values = await verifySharedColumn(reference, await readCompressed(sourceRoot, descriptor));
  assert.equal(values.length, descriptor.row_count, 'Source shared column row count');
  return values;
}

// Keep only tiny reference descriptors while determining exact cross-term reuse.
const payloads = new Map(), memberships = new Map();
for (const [term, descriptor] of terms) {
  if (!/^[a-zA-Z0-9_-]+$/.test(term)) throw new Error(`Unsafe term id: ${term}`);
  const payload = await readCompressed(sourceRoot, descriptor);
  assert.equal(payload.encoding, SHARED_SURVEY_ENCODING);
  assert.equal(payload.build_id, source.build_id, `Source build identity: ${term}`);
  assert.equal(payload.row_count, descriptor.row_count, `Source scalar row count: ${term}`);
  assert.ok(Number.isSafeInteger(payload.row_count) && payload.row_count >= 0, 'Safe source row count');
  assert.ok(payload.columns && typeof payload.columns === 'object' && !Array.isArray(payload.columns), 'Source columns');
  for (const column of Object.values(payload.columns)) {
    assert.deepEqual(Object.keys(column), ['reference'], 'Source reference descriptor');
    assert.match(column.reference, /^[a-f0-9]{64}$/);
    assert.ok(Object.hasOwn(source.survey.shared_columns, column.reference), 'Registered source reference');
    if (!memberships.has(column.reference)) memberships.set(column.reference, new Set());
    memberships.get(column.reference).add(term);
  }
  payloads.set(term, payload);
}
assert.deepEqual([...memberships.keys()].sort(), Object.keys(source.survey.shared_columns).sort(), 'No unreferenced source columns');

const groups = new Map();
for (const [reference, consumers] of memberships) {
  if (consumers.size === 1) continue;
  const key = JSON.stringify([...consumers].sort());
  if (!groups.has(key)) groups.set(key, []);
  groups.get(key).push(reference);
}
const bundles = Object.create(null), bundleForColumn = new Map();
for (const [consumers, references] of [...groups].sort(([a], [b]) => a.localeCompare(b))) {
  const columns = Object.create(null);
  for (const reference of references.sort()) columns[reference] = await readColumn(reference);
  const bundle = { encoding: 'rna-survey-column-bundle-1', columns };
  const reference = hash(JSON.stringify(bundle));
  bundles[reference] = {
    ...await writeCompressed(`survey/bundles/${reference}.json.gz`, bundle),
    content_sha256: reference, column_count: references.length, terms: JSON.parse(consumers),
  };
  for (const column of references) bundleForColumn.set(column, reference);
}

const scalarTerms = Object.create(null), validations = [];
let totalRows = 0, inlineColumns = 0;
for (const [term, descriptor] of terms) {
  const original = payloads.get(term);
  const columns = Object.create(null), inline = new Map();
  for (const [field, { reference }] of Object.entries(original.columns)) {
    if (bundleForColumn.has(reference)) columns[field] = { bundle: bundleForColumn.get(reference), column: reference };
    else {
      if (!inline.has(reference)) inline.set(reference, await readColumn(reference));
      columns[field] = inline.get(reference);
      inlineColumns++;
    }
  }
  const destination = {
    ...descriptor,
    ...await writeCompressed(`survey/scalars/${term}.json.gz`, { ...original, encoding: BUNDLED_SURVEY_ENCODING, columns }),
    encoding: BUNDLED_SURVEY_ENCODING,
  };
  scalarTerms[term] = destination;
  // Verify the bytes just written through the same strict codec used at runtime.
  const candidate = await expandBundledSurveyColumns(await readCompressed(output, destination), async reference => {
    assert.ok(Object.hasOwn(bundles, reference), 'Registered destination bundle');
    return (await verifySurveyBundle(reference, await readCompressed(output, bundles[reference]))).bundle;
  });
  const expected = await expandSharedSurveyColumns(original, readColumn);
  assert.deepEqual(candidate, expected, `All transport metadata and columns: ${term}`);
  const expectedRows = decodeSurveyRows(expected), candidateRows = decodeSurveyRows(candidate);
  assert.deepEqual(candidateRows, expectedRows, `Every decoded row and field: ${term}`);
  const semanticSha256 = hash(JSON.stringify(expectedRows));
  assert.equal(hash(JSON.stringify(candidateRows)), semanticSha256);
  const bundleCount = new Set(Object.values(columns).filter(value => !Array.isArray(value)).map(value => value.bundle)).size;
  validations.push({ term, row_count: candidateRows.length, decoded_equal: true, semantic_sha256: semanticSha256, cold_resource_count: 1 + bundleCount });
  totalRows += candidateRows.length;
  if (validations.length % 10 === 0 || validations.length === terms.length) console.error(`Verified ${validations.length}/${terms.length} bundled scalar terms`);
}

const candidate = {
  schema_version: 'rna-survey-bundled-candidate-1', scalar_only: true,
  build_id: source.build_id, source_manifest: source.source_manifest,
  source_candidate: { path: sourceFile, sha256: hash(sourceBytes) },
  generated_at: new Date().toISOString(),
  limitation: 'Scalar-only experiment. Preserves source build identity; not a complete release or an activation manifest.',
  survey: {
    terms: source.survey.terms, opening_bins: source.survey.opening_bins,
    scalars: { ...source.survey.scalars, terms: scalarTerms }, bundles,
  },
};
const candidateJson = JSON.stringify(candidate, null, 2) + '\n';
await writeFile(path.join(output, 'candidate.json'), candidateJson, { flag: 'wx' });
const bytes = descriptors => Object.values(descriptors).reduce((sum, value) => sum + value.bytes, 0);
const priorBytes = bytes(source.survey.scalars.terms) + bytes(source.survey.shared_columns);
const scalarBytes = bytes(scalarTerms), bundleBytes = bytes(bundles);
const report = {
  schema_version: 'rna-survey-bundled-candidate-validation-1',
  started_at: startedAt, completed_at: new Date().toISOString(),
  source_build_id: source.build_id, source_candidate: candidate.source_candidate,
  output_directory: output, scalar_only: true, all_terms_verified: true,
  term_count: terms.length, row_count: totalRows,
  source_shared_resource_bytes: priorBytes,
  candidate_scalar_bytes: scalarBytes, candidate_bundle_bytes: bundleBytes,
  candidate_scalar_resource_bytes: scalarBytes + bundleBytes,
  candidate_index_bytes: Buffer.byteLength(candidateJson),
  candidate_bytes_including_index: scalarBytes + bundleBytes + Buffer.byteLength(candidateJson),
  saved_from_shared_resource_bytes: priorBytes - scalarBytes - bundleBytes,
  source_resource_count: terms.length + memberships.size,
  candidate_resource_count: terms.length + Object.keys(bundles).length,
  bundle_count: Object.keys(bundles).length, inline_column_occurrences: inlineColumns,
  minimum_term_cold_resource_count: Math.min(...validations.map(term => term.cold_resource_count)),
  maximum_term_cold_resource_count: Math.max(...validations.map(term => term.cold_resource_count)),
  validation: 'Compressed checksums, all shared column and bundle content hashes, transport metadata, every decoded row and field, semantic SHA256 per term',
  limitation: 'Scalar resources and candidate index only; excludes this validation report. Browser request behavior requires separate validation. No active assets were modified.',
  terms: validations,
};
await writeFile(path.join(output, 'validation.json'), JSON.stringify(report, null, 2) + '\n', { flag: 'wx' });
const { terms: verifiedTerms, ...summary } = report;
console.log(JSON.stringify(summary, null, 2));
