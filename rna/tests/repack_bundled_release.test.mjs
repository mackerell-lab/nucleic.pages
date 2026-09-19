import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import path from 'node:path';
import os from 'node:os';
import {fileURLToPath} from 'node:url';
import {promisify} from 'node:util';
import {execFile} from 'node:child_process';
import {gzipSync, gunzipSync} from 'node:zlib';
import {sha256} from '../offline/output_scope.mjs';
import {validateRelease} from '../offline/assets.mjs';
import {verifyReleaseInventory} from '../offline/verify_release_inventory.mjs';
import {encodeFamilyRows, encodeSurveyRows, decodeFamilyRows, decodeSurveyRows} from '../core/survey-codec.js';
import {BUNDLED_FAMILY_ENCODING, FAMILY_BUNDLE_ENCODING, expandBundledFamilyColumns, verifyFamilyBundle} from '../core/bundled-family-codec.js';
import {BUNDLED_SURVEY_ENCODING, SURVEY_BUNDLE_ENCODING, expandBundledSurveyColumns, verifySurveyBundle} from '../core/bundled-survey-codec.js';

const execute = promisify(execFile);
const script = fileURLToPath(new URL('../offline/repack_bundled_release.mjs', import.meta.url));
const readJson = async file => JSON.parse(await fs.readFile(file));
const readPacked = async (root, descriptor) => JSON.parse(gunzipSync(await fs.readFile(path.join(root, descriptor.path))));
async function writePacked(root, relative, data) {
  const file = path.join(root, relative), raw = Buffer.from(JSON.stringify(data)), bytes = gzipSync(raw);
  await fs.mkdir(path.dirname(file), {recursive: true});
  await fs.writeFile(file, bytes);
  return {path: relative, sha256: sha256(bytes), bytes: bytes.length, uncompressed_bytes: raw.length};
}
async function bundleTransport(root, packed, kind, fields) {
  const family = kind === 'family', prefix = family ? 'families' : 'survey';
  const columns = {}, refs = new Map();
  for (const field of fields) {
    const column = packed.columns[field], reference = sha256(JSON.stringify(column));
    columns[reference] = column; refs.set(field, reference);
  }
  const bundle = {encoding: family ? FAMILY_BUNDLE_ENCODING : SURVEY_BUNDLE_ENCODING, columns};
  const reference = sha256(JSON.stringify(bundle));
  const descriptor = {...await writePacked(root, `${prefix}/bundles/${reference}.json.gz`, bundle),
    content_sha256: reference, column_count: Object.keys(columns).length};
  const result = structuredClone(packed);
  result.encoding = family ? BUNDLED_FAMILY_ENCODING : BUNDLED_SURVEY_ENCODING;
  for (const [field, column] of refs) result.columns[field] = {bundle: reference, column};
  return {packed: result, registry: {[reference]: descriptor}};
}
async function expand(root, packed, manifest, kind) {
  const family = kind === 'family', registry = family ? manifest.family_bundles : manifest.survey.bundles;
  if (packed.encoding !== (family ? BUNDLED_FAMILY_ENCODING : BUNDLED_SURVEY_ENCODING)) return packed;
  const expandColumns = family ? expandBundledFamilyColumns : expandBundledSurveyColumns;
  const verify = family ? verifyFamilyBundle : verifySurveyBundle;
  return expandColumns(packed, async reference => (await verify(reference, await readPacked(root, registry[reference]))).bundle);
}

async function fixture(t, {familyBundled = true} = {}) {
  const directory = await fs.mkdtemp(path.join(os.tmpdir(), 'rna-repack-bundles-'));
  t.after(() => fs.rm(directory, {recursive: true, force: true}));
  const sourceRoot = path.join(directory, 'source'); await fs.mkdir(sourceRoot);
  const rows = ['r1', 'r2'].map(id => ({id, pdb_id: 'TEST', entity_id: '1', model_id: '1',
    is_terminal_any: false, values: {chi: 1.2345678901234567}, statuses: {chi: 'available'}}));
  const scalars = rows.map(row => ({id: `angle:${row.id}`, term_id: 'angle', pdb_id: 'TEST', entity_id: '1',
    residue_id: row.id, model_id: '1', endpoint_entities: [{pdb_id: 'TEST', entity_id: '1'}],
    is_terminal_any: false, value: 7.123456789012345, status: 'ok'}));
  const familyPlain = encodeFamilyRows(rows, 'source-build');
  const family = familyBundled ? await bundleTransport(sourceRoot, familyPlain, 'family', ['id', 'pdb_id']) : {packed: familyPlain};
  const survey = await bundleTransport(sourceRoot, encodeSurveyRows(scalars, 'source-build'), 'scalar', ['id']);
  const manifest = {schema_version: 'rna-explorer-1', molecule_type: 'RNA', build_id: 'source-build', partial: false,
    counts: {entries: 1, entities: 1, residues: 2}, source: {candidate_count: 1},
    metadata: await writePacked(sourceRoot, 'metadata.json.gz', {entries: [{pdb_id: 'TEST', selected_model_id: '1'}], entities: [{pdb_id: 'TEST', entity_id: '1'}]}),
    families: [{...await writePacked(sourceRoot, 'families/backbone.json.gz', family.packed),
      id: 'backbone', label: 'Backbone', level: 'residue', parameters: [{id: 'chi', unit: 'degrees'}], row_count: 2, encoding: family.packed.encoding}],
    ...(familyBundled ? {family_bundles: family.registry} : {}), relations: {},
    survey: {terms: [{id: 'angle', label: 'Angle'}], opening_bins: [], bundles: survey.registry,
      scalars: {terms: {angle: {...await writePacked(sourceRoot, 'survey/scalars/angle.json.gz', survey.packed), row_count: 2, encoding: BUNDLED_SURVEY_ENCODING}}},
      coordinates: {groups: {}}},
    provenance: {decisions: await writePacked(sourceRoot, 'decisions.json.gz', [{pdb_id: 'TEST', accepted: true}]),
      build_stages: {geometry: {status: 'complete', signature: 'trusted-source'}},
      repack: {operation: 'previous-storage-transform', evidence: {path: 'historical-report.json'}},
      source_release: {build_id: 'original-science', repack: {operation: 'earlier-transform'}}},
  };
  const sourceFile = path.join(sourceRoot, 'manifest.json'); await fs.writeFile(sourceFile, JSON.stringify(manifest));
  return {directory, sourceRoot, sourceFile, manifest, rows, scalars};
}

async function candidateFor(sourceFile, directory, kind) {
  await fs.mkdir(directory);
  const sourceRoot = path.dirname(sourceFile), sourceBytes = await fs.readFile(sourceFile), source = JSON.parse(sourceBytes);
  const family = kind === 'family';
  const sourceDescriptor = family ? source.families[0] : source.survey.scalars.terms.angle;
  const original = await expand(sourceRoot, await readPacked(sourceRoot, sourceDescriptor), source, kind);
  const bundled = await bundleTransport(directory, original, kind, Object.keys(original.columns));
  const descriptor = {...sourceDescriptor, ...await writePacked(directory, sourceDescriptor.path, bundled.packed),
    encoding: family ? BUNDLED_FAMILY_ENCODING : BUNDLED_SURVEY_ENCODING};
  const candidate = {schema_version: family ? 'rna-family-bundled-candidate-1' : 'rna-survey-bundled-candidate-1',
    [family ? 'family_only' : 'scalar_only']: true, build_id: source.build_id,
    source_manifest: {path: sourceFile, sha256: sha256(sourceBytes)},
    ...(family ? {families: [descriptor], family_bundles: bundled.registry}
      : {survey: {...source.survey, scalars: {terms: {angle: descriptor}}, bundles: bundled.registry}})};
  const file = path.join(directory, 'candidate.json');
  await fs.writeFile(file, JSON.stringify(candidate));
  return {file, candidate, packed: bundled.packed, descriptor};
}

async function run(sourceFile, candidateFile, output, buildId) {
  const result = await execute(process.execPath, [script, sourceFile, candidateFile, output, buildId], {maxBuffer: 2 * 1024 * 1024});
  return JSON.parse(result.stdout);
}

test('family repack retains Survey bundles, complete rows, and nested provenance', async t => {
  const f = await fixture(t), c = await candidateFor(f.sourceFile, path.join(f.directory, 'family-candidate'), 'family');
  const output = path.join(f.directory, 'family-release');
  const report = await run(f.sourceFile, c.file, output, 'family-build');
  const file = path.join(output, 'manifest.json'), manifest = await readJson(file);
  assert.equal(report.all_decoded_equal, true); assert.equal(report.candidate_kind, 'family');
  assert.ok(report.family_candidate); assert.equal(report.scalar_candidate, undefined);
  assert.equal(manifest.provenance.repack.operation, 'lossless_bundled_family_transport');
  assert.equal(manifest.provenance.repack.family_encoding, BUNDLED_FAMILY_ENCODING);
  assert.deepEqual(manifest.provenance.source_release.repack, f.manifest.provenance.repack);
  assert.deepEqual(manifest.provenance.source_release.source_release, f.manifest.provenance.source_release);
  assert.deepEqual(manifest.provenance.source_release.build_stages, f.manifest.provenance.build_stages);
  assert.deepEqual(manifest.survey.bundles, f.manifest.survey.bundles);
  for (const descriptor of Object.values(manifest.survey.bundles)) {
    assert.deepEqual(await fs.readFile(path.join(output, descriptor.path)), await fs.readFile(path.join(f.sourceRoot, descriptor.path)));
  }
  assert.deepEqual(decodeFamilyRows(await expand(output, await readPacked(output, manifest.families[0]), manifest, 'family')), f.rows);
  assert.deepEqual(decodeSurveyRows(await expand(output, await readPacked(output, manifest.survey.scalars.terms.angle), manifest, 'scalar')), f.scalars);
  assert.equal((await validateRelease(file)).ok, true); assert.equal((await verifyReleaseInventory(file)).ok, true);
  await assert.rejects(run(f.sourceFile, c.file, output, 'family-build'), /Report already exists|EEXIST/);

  // A subsequent scalar repack starts from both bundled formats and preserves
  // the family registry plus the complete two-level storage history.
  const next = await candidateFor(file, path.join(f.directory, 'scalar-candidate'), 'scalar');
  const nextOutput = path.join(f.directory, 'scalar-release');
  const nextReport = await run(file, next.file, nextOutput, 'scalar-build');
  const nextFile = path.join(nextOutput, 'manifest.json'), nextManifest = await readJson(nextFile);
  assert.equal(nextReport.candidate_kind, 'scalar'); assert.ok(nextReport.scalar_candidate);
  assert.deepEqual(nextManifest.family_bundles, manifest.family_bundles);
  assert.deepEqual(nextManifest.provenance.source_release.repack, manifest.provenance.repack);
  assert.deepEqual(nextManifest.provenance.source_release.source_release, manifest.provenance.source_release);
  assert.equal((await validateRelease(nextFile)).ok, true); assert.equal((await verifyReleaseInventory(nextFile)).ok, true);
  assert.deepEqual(decodeFamilyRows(await expand(nextOutput, await readPacked(nextOutput, nextManifest.families[0]), nextManifest, 'family')), f.rows);
  assert.deepEqual(decodeSurveyRows(await expand(nextOutput, await readPacked(nextOutput, nextManifest.survey.scalars.terms.angle), nextManifest, 'scalar')), f.scalars);
});

test('family repack accepts plain families beside existing bundled Survey resources', async t => {
  const f = await fixture(t, {familyBundled: false});
  assert.equal((await validateRelease(f.sourceFile)).ok, true);
  assert.equal((await verifyReleaseInventory(f.sourceFile)).ok, true);
  const c = await candidateFor(f.sourceFile, path.join(f.directory, 'candidate'), 'family');
  const output = path.join(f.directory, 'release');
  await run(f.sourceFile, c.file, output, 'first-family-bundling');
  const file = path.join(output, 'manifest.json'), manifest = await readJson(file);
  assert.deepEqual(manifest.survey.bundles, f.manifest.survey.bundles);
  assert.equal(manifest.families[0].encoding, BUNDLED_FAMILY_ENCODING);
  assert.equal((await validateRelease(file)).ok, true);
  assert.equal((await verifyReleaseInventory(file)).ok, true);
  assert.deepEqual(decodeFamilyRows(await expand(output, await readPacked(output, manifest.families[0]), manifest, 'family')), f.rows);
});

test('repack rejects source identity, candidate scope, incomplete families and unsafe destinations', async t => {
  const f = await fixture(t), c = await candidateFor(f.sourceFile, path.join(f.directory, 'candidate'), 'family');
  const initial = structuredClone(c.candidate);
  const cases = [
    [candidate => {candidate.scalar_only = true;}, /Exactly one/],
    [candidate => {candidate.schema_version = 'unknown';}, /Candidate schema/],
    [candidate => {candidate.build_id = 'wrong';}, /Candidate source identity/],
    [candidate => {candidate.source_manifest.sha256 = '0'.repeat(64);}, /Candidate source manifest identity/],
    [candidate => {candidate.families = [];}, /Complete family registry/],
    [candidate => {candidate.families.push(candidate.families[0]);}, /Unique candidate family IDs/],
  ];
  for (let index = 0; index < cases.length; index++) {
    const [mutate, error] = cases[index], candidate = structuredClone(initial);
    mutate(candidate); await fs.writeFile(c.file, JSON.stringify(candidate));
    const output = path.join(f.directory, `rejected-${index}`);
    await assert.rejects(run(f.sourceFile, c.file, output, 'new-build'), error);
    await assert.rejects(fs.access(output), {code: 'ENOENT'});
  }
  await fs.writeFile(c.file, JSON.stringify(initial));
  await assert.rejects(run(f.sourceFile, c.file, path.join(f.directory, 'dotted-build'), 'new.build'), /Unsafe new build ID/);
  await assert.rejects(run(f.sourceFile, c.file, path.join(f.directory, 'same-build'), f.manifest.build_id), /New immutable build identity/);
  await assert.rejects(run(f.sourceFile, c.file, path.join(f.sourceRoot, 'nested'), 'new-build'), /must not overlap/);
  await assert.rejects(run(f.sourceFile, c.file, path.join(path.dirname(c.file), 'nested'), 'new-build'), /must not overlap/);
});

test('repack rejects rehashed numerical corruption and changed scientific descriptors', async t => {
  for (const kind of ['family', 'scalar']) {
    const f = await fixture(t), directory = path.join(f.directory, 'candidate');
    const c = await candidateFor(f.sourceFile, directory, kind);
    const key = kind === 'family' ? 'values' : 'value';
    c.packed.columns[key] = kind === 'family' ? [{chi: 99}, {chi: 99}] : [99, 99];
    Object.assign(c.descriptor, await writePacked(directory, c.descriptor.path, c.packed));
    await fs.writeFile(c.file, JSON.stringify(c.candidate));
    await assert.rejects(run(f.sourceFile, c.file, path.join(f.directory, 'corrupted'), 'new-build'), /All original .* transport metadata and columns/);
  }
  const f = await fixture(t), c = await candidateFor(f.sourceFile, path.join(f.directory, 'candidate'), 'family');
  c.descriptor.parameters = [{id: 'chi', unit: 'radians'}];
  await fs.writeFile(c.file, JSON.stringify(c.candidate));
  await assert.rejects(run(f.sourceFile, c.file, path.join(f.directory, 'changed-descriptor'), 'new-build'), /Unchanged candidate scientific descriptor/);
});
