import test from 'node:test';
import assert from 'node:assert/strict';
import { mkdtemp, mkdir, readFile, writeFile, rm, symlink, readdir } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { spawnSync } from 'node:child_process';
import { gzipSync, gunzipSync } from 'node:zlib';
import { createHash } from 'node:crypto';
import { expandPackedSurvey } from '../core/packed-survey-codec.js';
import { decodeSurveyRows } from '../core/survey-codec.js';
const tool = fileURLToPath(new URL('../offline/build_packed_survey_candidate.mjs', import.meta.url));
const hash = value => createHash('sha256').update(value).digest('hex');

async function fixture(t) {
  const root = await mkdtemp(path.join(tmpdir(), 'rna-packed-survey-candidate-'));
  t.after(() => rm(root, { recursive: true, force: true }));
  const source = path.join(root, 'source'); await mkdir(source);
  async function resource(relative, data, negativeZero = false) {
    let text = JSON.stringify(data);
    if (negativeZero) text = text.replace('"value":[0,', '"value":[-0,');
    const raw = Buffer.from(text), bytes = gzipSync(raw, { level: 9 });
    await mkdir(path.dirname(path.join(source, relative)), { recursive: true });
    await writeFile(path.join(source, relative), bytes);
    return { path: relative, bytes: bytes.length, uncompressed_bytes: raw.length, sha256: hash(bytes) };
  }
  const observations = ['entry:A:1', 'entry:A:2', 'entry:A:3', 'entry:A:4'];
  const column = hash(JSON.stringify(observations));
  const bundle = { encoding: 'rna-survey-column-bundle-1', columns: { [column]: observations } };
  const reference = hash(JSON.stringify(bundle));
  const bundleDescriptor = { ...await resource(`survey/bundles/${reference}.json.gz`, bundle), content_sha256: reference,
    column_count: 1, terms: ['a_term', 'b_term'] };
  const scalars = {}, expected = {};
  for (const term of ['a_term', 'b_term']) {
    const ids = observations.map((value, index) => term === 'a_term' ? `${value}|survey|${term}` : `provider:${index}`);
    const inline = { id: ids, term_id: observations.map(() => term), value: [-0, null, -180.00000000000003, 1e-300],
      status: ['ok', 'missing_atom', 'ok', 'ok'], observation_id: observations, opening: [null, -10, 2, 180], source_atom_pattern: observations.map(() => ["C2'", "O2'"]) };
    const payload = { encoding: 'rna-survey-bundled-columns-1', build_id: 'fixture', row_count: 4,
      columns: { ...inline, observation_id: { bundle: reference, column } }, provenance: { provider: 'independent-fixture' } };
    scalars[term] = { ...await resource(`survey/scalars/${term}.json.gz`, payload, true), row_count: 4,
      encoding: payload.encoding, fixture_descriptor: { retained: true } };
    expected[term] = observations.map((_value, index) => Object.fromEntries(Object.entries(inline).map(([key, values]) => [key, values[index]])));
  }
  const manifest = { schema_version: 'rna-explorer-1', molecule_type: 'RNA', build_id: 'fixture', survey: {
    terms: [{ term_id: 'a_term', unit: 'deg' }, { term_id: 'b_term', unit: 'angstrom' }], opening_bins: [{ id: 'middle', min: -8, max: 2 }],
    scalars: { terms: scalars, policy: 'one-term' }, bundles: { [reference]: bundleDescriptor }, coordinates: { groups: {} } } };
  const manifestFile = path.join(source, 'manifest.json');
  const save = () => writeFile(manifestFile, JSON.stringify(manifest)); await save();
  const run = (output = path.join(root, 'candidate')) => spawnSync(process.execPath, [tool, manifestFile, output], { encoding: 'utf8', timeout: 30000 });
  return { root, source, manifestFile, manifest, expected, run, save, resource, reference };
}

test('candidate authenticates complete resources, exact values and metadata, and preserves fallback IDs', async t => {
  const f = await fixture(t), before = await readFile(f.manifestFile);
  const result = f.run(); assert.equal(result.status, 0, result.stderr);
  const destination = path.join(f.root, 'candidate');
  const index = JSON.parse(await readFile(path.join(destination, 'candidate.json')));
  const report = JSON.parse(await readFile(path.join(destination, 'validation.json')));
  assert.equal(index.survey_only, true); assert.equal(index.schema_version, 'rna-survey-packed-candidate-1');
  assert.deepEqual(index.survey.terms, f.manifest.survey.terms);
  assert.deepEqual(index.survey.opening_bins, f.manifest.survey.opening_bins);
  assert.deepEqual(index.survey.coordinates, f.manifest.survey.coordinates);
  assert.deepEqual(index.survey.bundles, f.manifest.survey.bundles);
  assert.equal(index.survey.scalars.policy, 'one-term');
  assert.equal(report.term_count, 2); assert.equal(report.bundle_count, 1); assert.equal(report.row_count, 8);
  assert.equal(report.null_count, 2); assert.equal(report.totals.factored_id_terms, 1);
  assert.equal(report.source_resource_count, 3); assert.equal(report.candidate_resource_count, 3);
  assert.equal(report.minimum_term_cold_resource_count, 2); assert.equal(report.maximum_term_cold_resource_count, 2);
  assert.equal(report.source_manifest_stable, true); assert.ok((await readFile(f.manifestFile)).equals(before));
  for (const [term, descriptor] of Object.entries(index.survey.scalars.terms)) {
    assert.deepEqual(descriptor.fixture_descriptor, { retained: true });
    const packed = JSON.parse(gunzipSync(await readFile(path.join(destination, descriptor.path))));
    const decoded = decodeSurveyRows(await expandPackedSurvey(packed, async reference => JSON.parse(gunzipSync(await readFile(path.join(destination, index.survey.bundles[reference].path))))));
    assert.deepEqual(decoded, f.expected[term]); assert.ok(Object.is(decoded[0].value, -0));
  }
  assert.equal((await readdir(destination)).includes('manifest.json'), false, 'candidate cannot masquerade as activation manifest');
});

test('candidate refuses existing output without modifying its sentinel', async t => {
  const f = await fixture(t), output = path.join(f.root, 'candidate'); await mkdir(output);
  await writeFile(path.join(output, 'sentinel'), 'keep');
  const result = f.run(output); assert.notEqual(result.status, 0); assert.match(result.stderr, /EEXIST/);
  assert.equal(await readFile(path.join(output, 'sentinel'), 'utf8'), 'keep');
});

test('candidate rejects corrupted scalar and bundle bytes before success', async t => {
  for (const kind of ['scalar', 'bundle']) {
    const f = await fixture(t), descriptor = kind === 'scalar' ? f.manifest.survey.scalars.terms.a_term : f.manifest.survey.bundles[f.reference];
    const file = path.join(f.source, descriptor.path), bytes = await readFile(file); bytes[bytes.length - 1] ^= 1; await writeFile(file, bytes);
    const result = f.run(); assert.notEqual(result.status, 0); assert.match(result.stderr, /Compressed checksum/);
    await assert.rejects(readFile(path.join(f.root, 'candidate', 'validation.json')), { code: 'ENOENT' });
  }
});

test('candidate rejects descriptor rows assigned to the wrong declared term', async t => {
  const f = await fixture(t), descriptor = f.manifest.survey.scalars.terms.a_term;
  const payload = JSON.parse(gunzipSync(await readFile(path.join(f.source, descriptor.path)))); payload.columns.term_id[1] = 'wrong_term';
  Object.assign(descriptor, await f.resource(descriptor.path, payload)); await f.save();
  const result = f.run(); assert.notEqual(result.status, 0); assert.match(result.stderr, /Every source row belongs to its declared term/);
});

test('candidate rejects path traversal, source overlap and symlink resources', async t => {
  const f = await fixture(t);
  assert.match(f.run(path.join(f.source, 'candidate')).stderr, /must not overlap/);
  const descriptor = f.manifest.survey.scalars.terms.a_term, original = descriptor.path;
  descriptor.path = 'survey/scalars/../escape.json.gz'; await f.save();
  assert.match(f.run().stderr, /Unsafe Survey resource path/);
  descriptor.path = original; await f.save();
  const file = path.join(f.source, original), replacement = path.join(f.root, 'outside.json.gz');
  await writeFile(replacement, await readFile(file)); await rm(file); await symlink(replacement, file);
  assert.match(f.run().stderr, /Resource path or symlink escaped/);
});
