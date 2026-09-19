import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import os from 'node:os';
import path from 'node:path';
import {gzipSync} from 'node:zlib';
import {createHash} from 'node:crypto';
import {validateRelease} from '../offline/assets.mjs';
import {encodeSurveyRows} from '../core/survey-codec.js';
import {BUNDLED_SURVEY_ENCODING, SURVEY_BUNDLE_ENCODING} from '../core/bundled-survey-codec.js';

const hash = bytes => createHash('sha256').update(bytes).digest('hex');

async function fixture(t) {
  const directory = await fs.mkdtemp(path.join(os.tmpdir(), 'rna-release-integrity-'));
  const root = path.join(directory, 'release');
  await fs.mkdir(root);
  t.after(() => fs.rm(directory, {recursive: true, force: true}));
  const write = async (relative, value) => {
    const file = path.join(root, relative), raw = Buffer.from(JSON.stringify(value)), bytes = gzipSync(raw);
    await fs.mkdir(path.dirname(file), {recursive: true});
    await fs.writeFile(file, bytes);
    return {path: relative, sha256: hash(bytes), bytes: bytes.length, uncompressed_bytes: raw.length};
  };
  const scalar = encodeSurveyRows([{
    id: 'angle:residue', term_id: 'angle', pdb_id: 'TEST', entity_id: '1',
    residue_id: 'residue', endpoint_entities: [{pdb_id: 'TEST', entity_id: '1'}],
    is_terminal_any: false, value: 1.2345678901234567, status: 'ok',
  }], 'integrity-release');
  scalar.encoding = BUNDLED_SURVEY_ENCODING;
  const values = scalar.columns.value, column = hash(JSON.stringify(values));
  const bundle = {encoding: SURVEY_BUNDLE_ENCODING, columns: {[column]: values}};
  const reference = hash(JSON.stringify(bundle)), bundlePath = `survey/bundles/${reference}.json.gz`;
  scalar.columns.value = {bundle: reference, column};
  const manifest = {
    build_id: 'integrity-release', partial: false,
    counts: {entries: 1, entities: 1, residues: 1}, source: {candidate_count: 1},
    metadata: await write('metadata.json.gz', {
      entries: [{pdb_id: 'TEST'}], entities: [{pdb_id: 'TEST', entity_id: '1'}],
    }),
    families: [{...await write('backbone.json.gz', [{id: 'residue', pdb_id: 'TEST',
      entity_id: '1', is_terminal_any: false, values: {}, statuses: {}}]),
    id: 'backbone', level: 'residue', parameters: [], row_count: 1}],
    relations: {}, survey: {scalars: {terms: {}}, coordinates: {groups: {}},
      bundles: {[reference]: {...await write(bundlePath, bundle), content_sha256: reference}}},
    provenance: {decisions: await write('decisions.json.gz', [{pdb_id: 'TEST', accepted: true}])},
  };
  const writeScalar = async () => {
    manifest.survey.scalars.terms.angle = {...await write('survey/scalars/angle.json.gz', scalar),
      encoding: BUNDLED_SURVEY_ENCODING, row_count: 1};
  };
  const save = () => fs.writeFile(path.join(root, 'manifest.json'), JSON.stringify(manifest));
  await writeScalar();
  await save();
  return {root, directory, scalar, bundle, reference, bundlePath, manifest, write, writeScalar, save,
    validate: () => validateRelease(path.join(root, 'manifest.json'))};
}

test('full release validates authenticated bundle inventory and exact transport sizes', async t => {
  const f = await fixture(t);
  assert.equal((await f.validate()).ok, true);
});

test('release rejects rehashed wrong or missing build identities', async t => {
  const f = await fixture(t);
  f.scalar.build_id = 'different-release';
  await f.writeScalar(); await f.save();
  await assert.rejects(f.validate(), /build ID mismatch/);
  delete f.scalar.build_id;
  await f.writeScalar(); await f.save();
  await assert.rejects(f.validate(), /build ID mismatch/);
  f.manifest.partial = true;
  await f.save();
  assert.equal((await f.validate()).ok, true, 'Legacy partial encoded payloads may omit build identity');
  f.scalar.build_id = 'different-release';
  await f.writeScalar(); await f.save();
  await assert.rejects(f.validate(), /build ID mismatch/, 'A declared build ID must match even in partial releases');
});

test('release rejects descriptor encoding mismatch', async t => {
  const f = await fixture(t);
  f.manifest.survey.scalars.terms.angle.encoding = 'rna-survey-columnar-1';
  await f.save();
  await assert.rejects(f.validate(), /encoding mismatch/);
  delete f.manifest.survey.scalars.terms.angle.encoding;
  await f.save();
  await assert.rejects(f.validate(), /encoding mismatch/);
});

test('release rejects compressed and uncompressed asset byte mismatches', async t => {
  const f = await fixture(t), descriptor = f.manifest.metadata;
  for (const field of ['bytes', 'uncompressed_bytes']) {
    const value = descriptor[field];
    descriptor[field]++;
    await f.save();
    await assert.rejects(f.validate(), /byte size mismatch/);
    descriptor[field] = value;
  }
});

test('full release requires every referenced bundle in its registry', async t => {
  const f = await fixture(t);
  delete f.manifest.survey.bundles[f.reference];
  await f.save();
  await assert.rejects(f.validate(), /Unregistered Survey bundle/);
});

test('release checks bundle registry path, hashes, and both byte sizes', async t => {
  const f = await fixture(t), descriptor = f.manifest.survey.bundles[f.reference];
  for (const [field, invalid, error] of [
    ['path', 'elsewhere.json.gz', /registry descriptor/],
    ['content_sha256', '0'.repeat(64), /registry descriptor/],
    ['sha256', '0'.repeat(64), /Asset hash mismatch/],
    ['bytes', descriptor.bytes + 1, /compressed byte size mismatch/],
    ['uncompressed_bytes', descriptor.uncompressed_bytes + 1, /uncompressed byte size mismatch/],
    ['uncompressed_bytes', undefined, /registry descriptor/],
  ]) {
    const value = descriptor[field];
    descriptor[field] = invalid;
    await f.save();
    await assert.rejects(f.validate(), error, field);
    descriptor[field] = value;
  }
});

test('bundle semantic hash survives rehashed registry corruption', async t => {
  const f = await fixture(t);
  f.bundle.columns[f.scalar.columns.value.column] = [99];
  f.manifest.survey.bundles[f.reference] = {
    ...await f.write(f.bundlePath, f.bundle), content_sha256: f.reference,
  };
  await f.save();
  await assert.rejects(f.validate(), /Survey bundle hash mismatch/);
});

test('release validates unused registry bundles before rejecting their presence', async t => {
  const f = await fixture(t), descriptor = f.manifest.survey.bundles[f.reference];
  f.scalar.columns.value = [1.2345678901234567];
  await f.writeScalar(); await f.save();
  await assert.rejects(f.validate(), /Unreferenced Survey bundle/);
  descriptor.sha256 = '0'.repeat(64);
  await f.save();
  await assert.rejects(f.validate(), /Asset hash mismatch/);
});

test('release rejects missing bundle files and retries after restoration', async t => {
  const f = await fixture(t), file = path.join(f.root, f.bundlePath);
  await fs.rename(file, `${file}.saved`);
  await assert.rejects(f.validate(), {code: 'ENOENT'});
  await fs.rename(`${file}.saved`, file);
  assert.equal((await f.validate()).ok, true);
});

test('release rejects symlinks escaping its root despite valid bytes', async t => {
  const f = await fixture(t), file = path.join(f.root, f.bundlePath);
  const outside = path.join(f.directory, 'outside.json.gz');
  await fs.rename(file, outside);
  await fs.symlink(outside, file);
  await assert.rejects(f.validate(), /escapes release root/);
});
