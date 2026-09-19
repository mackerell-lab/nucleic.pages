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
  const root = await fs.mkdtemp(path.join(os.tmpdir(), 'rna-bundled-validation-'));
  t.after(() => fs.rm(root, {recursive: true, force: true}));
  const write = async (relative, value) => {
    const file = path.join(root, relative), bytes = gzipSync(JSON.stringify(value));
    await fs.mkdir(path.dirname(file), {recursive: true});
    await fs.writeFile(file, bytes);
    return {path: relative, sha256: hash(bytes)};
  };
  const bundle = async values => {
    const column = hash(JSON.stringify(values));
    const payload = {encoding: SURVEY_BUNDLE_ENCODING, columns: {[column]: values}};
    const reference = hash(JSON.stringify(payload));
    await write(`survey/bundles/${reference}.json.gz`, payload);
    return {bundle: reference, column};
  };
  const scalar = encodeSurveyRows([{
    id: 'angle:residue', term_id: 'angle', pdb_id: 'TEST', entity_id: '1',
    residue_id: 'residue', endpoint_entities: [{pdb_id: 'TEST', entity_id: '1'}],
    is_terminal_any: false, value: 1.2345678901234567, status: 'ok',
  }], 'bundled-validation');
  scalar.encoding = BUNDLED_SURVEY_ENCODING;
  scalar.columns.value = await bundle(scalar.columns.value);
  scalar.columns.residue_id = await bundle(scalar.columns.residue_id);
  const manifest = {
    build_id: 'bundled-validation', partial: true,
    counts: {entries: 1, entities: 1, residues: 1}, source: {candidate_count: 1},
    metadata: await write('metadata.json.gz', {
      entries: [{pdb_id: 'TEST'}], entities: [{pdb_id: 'TEST', entity_id: '1'}],
    }),
    families: [{...await write('backbone.json.gz', [{id: 'residue', pdb_id: 'TEST',
      entity_id: '1', is_terminal_any: false, values: {}, statuses: {}}]),
    id: 'backbone', level: 'residue', parameters: [], row_count: 1}],
    relations: {}, survey: {scalars: {terms: {}}, coordinates: {groups: {}}},
    provenance: {decisions: await write('decisions.json.gz', [{pdb_id: 'TEST', accepted: true}])},
  };
  const save = async () => {
    manifest.survey.scalars.terms.angle = {...await write('survey/scalars/angle.json.gz', scalar),
      encoding: BUNDLED_SURVEY_ENCODING, row_count: 1};
    await fs.writeFile(path.join(root, 'manifest.json'), JSON.stringify(manifest));
  };
  await save();
  return {root, scalar, write, bundle, save, validate: () => validateRelease(path.join(root, 'manifest.json'))};
}

test('release validator accepts mixed inline and bundled Survey columns', async t => {
  const f = await fixture(t);
  assert.equal((await f.validate()).ok, true);
});

test('bundled release validation retains scientific value and status checks', async t => {
  const f = await fixture(t);
  f.scalar.columns.value = await f.bundle([null]);
  await f.save();
  const invalid = await f.validate();
  assert.equal(invalid.ok, false);
  assert.ok(invalid.errors.includes('Survey value/status: angle'));
  f.scalar.columns.status = ['missing_atoms'];
  await f.save();
  assert.equal((await f.validate()).ok, true);
});

test('release validator rejects corrupted Survey bundle content', async t => {
  const f = await fixture(t), reference = f.scalar.columns.value;
  await f.write(`survey/bundles/${reference.bundle}.json.gz`, {
    encoding: SURVEY_BUNDLE_ENCODING, columns: {[reference.column]: [99]},
  });
  await assert.rejects(f.validate(), /hash mismatch/i);
});

test('release validator rejects missing Survey bundles and allows retry', async t => {
  const f = await fixture(t), reference = f.scalar.columns.value;
  const file = path.join(f.root, `survey/bundles/${reference.bundle}.json.gz`);
  await fs.rename(file, `${file}.saved`);
  await assert.rejects(f.validate(), {code: 'ENOENT'});
  await fs.rename(`${file}.saved`, file);
  assert.equal((await f.validate()).ok, true);
});

test('release validator rejects invalid bundle references and missing columns', async t => {
  const f = await fixture(t), original = {...f.scalar.columns.value};
  f.scalar.columns.value.bundle = '../escape';
  await f.save();
  await assert.rejects(f.validate(), /invalid/i);
  f.scalar.columns.value = {...original, column: '0'.repeat(64)};
  await f.save();
  await assert.rejects(f.validate(), /column/i);
});

test('release validator rejects bundled and inline row length mismatches', async t => {
  const f = await fixture(t), original = {...f.scalar.columns.value};
  f.scalar.columns.value = await f.bundle([]);
  await f.save();
  await assert.rejects(f.validate(), /length mismatch/i);
  f.scalar.columns.value = original;
  f.scalar.columns.status = [];
  await f.save();
  await assert.rejects(f.validate(), /length mismatch/i);
});
