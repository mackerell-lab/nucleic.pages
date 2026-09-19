import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import os from 'node:os';
import path from 'node:path';
import {gzipSync} from 'node:zlib';
import {createHash} from 'node:crypto';
import {validateRelease} from '../offline/assets.mjs';
import {encodeSurveyRows} from '../core/survey-codec.js';
import {shareSurveyColumns, SHARED_SURVEY_ENCODING} from '../core/shared-survey-codec.js';

const hash = bytes => createHash('sha256').update(bytes).digest('hex');

async function fixture(t) {
  const root = await fs.mkdtemp(path.join(os.tmpdir(), 'rna-shared-validation-'));
  t.after(() => fs.rm(root, {recursive: true, force: true}));
  const write = async (relative, value) => {
    const file = path.join(root, relative), bytes = gzipSync(JSON.stringify(value));
    await fs.mkdir(path.dirname(file), {recursive: true});
    await fs.writeFile(file, bytes);
    return {path: relative, sha256: hash(bytes)};
  };
  const column = async values => {
    const reference = hash(JSON.stringify(values));
    await write(`survey/columns/${reference}.json.gz`, values);
    return reference;
  };
  const scalar = await shareSurveyColumns(encodeSurveyRows([{
    id: 'angle:residue', term_id: 'angle', pdb_id: 'TEST', entity_id: '1',
    residue_id: 'residue', endpoint_entities: [{pdb_id: 'TEST', entity_id: '1'}],
    is_terminal_any: false, value: 1.2345678901234567, status: 'ok',
  }], 'shared-validation'), column);
  const manifest = {
    build_id: 'shared-validation', partial: true,
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
      encoding: SHARED_SURVEY_ENCODING, row_count: 1};
    await fs.writeFile(path.join(root, 'manifest.json'), JSON.stringify(manifest));
  };
  await save();
  return {root, scalar, write, column, save, validate: () => validateRelease(path.join(root, 'manifest.json'))};
}

test('release validator expands shared columns and retains scientific checks', async t => {
  const f = await fixture(t);
  assert.equal((await f.validate()).ok, true);
  f.scalar.columns.value.reference = await f.column([null]);
  await f.save();
  const invalid = await f.validate();
  assert.equal(invalid.ok, false);
  assert.ok(invalid.errors.includes('Survey value/status: angle'));
});

test('release validator rejects corrupted shared column contents', async t => {
  const f = await fixture(t);
  await f.write(`survey/columns/${f.scalar.columns.value.reference}.json.gz`, [99]);
  await assert.rejects(f.validate(), /Shared Survey column hash mismatch/);
});

test('release validator rejects missing shared columns and can retry', async t => {
  const f = await fixture(t);
  const reference = f.scalar.columns.value.reference;
  const file = path.join(f.root, `survey/columns/${reference}.json.gz`);
  await fs.rename(file, `${file}.saved`);
  await assert.rejects(f.validate(), {code: 'ENOENT'});
  await fs.rename(`${file}.saved`, file);
  assert.equal((await f.validate()).ok, true);
});

test('release validator rejects shared column traversal and row mismatch', async t => {
  const f = await fixture(t);
  f.scalar.columns.value.reference = '../escape';
  await f.save();
  await assert.rejects(f.validate(), /Invalid shared Survey content hash/);
  f.scalar.columns.value.reference = await f.column([]);
  await f.save();
  await assert.rejects(f.validate(), /Shared Survey column length mismatch/);
});
