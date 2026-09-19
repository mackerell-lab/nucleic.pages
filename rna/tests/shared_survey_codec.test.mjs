import test from 'node:test';
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { encodeSurveyRows, decodeSurveyRows } from '../core/survey-codec.js';
import { shareSurveyColumns, expandSharedSurveyColumns } from '../core/shared-survey-codec.js';
import { RnaDataRepository } from '../core/repository.js';

test('Shared Survey columns preserve complete values and deduplicate identities', async () => {
  const assets = new Map();
  const write = async values => {
    const json = JSON.stringify(values), hash = createHash('sha256').update(json).digest('hex');
    assets.set(hash, JSON.parse(json)); return hash;
  };
  const rows = [{ id: 'a', residue_id: 'a', value: 1.2345678901234567, status: 'ok' },
    { id: 'b', residue_id: 'b', value: null, status: 'missing_atoms' }];
  const source = encodeSurveyRows(rows, 'test');
  const shared = await shareSurveyColumns(source, write);
  assert.equal(assets.size, 3);
  let reads = 0;
  const expanded = await expandSharedSurveyColumns(JSON.parse(JSON.stringify(shared)), async ref => { reads++; return assets.get(ref); });
  assert.equal(reads, 3);
  assert.deepEqual(decodeSurveyRows(expanded), rows);
  assert.deepEqual(decodeSurveyRows(expanded, ['id']), [{ id: 'a' }, { id: 'b' }]);
  assets.set(shared.columns.value.reference, [0]);
  await assert.rejects(expandSharedSurveyColumns(shared, async ref => assets.get(ref)), /length mismatch/);
});

test('Repository verifies shared content, validates projected-away values and retries', async () => {
  const assets = new Map();
  const source = encodeSurveyRows([{ id: 'a', value: 2 }, { id: 'b', value: null }], 'test');
  const shared = await shareSurveyColumns(source, async values => {
    const hash = createHash('sha256').update(JSON.stringify(values)).digest('hex');
    assets.set(`/survey/columns/${hash}.json.gz`, values); return hash;
  });
  assets.set('/manifest.json', { molecule_type: 'RNA', schema_version: 'rna-explorer-1', build_id: 'test',
    survey: { scalars: { terms: { angle: { path: 'angle.json', row_count: 2 } } } } });
  assets.set('/angle.json', shared);
  const calls = [];
  const repository = new RnaDataRepository({ manifestUrl: 'https://rna.test/manifest.json', fetchImpl: async url => {
    const key = new URL(url).pathname; calls.push(key);
    return new Response(JSON.stringify(assets.get(key)));
  } });
  const valuePath = `/survey/columns/${shared.columns.value.reference}.json.gz`;
  const original = assets.get(valuePath);
  assets.set(valuePath, [9, null]);
  await assert.rejects(repository.loadSurveyScalars('angle', { fields: ['id'] }), /hash mismatch/);
  assets.set(valuePath, original);
  const loaded = await repository.loadSurveyScalars('angle', { fields: ['id'] });
  assert.deepEqual(loaded.rows, [{ id: 'a' }, { id: 'b' }]);
  assert(!Object.hasOwn(loaded, 'columns'));
  assert.equal(calls.filter(key => key === valuePath).length, 2);
  repository.releaseSurvey('scalars', 'angle');
  shared.columns.id.reference = '../escape';
  await assert.rejects(repository.loadSurveyScalars('angle'), /content hash/);
  assert(!calls.some(key => key.includes('escape')));
});

test('Invalid references and counts fail closed; rejected loads can be retried', async () => {
  const source = encodeSurveyRows([{ id: 'a' }]);
  const shared = await shareSurveyColumns(source, async () => 'column');
  await assert.rejects(expandSharedSurveyColumns({ ...shared, row_count: -1 }, async () => []), /row count/);
  await assert.rejects(expandSharedSurveyColumns({ ...shared, columns: { id: { reference: '' } } }, async () => []), /reference/);
  await assert.rejects(expandSharedSurveyColumns(shared, async () => { throw new Error('HTTP failure'); }), /HTTP failure/);
  assert.deepEqual(decodeSurveyRows(await expandSharedSurveyColumns(shared, async () => ['a'])), [{ id: 'a' }]);
});

test('Shared columns overlap at most four transfers and await the full table', async () => {
  const source = encodeSurveyRows([Object.fromEntries(Array.from({ length: 9 }, (_, i) => [`field${i}`, i]))]);
  const shared = await shareSurveyColumns(source, async values => String(values[0]));
  let active = 0, maximum = 0, reads = 0;
  const expanded = await expandSharedSurveyColumns(shared, async ref => {
    active++; reads++; maximum = Math.max(maximum, active);
    await new Promise(resolve => setTimeout(resolve, 1));
    active--; return [Number(ref)];
  });
  assert.equal(maximum, 4);
  assert.equal(active, 0);
  assert.equal(reads, 9);
  assert.deepEqual(expanded, source);
});
