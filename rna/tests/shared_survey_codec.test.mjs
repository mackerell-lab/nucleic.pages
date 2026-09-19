import test from 'node:test';
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { encodeSurveyRows, decodeSurveyRows } from '../core/survey-codec.js';
import { shareSurveyColumns, expandSharedSurveyColumns } from '../core/shared-survey-codec.js';

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

test('Invalid references and counts fail closed; rejected loads can be retried', async () => {
  const source = encodeSurveyRows([{ id: 'a' }]);
  const shared = await shareSurveyColumns(source, async () => 'column');
  await assert.rejects(expandSharedSurveyColumns({ ...shared, row_count: -1 }, async () => []), /row count/);
  await assert.rejects(expandSharedSurveyColumns({ ...shared, columns: { id: { reference: '' } } }, async () => []), /reference/);
  await assert.rejects(expandSharedSurveyColumns(shared, async () => { throw new Error('HTTP failure'); }), /HTTP failure/);
  assert.deepEqual(decodeSurveyRows(await expandSharedSurveyColumns(shared, async () => ['a'])), [{ id: 'a' }]);
});
