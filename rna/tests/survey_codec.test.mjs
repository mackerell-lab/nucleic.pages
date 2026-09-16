import test from 'node:test';
import assert from 'node:assert/strict';
import { decodeSurveyRows, encodeSurveyRows, SURVEY_COLUMNAR_ENCODING } from '../core/survey-codec.js';

test('Columnar Survey encoding round-trips raw identities, nulls and nested values', () => {
  const rows = [
    { id: 'r1', value: 12.3456789012345, status: 'ok', pair_id: null, residue_ids: ['r1'], endpoint_entities: [{ pdb_id: '1RNA', entity_id: '1' }] },
    { id: 'r2', value: null, status: 'missing_atoms', pair_id: 'p2', residue_ids: ['r2', 'r3'], endpoint_entities: [{ pdb_id: '2RNA', entity_id: '1' }] },
  ];
  const encoded = encodeSurveyRows(rows, 'full_test');
  assert.equal(encoded.encoding, SURVEY_COLUMNAR_ENCODING); assert.equal(encoded.row_count, rows.length);
  assert.deepEqual(decodeSurveyRows(encoded), rows);
  assert.deepEqual(decodeSurveyRows(rows), rows);
});

test('Columnar Survey decoding rejects inconsistent columns and unknown encodings', () => {
  assert.throws(() => decodeSurveyRows({ encoding: SURVEY_COLUMNAR_ENCODING, row_count: 2, columns: { id: ['r1'] } }), /column length/);
  assert.throws(() => decodeSurveyRows({ encoding: 'other', columns: {} }), /Unsupported/);
});
