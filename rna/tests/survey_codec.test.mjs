import test from 'node:test';
import assert from 'node:assert/strict';
import { decodeCoordinateRows, decodeFamilyRows, decodeInteractionRows, decodeSurveyRows, encodeCoordinateRows, encodeFamilyRows, encodeInteractionRows, encodeSurveyRows, COORDINATE_COLUMNAR_ENCODING, FAMILY_COLUMNAR_ENCODING, INTERACTION_COLUMNAR_ENCODING, SURVEY_COLUMNAR_ENCODING } from '../core/survey-codec.js';

test('Columnar Survey encoding round-trips raw identities, nulls and nested values', () => {
  const rows = [
    { id: 'r1', value: 12.3456789012345, status: 'ok', pair_id: null, residue_ids: ['r1'], endpoint_entities: [{ pdb_id: '1RNA', entity_id: '1' }] },
    { id: 'r2', value: null, status: 'missing_atoms', pair_id: 'p2', residue_ids: ['r2', 'r3'], endpoint_entities: [{ pdb_id: '2RNA', entity_id: '1' }] },
  ];
  const encoded = encodeSurveyRows(rows, 'full_test');
  assert.equal(encoded.encoding, SURVEY_COLUMNAR_ENCODING); assert.equal(encoded.row_count, rows.length);
  assert.deepEqual(decodeSurveyRows(encoded), rows);
  assert.deepEqual(decodeSurveyRows(encoded, ['id', 'value']), [{ id: 'r1', value: 12.3456789012345 }, { id: 'r2', value: null }]);
  assert.deepEqual(decodeSurveyRows(rows), rows);
});

test('Columnar Survey decoding rejects inconsistent columns and unknown encodings', () => {
  assert.throws(() => decodeSurveyRows({ encoding: SURVEY_COLUMNAR_ENCODING, row_count: 2, columns: { id: ['r1'] } }), /column length/);
  assert.throws(() => decodeSurveyRows({ encoding: 'other', columns: {} }), /Unsupported/);
});

test('Columnar coordinate encoding round-trips geometry rows and rejects bad columns', () => {
  const rows = [{ id: 'a', atom_label: 'A.N1', x: 1.25, y: null, residue_ids: ['a'] }, { id: 'b', atom_label: 'A.C2', x: 2.5, y: 3 }];
  const encoded = encodeCoordinateRows(rows, 'full_test');
  assert.equal(encoded.encoding, COORDINATE_COLUMNAR_ENCODING);
  assert.deepEqual(decodeCoordinateRows(JSON.parse(JSON.stringify(encoded))), rows);
  assert.deepEqual(decodeCoordinateRows(rows), rows);
  assert.throws(() => decodeCoordinateRows({ encoding: COORDINATE_COLUMNAR_ENCODING, row_count: 2, columns: { id: ['a'] } }), /column length/);
  assert.throws(() => decodeCoordinateRows({ ...encoded, row_count: -1 }), /row count/);
  assert.throws(() => decodeCoordinateRows({ ...encoded, missing: { id: [2] } }), /missing-field/);
});

test('Dictionary family encoding preserves nulls, absent fields and numeric values', () => {
  const rows = [{ id: 'r1', base: 'U', value: 1.25 }, { id: 'r2', base: 'U', value: null, optional: false }];
  const encoded = encodeFamilyRows(rows, 'full_test');
  assert.equal(encoded.encoding, FAMILY_COLUMNAR_ENCODING);
  assert.deepEqual(decodeFamilyRows(JSON.parse(JSON.stringify(encoded))), rows);
  assert.throws(() => decodeFamilyRows({ ...encoded, columns: { id: { dictionary: ['r1'], indices: [1, 0] } } }), /dictionary/);
});

test('Dictionary interaction encoding preserves directed graph identities and categories', () => {
  const rows = [{ id: 'e1', family: 'cWW', categories: ['basepair', 'basepair_detail'], near: false, residue1_id: 'r1', residue2_id: 'r2' },
    { id: 'e2', family: 'cWW', categories: ['basepair', 'basepair_detail'], near: true, residue1_id: 'r2', residue2_id: 'r1' }];
  const encoded = encodeInteractionRows(rows, 'full_test');
  assert.equal(encoded.encoding, INTERACTION_COLUMNAR_ENCODING);
  assert.deepEqual(decodeInteractionRows(JSON.parse(JSON.stringify(encoded))), rows);
});
