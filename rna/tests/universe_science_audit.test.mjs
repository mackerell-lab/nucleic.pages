import test from 'node:test';
import assert from 'node:assert/strict';
import { universeSummary } from '../views/universe-summary.js';

const rna = (pdb_id, entity_id, extra = {}) => ({ pdb_id, entity_id,
  type: 'polymer', polymer_type: 'polyribonucleotide', ...extra });
const axis = (result, key) => result.annotations.find(item => item.id === key);

test('Universe PDB and RNA entity denominators exclude nonpolymers without counting chain copies', () => {
  const entries = [{ pdb_id: '1AAA' }, { pdb_id: '2BBB' }, { pdb_id: '3CCC' }];
  const entities = [rna('1aaa', '1', { functions: ['ribozyme', 'ribozyme', 'aptamer'], chains: ['A', 'B'] }),
    rna('1AAA', '2', { functions: ['ribozyme'] }), rna('2BBB', '1', { functions: [] }),
    { pdb_id: '1AAA', entity_id: '3', type: 'non-polymer', functions: ['ligand'] },
    { pdb_id: '2BBB', entity_id: '2', type: 'water', functions: ['water'] },
    { pdb_id: '2BBB', entity_id: '3', type: 'polymer', polymer_type: 'polypeptide(L)', functions: ['protein'] },
    rna('OTHER', '1', { functions: ['outside'] })];
  const result = universeSummary(entries, entities, []), functions = axis(result, 'functions');
  assert.equal(result.entries, 3); assert.equal(result.entities, 3);
  assert.equal(result.entriesWithoutEntityMetadata, 1);
  assert.equal(result.annotatedEntries, 1); assert.equal(result.annotatedEntities, 2);
  assert.equal(functions.unknownEntities, 1); assert.equal(functions.entriesWithoutAnnotatedEntity, 2);
  assert.deepEqual(functions.tags.map(tag => [tag.id, tag.entities, tag.entries]),
    [['aptamer', 1, 1], ['ribozyme', 2, 1]]);
  assert.equal(functions.tags.reduce((sum, tag) => sum + tag.entities, 0), 3,
    'Overlapping tag memberships can exceed annotated entities');
});

test('Recorded entity inventory preserves explicit emptiness and does not relabel entry-scoped tags', () => {
  const result = universeSummary([{ pdb_id: '1AAA', annotation_scope: 'entry', functions: ['entry-function'] }], [
    rna('1AAA', '1', { functions: [], function_tags: ['hidden'], annotations: { functions: ['hidden'] } }),
    rna('1AAA', '2', { functions: null, function_tags: ['hidden'], structures: ['stem'] }),
    rna('1AAA', '3', { annotations: { functions: ['nested'] } }),
    rna('1AAA', '4', { function_tags: ['alias'] }),
  ], []);
  assert.deepEqual(axis(result, 'functions').tags.map(tag => tag.id), ['alias', 'nested']);
  assert.equal(axis(result, 'functions').unknownEntities, 2);
  assert.equal(result.annotatedEntities, 3);
  assert.equal(result.annotatedEntries, 1);
});

test('Component profile counts distinguish booleans, object verdicts and unknown without inferring HET mass', () => {
  const entries = [
    { pdb_id: '1AAA', profiles: { conservative: true, relaxed: { passed: true }, mw100: { eligible: true } } },
    { pdb_id: '2BBB', component_profiles: { conservative: false, relaxed: { passed: false }, mw100: { eligible: false } } },
    { pdb_id: '3CCC', profiles: { conservative: 'yes', relaxed: { passed: 'true' }, mw100: {} }, het_name: 'magnesium' },
    { pdb_id: '4DDD', profiles: ['relaxed', 'mw100'] },
    { pdb_id: '5EEE', eligibility: { profiles: { conservative: { passed: false, eligible: true }, relaxed: null } } },
  ];
  const result = universeSummary(entries, [], []);
  assert.deepEqual(result.profiles.map(p => [p.id, p.entries, p.unknown]),
    [['conservative', 2, 1], ['relaxed', 2, 2], ['mw100', 2, 2]]);
});

test('Stored family inventory does not infer residues, coverage, finite values or missing metadata', () => {
  const result = universeSummary([{ pdb_id: '1AAA', residue_count: 10 }], [rna('1AAA', '1')],
    [{ id: 'backbone', row_count: 10 }, { id: 'sugar', row_count: 10 }, { id: 'base_pair', row_count: 3 }]);
  assert.equal(result.familyRows, 23);
  assert.equal(result.entities, 1);
  assert(!Object.hasOwn(result, 'residues'));
  assert(!Object.hasOwn(result, 'coverage'));
  const unavailable = universeSummary([{ pdb_id: '1AAA' }], [], [{ id: 'backbone' }]);
  assert.equal(unavailable.familyRows, null);
  assert.equal(unavailable.families[0].rows, null);
  assert.equal(unavailable.entriesWithoutEntityMetadata, 1);
});
