import test from 'node:test';
import assert from 'node:assert/strict';
import { enrichEntryMetadata } from '../offline/metadata_enrichment.mjs';

function fixture() {
  const source_sha256 = 'a'.repeat(64);
  return { entry: { source_sha256, entities: [
    { entity_id: '1', type: 'polymer', polymer_type: 'polyribonucleotide', sequence: [
      { seq_id: '1', mon_id: 'U' }, { seq_id: '2', mon_id: 'G' },
      { seq_id: '2', mon_id: 'A', hetero: 'y' },
    ] },
  ], residues: [{ id: 'u', entity_id: '1', label_asym_id: 'A', label_seq_id: '1',
    comp_id: 'U', model_id: '1', atoms: { O2: [1, 2, 3], "O2'": [2, 3, 4] } }],
  components: [], chemical_components: [], declared_components: { nonpolymer: [] }, explicit_connections: [],
  coordinate_policy: { id: 'single_deposited_model_v1', scope: 'deposited_asymmetric_unit', model_id: '1', model_count: 10 } },
  registry: { source_sha256, struct_asym: [{ id: 'A', entity_id: '1' }, { id: 'B', entity_id: '1' }] } };
}

test('Coverage counts positions, including fully unmodeled chain copies, without multiplying NMR models', () => {
  const f = fixture();
  f.entry.residues.push({ ...f.entry.residues[0], id: 'alternative', altloc: 'B' });
  const before = structuredClone(f), result = enrichEntryMetadata(f.entry, f.registry);
  assert.equal(result.coverage.status, 'available');
  assert.equal(result.coverage.observed_position_count, 1);
  assert.equal(result.coverage.declared_position_count, 4);
  assert.equal(result.coverage.fraction, 0.25);
  assert.equal(result.coverage.selected_model_id, '1');
  assert.deepEqual(f, before);
});

test('Empty and nonfinite retained-atom groups yield known zero coverage, not full residue coverage', () => {
  const f = fixture();
  f.entry.residues[0].atoms = {};
  f.entry.residues.push({ ...f.entry.residues[0], id: 'bad', label_seq_id: '2', comp_id: 'G',
    atoms: { P: [NaN, 1, 2], O2: [1, null, 3], "O2'": [] } });
  const result = enrichEntryMetadata(f.entry, f.registry);
  assert.equal(result.coverage.status, 'available');
  assert.equal(result.coverage.observed_position_count, 0);
  assert.equal(result.coverage.declared_position_count, 4);
  assert.equal(result.coverage.fraction, 0);
  assert.equal(f.entry.residues.length, 2, 'Legacy residue-group inventory is not overwritten');
});

test('Ambiguous model, chain, entity and declared identities never produce a plausible coverage ratio', () => {
  for (const mutate of [
    f => { f.entry.residues[0].model_id = '2'; },
    f => { f.entry.coordinate_policy.model_id = null; },
    f => { f.entry.residues[0].entity_id = 'missing'; },
    f => { f.entry.residues[0].label_asym_id = 'missing'; },
    f => { f.entry.residues[0].label_seq_id = '9'; },
    f => { f.entry.residues[0].comp_id = 'C'; },
    f => { f.registry.struct_asym.push({ id: 'A', entity_id: '1' }); },
    f => { f.entry.entities.push(structuredClone(f.entry.entities[0])); },
  ]) {
    const f = fixture(); mutate(f);
    const result = enrichEntryMetadata(f.entry, f.registry);
    assert.equal(result.coverage.status, 'unknown');
    assert.equal(result.coverage.observed_position_count, null);
    assert.equal(result.coverage.declared_position_count, null);
    assert.equal(result.coverage.fraction, null);
    assert(result.coverage.reasons.length);
  }
});

test('Pure helper checks digest consistency but cannot authenticate a self-reported registry', () => {
  const f = fixture();
  assert.throws(() => enrichEntryMetadata(f.entry, { ...f.registry, source_sha256: 'b'.repeat(64) }), /source identity/);
  // Deliberately demonstrate the caller trust boundary. No source bytes are
  // supplied to this pure helper, so a copied SHA cannot prove chain provenance.
  f.registry.struct_asym.push({ id: 'C', entity_id: '1' });
  assert.equal(enrichEntryMetadata(f.entry, f.registry).coverage.declared_position_count, 6);
});

test('Associated water and bromide facts do not become ion taxonomy or selected physical-bond claims', () => {
  const f = fixture();
  for (const [entity_id, type, chain] of [['2', 'non-polymer', 'C'], ['3', 'water', 'D'], ['4', 'non-polymer', 'E']]) {
    f.entry.entities.push({ entity_id, type }); f.registry.struct_asym.push({ id: chain, entity_id });
  }
  f.entry.chemical_components = [{ id: 'BR', name: 'BROMIDE ION', type: 'non-polymer' }];
  f.entry.declared_components.nonpolymer = [{ entity_id: '2', comp_id: 'BR' },
    { entity_id: '3', comp_id: 'HOH', name: 'water' }, { entity_id: '4', comp_id: 'X' }];
  const component = (id, comp_id, entity_id, label_asym_id, atoms, is_water = false) => ({
    id, comp_id, entity_id, label_asym_id, model_id: '1', atoms, is_water });
  f.entry.components = [component('br', 'BR', '2', 'C', { BR: [1, 2, 3] }),
    component('water1', 'HOH', '3', 'D', { O: [0, 0, 0], H1: [1, 0, 0], H2: [-1, 0, 0] }, true),
    component('water2', 'HOH', '3', 'D', { O: [0, 0, 1] }, true),
    component('water-empty', 'HOH', '3', 'D', {}, true), component('empty-x', 'X', '4', 'E', {})];
  const connection = (conn_type_id, firstChain, firstComp, secondChain, secondComp) => ({ conn_type_id,
    ptnr1_label_asym_id: firstChain, ptnr1_label_comp_id: firstComp,
    ptnr2_label_asym_id: secondChain, ptnr2_label_comp_id: secondComp });
  f.entry.explicit_connections = [connection('covale', 'A', 'U', 'C', 'BR'),
    connection('covale', 'C', 'BR', 'E', 'X'), connection('metalc', 'A', 'U', 'C', 'BR'),
    connection('covale', 'A', 'U', 'UNDECLARED', 'BR'),
    connection('hydrog', 'A', 'U', 'D', 'HOH')];
  const result = enrichEntryMetadata(f.entry, f.registry).associated_components;
  const br = result.observed.find(row => row.comp_id === 'BR'), water = result.observed.find(row => row.comp_id === 'HOH');
  assert.equal(br.observed_instance_count, 1); assert.equal(br.observed_atom_count, 1);
  assert.equal(br.declared_covalent_connection_count, 2); assert.equal(br.declared_rna_covalent_connection_count, 1);
  assert.equal(br.name, 'BROMIDE ION'); assert.equal(br.ion, undefined); assert.equal(br.mass, undefined);
  assert.equal(water.observed_instance_count, 2); assert.equal(water.water_instance_count, 2);
  assert.equal(water.observed_atom_count, 4); assert.equal(water.declared_covalent_connection_count, 0);
  assert.deepEqual(result.declared_only, [{ comp_id: 'X', name: null }]);
  assert.equal(result.connection_scope, 'deposited_annotations_not_verified_selected_atom_bonds');
});

test('Nonpolymer observations from wrong models or conflicting chain ownership fail closed', () => {
  for (const mutate of [
    row => { row.model_id = '2'; },
    row => { row.label_asym_id = 'unknown'; },
    row => { row.entity_id = '1'; },
  ]) {
    const f = fixture();
    f.entry.entities.push({ entity_id: '2', type: 'non-polymer' });
    f.registry.struct_asym.push({ id: 'C', entity_id: '2' });
    const row = { id: 'component', comp_id: 'BR', entity_id: '2', label_asym_id: 'C',
      model_id: '1', atoms: { BR: [0, 0, 0] } };
    mutate(row); f.entry.components.push(row);
    assert.throws(() => enrichEntryMetadata(f.entry, f.registry), /component identity/);
  }
});
