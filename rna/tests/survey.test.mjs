import test from 'node:test';
import assert from 'node:assert/strict';
import { computeSurvey } from '../offline/survey.mjs';
import { TERM_REGISTRY, BASE_ATOMS, publicBaseGeometryConfig, openingBinForValue } from '../offline/survey_terms.mjs';

const residue = (id, base, atoms = {}) => ({ id, comp_id: base, atoms, label_asym_id: 'RNA' });
const entry = residues => ({ id: 'TEST', residues });
const termRow = (result, termId, residueId) => result.scalars.find(row => row.term_id === termId && (!residueId || row.residue_id === residueId));
const near = (actual, expected, epsilon = 1e-7) => assert.ok(Math.abs(actual - expected) < epsilon, `${actual} != ${expected}`);

test('RNA registry has 48 angles, 48 torsions and two exact pair distances', () => {
  assert.equal(TERM_REGISTRY.length, 98);
  assert.equal(new Set(TERM_REGISTRY.map(term => term.term_id)).size, 98);
  assert.equal(TERM_REGISTRY.filter(term => term.value_type === 'linear_angle').length, 48);
  assert.equal(TERM_REGISTRY.filter(term => term.value_type === 'dihedral').length, 48);
  assert.equal(TERM_REGISTRY.filter(term => term.value_type === 'distance').length, 2);
  assert.ok(!JSON.stringify(TERM_REGISTRY).includes('C7'));
  assert.ok(!BASE_ATOMS.U.includes('C7'));
  assert.ok(BASE_ATOMS.U.includes('O2'));
  assert.ok(!BASE_ATOMS.U.includes("O2'"));
  assert.deepEqual(TERM_REGISTRY.find(term => term.term_id === 'same_pair_g_c6__c_n4').source_atom_pattern, ['G.C6', 'C.N4']);
  assert.deepEqual(TERM_REGISTRY.find(term => term.term_id === 'same_pair_a_n6__u_o4').source_atom_pattern, ['A.N6', 'U.O4']);
  assert.equal(TERM_REGISTRY.filter(term => !term.requires_pair).length, 96);
  const config = publicBaseGeometryConfig();
  config.terms[0].source_atom_pattern[0] = 'BAD';
  assert.equal(TERM_REGISTRY[0].source_atom_pattern[0], 'N1');
});

test('unpaired residues retain all 96 applicable observations with explicit missing status', () => {
  const result = computeSurvey(entry(['A', 'C', 'G', 'U'].map(base => residue(base, base))));
  assert.equal(result.scalars.length, 96);
  assert.ok(result.scalars.every(row => row.status === 'missing_atoms' && row.value === null));
  assert.equal(result.diagnostics.length, 4);
  assert.equal(result.coordinates.length, 0);
});

test('analytic 90-degree angle and signed positive/negative 90-degree torsions', () => {
  const atoms = { N1: [1, 0, 0], C2: [0, 0, 0], N3: [0, 1, 0], C4: [0, 1, 1] };
  const positive = computeSurvey(entry([residue('u', 'U', atoms)]));
  near(termRow(positive, 'u_n1_c2_n3').value, 90);
  // Looking along C2->N3, N1->C4 has the standard signed torsion -90.
  near(termRow(positive, 'u_n1_c2_n3_c4').value, -90);
  const negative = computeSurvey(entry([residue('u', 'U', { ...atoms, C4: [0, 1, -1] })]));
  near(termRow(negative, 'u_n1_c2_n3_c4').value, 90);
});

test('O2 and O2-prime are never interchangeable; missing versus degenerate stays separate', () => {
  const result = computeSurvey(entry([residue('u', 'U', {
    N1: [1, 0, 0], C2: [0, 0, 0], N3: [0, 1, 0], "O2'": [0, 0, 1], C4: [0, 2, 0],
  })]));
  assert.equal(termRow(result, 'u_n1_c2_o2').status, 'missing_atoms');
  assert.equal(termRow(result, 'u_n1_c2_n3_c4').status, 'degenerate_geometry');
  assert.equal(termRow(result, 'u_n1_c2_n3').status, 'ok');
});

test('GC and reversed AU distances use specified atoms, retain isolated pairs and exclude near pairs', () => {
  const residues = [residue('g', 'G', { C6: [0, 0, 0], O6: [100, 0, 0] }),
    residue('c', 'C', { N4: [3, 4, 0] }), residue('a', 'A', { N6: [0, 0, 0] }),
    residue('u', 'U', { O4: [0, 0, 7], "O2'": [200, 0, 0] })];
  const pairs = [
    { id: 'gc', residue1_id: 'g', residue2_id: 'c', family: 'cWW', stem_eligible: false },
    { id: 'ua', residue1_id: 'u', residue2_id: 'a', family: 'cWW' },
    { id: 'near', residue1_id: 'g', residue2_id: 'c', family: 'ncWW' },
    { id: 'near-flag', residue1_id: 'g', residue2_id: 'c', family: 'cWW', near: true },
    { id: 'alternative-flag', residue1_id: 'u', residue2_id: 'a', family: 'cWW', alternative: true },
    { id: 'alternative-family', residue1_id: 'u', residue2_id: 'a', family: 'cWWa' },
  ];
  const result = computeSurvey(entry(residues), { pairs });
  near(termRow(result, 'same_pair_g_c6__c_n4').value, 5);
  near(termRow(result, 'same_pair_a_n6__u_o4').value, 7);
  assert.equal(result.scalars.filter(row => row.observation_level === 'pair').length, 2);
  assert.equal(result.scalars.filter(row => row.observation_level === 'residue').length, 96);
  assert.equal(result.interaction_links.length, 12);
});

test('standard U coordinates recover x3DNA reference under rigid rotation and translation', () => {
  // Independently transcribed from x3dna-v2.4/config/Atomic_U.pdb. Actual U,
  // including O4 out-of-plane coordinate; no thymine reference or C7 atom.
  const reference = { N1: [-1.284, 4.500, 0], C2: [-1.462, 3.131, 0],
    N3: [-0.302, 2.397, 0], C4: [0.989, 2.884, 0], C5: [1.089, 4.311, 0],
    C6: [-0.024, 5.053, 0], O4: [1.935, 2.094, -0.001] };
  const atoms = Object.fromEntries(Object.entries(reference).map(([atom, [x, y, z]]) => [atom, [-y + 10, x - 4, z + 3]]));
  const result = computeSurvey(entry([residue('u', 'U', atoms)]));
  assert.equal(result.coordinates.length, 7);
  for (const row of result.coordinates) {
    assert.equal(row.anchor_frame, 'rna_standard_base');
    assert.equal(row.frame_reference, 'x3dna_2.4_Atomic_U');
    [row.x, row.y, row.z].forEach((value, axis) => near(value, reference[row.atom_name][axis], 1e-6));
  }
});

test('invalid stable identities and pair mappings fail rather than merging observations', () => {
  assert.throws(() => computeSurvey(entry([residue('x', 'A'), residue('x', 'U')])), /Duplicate/);
  assert.throws(() => computeSurvey(entry([residue('', 'A')])), /stable residue/);
  assert.throws(() => computeSurvey(entry([residue('a', 'A')]), { pairs: [{ id: 'p', residue1_id: 'a', residue2_id: 'missing' }] }), /Invalid/);
});

test('opening bins have explicit inclusive boundaries and missing state', () => {
  assert.equal(openingBinForValue(-16), 'small');
  assert.equal(openingBinForValue(-8), 'middle');
  assert.equal(openingBinForValue(2), 'large');
  assert.equal(openingBinForValue(10), 'large');
  assert.equal(openingBinForValue(11), 'outside');
  assert.equal(openingBinForValue(null), 'missing');
});
