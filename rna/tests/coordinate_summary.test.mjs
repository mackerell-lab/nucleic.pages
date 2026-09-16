import test from 'node:test';
import assert from 'node:assert/strict';
import { CoordinateSummary } from '../core/coordinates.js';

test('Coordinate populations preserve many-pair residues and model/entry identity', () => {
  const summary = new CoordinateSummary();
  const row = { pdb_id: '1RNA', model_id: '1', pair_id: 'p1', target_residue_id: 'r1', context: 'G-C', atom_label: 'paired_G.N1', x: 0, y: 0, z: 0 };
  summary.add(row);
  summary.add({ ...row, pair_id: 'p2', x: 2 });
  summary.add({ ...row, pdb_id: '2RNA', x: 4 });
  summary.add({ ...row, model_id: '30', x: 6 });
  const result = summary.results()[0];
  assert.equal(result.n, 4); assert.equal(result.residues, 3); assert.equal(result.pairs, 4); assert.equal(result.entries, 2);
  assert.deepEqual(result.mean, [3, 0, 0]); assert.equal(result.rms, Math.sqrt(5));
});

test('Unpaired frames, context separation, missing coordinates and stable RMS', () => {
  const summary = new CoordinateSummary();
  const row = { pdb_id: '1RNA', residue_id: 'U1', context: 'U', atom_label: 'anchor_U.O2', x: 1e12, y: 0, z: 0 };
  summary.add(row); summary.add({ ...row, residue_id: 'U2', x: 1e12 + 2 });
  summary.add({ ...row, atom_label: "anchor_U.O2'", x: 5 });
  summary.add({ ...row, context: 'A', atom_label: 'anchor_A.N1' });
  summary.add({ ...row, xyz: [] }); summary.add({ ...row, z: null });
  summary.add({ ...row, status: 'missing_atoms' });
  const result = summary.results().find(item => item.atom_label === 'anchor_U.O2');
  assert.equal(summary.results().length, 3); assert.equal(result.n, 2); assert.equal(result.residues, 2);
  assert.equal(result.pairs, null); assert.equal(result.mean[0], 1e12 + 1); assert.equal(result.rms, 1);
  assert.deepEqual(new CoordinateSummary().results(), []);
});
