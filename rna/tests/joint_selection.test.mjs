import test from 'node:test';
import assert from 'node:assert/strict';
import { jointSelectionSpecs } from '../core/joint-selection.js';
import { selectRows } from '../core/selection.js';
import { join } from '../core/joints.js';

test('Pair context and independent RNA endpoint pucker preserve multiplexed identities', () => {
  const selection = { contexts: ['AU'], puckerStates: [], includeEnds: true };
  const joint = { mode: 'relation', residueContexts: ['U'], residuePuckers: ["C3'-endo"] };
  const specs = jointSelectionSpecs(selection, joint, 'pair', 'residue');
  const pairs = [{ id: 'p1', level: 'pair', context: 'AU' }, { id: 'p2', level: 'pair', context: 'AU' }, { id: 'p3', level: 'pair', context: 'GC' }];
  const residues = [{ id: 'a', comp_id: 'A', pucker_class: "C3'-endo" }, { id: 'u', comp_id: 'U', pucker_class: "C3'-endo" }, { id: 'u2', comp_id: 'U', pucker_class: "C2'-endo" }];
  const relations = [['p1', 'a'], ['p1', 'u'], ['p2', 'u'], ['p2', 'u2'], ['p3', 'u']].map(([pair_id, residue_id]) => ({ kind: 'pair_residue', pair_id, residue_id, endpoint_role: residue_id === 'a' ? 'first' : 'second' }));
  const selectedPairs = selectRows(pairs, {}, specs.left).rows, selectedResidues = selectRows(residues, {}, specs.right).rows;
  const result = join(selectedPairs, selectedResidues, { type: 'relation', relations, leftKey: 'pair_id', rightKey: 'residue_id' });
  assert.deepEqual(result.points.map(point => [point.pair_id, point.residue_id]), [['p1', 'u'], ['p2', 'u']]);
  assert.deepEqual(selection.contexts, ['AU']);
  assert.deepEqual(specs.right.contexts, ['U']);
  assert.deepEqual(jointSelectionSpecs(selection, joint, 'residue', 'pair').left.contexts, ['U']);
});

test('Joint selection rejects incompatible levels and keeps same-observation filters', () => {
  const selection = { contexts: ['U'], puckerStates: ["C3'-endo"] };
  for (const levels of [['step', 'residue'], ['pair', 'pair'], ['residue', 'residue']]) {
    assert.equal(jointSelectionSpecs(selection, { mode: 'relation' }, ...levels).valid, false);
  }
  assert.equal(jointSelectionSpecs(selection, { mode: 'identity' }, 'pair', 'residue').valid, false);
  const identity = jointSelectionSpecs(selection, { mode: 'identity' }, 'residue', 'residue');
  assert.equal(identity.left, selection); assert.equal(identity.right, selection);
  const reverse = jointSelectionSpecs(selection, { mode: 'relation', residueContexts: ['A'] }, 'residue', 'pair');
  assert.deepEqual(reverse.left.contexts, ['A']); assert.deepEqual(reverse.right.contexts, []);
});
