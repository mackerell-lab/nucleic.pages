import test from 'node:test';
import assert from 'node:assert/strict';
import { join } from '../core/joints.js';

test('Hidden relation endpoint choices cannot filter same-observation joins', () => {
  const left = [{ id: 'u', values: { chi: -120 } }, { id: 'g', values: { chi: -90 } }];
  const right = [{ id: 'g', values: { delta: 82 } }, { id: 'u', values: { delta: 81 } }];
  for (const type of [undefined, 'identity', 'same_level', 'same_residue']) {
    const baseline = join(left, right, { type, endpoint: 'both' });
    for (const endpoint of ['first', 'second']) {
      const spec = { type, endpoint };
      const result = join(left, right, spec);
      assert.deepEqual(result.points, baseline.points, `${type ?? 'default'} rejected ${endpoint}`);
      assert.deepEqual(result.diagnostics, baseline.diagnostics);
      assert.equal(result.joinSpec.endpoint, 'both');
      assert.equal(spec.endpoint, endpoint, 'Stored relation-side preference was changed');
    }
  }
});

test('Relation endpoint filtering remains effective in either axis orientation', () => {
  const pairs = [{ id: 'p', values: { opening: 12 } }];
  const residues = [{ id: 'u', values: { chi: -120 } }, { id: 'g', values: { chi: -90 } }];
  const relations = [
    { id: 'r1', pair_id: 'p', residue_id: 'u', endpoint_role: 'first' },
    { id: 'r2', pair_id: 'p', residue_id: 'g', endpoint_role: 'second' },
  ];
  for (const reverse of [false, true]) {
    for (const [endpoint, expected] of [['first', ['u']], ['second', ['g']], ['both', ['u', 'g']]]) {
      const result = join(reverse ? residues : pairs, reverse ? pairs : residues, {
        type: 'relation', endpoint, relations,
        xParameter: { level: reverse ? 'residue' : 'pair' },
        yParameter: { level: reverse ? 'pair' : 'residue' },
      });
      assert.deepEqual(result.points.map(point => point.residue_id), expected);
      assert.equal(result.joinSpec.endpoint, endpoint);
      assert(result.points.every(point => point.pair_id === 'p'));
    }
  }
});
