import test from 'node:test';
import assert from 'node:assert/strict';
import { summary, correlation } from '../math/numeric.js';

test('Circular summaries do not allocate percentile arrays that are discarded', () => {
  const angles = [359, 1, NaN];
  // The circular result needs only its weighted moments, never a percentile copy.
  angles.filter = () => { throw new Error('Circular summary allocated a percentile array'); };
  Object.freeze(angles);
  for (const weights of [null, [2, 2, 1]]) {
    const result = summary(angles, { period: 360, weights });
    assert.equal(result.n, 2);
    assert.equal(result.totalWeight, weights ? 4 : 2);
    assert(Math.abs(result.mean) < 1e-10 || Math.abs(result.mean - 360) < 1e-10);
    assert(Math.abs(result.resultant - Math.cos(Math.PI / 180)) < 1e-14);
    assert.equal(result.p05, null); assert.equal(result.p95, null);
    assert.equal(result.quantilePolicy, 'not_defined_for_circle');
  }
});

test('Nonperiodic degree-valued angles retain type-seven percentiles', () => {
  const result = summary([150, 30, NaN, 90, Infinity]);
  assert.equal(result.mean, 90);
  assert.equal(result.p05, 36); assert.equal(result.p95, 144);
  assert.equal(result.quantilePolicy, 'raw_linear_interpolation_type_7');
  const weighted = summary([30, 90, 150], { weights: [1, 2, 1] });
  assert.equal(weighted.mean, 90);
  assert.equal(weighted.p05, null); assert.equal(weighted.p95, null);
  assert.equal(weighted.quantilePolicy, 'not_computed_for_weighted_series');
});

test('Circular empty, undefined direction and arbitrary periods retain semantics', () => {
  assert.equal(summary([], { period: 360 }).meanStatus, 'empty');
  const flat = summary([0, 90, 180, 270], { period: 360 });
  assert.equal(flat.mean, null); assert.equal(flat.std, null);
  assert.equal(flat.meanStatus, 'undefined_mean_direction');
  const half = summary([179, 1], { period: 180 });
  assert(Math.abs(half.mean) < 1e-10 || Math.abs(half.mean - 180) < 1e-10);
  assert.equal(half.quantilePolicy, 'not_defined_for_circle');
  const correlated = correlation([355, 0, 5, 10], [5, 10, 15, 20], 360, 360);
  assert(Math.abs(correlated.r - 1) < 1e-12);
  assert.equal(correlated.r2, null);
});
