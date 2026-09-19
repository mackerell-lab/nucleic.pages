import test from 'node:test';
import assert from 'node:assert/strict';
import { chooseLinearRange } from '../core/linear-axis.js';

test('advisory ranges preserve stable axes for empty and narrow subsets', () => {
  const defaults = Object.freeze([-10, 10]);
  for (const values of [[], [NaN, Infinity, -Infinity, null, undefined, '2'], [2], [2, 2.01], [-10, 10]]) {
    const result = chooseLinearRange(values, { defaultRange: defaults });
    assert.deepEqual(result, [-10, 10]);
    assert.notEqual(result, defaults);
  }
  assert.deepEqual(chooseLinearRange([]), [0, 1]);
  assert.deepEqual(chooseLinearRange([NaN, Infinity, null]), [0, 1]);
});

test('advisory defaults expand to include outliers with DNA padding', () => {
  assert.deepEqual(chooseLinearRange([-20, 4], { defaultRange: [-10, 10] }), [-21.2, 10]);
  assert.deepEqual(chooseLinearRange([-4, 20], { defaultRange: [-10, 10] }), [-10, 21.2]);
  assert.deepEqual(chooseLinearRange([-100, 100], { defaultRange: [-10, 10] }), [-110, 110]);
  assert.deepEqual(chooseLinearRange([100], { defaultRange: [-10, 10] }), [-10, 110]);
  assert.deepEqual(chooseLinearRange([-100], { defaultRange: [-10, 10] }), [-110, 10]);
});

test('default-free linear padding matches DNA for singleton and narrow values', () => {
  assert.deepEqual(chooseLinearRange([0]), [-0.5, 0.5]);
  assert.deepEqual(chooseLinearRange([20]), [18, 22]);
  assert.deepEqual(chooseLinearRange([-20]), [-22, -18]);
  assert.deepEqual(chooseLinearRange([0, 10]), [-0.5, 10.5]);
  assert.deepEqual(chooseLinearRange([0, 100]), [-5, 105]);
  assert.deepEqual(chooseLinearRange([3, 3.2]), [2.5, 3.7]);
  assert.deepEqual(chooseLinearRange([20, 20 + 1e-10]), [20 - (20 + 1e-10) * 0.1, (20 + 1e-10) * 1.1]);
});

test('explicit requested range remains a validated hard clipping boundary', () => {
  const requested = Object.freeze([0, 1]);
  const range = chooseLinearRange([-100, 100], { requested, defaultRange: [-200, 200] });
  assert.deepEqual(range, [0, 1]);
  assert.notEqual(range, requested);
  assert.deepEqual(chooseLinearRange([], { requested: new Float64Array([0, 1]) }), [0, 1]);
  assert.deepEqual(chooseLinearRange([10], { requested: null, defaultRange: [0, 20] }), [0, 20]);
  for (const invalid of [false, 0, '', '01', {}, [], [1], [0, 1, 2], [0, 0], [1, 0], [NaN, 1], [0, Infinity], [, 1]]) {
    assert.throws(() => chooseLinearRange([2], { requested: invalid }), /Invalid plot range/);
    assert.throws(() => chooseLinearRange([2], { defaultRange: invalid }), /Invalid default plot range/);
  }
});

test('linear ranges do not infer circularity or discard finite observations', () => {
  const values = [-5, 0, 15, 90, 180, 220];
  assert.deepEqual(chooseLinearRange(values, { defaultRange: [0, 180] }), [-16.25, 231.25]);
  for (const data of [values, [-1000, -2, 0, 4, 1200], [1e-12], [-1e10, 1e10], [-Number.MAX_VALUE, Number.MAX_VALUE]]) {
    const before = [...data];
    const range = chooseLinearRange(data, { defaultRange: [0, 180] });
    assert.ok(range.every(Number.isFinite));
    for (const value of data) assert.ok(value >= range[0] && value <= range[1]);
    assert.deepEqual(data, before);
  }
  assert.deepEqual(chooseLinearRange(new Float64Array([NaN, -20, 20])), [-22, 22]);
});
