import test from 'node:test';
import assert from 'node:assert/strict';
import { chooseCircularCut } from '../core/circular-axis.js';
import { distribution, histogram2D } from '../core/analysis.js';

const repeat = (value, count) => Array(count).fill(value);
const valuesFromCounts = (counts, period = 360) => counts.flatMap((count, index) => repeat((index + 0.5) * period / counts.length, count));

test('Auto keeps the canonical seam for empty and fewer than eight observations', () => {
  for (const values of [[], [10], [10, 11, 12, 13, 14, 15, 16]]) assert.equal(chooseCircularCut(values, 360, 72), 0);
  assert.equal(chooseCircularCut(repeat(10, 8), 360, 72), 195);
});

test('longest empty gap crosses zero and deterministic ties use the first gap', () => {
  // Occupied bins 3..5: the nine-bin gap starts at 6, midpoint rounds to 11.
  assert.equal(chooseCircularCut([95, 125, 155].flatMap(value => repeat(value, 3)), 360, 12), 330);
  // Two equal five-bin gaps start at 1 and 7. First midpoint rounds to bin 4.
  const bimodal = [...repeat(15, 8), ...repeat(195, 8)];
  assert.equal(chooseCircularCut(bimodal, 360, 12), 120);
  assert.equal(chooseCircularCut([...bimodal].reverse(), 360, 12), 120);
});

test('uniform and nearly uniform data retain zero instead of a noise-selected seam', () => {
  for (const bins of [1, 2, 5, 12, 72, 144]) {
    const values = valuesFromCounts(Array(bins).fill(20));
    assert.equal(chooseCircularCut(values, 360, bins), 0);
    assert.equal(chooseCircularCut([...values, 181], 360, bins), 0);
  }
});

test('fully occupied low-density valley chooses its center after smoothing', () => {
  // Seven-bin window plus nine-bin Gaussian support exactly spans bins 20..34
  // at center 27. Every other center includes at least one high-density bin.
  const counts = Array.from({ length: 72 }, (_, index) => index >= 20 && index <= 34 ? 1 : 20);
  assert.equal(chooseCircularCut(valuesFromCounts(counts), 360, 72), 135);
  // Repeat the same valley half a turn later: the earlier equal minimum wins.
  for (let index = 56; index <= 70; index++) counts[index] = 1;
  assert.equal(chooseCircularCut(valuesFromCounts(counts), 360, 72), 135);
});

test('arbitrary periods scale seams while explicit signed and wrap modes stay fixed', () => {
  for (const period of [1, 180, 360, 2 * Math.PI, 720]) {
    const counts = Array.from({ length: 72 }, (_, index) => index >= 20 && index <= 34 ? 1 : 20);
    assert.equal(chooseCircularCut(valuesFromCounts(counts, period), period, 72), 27 * period / 72);
    for (const mode of ['signed', 'signed_180']) assert.equal(chooseCircularCut([period / 4], period, 72, mode), -period / 2);
    for (const mode of ['wrap', 'wrap_360']) assert.equal(chooseCircularCut([period / 4], period, 72, mode), 0);
    assert.equal(chooseCircularCut([0, period, -period, 2 * period], period, 72), 0);
  }
});

test('invalid periods, bin counts and nonfinite Auto values are rejected', () => {
  for (const period of [0, -1, NaN, Infinity]) assert.throws(() => chooseCircularCut([], period, 72), /period/);
  for (const bins of [0, 1.5, 2049, NaN, Infinity]) assert.throws(() => chooseCircularCut([], 360, bins), /Bin count/);
  for (const value of [NaN, Infinity, null, '1']) assert.throws(() => chooseCircularCut([value], 360, 72), /finite/);
});

test('grouped circular curves share one seam and retain physical moments across modes', () => {
  const rows = Array.from({ length: 8 }, (_, index) => ({ id: `r${index}`, pdb_id: 'TEST', comp_id: index % 2 ? 'A' : 'U', angle: index % 2 ? 355 : 5 }));
  const parameter = { id: 'angle', unit: 'degree', period: 360 };
  const modes = ['auto', 'signed_180', 'wrap_360'].map(circularMode => distribution(rows, parameter, { circularMode, groupBy: 'base', bins: 72, sigma: 0 }));
  // 5 and 355 occupy bins 1 and 71: gap 2..70 rounds its midpoint to 37.
  assert.equal(modes[0].displayCut, 185);
  for (const result of modes) {
    assert.equal(result.coverage.finiteRows, 8);
    assert.equal(result.coverage.plottedRows, 8);
    for (const series of result.series) {
      assert.equal(series.displayCut, result.displayCut);
      assert.equal(series.counts.reduce((a, b) => a + b, 0), 4);
      const baseline = modes[0].series.find(item => item.key === series.key);
      assert.deepEqual(series.values, baseline.values);
      assert.equal(series.statistics.mean, series.key === 'A' ? 355 : 5);
      assert.equal(series.statistics.resultant, 1);
      assert.equal(Math.abs(series.statistics.std), 0);
      assert.equal(series.statistics.std, baseline.statistics.std);
    }
  }
});

test('joint circular axes retain every point with independent declared periods', () => {
  const points = Array.from({ length: 8 }, (_, index) => ({ id: `p${index}`, x: index % 2 ? 355 : 5, y: index % 2 ? 177.5 : 2.5 }));
  const xParameter = { id: 'x', period: 360 }, yParameter = { id: 'y', period: 180 };
  const auto = histogram2D(points, xParameter, yParameter, { bins: 72, sigma: 0 });
  assert.equal(auto.x[0], 187.5);
  assert.equal(auto.y[0], 93.75);
  for (const circularMode of ['auto', 'signed_180', 'wrap_360']) {
    const result = histogram2D(points, xParameter, yParameter, { circularMode, bins: 72, sigma: 0 });
    assert.deepEqual(result.points.map(point => [point.x, point.y]), points.map(point => [point.x, point.y]));
    assert.equal(result.z.flat().reduce((a, b) => a + b, 0), 1);
  }
});
