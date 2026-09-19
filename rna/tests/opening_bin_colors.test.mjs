import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import vm from 'node:vm';
import { distributionTraces } from '../views/panels.js';
import { distribution } from '../core/analysis.js';
import { createPlotSnapshot, csv } from '../core/export.js';

const dnaSource = readFileSync(new URL('../../js/pure-dna.js', import.meta.url), 'utf8');
const dnaBins = dnaSource.match(/const BASE_GEOMETRY_OPENING_BINS = ([^;]+);/)[1];
const dnaMeta = dnaSource.match(/const BASE_GEOMETRY_BIN_META = (\{[\s\S]*?\n\});/)[1];
const oracle = vm.runInNewContext(`({ bins: ${dnaBins}, meta: ${dnaMeta} })`);
const rowsFor = keys => keys.map((opening_bin, index) => ({ id: `r${index}`, pdb_id: 'TEST', opening_bin, angle: index * 10 }));
const resultFor = keys => distribution(rowsFor(keys), { id: 'angle', unit: 'deg' }, { groupBy: 'opening_bin', sigma: 0 });

test('opening-bin colors and order match DNA across missing bins and reordering', () => {
  for (const keys of [['large'], ['middle'], ['small'], ['large', 'small'], ['large', 'middle'], ['middle', 'small'], ['large', 'middle', 'small'], ['middle', 'large', 'small'], ['small', 'middle', 'large'], []]) {
    const result = resultFor(keys), before = structuredClone(result);
    const expectedKeys = Array.from(oracle.bins).filter(key => keys.includes(key));
    const traces = distributionTraces(result);
    assert.deepEqual(traces.map(trace => trace.name.split(' ')[0]), expectedKeys);
    assert.deepEqual(traces.map(trace => trace.line.color), expectedKeys.map(key => oracle.meta[key].color));
    assert.deepEqual(traces.map(trace => trace.fillcolor), expectedKeys.map(key => `${oracle.meta[key].color}18`));
    assert.deepEqual(result, before);
  }
});

test('actual filtered populations keep stable colors and exact scientific exports', () => {
  const rows = rowsFor(['large', 'small', 'middle', 'large', 'middle']);
  for (const selected of [rows, rows.filter(row => row.opening_bin !== 'small'), rows.filter(row => row.opening_bin === 'large')]) {
    const result = distribution(selected, { id: 'angle' }, { groupBy: 'opening_bin' });
    const snapshot = createPlotSnapshot({ result, buildId: 'test', parameter: result.parameter });
    const before = csv(snapshot), original = structuredClone(result);
    const traces = distributionTraces(snapshot.result, { traceStyle: 'line', groupBy: 'none' });
    for (const trace of traces) {
      const key = trace.name.split(' ')[0];
      assert.equal(trace.line.color, oracle.meta[key].color);
      assert.equal(trace.fill, 'none');
      const series = result.series.find(item => item.key === key);
      assert.deepEqual(trace.x, series.x);
      assert.deepEqual(trace.y, series.y);
    }
    assert.deepEqual(result, original);
    assert.equal(csv(snapshot), before);
  }
});

test('ordinary groups retain original index palette even with opening-like names', () => {
  const result = resultFor(['large', 'small', 'middle']);
  result.displaySpec.groupBy = 'base';
  const traces = distributionTraces(result);
  assert.deepEqual(traces.map(trace => trace.name.split(' ')[0]), ['large', 'small', 'middle']);
  assert.deepEqual(traces.map(trace => trace.line.color), ['#174a7e', '#8c3b2a', '#146c43']);
});

test('unknown opening groups remain visible after canonical bins', () => {
  const result = resultFor(['outside', 'large', 'missing', 'small']);
  const traces = distributionTraces(result);
  assert.deepEqual(traces.map(trace => trace.name.split(' ')[0]), ['small', 'large', 'outside', 'missing']);
  assert.equal(traces[0].line.color, oracle.meta.small.color);
  assert.equal(traces[1].line.color, oracle.meta.large.color);
  assert.equal(traces.length, result.series.length);
});
