import test from 'node:test';
import assert from 'node:assert/strict';
import { distribution } from '../core/analysis.js';
import { distributionTraces } from '../views/panels.js';
import { createPlotSnapshot, csv } from '../core/export.js';
import { seriesColor } from '../views/series-colors.js';

const parameter = { id: 'angle', period: 360 };
const rows = ['A', 'C', 'G', 'U'].map((comp_id, index) => ({ id: `r${index}`, pdb_id: 'TEST', comp_id, angle: index * 35,
  method: ['X-RAY DIFFRACTION', 'SOLUTION NMR', 'ELECTRON MICROSCOPY', 'OTHER'][index],
  functions: [`tag-${index}`, 'shared'], structures: [`structure-${index}`] }));
function resultFor(selected, groupBy) { return distribution(selected, parameter, { groupBy, sigma: 0 }); }
function traceColors(result) {
  const traces = distributionTraces(result);
  return new Map(result.series.map((series, index) => [series.key, traces[index].line.color]));
}

test('base and method color identities stay fixed under true subsets and reordering', () => {
  for (const [groupBy, expected] of [
    ['base', { A: '#174a7e', C: '#8c3b2a', G: '#146c43', U: '#8659a1' }],
    ['method', { xray: '#174a7e', nmr: '#8c3b2a', em: '#146c43', other: '#be882e' }],
  ]) for (const selected of [rows, [...rows].reverse(), rows.slice(1), [rows[3]], [rows[2], rows[0]]]) {
    const result = resultFor(selected, groupBy);
    for (const [key, color] of traceColors(result)) assert.equal(color, expected[key]);
  }
});

test('arbitrary annotation and ordered context keys retain colors without sorting science', () => {
  const contextRows = rows.map((row, index) => ({ ...row, comp_id: undefined, context: ['AU', 'UA', 'CG', 'GC'][index] }));
  for (const [groupBy, population] of [['function', rows], ['structure', rows], ['base', contextRows]]) {
    const baseline = traceColors(resultFor(population, groupBy));
    for (const selected of [[...population].reverse(), population.slice(1), [population[3]]]) {
      const result = resultFor(selected, groupBy), before = structuredClone(result);
      const snapshot = createPlotSnapshot({ result, buildId: 'fixture', parameter });
      const beforeCsv = csv(snapshot);
      const colors = traceColors(result);
      for (const [key, color] of colors) assert.equal(color, baseline.get(key));
      const traces = distributionTraces(result, { groupBy: 'none' });
      assert.deepEqual(traces.map(trace => trace.name.split(' (n=')[0]), result.series.map(series => series.label));
      result.series.forEach((series, index) => {
        assert.deepEqual(traces[index].x, series.x); assert.deepEqual(traces[index].y, series.y);
      });
      assert.deepEqual(result, before);
      assert.equal(csv(snapshot), beforeCsv);
    }
  }
});

test('colors use raw keys instead of repeated presentation labels', () => {
  const result = resultFor(rows, 'base');
  for (const series of result.series) series.label = 'Same label';
  assert.deepEqual(distributionTraces(result).map(trace => trace.line.color), ['#174a7e', '#8c3b2a', '#146c43', '#8659a1']);
  result.series.reverse();
  assert.deepEqual(distributionTraces(result).map(trace => trace.line.color), ['#8659a1', '#146c43', '#8c3b2a', '#174a7e']);
});

test('unknown and all-observation groups have explicit neutral and default colors', () => {
  const unknownRows = [{ id: 'unknown', angle: 5 }];
  for (const groupBy of ['base', 'method', 'function', 'structure']) {
    assert.equal(distributionTraces(resultFor(unknownRows, groupBy))[0].line.color, '#6a6256');
  }
  assert.equal(distributionTraces(resultFor(rows, 'none'))[0].line.color, '#174a7e');
});

test('arbitrary key mapping is bounded, namespaced and independent of previous calls', () => {
  const before = new Map(['rRNA', 'AU', 'UA', '__proto__', 'constructor', 'タグ', ''].map(key => [key, seriesColor(key, 'function')]));
  const palette = new Set(['#174a7e', '#8c3b2a', '#146c43', '#8659a1', '#be882e', '#32898c', '#ae567e', '#6a6256']);
  for (let index = 0; index < 10000; index++) assert(palette.has(seriesColor(`arbitrary-${index}`, 'function')));
  for (const [key, color] of before) assert.equal(seriesColor(key, 'function'), color);
  assert.equal(seriesColor('rRNA', 'function'), '#32898c');
  assert.equal(seriesColor('rRNA', 'structure'), '#146c43');
  assert.notEqual(seriesColor('AU', 'base'), seriesColor('UA', 'base'));
  assert.equal(seriesColor('cWW', 'interaction'), seriesColor('cWW', 'interactionFamily'));
  assert.match(seriesColor('constructor', 'method'), /^#[0-9a-f]{6}$/);
  assert.match(seriesColor('__proto__', 'base'), /^#[0-9a-f]{6}$/);
});
