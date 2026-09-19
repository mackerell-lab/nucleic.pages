import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import vm from 'node:vm';
import { summaryCards, number } from '../views/panels.js';
import { displayStatisticValue } from '../views/statistic-display.js';
import { csv } from '../core/export.js';

const dna = readFileSync(new URL('../../js/pure-dna.js', import.meta.url), 'utf8');
const actualFunction = name => { const start = dna.indexOf(`function ${name}(`); assert(start >= 0); return dna.slice(start, dna.indexOf('\nfunction ', start + 1)); };
const oracle = vm.createContext({});
vm.runInContext(['wrapCircular', 'circularDisplayValue'].map(actualFunction).join('\n'), oracle);

function renderedMetrics(result) {
  const previous = globalThis.document;
  const node = () => ({ children: [], setAttribute() {}, append(...nodes) { this.children.push(...nodes); }, replaceChildren(...nodes) { this.children = nodes; } });
  globalThis.document = { createElement: node };
  try {
    const root = node(); summaryCards(root, result);
    const metrics = root.children[0].children.find(child => child.className === 'metric-list');
    return Object.fromEntries(metrics.children.map(metric => metric.children.map(child => child.textContent)));
  } finally { globalThis.document = previous; }
}

const resultFor = (parameter, circularMode, cut, mean, peak = mean) => ({
  parameter, displaySpec: { circularMode }, displayCut: cut,
  series: [{ key: 'A', label: 'A', displayCut: cut, values: [mean], rows: [{ id: 'r1', pdb_id: 'TEST' }], statistics: { n: 1, pdbCount: 1, mean, peak, std: 12, resultant: 0.75, p05: mean, p95: mean } }],
});

test('Auto card mean and peak remain canonical below a shifted histogram seam', () => {
  const result = resultFor({ id: 'angle', period: 360, unit: 'deg' }, 'auto', 185, 5, 2.5);
  const raw = JSON.stringify(result); const exported = csv({ result });
  const shown = renderedMetrics(result);
  assert.equal(shown.Mean, '5 deg'); assert.equal(shown['Smoothed peak'], '2.5 deg');
  assert.equal(shown['Circular std. deviation'], '12 deg'); assert.equal(shown['Resultant length'], '0.75');
  assert.equal(JSON.stringify(result), raw); assert.equal(csv({ result }), exported);
});

test('Actual DNA oracle matches rendered cards across periods, modes and seams', () => {
  let comparisons = 0;
  for (const period of [180, 360, 720, 2 * Math.PI]) for (const mode of ['auto', 'wrap_360', 'signed_180']) {
    const parameter = { id: 'angle', period, unit: period === 2 * Math.PI ? 'rad' : 'deg' };
    for (const fraction of [-2, -1, -0.75, -0.5, -0.01, 0, 0.01, 0.5, 0.75, 1, 2]) {
      const value = fraction * period, cut = mode === 'signed_180' ? -period / 2 : 0.51 * period;
      const result = resultFor(parameter, mode, cut, value);
      const expected = oracle.circularDisplayValue(value, mode, period);
      const shown = renderedMetrics(result);
      assert.equal(displayStatisticValue(value, parameter, mode), expected);
      assert.equal(shown.Mean, `${number(expected)} ${parameter.unit}`);
      assert.equal(shown['Smoothed peak'], `${number(expected)} ${parameter.unit}`);
      comparisons++;
    }
  }
  assert.equal(comparisons, 132);
});

test('Signed positive half endpoint and undefined means match scientific policy', () => {
  for (const mean of [180, -180, 540]) assert.equal(renderedMetrics(resultFor({ period: 360, unit: 'deg' }, 'signed_180', -180, mean)).Mean, '180 deg');
  for (const mean of [null, NaN, Infinity]) {
    const shown = renderedMetrics(resultFor({ period: 360, unit: 'deg' }, 'auto', 185, mean));
    assert.equal(shown.Mean, 'Undefined'); assert.equal(shown['Smoothed peak'], 'Undefined');
  }
});

test('Degree-valued linear terms and distances never acquire angular wrapping', () => {
  for (const unit of ['deg', 'Å']) for (const mode of ['auto', 'wrap_360', 'signed_180']) {
    const result = resultFor({ id: 'linear', unit, period: null }, mode, 185, 370);
    const shown = renderedMetrics(result);
    assert.equal(shown.Mean, `370 ${unit}`); assert.equal(shown.P05, `370 ${unit}`); assert.equal(shown.P95, `370 ${unit}`);
    assert.equal(shown['Std. deviation'], `12 ${unit}`);
  }
});
