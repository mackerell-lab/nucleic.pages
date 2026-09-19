import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import vm from 'node:vm';
import { surveyRankingDisplay, surveyRankingMean } from '../views/survey-ranking.js';

const source = readFileSync(new URL('../../js/pure-dna.js', import.meta.url), 'utf8');
function dnaFunction(name) {
  const start = source.indexOf(`function ${name}(`);
  assert(start >= 0, `Missing actual DNA oracle ${name}`);
  return source.slice(start, source.indexOf('\nfunction ', start + 1));
}
const oracle = vm.createContext({ state: {} });
vm.runInContext(['wrapCircular', 'circularDisplayValue', 'displayBaseGeometryValue', 'unitLabel'].map(dnaFunction).join('\n'), oracle);

test('Survey means match actual DNA display for declared periods and all modes', () => {
  let comparisons = 0;
  for (const period of [null, 180, 360, 720, 2 * Math.PI]) {
    const parameter = { period, unit: period === 2 * Math.PI ? 'rad' : 'deg' };
    const metadata = { isCircular: period !== null, period };
    for (const circularMode of ['wrap_360', 'signed_180', 'auto']) {
      oracle.state.circularMode = circularMode;
      const scale = period || 360;
      const values = [null, undefined, NaN, Infinity, -Infinity, -0, 0, scale / 2, -scale / 2, scale, -scale, 0.00049, -0.00049];
      for (let index = -100; index <= 100; index++) values.push(index * scale / 37);
      for (const value of values) {
        assert.equal(surveyRankingMean(value, parameter, circularMode), oracle.displayBaseGeometryValue(value, metadata, 3));
        comparisons++;
      }
    }
  }
  assert.equal(comparisons, 3210);
});

test('Bond angles stay linear; circular means retain DNA positive half endpoint', () => {
  assert.equal(surveyRankingMean(350, { period: 360 }, 'signed_180'), '-10.000');
  assert.equal(surveyRankingMean(350, { unit: 'deg', period: null }, 'signed_180'), '350.000');
  assert.equal(surveyRankingMean(-180, { period: 360 }, 'signed_180'), '180.000');
  assert.equal(surveyRankingMean(270, { period: 180 }, 'signed_180'), '90.000');
  assert.equal(surveyRankingMean(350, { period: 360 }, 'auto'), '350.000');
});

test('Units, precision and signed separation follow DNA without changing cached ranks', () => {
  for (const unit of ['A', 'Å', 'deg', 'degree', 'rad', '', undefined]) {
    const rank = Object.freeze({ term: Object.freeze({ id: 'torsion', label: 'Torsion', period: 360, unit }), context: 'AU', means: Object.freeze([350, 0, 10]), difference: -20.12567 });
    const shown = surveyRankingDisplay(rank, { circularMode: 'signed_180' });
    assert.equal(shown.unit, oracle.unitLabel(unit));
    assert.deepEqual(shown.means, ['-10.000', '0.000', '10.000']);
    assert.equal(shown.difference, '-20.1257');
    assert.equal(rank.means[0], 350); assert.equal(rank.difference, -20.12567);
  }
  for (const difference of [null, NaN, Infinity]) assert.equal(surveyRankingDisplay({ term: {}, means: [], difference }).difference, '-');
});

test('Current ranking row means exactly one selected term and context', () => {
  const rank = { term: { id: 'bond' }, context: 'G', means: [1, 2, 3], difference: 2 };
  for (const [termId, contexts, expected] of [['bond', ['G'], true], ['bond', [], false], ['bond', ['A', 'G'], false], ['bond', ['A'], false], ['other', ['G'], false]]) {
    const shown = surveyRankingDisplay(rank, { termId, contexts });
    assert.equal(shown.active, expected);
    assert.deepEqual(shown.means, ['1.000', '2.000', '3.000']);
    assert.equal(shown.difference, '2.0000');
  }
});
