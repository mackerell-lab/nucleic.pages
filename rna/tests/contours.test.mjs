import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import vm from 'node:vm';
import { contourSpacing, jointContourConfig } from '../core/contours.js';

const dnaSource = readFileSync(new URL('../../js/pure-dna.js', import.meta.url), 'utf8');
// Execute the independent, unmodified DNA implementation as the oracle.
const source = dnaSource.slice(dnaSource.indexOf('function nextNiceStepAtLeast('), dnaSource.indexOf('function buildJointPlotTraces('));
const options = dnaSource.match(/const JOINT_CONTOUR_WIDTH_OPTIONS = (\[[\s\S]*?\]);/)[1];
const targetSource = dnaSource.slice(dnaSource.indexOf('function currentJointContourTargetLevels('), dnaSource.indexOf('function nextNiceStepAtLeast('));
const oracle = vm.createContext({ state: {} });
vm.runInContext(`const JOINT_CONTOUR_WIDTH_OPTIONS = ${options};\n${targetSource}`, oracle);
vm.runInContext(source, oracle);
const plain = value => JSON.parse(JSON.stringify(value));

test('RNA numeric spacing matches actual DNA code across scales and presets', () => {
  let cases = 0;
  for (let exponent = -300; exponent <= 300; exponent += 3) {
    for (const factor of [1, 1.4, 2.5, 6.7, 9.99]) {
      const range = factor * 10 ** exponent;
      const spacing = contourSpacing(range);
      assert.deepEqual(spacing, plain(oracle.buildContourSpacingSet(range)));
      assert(spacing.wide > spacing.standard && spacing.standard > spacing.tight && spacing.tight > 0);
      for (const [contourCount, preset] of [[6, 'wide'], [12, 'standard'], [24, 'tight']]) {
        oracle.state = { jointContourWidth: preset, jointContourLabels: 'on' };
        const low = exponent > 0 ? -range / 2 : 0;
        assert.deepEqual(jointContourConfig(low, low + range, { contourCount, labels: true }), plain(oracle.buildJointContourConfig(low, low + range)));
        cases++;
      }
    }
  }
  assert.equal(cases, 3015);
});

test('known nice intervals are explicit in linear and logarithmic display units', () => {
  assert.deepEqual(contourSpacing(1), { wide: 0.2, standard: 0.1, tight: 0.05 });
  for (const [count, size] of [[6, 2], [12, 1], [24, 0.5]]) {
    assert.deepEqual(jointContourConfig(-8, -1, { contourCount: count }), { autocontour: false, contours: { showlabels: false, start: -8, end: -1, size } });
  }
});

test('empty, flat, reversed and nonfinite ranges use bounded automatic fallback', () => {
  for (const [zmin, zmax] of [[0, undefined], [0, 0], [-8, -8], [1, 0], [NaN, 1], [0, Infinity], [-Number.MAX_VALUE, Number.MAX_VALUE], [0, Number.MIN_VALUE]]) {
    for (const [contourCount, preset, levels] of [[6, 'wide', 7], [12, 'standard', 10], [24, 'tight', 14]]) {
      const actual = jointContourConfig(zmin, zmax, { contourCount, labels: true });
      assert.deepEqual(actual, { autocontour: true, ncontours: levels, contours: { showlabels: true } });
      // DNA has no finite-range overflow guard; other fallback cases match it.
      if (!(Number.isFinite(zmin) && Number.isFinite(zmax) && !Number.isFinite(zmax - zmin))) {
        oracle.state = { jointContourWidth: preset, jointContourLabels: 'on' };
        assert.deepEqual(actual, plain(oracle.buildJointContourConfig(zmin, zmax)));
      }
    }
  }
  assert.deepEqual(jointContourConfig(0, 0, { contourCount: -1 }), { autocontour: true, ncontours: 10, contours: { showlabels: false } });
});

test('configuration does not mutate caller preferences or read observation data', () => {
  const preferences = Object.freeze({ contourCount: 24, labels: true });
  assert.equal(jointContourConfig(0, 1e-20, preferences).contours.size, contourSpacing(1e-20).tight);
  assert.deepEqual(preferences, { contourCount: 24, labels: true });
});
