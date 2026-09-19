import test from 'node:test';
import assert from 'node:assert/strict';
import { distributionTraces } from '../views/panels.js';
import { distribution } from '../core/analysis.js';
import { createPlotSnapshot, csv } from '../core/export.js';

test('Distribution hover labels canonical and displayed angles without changing source values', () => {
  for (const period of [360, 180]) for (const circularMode of ['signed_180', 'wrap_360', 'auto']) {
    const rows = [period - 10, 5, 20].map((value, index) => ({ id: `r${index}`, pdb_id: 'TEST', values: { angle: value } }));
    const result = distribution(rows, { id: 'angle', label: 'RNA torsion', unit: 'deg', period }, { circularMode, normalization: 'density', groupBy: 'none', sigma: 0 });
    const snapshot = createPlotSnapshot({ result, buildId: 'test', parameter: result.parameter });
    const before = csv(snapshot), original = structuredClone(result);
    const trace = distributionTraces(snapshot.result, { traceStyle: 'line' })[0];
    assert.match(trace.hovertemplate, /RNA torsion \(deg\)/);
    assert.match(trace.hovertemplate, /View/);
    assert.match(trace.hovertemplate, /Angle/);
    assert.match(trace.hovertemplate, /Probability density \(smoothed\)/);
    assert.equal(trace.customdata.length, trace.x.length);
    trace.x.forEach((x, index) => assert.deepEqual(trace.customdata[index], [x, ((x % period) + period) % period]));
    assert.deepEqual(result, original);
    assert.equal(csv(snapshot), before);
    assert.notEqual(trace.x, snapshot.result.series[0].x);
    assert.notEqual(trace.y, snapshot.result.series[0].y);
  }
});

test('Linear hover retains units and normalization without inventing periodic angles', () => {
  for (const unit of ['Å', 'deg', '']) for (const normalization of ['probability', 'density']) {
    const result = { parameter: { id: 'linear', label: 'RNA measure', unit, period: null }, displaySpec: { normalization }, series: [{ key: 'U', x: [10, 20], y: [0, 1e-8] }] };
    const trace = distributionTraces(result, { normalization: normalization === 'density' ? 'probability' : 'density' })[0];
    assert(trace.hovertemplate.includes(`RNA measure${unit ? ` (${unit})` : ''}`));
    assert(trace.hovertemplate.includes(normalization === 'density' ? 'Probability density (smoothed)' : 'Probability (smoothed)'));
    assert(!trace.hovertemplate.includes('Angle'));
    assert(!trace.hovertemplate.includes('View'));
    assert.equal(trace.customdata, undefined);
    assert.deepEqual(trace.y, [0, 1e-8]);
  }
});

test('Hover defaults remain usable for minimal callers and labels are escaped', () => {
  const trace = distributionTraces({ parameter: { id: 'length' }, series: [{ key: 'All', x: [1], y: [1] }] }, { normalization: 'density' })[0];
  assert.match(trace.hovertemplate, /length/);
  assert.match(trace.hovertemplate, /Probability density/);
  const escaped = distributionTraces({ parameter: { id: 'x', label: 'A < B & C' }, series: [{ key: 'All', x: [1], y: [1] }] })[0];
  assert.match(escaped.hovertemplate, /A &lt; B &amp; C/);
});
