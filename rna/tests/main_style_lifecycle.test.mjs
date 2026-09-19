import test from 'node:test';
import assert from 'node:assert/strict';
import { PureRnaExplorer } from '../app/PureRnaExplorer.js';
import { createPlotSnapshot, restyleTraceSnapshot, csv } from '../core/export.js';
import { distribution } from '../core/analysis.js';

const deferred = () => { let resolve; const promise = new Promise(done => { resolve = done; }); return { promise, resolve }; };
function setup() {
  const nodes = new Map(), plots = [];
  const app = new PureRnaExplorer({ root: { querySelector: id => {
    if (!nodes.has(id)) nodes.set(id, { dataset: {}, disabled: false }); return nodes.get(id);
  } }, repository: {} });
  app.manifest = { build_id: 'trace-test' };
  let renders = 0;
  app.render = async () => {
    renders++;
    const result = distribution([{ id: 'u', values: { x: -0 } }, { id: 'a', values: { x: 1.2345678901234567 } }],
      { id: 'x', range: [-1, 2] }, { ...app.state.display, bins: 64 });
    app.snapshots.distribution = createPlotSnapshot({ result });
    if (app.state.survey.loaded) app.snapshots.survey = createPlotSnapshot({ result });
    app.status('', 'ready');
  };
  app.plot = async (node, traces) => { plots.push({ node, traces }); };
  return { app, plots, renders: () => renders };
}

test('Trace snapshots retain frozen numerical arrays and replace only style metadata', () => {
  const result = distribution([{ id: 'a', values: { x: -0 } }], { id: 'x', range: [-1, 1] }, { bins: 64, traceStyle: 'filled' });
  const before = createPlotSnapshot({ result, displaySpec: result.displaySpec });
  const after = restyleTraceSnapshot(before, 'line');
  assert.equal(after.result.series, before.result.series);
  assert.equal(after.result.coverage, before.result.coverage);
  assert.equal(after.selection_spec, before.selection_spec);
  assert.equal(after.display_spec.traceStyle, 'line'); assert.equal(after.result.displaySpec.traceStyle, 'line');
  assert.equal(before.display_spec.traceStyle, 'filled');
  assert.equal(csv(after).replaceAll(after.snapshot_id, before.snapshot_id), csv(before));
  assert.throws(() => { after.result.series[0].values[0] = 10; }, TypeError);
  assert.throws(() => restyleTraceSnapshot(Object.freeze({ ...before }), 'line'), /owned frozen/);
  assert.throws(() => restyleTraceSnapshot(before, 'bogus'), /valid style/);
});

test('Trace-only work retains main, Survey and joint science and refreshes joint cache key', async () => {
  const { app, plots, renders } = setup(); app.state.survey.loaded = true;
  await app.requestRender();
  app.snapshots.joint = createPlotSnapshot({ result: { kind: 'joint', points: [], z: [[1]] }, displaySpec: app.state.display, joinSpec: app.state.joint });
  const previous = { ...app.snapshots };
  await app.setDisplay({ traceStyle: 'line' });
  assert.equal(renders(), 1); assert.equal(plots.length, 2);
  for (const key of ['distribution', 'survey', 'joint']) {
    assert.notEqual(app.snapshots[key].snapshot_id, previous[key].snapshot_id);
    assert.equal(app.snapshots[key].display_spec.traceStyle, 'line');
    assert.equal(app.snapshots[key].result.displaySpec.traceStyle, 'line');
  }
  assert.equal(app.snapshots.distribution.result.series, previous.distribution.result.series);
  assert.equal(app.snapshots.survey.result.series, previous.survey.result.series);
  assert.equal(app.snapshots.joint.result.points, previous.joint.result.points);
  assert.equal(app.completedJointKey, app.jointAnalysisKey(app.state));
  assert.equal(app.$('appStatus').dataset.state, 'ready');
});

test('All nontrace state changes and absent completed panels force full rendering', async t => {
  for (const mutate of [a => { a.state.selection.contexts = ['U']; }, a => { a.state.display.sigma = 0; },
    a => { a.state.family2Id = 'new'; }, a => { a.state.joint.endpoint = 'nt1'; },
    a => { a.state.survey.termId = 'other'; }, a => { a.state.survey.loaded = true; },
    a => { a.state.futureField = 1; }, a => { a.manifest.build_id = 'new'; },
    a => { a.repository.releaseUrl = 'new'; }, a => { a.snapshots.distribution = null; },
    a => { a.fullRenderComplete = false; }]) {
    const { app, renders } = setup(); await app.requestRender();
    // Build/release identity is captured in the completed key, not mutable state.
    mutate(app); await app.setDisplay({ traceStyle: 'line' });
    assert.equal(renders(), 2);
  }
});

test('Failed style plots publish no snapshot and following style repairs the entire page', async t => {
  t.mock.method(console, 'error', () => {});
  const { app, renders } = setup(); app.state.survey.loaded = true; await app.requestRender();
  const previous = { ...app.snapshots }; let calls = 0;
  app.plot = async () => { if (++calls === 2) throw Error('Injected Survey plot failure'); };
  await app.setDisplay({ traceStyle: 'line' });
  assert.equal(app.$('appStatus').dataset.state, 'error'); assert.equal(app.fullRenderComplete, false);
  assert.equal(app.snapshots.distribution, previous.distribution); assert.equal(app.snapshots.survey, previous.survey);
  assert.equal(app.$('filteredCsvDownload').disabled, true);
  await app.setDisplay({ traceStyle: 'filled' }); assert.equal(renders(), 2); assert.equal(app.fullRenderComplete, true);
});

test('Rapid style requests cancel stale publication and repair all panels', async () => {
  const { app, renders } = setup(); await app.requestRender();
  const entered = deferred(), finish = deferred();
  app.plot = async () => { entered.resolve(); await finish.promise; };
  const first = app.setDisplay({ traceStyle: 'line' }); await entered.promise;
  const second = app.setDisplay({ traceStyle: 'filled' }); await second;
  const current = app.snapshots.distribution; finish.resolve(); await first;
  assert.equal(renders(), 2); assert.equal(app.snapshots.distribution, current); assert.equal(app.state.display.traceStyle, 'filled');
  assert.equal(app.fullRenderComplete, true);
});

test('Style rendering drains coordinate labels changed while Plotly is pending', async () => {
  const { app } = setup(); await app.requestRender();
  app.completedCoordinateKey = 'coords'; app.completedCoordinateLabels = 'all';
  const entered = deferred(), finish = deferred();
  app.plot = async () => { entered.resolve(); await finish.promise; };
  app.renderCoordinatePlot = async () => { app.completedCoordinateLabels = app.state.survey.coordinateLabels; };
  const pending = app.setDisplay({ traceStyle: 'line' }); await entered.promise;
  await app.setCoordinateLabels('none'); finish.resolve(); await pending;
  assert.equal(app.completedCoordinateLabels, 'none'); assert.equal(app.$('appStatus').dataset.state, 'ready');
});
