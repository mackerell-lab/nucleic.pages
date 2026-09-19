import test from 'node:test';
import assert from 'node:assert/strict';
import { PureRnaExplorer } from '../app/PureRnaExplorer.js';

const deferred = () => { let resolve; const promise = new Promise(done => { resolve = done; }); return { promise, resolve }; };
function setup() {
  const status = { dataset: { state: 'ready' } }, plot = {}, labels = { value: 'all' }, button = {};
  const nodes = { '#appStatus': status, '#coordinatePlot': plot, '#coordinateLabelsSelect': labels, '#jointCsvDownload': button };
  const app = new PureRnaExplorer({ root: { querySelector: id => nodes[id] }, repository: { loadFamily: async () => [] } });
  app.metadata = { entries: [] }; app.fullRenderComplete = true;
  app.completedCoordinateKey = 'coordinate-population'; app.completedCoordinateLabels = 'all'; app.coordinateSummary = [];
  app.snapshots.distribution = { snapshot_id: 'unchanged' };
  app.plot = async (_node, traces) => { plot.mode = traces[0].mode; };
  app.renderJoint = async () => {};
  return { app, status, plot, labels };
}

test('Label changes during joint-only loading are drained before ready', async () => {
  const { app, status, plot } = setup(), gate = deferred();
  app.renderJoint = async () => gate.promise;
  const joint = app.requestJointOnly();
  assert.equal(status.dataset.state, 'loading');
  await app.setCoordinateLabels('none');
  gate.resolve(); await joint;
  assert.equal(status.dataset.state, 'ready'); assert.equal(plot.mode, 'markers');
  assert.equal(app.completedCoordinateLabels, 'none');
  assert.equal(app.snapshots.distribution.snapshot_id, 'unchanged');
});

test('A joint revision drains a label commit superseded while queued', async () => {
  const { app, plot, status } = setup(), gate = deferred();
  app.commitQueue = gate.promise;
  const labels = app.setCoordinateLabels('none');
  const joint = app.requestJointOnly();
  gate.resolve(); await Promise.all([labels, joint]);
  assert.equal(status.dataset.state, 'ready'); assert.equal(plot.mode, 'markers');
  assert.equal(app.completedCoordinateLabels, 'none');
});

test('Latest label choice wins during awaited plotting without a new revision', async () => {
  const { app, plot } = setup(), gate = deferred(), entered = deferred();
  const original = app.plot; let calls = 0;
  app.plot = async (...args) => { if (++calls === 1) { entered.resolve(); await gate.promise; } return original(...args); };
  const first = app.setCoordinateLabels('none'); await entered.promise;
  const second = app.setCoordinateLabels('all'); gate.resolve(); await Promise.all([first, second]);
  assert.equal(plot.mode, 'markers+text'); assert.equal(app.completedCoordinateLabels, 'all');
  assert.equal(app.revision, 0);
});

test('A label click repairs failed full rendering through the full path', async () => {
  const { app, status, plot } = setup(); let repaired = 0;
  status.dataset.state = 'error'; app.fullRenderComplete = false;
  app.render = async ({ revision }) => { repaired++; await app.renderCoordinatePlot(app.coordinateSummary, revision, app.completedCoordinateKey); app.status('', 'ready'); };
  await app.setCoordinateLabels('none');
  assert.equal(repaired, 1); assert.equal(status.dataset.state, 'ready');
  assert.equal(plot.mode, 'markers'); assert.equal(app.fullRenderComplete, true);
});
