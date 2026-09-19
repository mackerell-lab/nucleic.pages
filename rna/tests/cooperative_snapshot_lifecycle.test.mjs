import test from 'node:test';
import assert from 'node:assert/strict';
import { PureRnaExplorer } from '../app/PureRnaExplorer.js';
import { histogram2D } from '../core/analysis.js';

const deferred = () => {
  let resolve, reject;
  const promise = new Promise((done, fail) => { resolve = done; reject = fail; });
  return { promise, resolve, reject };
};
const parameter = { id: 'x', level: 'residue' };

function setup(panel) {
  const nodes = new Map(), snapshots = [], entered = deferred(), pending = deferred();
  const app = new PureRnaExplorer({ root: { querySelector: id => {
    if (!nodes.has(id)) nodes.set(id, { dataset: {}, disabled: false, textContent: '' });
    return nodes.get(id);
  }, querySelectorAll: () => [] }, repository: {
    loadFamily: async () => ({ rows: [] }), loadSurveyScalars: async () => ({ rows: [] }),
  } });
  app.manifest = { build_id: 'snapshot-lifecycle', survey: { terms: [parameter], opening_bins: [] } };
  app.metadata = { entries: [] }; app.families = [];
  app.parameter = () => parameter;
  app.state.familyId = app.state.family2Id = 'backbone';
  app.state.parameterId = app.state.parameter2Id = 'x';
  app.checkpoint = async revision => app.current(revision);
  app.snapshot = options => { snapshots.push(options); entered.resolve(); return pending.promise; };
  app.commit = async () => assert.fail('Canceled snapshot reached a DOM commit');
  app.plot = async () => assert.fail('Canceled snapshot reached Plotly');
  for (const key of ['distribution', 'joint', 'survey']) app.snapshots[key] = { snapshot_id: `previous-${key}` };
  app.completedJointKey = 'previous-joint-key';
  if (panel === 'joint') {
    app.fullRenderComplete = true;
    app.renderJoint = (state, revision) => app.renderJointResult(histogram2D([], parameter, parameter, { bins: 4 }), state, revision, {});
  } else if (panel === 'survey') {
    app.render = request => {
      app.status('Updating Survey'); app.$('surveyCsvDownload').disabled = true;
      return app.renderSurvey(request.state, request.revision);
    };
  }
  const run = () => panel === 'joint' ? app.requestJointOnly() : app.requestRender();
  return { app, snapshots, entered, pending, run };
}

for (const panel of ['main', 'joint', 'survey']) {
  for (const completion of ['resolve', 'reject']) {
    test(`${panel} snapshot late ${completion} cannot publish over newer failure`, async t => {
      const logs = []; t.mock.method(console, 'error', error => logs.push(error));
      const { app, snapshots, entered, pending, run } = setup(panel);
      const previous = { ...app.snapshots }, previousJointKey = app.completedJointKey;
      const operation = run(); await entered.promise;
      assert.equal(snapshots[0].revision, app.revision, 'Caller omitted its captured revision');
      app.capture(); app.fullRenderComplete = false;
      app.status('Newer render failed', 'error');
      for (const id of ['filteredCsvDownload', 'jointCsvDownload', 'surveyCsvDownload', 'plotProvenanceDownload']) app.$(id).disabled = true;
      const newerResult = { series: [{ values: [42] }] };
      if (completion === 'resolve') pending.resolve({ snapshot_id: 'superseded', result: snapshots[0].result });
      else pending.reject(Object.assign(new Error('Old snapshot failed'), { name: 'AbortError' }));
      await operation;
      assert.deepEqual(app.snapshots, previous);
      assert.equal(app.completedJointKey, previousJointKey);
      assert.equal(app.$('appStatus').dataset.state, 'error');
      assert.equal(app.$('appStatus').textContent, 'Newer render failed');
      assert.equal(app.fullRenderComplete, false);
      for (const id of ['filteredCsvDownload', 'jointCsvDownload', 'surveyCsvDownload', 'plotProvenanceDownload']) assert.equal(app.$(id).disabled, true);
      newerResult.series[0].values.push(43);
      assert.deepEqual(newerResult.series[0].values, [42, 43]);
      assert.deepEqual(logs, []);
    });
  }
}

test('Joint styling while a full snapshot is pending repairs full panels instead of reusing', async t => {
  t.mock.method(console, 'error', () => {});
  const { app, entered, pending, run } = setup('main');
  const first = run(); await entered.promise;
  const repairs = [];
  app.render = async request => { repairs.push(request); app.status('Repaired', 'ready'); };
  app.renderJointResult = async () => assert.fail('Style cache bypassed pending full snapshot');
  app.state.joint.palette = 'ocean';
  await app.requestJointOnly();
  assert.equal(repairs.length, 1);
  assert.equal(repairs[0].state.joint.palette, 'ocean');
  assert.equal(app.fullRenderComplete, true);
  pending.reject(Object.assign(new Error('Superseded full snapshot'), { name: 'AbortError' }));
  await first;
  assert.equal(app.fullRenderComplete, true);
  assert.equal(app.fullRenderOwner, null);
  assert.equal(app.$('appStatus').textContent, 'Repaired');
});

test('Current snapshot failure keeps full repair required for the next joint style', async t => {
  t.mock.method(console, 'error', () => {});
  const { app, entered, pending, run } = setup('main');
  const operation = run(); await entered.promise;
  pending.reject(new Error('Snapshot scheduling failed')); await operation;
  assert.equal(app.$('appStatus').dataset.state, 'error');
  assert.equal(app.fullRenderComplete, false);
  assert.equal(app.$('filteredCsvDownload').disabled, true);
  let repairs = 0;
  app.render = async () => { repairs++; app.status('Repaired', 'ready'); };
  await app.requestJointOnly();
  assert.equal(repairs, 1); assert.equal(app.fullRenderComplete, true);
});

test('Real app snapshot constructor cancels using its captured revision and leaves newer graph mutable', async () => {
  const { app } = setup('main');
  app.snapshot = PureRnaExplorer.prototype.snapshot;
  const entered = deferred(), pending = deferred();
  app.checkpoint = async revision => { entered.resolve(revision); await pending.promise; return app.current(revision); };
  const request = app.capture();
  // Exceed the default 262,144-operation budget regardless of machine speed.
  const oldResult = { kind: 'distribution', parameter, series: [], payload: Array(300000).fill(1) };
  const operation = app.snapshot({ result: oldResult, revision: request.revision, buildId: app.manifest.build_id });
  const rejected = assert.rejects(operation, { name: 'AbortError' });
  assert.equal(await entered.promise, request.revision);
  const newer = app.capture();
  const newResult = { kind: 'distribution', parameter: { ...parameter }, series: [{ values: [42] }] };
  pending.resolve(); await rejected;
  assert.equal(app.revision, newer.revision);
  assert.equal(app.snapshots.distribution.snapshot_id, 'previous-distribution');
  assert.equal(Object.isFrozen(newResult), false);
  newResult.series[0].values.push(43);
  assert.deepEqual(newResult.series[0].values, [42, 43]);
});

test('Real main rendering cancels inside cooperative snapshot before exporting or plotting', async t => {
  const logs = []; t.mock.method(console, 'error', error => logs.push(error));
  const { app } = setup('main');
  app.snapshot = PureRnaExplorer.prototype.snapshot;
  app.metadata = { entries: [{ pdb_id: 'TEST' }] };
  app.state.selection = { components: 'all', methods: [], includeEnds: true };
  app.state.display = { groupBy: 'none', sigma: 0 };
  app.repository.loadFamily = async () => ({ rows: Array.from({ length: 10000 }, (_, index) => ({
    id: `r${index}`, pdb_id: 'TEST', comp_id: 'U', values: { x: index / 10 },
  })) });
  const entered = deferred(), pending = deferred(); let checkpoints = 0;
  app.checkpoint = async revision => {
    if (++checkpoints > 3) { entered.resolve(); await pending.promise; }
    return app.current(revision);
  };
  const operation = app.requestRender(); await entered.promise;
  app.capture(); app.status('Newer request failed', 'error');
  pending.resolve(); await operation;
  assert.equal(checkpoints, 4, 'Expected three analysis checkpoints then snapshot work');
  assert.equal(app.snapshots.distribution.snapshot_id, 'previous-distribution');
  assert.equal(app.$('filteredCsvDownload').disabled, true);
  assert.equal(app.$('appStatus').textContent, 'Newer request failed');
  assert.equal(app.fullRenderComplete, false);
  assert.deepEqual(logs, []);
});
