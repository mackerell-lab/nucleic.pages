import test from 'node:test';
import assert from 'node:assert/strict';
import { PureRnaExplorer } from '../app/PureRnaExplorer.js';

const deferred = () => { let resolve, reject; const promise = new Promise((done, fail) => { resolve = done; reject = fail; }); return { promise, resolve, reject }; };
function setup() {
  const status = { dataset: { state: 'ready' }, textContent: '' }, button = { disabled: false };
  const pending = [];
  const app = new PureRnaExplorer({ root: { querySelector: id => id === '#appStatus' ? status : button },
    repository: { loadFamily: () => { const gate = deferred(); pending.push(gate); return gate.promise; } } });
  app.metadata = { entries: [] };
  // Simulate an already completed initial full render for joint-only tests.
  app.fullRenderComplete = true;
  app.snapshots.distribution = { snapshot_id: 'keep-1d' };
  app.renderJoint = async (_state, revision) => {
    if (app.current(revision)) { app.snapshots.joint = { revision }; button.disabled = false; }
  };
  return { app, status, button, pending };
}

test('Joint refresh disables stale CSV until the current request completes', async () => {
  const { app, status, button, pending } = setup();
  const first = app.requestJointOnly();
  assert.equal(status.dataset.state, 'loading');
  assert.equal(button.disabled, true);
  const second = app.requestJointOnly();
  pending[0].resolve([]); await first;
  assert.equal(status.dataset.state, 'loading', 'Superseded request reported ready');
  assert.equal(button.disabled, true);
  pending[1].resolve([]); await second;
  assert.equal(status.dataset.state, 'ready');
  assert.equal(button.disabled, false);
  assert.equal(app.snapshots.joint.revision, 2);
  assert.equal(app.snapshots.distribution.snapshot_id, 'keep-1d');
});

test('Joint controls during a full refresh preserve all pending panels', async () => {
  const { app } = setup();
  const full = [];
  app.render = async request => {
    const gate = deferred(); full.push({ gate, request });
    await gate.promise;
    if (!app.current(request.revision)) return;
    for (const panel of ['distribution', 'survey', 'coordinates', 'joint']) {
      app.snapshots[panel] = { contexts: request.state.selection.contexts, palette: request.state.joint.palette };
    }
  };
  const first = app.setSelection({ contexts: ['U'] });
  app.state.joint.palette = 'viridis';
  const second = app.requestJointOnly();
  assert.equal(full.length, 2, 'Joint-only refresh abandoned unfinished full panels');
  full[0].gate.resolve(); await first;
  app.state.joint.palette = 'ocean';
  const third = app.requestJointOnly();
  assert.equal(full.length, 3, 'Old completion cleared the newer full-render owner');
  full[1].gate.resolve(); await second;
  full[2].gate.resolve(); await third;
  for (const panel of ['distribution', 'survey', 'coordinates', 'joint']) {
    assert.deepEqual(app.snapshots[panel], { contexts: ['U'], palette: 'ocean' });
  }
  assert.equal(app.fullRenderOwner, null);
});


test('Joint controls repair a failed full render and then resume independent updates', async t => {
  t.mock.method(console, 'error', () => {});
  const { app, status, pending } = setup();
  let fullCalls = 0;
  app.render = async request => {
    fullCalls++;
    if (fullCalls === 1) throw new Error('Injected full-panel failure');
    app.snapshots.distribution = { snapshot_id: 'repaired', contexts: request.state.selection.contexts };
    app.status('', 'ready');
  };
  await app.setSelection({ contexts: ['U'] });
  assert.equal(status.dataset.state, 'error');
  assert.equal(app.fullRenderOwner, null);
  const repair = app.requestJointOnly();
  assert.equal(fullCalls, 2, 'A failed full render was mistaken for completed panels');
  await repair;
  assert.equal(status.dataset.state, 'ready');
  assert.deepEqual(app.snapshots.distribution.contexts, ['U']);
  assert.equal(app.fullRenderComplete, true);
  const independent = app.requestJointOnly();
  assert.equal(fullCalls, 2, 'Successful recovery disabled the joint-only optimization');
  pending[0].resolve([]); await independent;
  assert.equal(app.snapshots.distribution.snapshot_id, 'repaired');
});

test('An older successful render cannot mark a newer failed render complete', async t => {
  t.mock.method(console, 'error', () => {});
  const { app, status } = setup();
  const full = [];
  app.render = async request => {
    const gate = deferred(); full.push({ gate, request });
    await gate.promise;
    if (app.current(request.revision)) {
      app.snapshots.distribution = { contexts: request.state.selection.contexts };
      app.status('', 'ready');
    }
  };
  const first = app.setSelection({ contexts: ['A'] });
  const second = app.setSelection({ contexts: ['U'] });
  full[1].gate.reject(new Error('Newer full render failed')); await second;
  full[0].gate.resolve(); await first;
  assert.equal(status.dataset.state, 'error');
  assert.equal(app.fullRenderComplete, false, 'Stale success marked failed current panels complete');
  const repair = app.requestJointOnly();
  assert.equal(full.length, 3, 'Joint click did not repair the newer failed selection');
  full[2].gate.resolve(); await repair;
  assert.equal(app.fullRenderComplete, true);
  assert.deepEqual(app.snapshots.distribution.contexts, ['U']);
});
