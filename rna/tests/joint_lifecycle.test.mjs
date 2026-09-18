import test from 'node:test';
import assert from 'node:assert/strict';
import { PureRnaExplorer } from '../app/PureRnaExplorer.js';

const deferred = () => { let resolve; const promise = new Promise(done => { resolve = done; }); return { promise, resolve }; };
function setup() {
  const status = { dataset: { state: 'ready' }, textContent: '' }, button = { disabled: false };
  const pending = [];
  const app = new PureRnaExplorer({ root: { querySelector: id => id === '#appStatus' ? status : button },
    repository: { loadFamily: () => { const gate = deferred(); pending.push(gate); return gate.promise; } } });
  app.metadata = { entries: [] };
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
