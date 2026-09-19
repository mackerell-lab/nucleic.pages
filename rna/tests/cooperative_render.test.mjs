import test from 'node:test';
import assert from 'node:assert/strict';
import { NucleicAcidExplorer } from '../app/NucleicAcidExplorer.js';
import { PureRnaExplorer } from '../app/PureRnaExplorer.js';
const deferred = () => { let resolve; const promise = new Promise(done => { resolve = done; }); return { promise, resolve }; };
const root = () => ({ querySelector: () => ({ dataset: {} }), querySelectorAll: () => [] });

test('Checkpoint gives queued ordinary input a task and skips stale scheduling', async t => {
  const previous = Object.getOwnPropertyDescriptor(globalThis, 'scheduler');
  t.after(() => previous ? Object.defineProperty(globalThis, 'scheduler', previous) : delete globalThis.scheduler);
  Object.defineProperty(globalThis, 'scheduler', { configurable: true, value: { yield() { throw Error('Boosted scheduler continuation must not be used'); } } });
  const app = new NucleicAcidExplorer({ root: root(), repository: {} });
  const events = [];
  setTimeout(() => { events.push('input'); app.capture(); }, 0);
  const pending = app.checkpoint(0); events.push('scheduled');
  assert.equal(await pending, false); assert.deepEqual(events, ['scheduled', 'input']);
  const timer = t.mock.method(globalThis, 'setTimeout', () => { throw Error('Stale work scheduled another task'); });
  assert.equal(await app.checkpoint(0), false); timer.mock.restore();
  assert.equal(await app.checkpoint(1), true);
});

test('Timer checkpoint permits task input and disposal without relying on animation frames', async t => {
  const previous = Object.getOwnPropertyDescriptor(globalThis, 'scheduler');
  t.after(() => previous ? Object.defineProperty(globalThis, 'scheduler', previous) : delete globalThis.scheduler);
  Object.defineProperty(globalThis, 'scheduler', { configurable: true, value: undefined });
  const app = new NucleicAcidExplorer({ root: root(), repository: {} });
  const input = new Promise(resolve => setTimeout(() => { app.dispose(); resolve(); }, 0));
  const checkpoint = app.checkpoint(0); await input;
  assert.equal(await checkpoint, false);
});

test('Overview checkpoint cancels before constructing another card and releases queued commit', async t => {
  const previous = globalThis.document;
  t.after(() => { globalThis.document = previous; });
  let created = 0;
  globalThis.document = { createElement: () => { created++; return { setAttribute() {}, append() {}, addEventListener() {} }; } };
  const container = { querySelectorAll: () => [], replaceChildren() {}, append() {} };
  const app = new PureRnaExplorer({ root: { querySelector: () => container }, repository: {} });
  app.parameters = () => [{ id: 'alpha', period: 360 }, { id: 'beta', period: 360 }];
  const gate = deferred(), entered = deferred(); let checkpoints = 0, plots = 0;
  app.checkpoint = async revision => { if (++checkpoints === 2) { entered.resolve(); await gate.promise; } return app.current(revision); };
  app.plot = async () => { plots++; };
  const state = { familyId: 'backbone', parameterId: 'alpha', selection: {}, display: {} };
  const old = app.commit(0, () => app.renderFamilyOverview([], state, 0));
  await entered.promise; const createdBefore = created;
  const request = app.capture(); let newestCommitted = false;
  const newest = app.commit(request.revision, () => { newestCommitted = true; });
  gate.resolve(); assert.equal(await old, false); assert.equal(await newest, true);
  assert(newestCommitted); assert.equal(plots, 1); assert.equal(created, createdBefore); assert.equal(app.overviewKey, null);
});

test('Main render cancellation at its first cached-data checkpoint performs no scientific or DOM work', async () => {
  const nodes = new Map(); const app = new PureRnaExplorer({ root: { querySelector: id => {
    if (!nodes.has(id)) nodes.set(id, { dataset: {}, disabled: false }); return nodes.get(id);
  } }, repository: { loadFamily: async () => [] } });
  app.parameter = () => { throw Error('Cancelled work reached parameter selection'); };
  const gate = deferred(), entered = deferred();
  app.checkpoint = async revision => { entered.resolve(); await gate.promise; return app.current(revision); };
  const pending = app.requestRender(); await entered.promise;
  app.capture(); app.status('New request failed', 'error'); gate.resolve(); await pending;
  assert.equal(app.$('appStatus').dataset.state, 'error'); assert.equal(app.$('filteredCsvDownload').disabled, true);
  assert.equal(app.fullRenderComplete, false); assert.deepEqual(app.snapshots, {});
});
