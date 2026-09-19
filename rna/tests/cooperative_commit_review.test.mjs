import test from 'node:test';
import assert from 'node:assert/strict';
import { PureRnaExplorer } from '../app/PureRnaExplorer.js';

const deferred = () => { let resolve; const promise = new Promise(done => { resolve = done; }); return { promise, resolve }; };
function explorer(repository = {}) {
  const nodes = new Map();
  return new PureRnaExplorer({ root: { querySelector: id => {
    if (!nodes.has(id)) nodes.set(id, { dataset: {}, disabled: false });
    return nodes.get(id);
  } }, repository });
}

test('Yielding inside a commit cannot overwrite a newer failed render or prevent repair', async t => {
  const app = explorer(), entered = deferred(), resume = deferred();
  t.mock.method(console, 'error', () => {});
  t.mock.method(globalThis, 'setTimeout', callback => {
    entered.resolve(); resume.promise.then(callback); return 0;
  });
  let calls = 0;
  app.render = async ({ revision }) => {
    const call = ++calls;
    app.status('loading'); app.$('filteredCsvDownload').disabled = true;
    if (call === 2) throw Error('Newer render failed before committing');
    await app.commit(revision, async () => {
      if (call === 1 && !await app.checkpoint(revision)) return;
      app.snapshots.distribution = { revision };
      app.$('filteredCsvDownload').disabled = false;
      app.status('', 'ready');
    });
  };
  const older = app.requestRender(); await entered.promise;
  await app.requestRender();
  assert.equal(app.$('appStatus').dataset.state, 'error');
  resume.resolve(); await older;
  assert.equal(app.$('appStatus').dataset.state, 'error');
  assert.equal(app.$('filteredCsvDownload').disabled, true);
  assert.equal(app.snapshots.distribution, undefined);
  assert.equal(app.fullRenderComplete, false);
  await app.requestJointOnly();
  assert.equal(calls, 3);
  assert.equal(app.snapshots.distribution.revision, app.revision);
  assert.equal(app.fullRenderComplete, true);
  assert.equal(app.$('filteredCsvDownload').disabled, false);
});

test('A superseded residue lookup cannot replace current joint controls', async () => {
  const entered = deferred(), resume = deferred(), controls = [];
  let primaryLoads = 0;
  const app = explorer({ loadFamily: async family => {
    if (family === 'residue' && ++primaryLoads === 2) { entered.resolve(); await resume.promise; }
    return { rows: [] };
  }, loadRelations: async () => ({ rows: [] }) });
  app.manifest = { relations: { observations: {} } }; app.metadata = { entries: [] };
  app.parameter = family => ({ id: 'value', level: family === 'residue' ? 'residue' : 'pair' });
  app.checkpoint = async revision => app.current(revision);
  app.updateJointResidueControls = table => { controls.push(table); };
  app.renderJointResult = async () => { throw Error('Stale joint reached result rendering'); };
  app.state.familyId = 'residue'; app.state.family2Id = 'pair';
  app.state.parameterId = app.state.parameter2Id = 'value'; app.state.joint.mode = 'relation';
  const request = app.capture();
  const pending = app.renderJoint(request.state, request.revision, { rows: [] });
  await entered.promise;
  app.capture();
  resume.resolve(); await pending;
  assert.deepEqual(controls, [null], 'A stale primary lookup updated newer residue controls');
});
