import test from 'node:test';
import assert from 'node:assert/strict';
import { PureRnaExplorer } from '../app/PureRnaExplorer.js';

test('Overview releases Plotly instances before detaching their nodes', async () => {
  const plots = [{ id: 'old-alpha' }, { id: 'old-beta' }];
  const events = [];
  const container = {
    querySelectorAll(selector) { assert.equal(selector, '.rna-mini-plot'); return plots; },
    replaceChildren() { events.push('detach'); },
  };
  const app = new PureRnaExplorer({ root: { querySelector: () => container }, repository: {} });
  app.plotly = { purge: plot => events.push(plot.id) };
  app.parameters = () => [];
  await app.renderFamilyOverview([], { familyId: 'backbone' }, 0);
  assert.deepEqual(events, ['old-alpha', 'old-beta', 'detach']);
});

test('Cancelled overview is never reused as a completed overview', async () => {
  const previousDocument = globalThis.document;
  const node = () => ({ setAttribute() {}, append() {}, addEventListener() {} });
  globalThis.document = { createElement: node };
  try {
    let rebuilds = 0, finish, entered;
    const plotEntered = new Promise(resolve => { entered = resolve; });
    const container = { querySelectorAll: () => [], replaceChildren() { rebuilds++; }, append() {} };
    const app = new PureRnaExplorer({ root: { querySelector: () => container }, repository: {} });
    app.parameters = () => [{ id: 'alpha', period: 360 }];
    app.plot = () => new Promise(resolve => { finish = resolve; entered(); });
    const state = { familyId: 'backbone', parameterId: 'alpha', selection: {}, display: {} };
    const pending = app.renderFamilyOverview([], state, 0);
    await plotEntered;
    app.capture(); finish(); await pending;
    assert.equal(app.overviewKey, null);
    app.plot = async () => {};
    await app.renderFamilyOverview([], state, 1);
    assert.equal(rebuilds, 2, 'Cancelled partial overview was reused');
    assert.equal(typeof app.overviewKey, 'string');
    await app.renderFamilyOverview([], state, 1);
    assert.equal(rebuilds, 2, 'Completed identical overview was rebuilt');
  } finally { globalThis.document = previousDocument; }
});
