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
