import test from 'node:test';
import assert from 'node:assert/strict';
import { PureRnaExplorer } from '../app/PureRnaExplorer.js';

function setup(t, terms) {
  const previousDocument = globalThis.document;
  const node = () => ({ dataset: {}, children: [], textContent: '', disabled: false,
    setAttribute() {}, append(...children) { this.children.push(...children); },
    replaceChildren(...children) { this.children = children; } });
  globalThis.document = { createElement: node };
  t.after(() => { globalThis.document = previousDocument; });
  const nodes = new Map();
  const app = new PureRnaExplorer({ root: { querySelector: id => {
    if (!nodes.has(id)) nodes.set(id, node()); return nodes.get(id);
  } }, repository: {} });
  app.manifest = { survey: { terms } };
  app.rankingCountsRendered = true;
  app.state.survey.loaded = true; app.state.survey.ranking = true;
  app.setRankingCounts('computed', [{ sufficient: true }]);
  app.surveyRanks = [{ term: { id: 'old' }, context: 'U', sufficient: true }];
  app.$('baseGeometryRankingBody').append(node());
  return app;
}

for (const [name, terms, group] of [
  ['empty registry', [], 'all'],
  ['group with no available terms', [{ id: 't', group: 'current' }], 'previous-group'],
]) test(`Ranking counts finish without stale rows for ${name}`, async t => {
  const app = setup(t, terms);
  app.state.survey.group = group;
  const request = app.capture();
  await app.renderSurvey(request.state, request.revision);
  assert.equal(app.rankingCounts.status, 'computed', 'Completed empty scope still reports pending ranking');
  assert.equal(app.rankingCounts.rows, 0);
  assert.equal(app.rankingCounts.sufficient, 0);
  assert.deepEqual(app.surveyRanks, []);
  assert.equal(app.$('baseGeometryRankRows').textContent, '0');
  assert.equal(app.$('baseGeometrySufficientRankRows').textContent, '0');
});
