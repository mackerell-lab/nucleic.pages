import test from 'node:test';
import assert from 'node:assert/strict';
import { PureRnaExplorer } from '../app/PureRnaExplorer.js';

const deferred = () => { let resolve, reject; const promise = new Promise((done, fail) => { resolve = done; reject = fail; }); return { promise, resolve, reject }; };
const term = { id: 't', label: 'Angle', level: 'residue' };
const observations = () => ['small', 'middle', 'large'].flatMap((opening_bin, bin) => [0, 1].map(index => ({
  id: `u-${bin}-${index}`, pdb_id: 'TEST', context: 'U', term_id: 't', opening_bin, value: bin * 10 + index,
}))).concat(['small', 'large'].map((opening_bin, index) => ({
  id: `c-${index}`, pdb_id: 'TEST', context: 'C', term_id: 't', opening_bin, value: 3 + index,
})));

function setup(t) {
  const previousDocument = globalThis.document;
  const node = () => ({ children: [], textContent: '', dataset: {}, disabled: false,
    setAttribute(key, value) { this[key] = value; }, addEventListener() {}, append(...values) { this.children.push(...values); },
    replaceChildren(...values) { this.children = values; } });
  globalThis.document = { createElement: node };
  t.after(() => { globalThis.document = previousDocument; });
  const nodes = new Map(), loads = [];
  const app = new PureRnaExplorer({ root: { querySelector: id => {
    if (!nodes.has(id)) nodes.set(id, node()); return nodes.get(id);
  } }, repository: { loadSurveyScalars: async id => { loads.push(id); return { rows: observations() }; } } });
  app.metadata = { entries: [{ pdb_id: 'TEST' }] };
  app.state.selection = { components: 'all', methods: [], contexts: [], includeEnds: true };
  app.state.survey.minimum = 1;
  app.checkpoint = async revision => app.current(revision);
  app.openingIncidences = rows => rows;
  app.rankingCountsRendered = true;
  const run = (terms = [term]) => { const request = app.capture(); return app.renderOpeningRanking(terms, request.state, request.revision, {}); };
  return { app, loads, run };
}

test('Computed counts use exact term/context ranks; minimum only changes sufficient count', async t => {
  const { app, loads, run } = setup(t);
  assert.equal(app.rankingCountValue('rows'), 'Not computed');
  await run();
  assert.deepEqual(app.rankingCounts, { status: 'computed', rows: 2, sufficient: 1 });
  assert.equal(app.$('baseGeometryRankRows').textContent, '2');
  assert.equal(app.$('baseGeometrySufficientRankRows').textContent, '1');
  const science = app.surveyRanks.map(({ sufficient, ...rank }) => rank);
  app.state.survey.minimum = 20; await run();
  assert.deepEqual(app.rankingCounts, { status: 'computed', rows: 2, sufficient: 0 });
  assert.deepEqual(app.surveyRanks.map(({ sufficient, ...rank }) => rank).sort((a, b) => a.context.localeCompare(b.context)), science.sort((a, b) => a.context.localeCompare(b.context)));
  assert.deepEqual(loads, ['t'], 'Minimum-only change reloaded scientific tables');
});

test('An explanatory empty table row never counts as a computed result', async t => {
  const { app, run } = setup(t);
  await run([]);
  assert.equal(app.$('baseGeometryRankingBody').children.length, 1);
  assert.deepEqual(app.rankingCounts, { status: 'computed', rows: 0, sufficient: 0 });
  assert.equal(app.$('baseGeometryRankRows').textContent, '0');
});

test('Pending and current failure show lifecycle states; retry publishes fresh counts', async t => {
  const { app, run } = setup(t), pending = deferred();
  const original = app.repository.loadSurveyScalars;
  app.repository.loadSurveyScalars = () => pending.promise;
  const first = run(), failed = assert.rejects(first, /Ranking load failed/);
  assert.equal(app.$('baseGeometryRankRows').textContent, 'Computing…');
  pending.reject(new Error('Ranking load failed')); await failed;
  assert.equal(app.$('baseGeometryRankRows').textContent, 'Unavailable');
  assert.equal(app.$('baseGeometrySufficientRankRows').textContent, 'Unavailable');
  app.repository.loadSurveyScalars = original; await run();
  assert.deepEqual(app.rankingCounts, { status: 'computed', rows: 2, sufficient: 1 });
});

test('Superseded ranking rejection cannot overwrite newer computed counts', async t => {
  const { app, run } = setup(t), pending = deferred();
  const original = app.repository.loadSurveyScalars;
  app.repository.loadSurveyScalars = () => pending.promise;
  const first = run(), failed = assert.rejects(first, /Old failure/);
  app.repository.loadSurveyScalars = original; await run();
  pending.reject(new Error('Old failure')); await failed;
  assert.deepEqual(app.rankingCounts, { status: 'computed', rows: 2, sufficient: 1 });
  assert.equal(app.$('baseGeometryRankRows').textContent, '2');
});

test('Reset clears computed status and old completion cannot restore it', async t => {
  const { app, run } = setup(t), pending = deferred();
  app.repository.loadSurveyScalars = () => pending.promise;
  const first = run();
  app.updateSelectors = app.renderControls = app.renderSurveyRankingControls = () => {};
  app.requestRender = async () => { app.capture(); };
  await app.resetFilters();
  assert.equal(app.$('baseGeometryRankRows').textContent, 'Not computed');
  assert.equal(app.$('baseGeometrySufficientRankRows').textContent, 'Not computed');
  pending.resolve({ rows: observations() }); await first;
  assert.deepEqual(app.surveyRanks, []);
  assert.equal(app.rankingCountValue('rows'), 'Not computed');
});

test('Survey failure before ranking begins marks existing computed counts unavailable', async t => {
  const { app, run } = setup(t);
  await run();
  app.surveyTerms = () => [term];
  app.state.survey.ranking = true;
  app.repository.loadSurveyScalars = async () => { throw new Error('Current scalar unavailable'); };
  const request = app.capture();
  await assert.rejects(app.renderSurvey(request.state, request.revision), /Current scalar unavailable/);
  assert.equal(app.$('baseGeometryRankRows').textContent, 'Unavailable');
});

test('An earlier main-panel failure cannot leave old ranking counts presented as current', async t => {
  t.mock.method(console, 'error', () => {});
  const { app, run } = setup(t); await run();
  app.state.survey.loaded = app.state.survey.ranking = true;
  const pending = deferred();
  app.render = () => pending.promise;
  const render = app.requestRender();
  assert.equal(app.$('baseGeometryRankRows').textContent, 'Computing…');
  pending.reject(new Error('Main family unavailable')); await render;
  assert.equal(app.$('baseGeometryRankRows').textContent, 'Unavailable');
  assert.equal(app.$('baseGeometrySufficientRankRows').textContent, 'Unavailable');
});

test('Superseded full-panel failure cannot overwrite newer committed ranking counts', async t => {
  t.mock.method(console, 'error', () => {});
  const { app } = setup(t), pending = deferred();
  app.state.survey.loaded = app.state.survey.ranking = true;
  app.render = () => pending.promise;
  const older = app.requestRender();
  app.render = request => app.renderOpeningRanking([term], request.state, request.revision, {});
  await app.requestRender();
  pending.reject(new Error('Older main failure')); await older;
  assert.deepEqual(app.rankingCounts, { status: 'computed', rows: 2, sufficient: 1 });
  assert.equal(app.$('baseGeometryRankRows').textContent, '2');
});

for (const ranking of [false, true]) {
  test(`Empty registry creates first-render metrics with ranking ${ranking}`, async t => {
    const { app, loads } = setup(t);
    app.surveyTerms = () => [];
    app.rankingCountsRendered = false;
    app.state.survey.ranking = ranking;
    const request = app.capture(); await app.renderSurvey(request.state, request.revision);
    assert.equal(app.rankingCounts.status, ranking ? 'computed' : 'not_computed');
    const metrics = app.$('baseGeometryStats').children.flatMap(node => node.children);
    assert.equal(metrics.find(node => node.id === 'baseGeometrySurveyTerms').textContent, '0');
    assert.equal(metrics.find(node => node.id === 'baseGeometryRankRows').textContent, ranking ? '0' : 'Not computed');
    assert.equal(metrics.find(node => node.id === 'baseGeometrySufficientRankRows').textContent, ranking ? '0' : 'Not computed');
    assert.equal(app.rankingCountsRendered, true);
    assert.equal(app.$('surveyCsvDownload').disabled, true);
    assert.deepEqual(loads, []);
  });
}

test('Empty group clears previous Survey and ranking instead of exposing stale exports', async t => {
  const { app, run } = setup(t); await run();
  const released = [], purged = [];
  app.repository.releaseSurvey = (...args) => released.push(args);
  app.plotly = { purge: node => purged.push(node) };
  app.lastSurveyTerm = 'old'; app.snapshots.survey = { snapshot_id: 'old' };
  app.surveyTerms = () => [{ ...term, group: 'available' }];
  app.state.survey.group = 'empty'; app.state.survey.ranking = true;
  const request = app.capture(); await app.renderSurvey(request.state, request.revision);
  assert.deepEqual(app.rankingCounts, { status: 'computed', rows: 0, sufficient: 0 });
  assert.deepEqual(app.surveyRanks, []);
  assert.equal(app.snapshots.survey, null);
  assert.equal(app.$('surveyCsvDownload').disabled, true);
  assert.equal(app.$('baseGeometryTermSelect').disabled, true);
  assert.equal(app.state.survey.termId, '');
  assert.equal(app.lastSurveyTerm, null);
  assert.deepEqual(released, [['scalars', 'old']]);
  assert.deepEqual(purged, [app.$('baseGeometryPlot')]);
  assert.equal(app.$('baseGeometryRankingBody').children.length, 1);
  assert.equal(app.$('surveyCoverageBody').children.length, 0);
  app.manifest = { build_id: 'recovery', survey: { opening_bins: [] } };
  app.plot = async () => {};
  app.snapshot = async options => ({ result: options.result });
  app.state.survey.group = 'available'; app.state.survey.ranking = false;
  const recovery = app.capture(); await app.renderSurvey(recovery.state, recovery.revision);
  assert.equal(app.$('baseGeometryTermSelect').disabled, false);
  assert.equal(app.$('surveyCsvDownload').disabled, false);
  assert.equal(app.state.survey.termId, 't');
  assert(app.snapshots.survey.result);
  assert.equal(app.rankingCountValue('rows'), 'Not computed');
});
