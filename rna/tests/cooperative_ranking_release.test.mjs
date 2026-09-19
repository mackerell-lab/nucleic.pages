import test from 'node:test';
import assert from 'node:assert/strict';
import { PureRnaExplorer } from '../app/PureRnaExplorer.js';
import { RnaDataRepository } from '../core/repository.js';

function setup() {
  const repository = new RnaDataRepository({ manifestUrl: 'http://localhost/manifest.json' });
  repository.loadManifest = async () => ({ survey: { scalars: { terms: { transient: { path: 'transient.json' } } } } });
  repository.readJson = async () => ({ rows: [] });
  const button = { disabled: false, textContent: '' };
  const app = new PureRnaExplorer({ root: { querySelector: () => button, querySelectorAll: () => [] }, repository });
  app.lastSurveyTerm = 'active';
  app.metadata = { entries: [] };
  const ranking = () => {
    const request = app.capture();
    return app.renderOpeningRanking([{ id: 'transient' }], request.state, request.revision, {});
  };
  const retained = () => [...repository.promises.keys()].filter(key => key.startsWith('survey:scalars:transient'));
  return { app, repository, button, ranking, retained };
}

test('Canceled ranking releases its decoded projected scalar table', async () => {
  const { app, ranking, retained, button } = setup();
  app.checkpoint = async () => { app.capture(); return false; };
  await ranking();
  assert.deepEqual(retained(), []);
  assert.equal(app.rankingOwner, null);
  assert.equal(button.disabled, false);
});

test('Disposed ranking releases its decoded projected scalar table', async () => {
  const { app, ranking, retained } = setup();
  app.checkpoint = async () => { app.dispose(); return false; };
  await ranking();
  assert.deepEqual(retained(), []);
});

test('Scientific ranking failure releases decoded tables and preserves retry', async () => {
  const { app, ranking, retained, button } = setup();
  app.checkpoint = async revision => app.current(revision);
  app.surveyRows = () => { throw Error('Injected ranking processing failure'); };
  await assert.rejects(ranking(), /Injected ranking processing failure/);
  assert.deepEqual(retained(), []);
  assert.equal(app.rankingOwner, null);
  assert.equal(button.disabled, false);
});

test('A canceled ranking preserves a term newly owned by the active Survey', async () => {
  const { app, repository, ranking, retained } = setup();
  let active;
  app.checkpoint = async () => {
    active = await repository.loadSurveyScalars('transient');
    app.lastSurveyTerm = 'transient';
    app.capture();
    return false;
  };
  await ranking();
  assert(retained().includes('survey:scalars:transient'));
  assert.equal(await repository.loadSurveyScalars('transient'), active);
});

test('Failed ranking asset load evicts the rejected promise and retries', async () => {
  const { app, repository, ranking, retained, button } = setup();
  let loads = 0;
  repository.readJson = async () => {
    if (++loads === 1) throw Error('Injected scalar load failure');
    return { rows: [] };
  };
  await assert.rejects(ranking(), /Injected scalar load failure/);
  assert.deepEqual(retained(), []);
  assert.equal(button.disabled, false);
  app.checkpoint = async () => { app.capture(); return false; };
  await ranking();
  assert.equal(loads, 2);
  assert.deepEqual(retained(), []);
});
