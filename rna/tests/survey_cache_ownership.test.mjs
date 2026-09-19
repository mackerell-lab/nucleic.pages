import test from 'node:test';
import assert from 'node:assert/strict';
import { PureRnaExplorer } from '../app/PureRnaExplorer.js';
import { RnaDataRepository } from '../core/repository.js';

const deferred = () => { let resolve; const promise = new Promise(done => { resolve = done; }); return { promise, resolve }; };

function setup() {
  const repository = new RnaDataRepository({ manifestUrl: 'http://localhost/manifest.json' });
  repository.loadSurveyScalars = id => repository.cached(`survey:scalars:${id}`, async () => ({ rows: [] }));
  repository.iterateSurveyCoordinates = async function* () { yield { rows: [] }; };
  const released = [];
  const release = repository.releaseSurvey.bind(repository);
  repository.releaseSurvey = (...args) => { released.push(args); return release(...args); };
  const app = new PureRnaExplorer({ root: {
    querySelector: () => { throw Error('A canceled commit must not access plot DOM'); },
    querySelectorAll: () => [],
  }, repository });
  app.manifest = { build_id: 'ownership', survey: {
    terms: [{ id: 'old', level: 'residue', group: 'linear' }], opening_bins: [],
    coordinates: { groups: { old: {}, 'new-active': {} } },
  } };
  app.families = []; app.metadata = { entries: [] };
  app.state.survey.termId = 'old'; app.state.survey.coordinateGroup = 'old';
  app.checkpoint = async revision => app.current(revision);
  const gate = deferred(), entered = deferred();
  app.commitQueue = gate.promise;
  const commit = app.commit.bind(app);
  app.commit = (...args) => { const operation = commit(...args); entered.resolve(); return operation; };
  return { app, repository, released, gate, entered };
}

for (const cancel of ['superseded', 'disposed']) {
  test(`${cancel} scalar commit cannot release or replace the current Survey owner`, async () => {
    const { app, repository, released, gate, entered } = setup();
    const activeTable = await repository.loadSurveyScalars('new-active');
    const request = app.capture();
    const pending = app.renderSurvey(request.state, request.revision);
    await entered.promise;
    if (cancel === 'disposed') app.dispose(); else app.capture();
    app.lastSurveyTerm = 'new-active';
    const activeSnapshot = { owner: 'new-active' }; app.snapshots.survey = activeSnapshot;
    gate.resolve(); await pending;
    assert.equal(app.lastSurveyTerm, 'new-active');
    assert.equal(app.snapshots.survey, activeSnapshot);
    assert.deepEqual(released, []);
    assert.equal(await repository.loadSurveyScalars('new-active'), activeTable);
  });

  test(`${cancel} coordinate commit cannot release or replace the current group owner`, async () => {
    const { app, repository, released, gate, entered } = setup();
    const activeTable = await repository.cached('survey:coordinates:collected:new-active', async () => ({ rows: [] }));
    const request = app.capture();
    const pending = app.renderCoordinates(request.state, request.revision);
    await entered.promise;
    if (cancel === 'disposed') app.dispose(); else app.capture();
    app.lastCoordinateGroup = 'new-active'; app.completedCoordinateKey = 'new-key';
    const activeSummary = []; app.coordinateSummary = activeSummary;
    gate.resolve(); await pending;
    assert.equal(app.lastCoordinateGroup, 'new-active');
    assert.equal(app.completedCoordinateKey, 'new-key');
    assert.equal(app.coordinateSummary, activeSummary);
    assert.deepEqual(released, []);
    assert.equal(await repository.promises.get('survey:coordinates:collected:new-active'), activeTable);
  });
}
