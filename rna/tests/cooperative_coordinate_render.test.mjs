import test from 'node:test';
import assert from 'node:assert/strict';
import { PureRnaExplorer } from '../app/PureRnaExplorer.js';
import { CoordinateSummary } from '../core/coordinates.js';
const node = () => ({ children: [], value: '', setAttribute(key, value) { this[key] = value; }, append(...children) { this.children.push(...children); }, replaceChildren(...children) { this.children = children; } });
function setup(t, cancelAt) {
  const previous = globalThis.document; globalThis.document = { createElement: node };
  t.after(() => { globalThis.document = previous; });
  const nodes = new Map(), events = { yielded: [], accumulated: [], plots: 0, closed: 0, releases: [] };
  let app, resolveInput;
  const input = new Promise(resolve => { resolveInput = resolve; });
  const enqueueInput = () => setTimeout(() => { app.capture(); resolveInput(); }, 0);
  const rows = Array.from({ length: 6 }, (_, i) => ({ id: `a${i}`, pdb_id: 'TEST', residue_id: `r${i}`, context: 'G', atom_label: 'N1', xyz: [i, 2 * i, -i], status: 'available' }));
  app = new PureRnaExplorer({ root: { querySelector(id) { if (!nodes.has(id)) nodes.set(id, node()); return nodes.get(id); } }, repository: {
    async *iterateSurveyCoordinates() {
      try {
        for (let index = 0; index < 3; index++) {
          if (cancelAt === index) enqueueInput();
          events.yielded.push(index); yield { rows: rows.slice(2 * index, 2 * index + 2) };
        }
      } finally { events.closed++; if (cancelAt === 'final') enqueueInput(); }
    },
    releaseSurvey: (...args) => events.releases.push(args),
  } });
  app.manifest = { build_id: 'coordinate-checkpoint', survey: { coordinates: { groups: { rna_standard_base_G: {} } } } };
  app.metadata = { entries: [{ pdb_id: 'TEST', method: 'X-RAY DIFFRACTION', resolution: 2, profiles: { relaxed: true } }] };
  app.state.survey.coordinateGroup = 'rna_standard_base_G'; app.lastCoordinateGroup = 'retained-old-group';
  app.coordinateSummary = [{ sentinel: 'keep completed summary until replacement is complete' }];
  app.renderCoordinatePlot = async () => { events.plots++; };
  const originalAdd = CoordinateSummary.prototype.add;
  t.mock.method(CoordinateSummary.prototype, 'add', function (row) { events.accumulated.push(row.id); return originalAdd.call(this, row); });
  return { app, events, input, rows };
}

for (const index of [0, 1, 2]) test(`Queued ordinary input cancels coordinate partition ${index} before accumulation`, async t => {
  const { app, events, input, rows } = setup(t, index);
  const prior = app.coordinateSummary, request = app.capture();
  await app.renderCoordinates(request.state, request.revision); await input;
  assert.deepEqual(events.accumulated, rows.slice(0, index * 2).map(row => row.id));
  assert.deepEqual(events.yielded, Array.from({ length: index + 1 }, (_, i) => i));
  assert.equal(events.closed, 1, 'Early return must close the async iterator');
  assert.equal(events.plots, 0); assert.equal(app.coordinateSummary, prior); assert.equal(app.completedCoordinateKey, null);
  assert.equal(app.lastCoordinateGroup, 'retained-old-group'); assert.deepEqual(events.releases, []);
});

test('Queued input after final partition prevents stale result publication and cache release', async t => {
  const { app, events, input, rows } = setup(t, 'final');
  const prior = app.coordinateSummary, request = app.capture();
  await app.renderCoordinates(request.state, request.revision); await input;
  assert.deepEqual(events.accumulated, rows.map(row => row.id)); assert.equal(events.closed, 1);
  assert.equal(events.plots, 0); assert.equal(app.coordinateSummary, prior); assert.equal(app.completedCoordinateKey, null);
  assert.equal(app.lastCoordinateGroup, 'retained-old-group'); assert.deepEqual(events.releases, []);
});

test('Current streamed coordinates accumulate every row once and completed reuse adds no checkpoints', async t => {
  const { app, events, rows } = setup(t, null);
  const request = app.capture(); await app.renderCoordinates(request.state, request.revision);
  assert.deepEqual(events.accumulated, rows.map(row => row.id)); assert.equal(events.closed, 1); assert.equal(events.plots, 1);
  assert.deepEqual(app.coordinateSummary[0].mean, [2.5, 5, -2.5]); assert.equal(app.coordinateSummary[0].n, 6);
  assert.equal(app.coordinateSummary[0].residues, 6); assert.equal(app.coordinateSummary[0].rms, Math.sqrt(17.5));
  assert.deepEqual(events.releases, [['coordinates', 'retained-old-group']]);
  app.completedCoordinateLabels = app.state.survey.coordinateLabels;
  app.checkpoint = async () => { throw Error('Completed coordinate reuse must not add task delays'); };
  const next = app.capture(); await app.renderCoordinates(next.state, next.revision);
  assert.equal(events.plots, 1); assert.equal(events.closed, 1); assert.equal(events.accumulated.length, 6);
});
