import test from 'node:test';
import assert from 'node:assert/strict';
import { PureRnaExplorer } from '../app/PureRnaExplorer.js';

const deferred = () => {
  let resolve, reject;
  const promise = new Promise((done, fail) => { resolve = done; reject = fail; });
  return { promise, resolve, reject };
};

function setup(t) {
  const originalDocument = globalThis.document;
  const node = () => ({ dataset: {}, children: [], setAttribute() {}, append(...children) { this.children.push(...children); },
    replaceChildren(...children) { this.children = children; } });
  globalThis.document = { createElement: node };
  t.after(() => { globalThis.document = originalDocument; });
  const nodes = new Map(), loads = [], plots = [];
  const rows = [
    { id: 'u', pdb_id: '1AAA', comp_id: 'U', values: { x: 10, y: 20, z: 30 } },
    { id: 'a', pdb_id: '1AAA', comp_id: 'A', values: { x: 40, y: 60, z: 90 } },
  ];
  const app = new PureRnaExplorer({ root: { querySelector: id => {
    if (!nodes.has(id)) nodes.set(id, node());
    return nodes.get(id);
  } }, repository: { loadFamily: async id => {
    loads.push(id);
    return { rows: id === 'sugar' ? rows.map(row => ({ ...row, values: { x: row.values.x + 100, y: row.values.y + 100, z: row.values.z + 100 } })) : rows };
  } }, plotly: { purge() {} } });
  app.manifest = { build_id: 'joint-style-test' }; app.families = [];
  app.metadata = { entries: [{ pdb_id: '1AAA', method: 'X-RAY DIFFRACTION', resolution: 2 }] };
  app.parameter = (_family, id) => ({ id, label: id, level: 'residue', range: [0, 300] });
  app.updateJointResidueControls = () => {};
  app.plot = async (target, traces, layout) => { plots.push({ target, traces, layout }); };
  app.fullRenderComplete = true;
  app.state.familyId = 'backbone'; app.state.family2Id = 'backbone';
  app.state.parameterId = 'x'; app.state.parameter2Id = 'y';
  app.state.selection = { components: 'all', methods: [], includeEnds: true };
  app.state.display = { sigma: 0, fine: false, normalization: 'probability', circularMode: 'wrap_360' };
  app.snapshots.distribution = { snapshot_id: 'unchanged-main' };
  return { app, nodes, loads, plots, status: app.$('appStatus'), button: app.$('jointCsvDownload') };
}

async function seeded(t) {
  const state = setup(t);
  await state.app.requestJointOnly();
  assert.equal(state.status.dataset.state, 'ready');
  assert.equal(state.app.snapshots.joint.result.points.length, 2);
  assert.equal(typeof state.app.completedJointKey, 'string');
  state.loads.length = 0; state.plots.length = 0;
  return state;
}

test('Joint visual controls reuse one completed analysis with fresh export controls', async t => {
  const { app, loads, button } = await seeded(t);
  const initial = app.snapshots.joint, result = initial.result;
  for (const patch of [{ palette: 'ocean' }, { colorScale: 'log' }, { type: 'contour' },
    { labels: true }, { contourCount: 24 }, { type: 'filled_contour' }]) {
    Object.assign(app.state.joint, patch);
    const previous = app.snapshots.joint;
    await app.requestJointOnly();
    assert.equal(app.snapshots.joint.result, result, JSON.stringify(patch));
    assert.notEqual(app.snapshots.joint, previous, 'Export provenance retained an old visual setting');
    assert.deepEqual(app.snapshots.joint.join_spec, app.state.joint);
    assert.deepEqual(app.snapshots.joint.provenance.axis_selections, initial.provenance.axis_selections);
    assert.deepEqual(app.snapshots.joint.provenance.join_diagnostics, initial.provenance.join_diagnostics);
    assert.equal(button.disabled, false);
  }
  assert.deepEqual(loads, [], 'A visual control reloaded and selected a measurement family');
  assert.equal(app.snapshots.distribution.snapshot_id, 'unchanged-main');
});

test('Selection, axes and numerical controls invalidate a completed joint analysis', async t => {
  const changes = [
    ['primary family with shared parameter IDs', state => { state.familyId = 'sugar'; }],
    ['secondary family with shared parameter IDs', state => { state.family2Id = 'sugar'; }],
    ['primary parameter', state => { state.parameterId = 'z'; }],
    ['secondary parameter', state => { state.parameter2Id = 'z'; }],
    ...Object.entries({ components: 'strict', methods: ['xray'], resolutionMax: 1.5, contexts: ['U'],
      functions: ['unknown'], structures: ['unknown'], subtypes: ['unknown'], puckerStates: ["C3'-endo"],
      includeEnds: false, pairPolicy: 'all', interactionFamilies: ['cWW'], stemOnly: true,
      annotationEndpointPolicy: 'any', chiStates: ['anti'], resolutionKnown: true, search: '1AAA' })
      .map(([key, value]) => [`selection.${key}`, state => { state.selection[key] = value; }]),
    ...Object.entries({ sigma: 1.6, fine: true, normalization: 'density', circularMode: 'signed', range: [0, 200], groupBy: 'function', traceStyle: 'line' })
      .map(([key, value]) => [`display.${key}`, state => { state.display[key] = value; }]),
    ...Object.entries({ mode: 'relation', endpoint: 'first', residueContexts: ['U'], residuePuckers: ["C3'-endo"] })
      .map(([key, value]) => [`joint.${key}`, state => { state.joint[key] = value; }]),
  ];
  for (const [label, mutate] of changes) await t.test(label, async sub => {
    const { app, loads } = await seeded(sub), previous = app.snapshots.joint.result;
    mutate(app.state);
    await app.requestJointOnly();
    assert(loads.length > 0, `${label} reused the previous population`);
    assert.notEqual(app.snapshots.joint?.result, previous, `${label} reused the previous histogram`);
    if (label === 'primary family with shared parameter IDs') assert.deepEqual(app.snapshots.joint.result.points.map(point => point.x), [110, 140]);
    if (label === 'secondary family with shared parameter IDs') assert.deepEqual(app.snapshots.joint.result.points.map(point => point.y), [120, 160]);
  });
});

test('A failed visual update keeps CSV disabled and retries the valid analysis', async t => {
  t.mock.method(console, 'error', () => {});
  const { app, loads, status, button } = await seeded(t), previous = app.snapshots.joint;
  const originalPlot = app.plot;
  app.plot = async () => { throw new Error('Injected Plotly failure'); };
  app.state.joint.palette = 'ocean';
  await app.requestJointOnly();
  assert.equal(status.dataset.state, 'error'); assert.equal(button.disabled, true);
  assert.equal(app.snapshots.joint, previous, 'A failed plot published a successful export');
  app.plot = originalPlot;
  await app.requestJointOnly();
  assert.equal(status.dataset.state, 'ready'); assert.equal(button.disabled, false);
  assert.equal(app.snapshots.joint.result, previous.result);
  assert.equal(app.snapshots.joint.join_spec.palette, 'ocean');
  assert.deepEqual(loads, []);
});

test('Rapid style updates commit only the latest completed controls', async t => {
  const { app, loads, status, button } = await seeded(t), previous = app.snapshots.joint;
  const entered = deferred(), finish = deferred(), originalPlot = app.plot;
  let count = 0;
  app.plot = async (...args) => { if (++count === 1) { entered.resolve(); await finish.promise; } await originalPlot(...args); };
  app.state.joint.palette = 'ocean';
  const first = app.requestJointOnly(); await entered.promise;
  app.state.joint.palette = 'viridis'; app.state.joint.type = 'contour';
  const second = app.requestJointOnly();
  assert.equal(button.disabled, true); assert.equal(status.dataset.state, 'loading');
  assert.equal(app.snapshots.joint, previous);
  finish.resolve(); await Promise.all([first, second]);
  assert.equal(status.dataset.state, 'ready'); assert.equal(button.disabled, false);
  assert.equal(app.snapshots.joint.result, previous.result);
  assert.equal(app.snapshots.joint.join_spec.palette, 'viridis');
  assert.equal(app.snapshots.joint.join_spec.type, 'contour');
  assert.deepEqual(loads, []);
});

test('None cancels pending visual reuse and cannot resurrect a discarded analysis', async t => {
  const { app, loads, status, button } = await seeded(t), previous = app.snapshots.joint.result;
  const entered = deferred(), finish = deferred(), originalPlot = app.plot;
  app.plot = async (...args) => { entered.resolve(); await finish.promise; await originalPlot(...args); };
  app.state.joint.palette = 'ocean';
  const visual = app.requestJointOnly(); await entered.promise;
  app.state.family2Id = ''; app.state.parameter2Id = '';
  const clear = app.requestJointOnly(); finish.resolve(); await Promise.all([visual, clear]);
  assert.equal(app.snapshots.joint, null); assert.equal(app.completedJointKey, null);
  assert.equal(button.disabled, true); assert.equal(status.dataset.state, 'ready');
  loads.length = 0; app.plot = originalPlot;
  app.state.family2Id = 'backbone'; app.state.parameter2Id = 'y';
  await app.requestJointOnly();
  assert(loads.length > 0, 'Restoring axes reused a discarded result');
  assert.notEqual(app.snapshots.joint.result, previous);
  assert.equal(button.disabled, false);
});

test('A matching joint key never bypasses incomplete full-panel repair', async t => {
  const { app } = await seeded(t), revisions = [];
  const pending = [];
  app.render = async request => { revisions.push(request.revision); const gate = deferred(); pending.push(gate); await gate.promise; };
  app.fullRenderComplete = false;
  const first = app.requestJointOnly();
  assert.equal(pending.length, 1, 'Failed full panels were bypassed by cached joint data');
  app.state.joint.palette = 'ocean';
  const second = app.requestJointOnly();
  assert.equal(pending.length, 2, 'An active full render was bypassed by cached joint data');
  pending[0].resolve(); await first;
  assert.equal(app.fullRenderComplete, false);
  assert.notEqual(app.fullRenderOwner, null);
  pending[1].resolve(); await second;
  assert.equal(app.fullRenderComplete, true); assert.equal(app.fullRenderOwner, null);
  assert(revisions[1] > revisions[0]);
});

test('Reused joint styling drains pending coordinate labels before reporting ready', async t => {
  const { app, loads, plots, status } = await seeded(t);
  app.completedCoordinateKey = 'retained-coordinate-population'; app.completedCoordinateLabels = 'all'; app.coordinateSummary = [];
  const entered = deferred(), finish = deferred(), originalPlot = app.plot;
  app.plot = async (...args) => {
    if (args[0] === app.$('jointPlot')) { entered.resolve(); await finish.promise; }
    await originalPlot(...args);
  };
  app.state.joint.palette = 'ocean';
  const pending = app.requestJointOnly(); await entered.promise;
  await app.setCoordinateLabels('none');
  finish.resolve(); await pending;
  assert.equal(status.dataset.state, 'ready'); assert.equal(app.completedCoordinateLabels, 'none');
  assert.equal(plots.find(plot => plot.target === app.$('coordinatePlot')).traces[0].mode, 'markers');
  assert.equal(app.snapshots.distribution.snapshot_id, 'unchanged-main');
  assert.deepEqual(loads, []);
});


test('A superseded failed style commit does not prevent the newest style from completing', async t => {
  const { app, loads, status, button } = await seeded(t), result = app.snapshots.joint.result;
  const entered = deferred(), finish = deferred(), originalPlot = app.plot;
  let count = 0;
  app.plot = async (...args) => {
    if (++count === 1) { entered.resolve(); await finish.promise; throw new Error('Superseded Plotly failure'); }
    await originalPlot(...args);
  };
  app.state.joint.palette = 'ocean';
  const first = app.requestJointOnly(); await entered.promise;
  app.state.joint.palette = 'viridis';
  const second = app.requestJointOnly(); finish.resolve(); await Promise.all([first, second]);
  assert.equal(status.dataset.state, 'ready'); assert.equal(button.disabled, false);
  assert.equal(app.snapshots.joint.result, result);
  assert.equal(app.snapshots.joint.join_spec.palette, 'viridis');
  assert.deepEqual(loads, []);
});

test('An axis change during pending visual work commits the new population', async t => {
  const { app, loads, status } = await seeded(t), result = app.snapshots.joint.result;
  const entered = deferred(), finish = deferred(), originalPlot = app.plot;
  let count = 0;
  app.plot = async (...args) => { if (++count === 1) { entered.resolve(); await finish.promise; } await originalPlot(...args); };
  app.state.joint.palette = 'ocean';
  const visual = app.requestJointOnly(); await entered.promise;
  app.state.family2Id = 'sugar';
  const changed = app.requestJointOnly(); finish.resolve(); await Promise.all([visual, changed]);
  assert.equal(status.dataset.state, 'ready'); assert(loads.length > 0);
  assert.notEqual(app.snapshots.joint.result, result);
  assert.deepEqual(app.snapshots.joint.result.points.map(point => point.y), [120, 160]);
  assert.equal(app.completedJointKey, app.jointAnalysisKey(app.state));
});
