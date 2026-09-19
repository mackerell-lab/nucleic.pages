import test from 'node:test';
import assert from 'node:assert/strict';
import { PureRnaExplorer } from '../app/PureRnaExplorer.js';
import { selectRows } from '../core/selection.js';
import { csv, createPlotSnapshot } from '../core/export.js';

const deferred = () => { let resolve; const promise = new Promise(done => { resolve = done; }); return { promise, resolve }; };
const metadata = { entries: [{ pdb_id: 'TEST' }] };
const residues = [
  { id: 'u', pdb_id: 'TEST', comp_id: 'U', values: { torsion: -5 } },
  { id: 'a', pdb_id: 'TEST', comp_id: 'A', values: { torsion: 5 } },
];
const pairs = [{ id: 'p', pdb_id: 'TEST', residue1_id: 'u', residue2_id: 'a',
  pair_label: 'U-A', values: { opening: 10 } }];
const relations = [
  { id: 'first', pair_id: 'p', residue_id: 'u', endpoint_role: 'first' },
  { id: 'second', pair_id: 'p', residue_id: 'a', endpoint_role: 'second' },
];

function setup(reverse = false) {
  const nodes = new Map(), controls = [], outputs = [];
  const node = id => {
    if (!nodes.has(id)) nodes.set(id, { disabled: false, textContent: '', dataset: {} });
    return nodes.get(id);
  };
  const app = new PureRnaExplorer({ root: { querySelector: node }, repository: {
    loadFamily: async id => ({ rows: id === 'pairs' ? pairs : residues }),
    loadRelations: async () => ({ rows: relations }),
  } });
  app.metadata = metadata; app.manifest = { build_id: 'science-audit', relations: { observations: {} } };
  app.parameter = (family, id) => ({ id, level: family === 'pairs' ? 'pair' : 'residue', period: family === 'pairs' ? null : 360 });
  app.updateJointResidueControls = table => { controls.push(table); };
  app.checkpoint = async revision => app.current(revision);
  app.renderJointResult = async (result, state) => {
    outputs.push(createPlotSnapshot({ result, selectionSpec: state.selection, displaySpec: state.display,
      joinSpec: state.joint, buildId: 'science-audit', snapshot_id: 'audit-355-snapshot' }));
  };
  app.state.familyId = reverse ? 'residues' : 'pairs'; app.state.parameterId = reverse ? 'torsion' : 'opening';
  app.state.family2Id = reverse ? 'pairs' : 'residues'; app.state.parameter2Id = reverse ? 'opening' : 'torsion';
  app.state.selection = { components: 'all', methods: [], contexts: [], includeEnds: true };
  app.state.joint.mode = 'relation'; app.state.joint.endpoint = 'both';
  app.state.display = { sigma: 0, normalization: 'probability', fine: false, circularMode: 'wrap_360' };
  const state = structuredClone(app.state);
  return { app, state, nodes, controls, outputs,
    left: selectRows(reverse ? residues : pairs, metadata, state.selection) };
}

test('Cancelled relation load cannot construct an old joint snapshot or read its rows', async () => {
  const f = setup(), entered = deferred(), finish = deferred();
  f.app.repository.loadRelations = async () => { entered.resolve(); return finish.promise; };
  const pending = f.app.renderJoint(f.state, 0, f.left);
  await entered.promise;
  const controlCount = f.controls.length;
  f.app.capture();
  const poison = [];
  Object.defineProperty(poison, '0', { get() { assert.fail('Stale relation rows reached the scientific join'); } });
  finish.resolve({ rows: poison });
  await pending;
  assert.equal(f.controls.length, controlCount);
  assert.deepEqual(f.outputs, []);
});

test('Reverse joint awaiting its endpoint control table cannot mutate newer controls', async () => {
  const f = setup(true), entered = deferred(), finish = deferred();
  let primaryLoads = 0;
  f.app.repository.loadFamily = async id => {
    if (id === 'residues' && ++primaryLoads === 2) { entered.resolve(); return finish.promise; }
    return { rows: id === 'pairs' ? pairs : residues };
  };
  const pending = f.app.renderJoint(f.state, 0, f.left);
  await entered.promise;
  const controlCount = f.controls.length;
  f.app.capture(); finish.resolve({ rows: residues });
  await pending;
  assert.equal(f.controls.length, controlCount, 'Old endpoint controls were rebuilt after a newer selection');
  assert.deepEqual(f.outputs, []);
});

test('Cancellation after joining or histogramming cannot publish or freeze a stale result', async () => {
  for (const cancelAt of [3, 4]) {
    const f = setup(); let checkpointCount = 0;
    f.app.checkpoint = async revision => {
      if (++checkpointCount === cancelAt) f.app.capture();
      return f.app.current(revision);
    };
    await f.app.renderJoint(f.state, 0, f.left);
    assert.equal(checkpointCount, cancelAt);
    assert.deepEqual(f.outputs, []);
    assert.equal(Object.isFrozen(f.left.rows[0]), false);
    assert.deepEqual(f.left.rows[0].values, { opening: 10 });
  }
});

test('A cancelled ranking load neither reads old scalar rows nor releases a newer owner', async () => {
  const f = setup(), entered = deferred(), finish = deferred();
  f.app.repository.loadSurveyScalars = async () => { entered.resolve(); return finish.promise; };
  f.app.surveyRows = () => assert.fail('Cancelled ranking must not normalize old scalar rows');
  const state = structuredClone(f.app.state), term = { id: 'torsion', label: 'Torsion', period: 360 };
  const pending = f.app.renderOpeningRanking([term], state, 0, { residues: new Map(), pairs: new Map() });
  await entered.promise;
  f.app.capture();
  const newerOwner = {}; f.app.rankingOwner = newerOwner;
  const button = f.app.$('surveyRankingLoad'); button.disabled = true; button.textContent = 'New ranking';
  finish.resolve({ rows: [] }); await pending;
  assert.equal(f.app.rankingOwner, newerOwner);
  assert.equal(button.disabled, true); assert.equal(button.textContent, 'New ranking');
  assert([...f.app.rankingCache.values()].every(cache => cache.size === 0));
});

test('Already stale ranking entry cannot steal ownership or change controls', async () => {
  const f = setup(); f.app.capture();
  const owner = {}; f.app.rankingOwner = owner;
  const button = f.app.$('surveyRankingLoad'); button.disabled = true; button.textContent = 'Current ranking';
  await f.app.renderOpeningRanking([{ id: 'torsion' }], f.state, 0, {});
  assert.equal(f.app.rankingOwner, owner);
  assert.equal(button.disabled, true); assert.equal(button.textContent, 'Current ranking');
});

test('Uncancelled checkpoints retain hand-derived relation populations and raw CSV angles', async () => {
  for (const reverse of [false, true]) {
    const f = setup(reverse);
    f.state.joint.residueContexts = ['U'];
    await f.app.renderJoint(f.state, 0, f.left);
    assert.equal(f.outputs.length, 1);
    const snapshot = f.outputs[0], result = snapshot.result;
    assert.deepEqual(result.points.map(p => [p.left_id, p.right_id, p.x, p.y, p.endpoint_role]),
      [reverse ? ['u', 'p', -5, 10, 'first'] : ['p', 'u', 10, -5, 'first']]);
    assert.equal(result.z.flat().reduce((sum, value) => sum + value, 0), 1);
    assert.equal(result.coverage.plottedPoints, 1);
    // This fixture contains no quoted fields. Inspect measurement columns,
    // not identity text: a valid snapshot ID can itself contain "355".
    const lines = csv(snapshot).trimEnd().split('\r\n');
    assert.equal(lines.length, 2);
    const headers = lines[0].split(','), fields = lines[1].split(',');
    assert.equal(fields[headers.indexOf('snapshot_id')], 'audit-355-snapshot');
    assert.deepEqual(['x_value', 'y_value'].map(key => fields[headers.indexOf(key)]),
      reverse ? ['-5', '10'] : ['10', '-5']);
    assert.deepEqual(snapshot.selection_spec, f.state.selection);
    assert.deepEqual(snapshot.join_spec.residueContexts, ['U']);
  }
});
