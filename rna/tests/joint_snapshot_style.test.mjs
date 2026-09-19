import test from 'node:test';
import assert from 'node:assert/strict';
import { createPlotSnapshot, restyleJointSnapshot, csv } from '../core/export.js';

const result = () => ({ kind: 'joint', xParameter: { id: 'x' }, yParameter: { id: 'y' },
  points: [{ left: { id: 'a', pdb_id: '1AAA' }, right: { id: 'a', pdb_id: '1AAA' }, x: -0, y: 1.2345678901234567 }],
  z: [[0, 1]], statistics: { r: null, r2: null } });

test('Joint restyling shares only an already owned frozen graph and updates provenance', () => {
  const old = createPlotSnapshot({ result: result(), buildId: 'release', snapshot_id: 'original',
    joinSpec: { mode: 'identity', endpoint: 'both', palette: 'hotspots', labels: false },
    selectionSpec: { contexts: ['U'] }, displaySpec: { bins: 72 }, provenance: { axis_selections: { x: { contexts: ['U'] } } } });
  const next = restyleJointSnapshot(old, { mode: 'relation', endpoint: 'first', palette: 'viridis', labels: true });
  assert.notEqual(next.snapshot_id, old.snapshot_id);
  assert.equal(next.result, old.result);
  assert.equal(next.provenance, old.provenance);
  assert.equal(next.selection_spec, old.selection_spec);
  assert.equal(next.display_spec, old.display_spec);
  // Only presentation keys can change through this API, even with bad callers.
  assert.deepEqual(next.join_spec, { mode: 'identity', endpoint: 'both', palette: 'viridis', labels: true });
  assert.equal(old.join_spec.palette, 'hotspots');
  assert.equal(old.join_spec.labels, false);
  assert.throws(() => { next.result.points[0].y = 0; }, TypeError);
  assert.throws(() => { next.join_spec.palette = 'jet'; }, TypeError);
  assert.throws(() => { next.provenance.axis_selections.x.contexts.push('A'); }, TypeError);
  assert.equal(csv(next).replaceAll(next.snapshot_id, old.snapshot_id), csv(old));
  const again = restyleJointSnapshot(next, { palette: 'hotspots', labels: false });
  assert.equal(again.result, old.result);
  assert.notEqual(again.snapshot_id, next.snapshot_id);
});

test('Shallow-frozen lookalikes and nonjoint snapshots cannot bypass graph freezing', () => {
  const fake = Object.freeze({ result: result(), join_spec: {} });
  assert.throws(() => restyleJointSnapshot(fake, { palette: 'viridis' }), /owned frozen joint snapshot/);
  const linear = createPlotSnapshot({ result: { kind: 'distribution', series: [] } });
  assert.throws(() => restyleJointSnapshot(linear, {}), /owned frozen joint snapshot/);
  assert.throws(() => restyleJointSnapshot(structuredClone(createPlotSnapshot({ result: result() })), {}), /owned frozen joint snapshot/);
});

test('Repeated fast style snapshots have distinct identities', () => {
  let snapshot = createPlotSnapshot({ result: result() });
  const ids = new Set([snapshot.snapshot_id]);
  for (let i = 0; i < 100; i++) {
    snapshot = restyleJointSnapshot(snapshot, { labels: Boolean(i % 2) });
    assert(!ids.has(snapshot.snapshot_id));
    ids.add(snapshot.snapshot_id);
  }
});
