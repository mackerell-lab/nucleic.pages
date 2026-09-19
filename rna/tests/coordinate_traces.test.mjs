import test from 'node:test';
import assert from 'node:assert/strict';
import { coordinateTraces } from '../views/coordinate-traces.js';

const averages = [
  { atom: 'G-C anchor_C.N1', context: 'G-C', atom_label: 'anchor_C.N1', mean: [0, -0, 1.123456789], n: 7, residues: 5, pairs: 6, entries: 3, rms: 0.28 },
  { atom: "U anchor_U.O2'", context: 'U', atom_label: "anchor_U.O2'", mean: [1, 2, 3], n: 11, residues: null, pairs: null, entries: 2, rms: 0.19 },
];

test('Hiding coordinate labels preserves every mean, identity and population count', () => {
  const before = structuredClone(averages);
  const labelled = coordinateTraces(averages)[0], markers = coordinateTraces(averages, 'none')[0];
  assert.equal(labelled.mode, 'markers+text'); assert.equal(markers.mode, 'markers');
  assert.deepEqual({ ...labelled, mode: 'markers' }, markers);
  assert.deepEqual(markers.x, [0, 1]); assert.deepEqual(markers.y, [-0, 2]);
  assert.deepEqual(markers.z, [1.123456789, 3]);
  assert.deepEqual(markers.customdata[0], ['G-C', 'anchor_C.N1', 7, 5, 6, 3, 0.28]);
  assert.deepEqual(markers.customdata[1], ['U', "anchor_U.O2'", 11, 'Unavailable', 'Not applicable', 2, 0.19]);
  assert.deepEqual(averages, before);
});

test('Coordinate text escapes markup and empty plots retain a valid trace', () => {
  const trace = coordinateTraces([{ ...averages[0], context: '<G&>', atom_label: '<N1>', atom: '<G&> <N1>' }])[0];
  assert.equal(trace.text[0], '&lt;G&amp;&gt; &lt;N1&gt;');
  assert.deepEqual(trace.customdata[0].slice(0, 2), ['&lt;G&amp;&gt;', '&lt;N1&gt;']);
  assert.deepEqual(coordinateTraces([])[0].x, []);
});
