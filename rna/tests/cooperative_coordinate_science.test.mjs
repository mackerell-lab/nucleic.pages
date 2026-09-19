import test from 'node:test';
import assert from 'node:assert/strict';
import { PureRnaExplorer } from '../app/PureRnaExplorer.js';
import { RnaDataRepository } from '../core/repository.js';

const atom = (id, context, label, xyz, extra = {}) => ({ id, pdb_id: '1AAA', model_id: '1',
  residue_id: id, context, atom_label: label, x: xyz[0], y: xyz[1], z: xyz[2], status: 'ok', ...extra });
const baseChunks = [
  [atom('r1', 'U', 'anchor_U.O2', [0, 0, 0]), atom('r1', 'U', "anchor_U.O2'", [10, 0, 0]),
    atom('bad', 'U', 'anchor_U.O2', [900, 0, 0], { status: 'missing_atoms' })],
  [atom('r1', 'U', 'anchor_U.O2', [2, 4, 0], { model_id: '2' }),
    atom('r2', 'U', "anchor_U.O2'", [14, 0, 0]), atom('bad-null', 'U', 'anchor_U.O2', [900, 0, null])],
  [atom('r1', 'U', 'anchor_U.O2', [4, 8, 0], { pdb_id: '2BBB' }),
    atom('c1', 'C', 'anchor_C.O2', [99, 0, 0])],
];
const pairChunks = [
  [atom('c1', 'C-G', 'anchor_C.O2', [0, 0, 0], { pair_id: 'p1', opening_bin: 'small' }),
    atom('g1', 'C-G', 'paired_G.O6', [5, 0, 0], { pair_id: 'p1', opening_bin: 'small' })],
  [atom('c1', 'C-G', 'anchor_C.O2', [2, 0, 0], { pair_id: 'p2', opening_bin: 'middle' }),
    atom('c1', 'C-G', 'anchor_C.O2', [4, 0, 0], { pdb_id: '2BBB', pair_id: 'p1', opening_bin: 'small' })],
  [atom('c2', 'G-C', 'anchor_C.O2', [8, 0, 0], { pair_id: 'p3', opening_bin: 'small' }),
    atom('c3', 'C-G', 'anchor_C.O2', [99, 0, 0], { pair_id: 'excluded', opening_bin: 'small' })],
];

function setup(t, group, chunks) {
  const oldDocument = globalThis.document;
  const node = tag => ({ tag, children: [], value: '', textContent: '', dataset: {},
    setAttribute(key, value) { this[key] = value; }, append(...values) { this.children.push(...values); },
    replaceChildren(...values) { this.children = values; } });
  globalThis.document = { createElement: node };
  t.after(() => { globalThis.document = oldDocument; });
  const nodes = new Map(), plots = [];
  const app = new PureRnaExplorer({ root: { querySelector: id => {
    if (!nodes.has(id)) nodes.set(id, node(id)); return nodes.get(id);
  } }, repository: { async *iterateSurveyCoordinates() { for (const rows of chunks) yield { rows }; } } });
  app.metadata = { entries: [{ pdb_id: '1AAA' }, { pdb_id: '2BBB' }] };
  app.manifest = { build_id: 'coordinate-science', survey: { coordinates: { groups: { [group]: {} } } } };
  app.state.selection = { methods: [], components: 'all', includeEnds: true, contexts: [] };
  app.state.survey.coordinateGroup = group;
  app.renderCoordinatePlot = async averages => { plots.push(structuredClone(averages)); };
  app.openingIndex = async () => ({ pairs: new Map(['p1', 'p2', 'p3'].map(id => [id, {}])) });
  return { app, plots };
}

const find = (rows, context, atomLabel) => {
  const result = rows.find(row => row.context === context && row.atom_label === atomLabel);
  assert(result, `Missing ${context}/${atomLabel}`); return result;
};
function check(row, expected) {
  for (const field of ['n', 'entries', 'pairs', 'residues']) assert.equal(row[field], expected[field], field);
  assert.deepEqual(row.mean, expected.mean);
  assert(Math.abs(row.rms - expected.rms) < 1e-12, `RMS ${row.rms} != ${expected.rms}`);
}

test('Coordinate checkpoints preserve per-atom RNA identities and hand-derived moments across chunks', async t => {
  // O2' is synthetic identity coverage; current released coordinate groups are
  // base heavy atoms and C1', not an assertion of released O2' observations.
  const f = setup(t, 'rna_standard_base_test', baseChunks);
  const request = f.app.capture(); await f.app.renderCoordinates(request.state, request.revision);
  assert.equal(f.plots.length, 1); assert.equal(f.app.coordinateSummary.length, 3);
  check(find(f.app.coordinateSummary, 'U', 'anchor_U.O2'),
    { n: 3, entries: 2, pairs: null, residues: 3, mean: [2, 4, 0], rms: Math.sqrt(40 / 3) });
  check(find(f.app.coordinateSummary, 'U', "anchor_U.O2'"),
    { n: 2, entries: 1, pairs: null, residues: 2, mean: [12, 0, 0], rms: 2 });
  check(find(f.app.coordinateSummary, 'C', 'anchor_C.O2'),
    { n: 1, entries: 1, pairs: null, residues: 1, mean: [99, 0, 0], rms: 0 });
  assert.deepEqual(f.plots[0], f.app.coordinateSummary);
  assert.equal(baseChunks[0][0].x, 0);
  assert.equal(baseChunks[1][2].z, null);
});

test('Pair coordinate checkpoints retain incidence counts, missing atom subsets and opening/context filters', async t => {
  const f = setup(t, 'cytosine_standard_pair_test', pairChunks);
  let request = f.app.capture(); await f.app.renderCoordinates(request.state, request.revision);
  assert.equal(f.app.coordinateSummary.length, 3);
  check(find(f.app.coordinateSummary, 'C-G', 'anchor_C.O2'),
    { n: 3, entries: 2, pairs: 3, residues: 2, mean: [2, 0, 0], rms: Math.sqrt(8 / 3) });
  check(find(f.app.coordinateSummary, 'C-G', 'paired_G.O6'),
    { n: 1, entries: 1, pairs: 1, residues: 1, mean: [5, 0, 0], rms: 0 });
  check(find(f.app.coordinateSummary, 'G-C', 'anchor_C.O2'),
    { n: 1, entries: 1, pairs: 1, residues: 1, mean: [8, 0, 0], rms: 0 });
  f.app.state.survey.coordinateContext = 'C-G'; f.app.state.survey.coordinateOpening = 'small';
  request = f.app.capture(); await f.app.renderCoordinates(request.state, request.revision);
  assert.equal(f.app.coordinateSummary.length, 2);
  check(find(f.app.coordinateSummary, 'C-G', 'anchor_C.O2'),
    { n: 2, entries: 2, pairs: 2, residues: 2, mean: [2, 0, 0], rms: 2 });
});

test('Cancelled coordinate rendering closes the real repository iterator without retaining partitions', async t => {
  const group = 'rna_standard_base_test', f = setup(t, group, baseChunks), calls = [];
  const manifest = { schema_version: 'rna-explorer-1', molecule_type: 'RNA', build_id: 'coordinate-science',
    survey: { coordinates: { groups: { [group]: { partitions: baseChunks.map((rows, i) => ({
      path: `part${i}.json`, row_count: rows.length, entry_ids: ['1AAA', '2BBB'],
    })) } } } } };
  const repository = new RnaDataRepository({ manifestUrl: 'https://coordinates.test/manifest.json', fetchImpl: async url => {
    const filename = new URL(url).pathname.slice(1); calls.push(filename);
    return new Response(JSON.stringify(filename === 'manifest.json' ? manifest
      : { build_id: 'coordinate-science', rows: baseChunks[Number(filename.match(/part(\d)/)[1])] }));
  } });
  const original = repository.iterateSurveyCoordinates.bind(repository); let closed = 0;
  repository.iterateSurveyCoordinates = async function* (...args) {
    try { yield* original(...args); } finally { closed++; }
  };
  f.app.repository = repository;
  const previousSummary = [{ previous: true }]; f.app.coordinateSummary = previousSummary;
  f.app.lastCoordinateGroup = 'previous';
  const releases = []; repository.releaseSurvey = (...args) => releases.push(args);
  const request = f.app.capture(); let checkpoints = 0;
  f.app.checkpoint = async () => { checkpoints++; f.app.capture(); return false; };
  await f.app.renderCoordinates(request.state, request.revision);
  assert.equal(checkpoints, 1); assert.equal(closed, 1);
  assert.deepEqual(calls, ['manifest.json', 'part0.json']);
  assert.equal(repository.coordinateRequests.size, 0);
  assert.deepEqual([...repository.promises.keys()], ['manifest']);
  assert.deepEqual(f.plots, []); assert.equal(f.app.coordinateSummary, previousSummary);
  assert.equal(f.app.lastCoordinateGroup, 'previous'); assert.deepEqual(releases, []);
  f.app.checkpoint = async revision => f.app.current(revision);
  const retry = f.app.capture(); await f.app.renderCoordinates(retry.state, retry.revision);
  assert.equal(closed, 2);
  assert.deepEqual(calls, ['manifest.json', 'part0.json', 'part0.json', 'part1.json', 'part2.json']);
  assert.equal(repository.coordinateRequests.size, 0);
  check(find(f.app.coordinateSummary, 'U', 'anchor_U.O2'),
    { n: 3, entries: 2, pairs: null, residues: 3, mean: [2, 4, 0], rms: Math.sqrt(40 / 3) });
});
