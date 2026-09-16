import test from 'node:test';
import assert from 'node:assert/strict';
import { gzipSync } from 'node:zlib';
import { RnaDataRepository } from '../core/repository.js';
import { selectRows } from '../core/selection.js';
import { distribution, histogram2D } from '../core/analysis.js';
import { join } from '../core/joints.js';
import { createPlotSnapshot, csv, provenance, escapeCsv } from '../core/export.js';
import { normalizeParameter } from '../core/registry.js';
import { smoothCounts, summary, correlation, wrapCircular } from '../math/numeric.js';
import { computeResidueObservables } from '../offline/residue_geometry.mjs';

const near = (actual, expected, tolerance = 1e-10) => assert.ok(Math.abs(actual - expected) < tolerance, `${actual} != ${expected}`);
const row = (id, value, extra = {}) => ({ id, pdb_id: '1ABC', entity_id: '1', comp_id: 'U', values: { chi: value, length: value }, statuses: {}, ...extra });
const torsion = { id: 'chi', label: 'Chi', unit: 'degrees', period: 360 };
const linear = { id: 'length', unit: 'angstrom', period: null };

test('Circular moments respect seams, arbitrary declared periods and undefined directions', () => {
  near(wrapCircular(-721), 359);
  near(summary([359, 1], { period: 360 }).mean, 0);
  near(summary([179, 1], { period: 180 }).mean, 0);
  assert.equal(summary([0, 90, 180, 270], { period: 360 }).mean, null);
  assert.equal(summary([], { period: 360 }).meanStatus, 'empty');
  near(summary([12, 12], { period: 360 }).std, 0);
  assert.equal(normalizeParameter({ id: 'angle', unit: 'degrees' }).isCircular, false);
});

test('Smoothing always returns floating values and supports sigma zero and tiny circular arrays', () => {
  const raw = new Uint32Array([0, 1, 0]);
  const output = smoothCounts(raw, 1.2, true);
  assert.ok(output instanceof Float64Array);
  assert.ok(output[0] > 0 && output[0] < 1);
  near(output.reduce((a, b) => a + b), 1);
  assert.deepEqual(Array.from(smoothCounts(raw, 0)), [0, 1, 0]);
  near(smoothCounts([1], 2, true)[0], 1);
  assert.deepEqual(Array.from(raw), [0, 1, 0]);
});

test('Stable raw linear moments and type-7 quantiles', () => {
  const stats = summary([1e12 + 1, 1e12 + 2, 1e12 + 3]);
  near(stats.mean, 1e12 + 2);
  near(stats.std, Math.sqrt(2 / 3));
  near(summary([0, 10]).p05, 0.5);
  assert.equal(summary([]).mean, null);
});

test('Linear and circular probability and density integrate to one after smoothing', () => {
  for (const parameter of [linear, torsion]) {
    for (const normalization of ['probability', 'density']) {
      const result = distribution([row('a', 0), row('b', 1), row('c', 359)], parameter,
        { groupBy: 'none', bins: 8, sigma: 1.2, normalization });
      const series = result.series[0];
      near(series.y.reduce((a, b) => a + b) * (normalization === 'density' ? series.binWidth : 1), 1);
      assert.equal(series.statistics.n, 3);
      near(series.counts.reduce((a, b) => a + b), 3);
    }
  }
  const empty = distribution([row('missing', null)], torsion);
  assert.equal(empty.series.length, 0);
  assert.equal(empty.coverage.unavailableRows, 1);
});

test('One frozen result preserves raw circular values and CSV against changed controls and rows', () => {
  const rows = [row('a', -181, { auth_asym_id: 'a\rb', label_seq_id: '8' }), row('b', 359)];
  const controls = { contexts: ['U'], includeEnds: true };
  const result = distribution(rows, torsion, { groupBy: 'none', circularMode: 'signed_180', sigma: 0 });
  const snapshot = createPlotSnapshot({ result, selectionSpec: controls, buildId: 'release-a' });
  const before = csv(snapshot);
  rows[0].values.chi = 12;
  result.series[0].values[0] = 44;
  controls.contexts.push('A');
  assert.equal(csv(snapshot), before);
  assert.match(before, /-181/);
  assert.match(before, /"a\rb"/);
  assert.deepEqual(snapshot.selection_spec.contexts, ['U']);
  assert.throws(() => { snapshot.result.series[0].values[0] = 7; }, TypeError);
  assert.equal(JSON.parse(provenance(snapshot)).build_id, 'release-a');
  assert.equal(escapeCsv('a,"b"\nc'), '"a,""b""\nc"');
});

test('Scoped RNA annotations do not propagate entry union to a different entity', () => {
  const metadata = { entries: [{ pdb_id: '1ABC', method: 'X-RAY DIFFRACTION', resolution: 2.5,
    profiles: { relaxed: true }, functions: ['riboswitch', 'rRNA'] }],
    entities: [{ pdb_id: '1ABC', entity_id: '1', functions: ['riboswitch'] }, { pdb_id: '1ABC', entity_id: '2', functions: ['rRNA'] }] };
  const rows = [row('one', 1), row('two', 2, { entity_id: '2' }), row('unknown', 3, { entity_id: '3' })];
  const selected = selectRows(rows, metadata, { components: 'relaxed', methods: ['xray'], resolutionMax: 3, functions: ['riboswitch'] });
  assert.deepEqual(selected.rows.map(item => item.id), ['one']);
  assert.equal(selected.coverage.eligibleEntries, 1);
  assert.deepEqual(selected.indices, [0]);
  assert.equal(selectRows(rows, metadata, { components: 'conservative' }).rows.length, 0);
});

test('NMR selection does not require inapplicable diffraction resolution', () => {
  const metadata = { entries: [{ pdb_id: '1ABC', method: 'SOLUTION NMR', profiles: ['relaxed'] }] };
  assert.equal(selectRows([row('a', 1)], metadata, { methods: ['nmr'], resolutionMax: 3, components: 'relaxed' }).rows.length, 1);
});

test('Actual RNA geometry output preserves U torsions, availability and terminal filtering', () => {
  const entry = { pdb_id: '1ABC', links: [], residues: [{ id: 'u1', comp_id: 'U', entity_id: '1',
    label_asym_id: 'A', label_seq_id: 1, atoms: { "O4'": [1, 0, 0], "C1'": [0, 0, 0], N1: [0, 1, 0], C2: [0, 1, 1] } }] };
  const computed = computeResidueObservables(entry);
  const metadata = { entries: [{ pdb_id: '1ABC', method: 'X-RAY DIFFRACTION', resolution: 2, profiles: { relaxed: true } }], entities: [] };
  const selected = selectRows(computed, metadata, { includeEnds: true, contexts: ['U'], functions: ['unknown'] });
  assert.equal(selected.rows.length, 1);
  const result = distribution(selected.rows, torsion, { groupBy: 'base', sigma: 0 });
  assert.equal(result.series[0].key, 'U');
  near(result.series[0].values[0], -90);
  assert.equal(selectRows(computed, metadata, { includeEnds: false }).rows.length, 0);
  assert.equal(distribution(selected.rows, { id: 'alpha' }).series.length, 0);
});

test('Pair and step class filters require all explicitly scoped endpoint entities', () => {
  const metadata = { entries: [{ pdb_id: '1ABC' }], entities: [
    { pdb_id: '1ABC', entity_id: '1', functions: ['riboswitch'] }, { pdb_id: '1ABC', entity_id: '2', functions: ['rRNA'] },
  ] };
  const pair = { id: 'pair1', pdb_id: '1ABC', residue_ids: ['r1', 'r2'], pair_label: 'G-U',
    endpoint_entities: [{ entity_id: '1' }, { entity_id: '2' }], values: { opening: 10 }, statuses: { opening: 'available' } };
  assert.equal(selectRows([pair], metadata, { functions: ['riboswitch'] }).rows.length, 0);
  assert.equal(selectRows([pair], metadata, { functions: ['riboswitch'], annotationEndpointPolicy: 'any' }).rows.length, 1);
  const all = selectRows([pair], metadata, { contexts: ['G-U'] });
  assert.equal(all.rows.length, 1);
  assert.equal(distribution(all.rows, { id: 'opening' }, { groupBy: 'function' }).series[0].key, 'Mixed / incomplete entity annotations');
});

test('Multi-tag groups retain overlapping membership counts; per-entry weights balance contributions', () => {
  const rows = [row('a', 0, { functions: ['A', 'B'] }), row('b', 10, { functions: ['A'] }), row('c', 20, { pdb_id: '2ABC', functions: ['B'] })];
  const grouped = distribution(rows, linear, { groupBy: 'function', sigma: 0 });
  assert.equal(grouped.coverage.plottedRows, 3);
  assert.equal(grouped.coverage.memberships, 4);
  const balanced = distribution(rows, linear, { groupBy: 'none', weighting: 'entry_equal', sigma: 0 });
  near(balanced.series[0].statistics.mean, 12.5);
  near(balanced.series[0].weights.reduce((a, b) => a + b), 2);
});

test('Same-level joins use stable identities across independently filtered context vocabularies', () => {
  const left = [row('step1', 1, { context: 'CG' }), row('step2', 2, { context: 'AU' })];
  const right = [row('step1', 3, { context: 'A-like' })];
  const result = join(left, right, { type: 'identity' });
  assert.equal(result.points.length, 1);
  assert.equal(result.points[0].left_id, 'step1');
  assert.equal(result.points[0].rightIndex, 0);
  assert.throws(() => join([row('x', 1), row('x', 2)], right), /duplicate/);
  assert.throws(() => join({ rows: left, build_id: 'a' }, { rows: right, build_id: 'b' }), /different RNA builds/);
});

test('Explicit endpoint joins preserve incidences and count equal-pair weights', () => {
  const pairs = [row('p1', 2)], residues = [row('r1', 4), row('r2', 6)];
  const relations = [
    { id: 'edge1', left_id: 'p1', right_id: 'r1', pair_id: 'p1', residue_id: 'r1', endpoint_role: 'first' },
    { id: 'edge2', left_id: 'p1', right_id: 'r2', pair_id: 'p1', residue_id: 'r2', endpoint_role: 'second' },
  ];
  const result = join(pairs, residues, { type: 'pair_residue', relations, endpoint: 'both', weighting: 'pair_equal' });
  assert.equal(result.diagnostics.emittedIncidences, 2);
  assert.equal(result.diagnostics.uniquePairs, 1);
  assert.equal(result.diagnostics.uniqueResidues, 2);
  near(result.points.reduce((a, b) => a + b.weight, 0), 1);
  assert.deepEqual(result.points.map(point => point.endpoint_role), ['first', 'second']);
});

test('Pipeline multiplexed relation schema joins typed endpoints in either axis direction', () => {
  const pairs = [row('p1', 2)], residues = [row('r1', 4), row('r2', 6)];
  const relations = [
    { id: 'p1/residue/1', kind: 'pair_residue', pair_id: 'p1', residue_id: 'r1', side: 1 },
    { id: 'p1/residue/2', kind: 'pair_residue', pair_id: 'p1', residue_id: 'r2', side: 2 },
    { id: 'step1/pair/1', kind: 'step_pair', step_id: 'step1', pair_id: 'p1', side: 1 },
  ];
  const forward = join(pairs, residues, { type: 'relation', relations, xParameter: { level: 'pair' }, yParameter: { level: 'residue' } });
  assert.equal(forward.points.length, 2);
  assert.deepEqual(forward.points.map(point => point.endpoint_role), ['first', 'second']);
  const reverse = join(residues, pairs, { type: 'relation', relations, xParameter: { level: 'residue' }, yParameter: { level: 'pair' } });
  assert.equal(reverse.points.length, 2);
  assert.equal(reverse.points[0].left_id, 'r1');
});

test('2D smoothed density has unit mass and circular correlation is seam invariant', () => {
  const input = [{ x: 355, y: 5 }, { x: 0, y: 10 }, { x: 5, y: 15 }];
  const result = histogram2D(input, torsion, torsion, { bins: 12, normalization: 'density', sigma: 1.2 });
  near(result.z.flat().reduce((a, b) => a + b) * result.binArea, 1);
  near(result.statistics.r, 1);
  assert.equal(result.statistics.r2, null);
  near(correlation([1, 2, 3], [3, 5, 7]).r, 1);
  assert.equal(correlation([1, 1, 1], [1, 2, 3]).status, 'constant_or_degenerate');
  assert.equal(correlation([1, 2, 3], [1, 2, 3], 360, null).status, 'circular_linear_not_supported');
});

test('Repository pins relative release URLs, deduplicates in-flight requests and retries rejection', async () => {
  const calls = new Map();
  let fail = true;
  const fetchImpl = async url => {
    calls.set(url, (calls.get(url) || 0) + 1);
    if (url.endsWith('/manifest.json')) return new Response(JSON.stringify({ build_id: 'test', manifest: 'releases/test/release.json' }));
    if (url.endsWith('/release.json')) return new Response(JSON.stringify({ schema_version: 'rna-explorer-1', molecule_type: 'RNA', build_id: 'test',
      families: [{ id: 'backbone', path: 'data.json.gz', row_count: 1 }] }));
    if (fail) { fail = false; return new Response('Temporary error', { status: 503 }); }
    return new Response(gzipSync(JSON.stringify({ build_id: 'test', rows: [row('a', 1)] })));
  };
  const repository = new RnaDataRepository({ manifestUrl: 'https://example.org/nucleic.pages/assets/pure_rna/manifest.json', fetchImpl });
  const first = repository.loadFamily('backbone'), second = repository.loadFamily('backbone');
  assert.equal(first, second);
  await assert.rejects(first, /503/);
  const table = await repository.loadFamily('backbone');
  assert.equal(table.rows.length, 1);
  assert.equal(calls.get('https://example.org/nucleic.pages/assets/pure_rna/releases/test/data.json.gz'), 2);
  assert.equal(calls.get('https://example.org/nucleic.pages/assets/pure_rna/manifest.json'), 1);
  assert.throws(() => { table.rows[0].values.chi = 3; }, TypeError);
});

test('Repository rejects incompatible schema, duplicate row IDs and cross-build assets', async () => {
  for (const fault of ['schema', 'duplicate', 'build']) {
    const repository = new RnaDataRepository({ manifestUrl: 'https://example.org/manifest.json', fetchImpl: async url => {
      if (url.endsWith('manifest.json')) return new Response(JSON.stringify({ schema_version: fault === 'schema' ? 'other' : 'rna-explorer-1', molecule_type: 'RNA', build_id: 'a', families: [{ id: 'test', path: 'rows.json' }] }));
      return new Response(JSON.stringify({ build_id: fault === 'build' ? 'b' : 'a', rows: [row('x', 1), row('x', 2)] }));
    } });
    await assert.rejects(repository.loadFamily('test'), fault === 'schema' ? /Unsupported/ : fault === 'build' ? /Cross-build/ : /duplicate/);
  }
});

test('Survey term loads stay separate from coordinates and can release cached partitions', async () => {
  const calls = [];
  const repository = new RnaDataRepository({ manifestUrl: 'https://example.org/manifest.json', fetchImpl: async url => {
    calls.push(url);
    if (url.endsWith('manifest.json')) return new Response(JSON.stringify({ schema_version: 'rna-explorer-1', molecule_type: 'RNA', build_id: 'a', families: [],
      survey: { scalars: { terms: { 'U.angle': { path: 'u.json', row_count: 1 } } }, coordinates: { groups: { U: { path: 'coords.json' } } } } }));
    return new Response(JSON.stringify({ build_id: 'a', rows: [row('x', 1)] }));
  } });
  assert.equal((await repository.loadSurveyScalars()).terms[0].id, 'U.angle');
  assert.equal(calls.length, 1);
  await repository.loadSurveyScalars('U.angle');
  await repository.loadSurveyScalars('U.angle');
  assert.equal(calls.length, 2);
  assert.ok(!calls.some(url => url.includes('coords')));
  repository.releaseSurvey('scalars', 'U.angle');
  await repository.loadSurveyScalars('U.angle');
  assert.equal(calls.length, 3);
  await repository.loadSurveyCoordinates('U');
  assert.equal(calls.length, 4);
});

test('Repository fetch preserves the browser global receiver', async () => {
  const repository = new RnaDataRepository({ manifestUrl: 'https://example.org/manifest.json', fetchImpl: function () {
    assert.equal(this, globalThis);
    return Promise.resolve(new Response(JSON.stringify({ schema_version: 'rna-explorer-1', molecule_type: 'RNA', families: [] })));
  } });
  assert.equal((await repository.loadManifest()).molecule_type, 'RNA');
});

function coordinateRepository() {
  const calls = [];
  const assets = {
    'manifest.json': { schema_version: 'rna-explorer-1', molecule_type: 'RNA', build_id: 'a', families: [],
      survey: { coordinates: { groups: { G: { label: 'Guanine', row_count: 3, partitions: [
        { path: 'ab.json', row_count: 2, entry_ids: ['1AAA', '1BBB'] },
        { path: 'c.json', row_count: 1, entry_ids: ['1CCC'] },
      ] } } } } },
    'ab.json': { build_id: 'a', rows: [row('a1', 1, { pdb_id: '1AAA' }), row('b1', 2, { pdb_id: '1BBB' })] },
    'c.json': { build_id: 'a', rows: [row('c1', 3, { pdb_id: '1CCC' })] },
  };
  const repository = new RnaDataRepository({ manifestUrl: 'https://example.org/manifest.json', fetchImpl: async url => {
    calls.push(url);
    return new Response(JSON.stringify(assets[url.split('/').at(-1)]));
  } });
  return { repository, calls, assets };
}

test('Coordinate streaming prunes entry partitions before fetch and retains no parsed chunks', async () => {
  const { repository, calls } = coordinateRepository();
  const received = [];
  for await (const chunk of repository.iterateSurveyCoordinates('G', { entryIds: ['1bbb'] })) {
    received.push(...chunk.rows.map(item => item.id));
    assert.equal(chunk.source_row_count, 2);
    assert.equal(chunk.selected_row_count, 1);
    assert.throws(() => { chunk.rows.push(row('bad', 1)); }, TypeError);
  }
  assert.deepEqual(received, ['b1']);
  assert.ok(!calls.some(url => url.endsWith('c.json')));
  assert.equal(repository.coordinateRequests.size, 0);
  assert.deepEqual([...repository.promises.keys()], ['manifest']);
  const beforeEmptySelection = calls.length;
  for await (const chunk of repository.iterateSurveyCoordinates('G', { entryIds: [] })) assert.fail('An empty selection must emit no chunks');
  assert.equal(calls.length, beforeEmptySelection);
});

test('Coordinate iteration can stop before fetching later chunks and honors cancellation', async () => {
  const { repository, calls } = coordinateRepository();
  for await (const chunk of repository.iterateSurveyCoordinates('G')) {
    assert.equal(chunk.rows.length, 2);
    break;
  }
  assert.equal(calls.filter(url => !url.endsWith('manifest.json')).length, 1);
  const controller = new AbortController();
  controller.abort(new Error('superseded'));
  await assert.rejects(async () => {
    for await (const chunk of repository.iterateSurveyCoordinates('G', { signal: controller.signal })) assert.fail('Aborted iteration cannot emit chunks');
  }, /superseded/);
  assert.equal(calls.length, 2);
});

test('Coordinate collection is bounded, uncached by default and retains at most one opt-in selection', async () => {
  const { repository, calls } = coordinateRepository();
  await assert.rejects(repository.loadSurveyCoordinates('G', { maxRows: 2 }), /use iterateSurveyCoordinates/);
  assert.equal(repository.coordinateRequests.size, 0);
  const selected = await repository.loadSurveyCoordinates('G', { entryIds: ['1AAA'], maxRows: 1 });
  assert.deepEqual(selected.rows.map(item => item.id), ['a1']);
  const previous = calls.length;
  await repository.loadSurveyCoordinates('G', { entryIds: ['1AAA'], maxRows: 1 });
  assert.equal(calls.length, previous + 1);
  await repository.loadSurveyCoordinates('G', { entryIds: ['1AAA'], maxRows: 1, retain: true });
  await repository.loadSurveyCoordinates('G', { entryIds: ['1BBB'], maxRows: 1, retain: true });
  assert.equal([...repository.promises.keys()].filter(key => key.startsWith('survey:coordinates:collected:')).length, 1);
  repository.releaseSurvey('coordinates', 'G');
  assert.deepEqual([...repository.promises.keys()], ['manifest']);
});

test('Coordinate partitions validate build identity and serialized row counts', async () => {
  for (const invalid of ['build', 'count']) {
    const { repository, assets } = coordinateRepository();
    if (invalid === 'build') assets['ab.json'].build_id = 'other';
    else assets['ab.json'].rows.pop();
    await assert.rejects(async () => {
      for await (const chunk of repository.iterateSurveyCoordinates('G')) assert.fail('Invalid partition cannot be delivered');
    }, invalid === 'build' ? /Cross-build/ : /row count mismatch/);
    assert.equal(repository.coordinateRequests.size, 0);
  }
});

test('Opening-conditioned CSV freezes source and endpoint identities with exact values', () => {
  const source = row('site1|pair|pair1', 12.125, { source_observation_id: 'site1', residue_id: 'residue1',
    pair_id: 'pair1', endpoint_role: 'first', opening: -2.375, opening_bin: '[-5, 0)' });
  const result = distribution([source], torsion, { groupBy: 'none', sigma: 0 });
  const snapshot = createPlotSnapshot({ result, buildId: 'opening-test', selectionSpec: { openingBin: '[-5, 0)' },
    dataHashes: { coordinates: 'hash-a' }, provenance: { annotation_endpoint_policy: 'all' } });
  const output = csv(snapshot), header = output.split('\r\n')[0].split(',');
  for (const column of ['source_observation_id', 'residue_id', 'pair_id', 'endpoint_role', 'opening', 'opening_bin']) assert.ok(header.includes(column));
  assert.match(output, /site1,residue1,pair1,first,-2\.375,"\[-5, 0\)"/);
  const record = snapshot.result.series[0].rows[0];
  source.opening = 99; source.opening_bin = 'changed'; source.pair_id = 'changed';
  assert.equal(csv(snapshot), output);
  assert.equal(record.pair_id, 'pair1');
  assert.equal(record.opening, -2.375);
  assert.equal(JSON.parse(provenance(snapshot)).data_hashes.coordinates, 'hash-a');
});

test('Family LRU cache retains three recent resolved tables and reloads evicted data', async () => {
  const calls = new Map();
  const repository = new RnaDataRepository({ manifestUrl: 'https://example.org/manifest.json', fetchImpl: async url => {
    const file = url.split('/').at(-1);
    calls.set(file, (calls.get(file) || 0) + 1);
    return new Response(JSON.stringify(file === 'manifest.json'
      ? { schema_version: 'rna-explorer-1', molecule_type: 'RNA', build_id: 'a', families: ['a', 'b', 'c', 'd'].map(id => ({ id, path: `${id}.json`, row_count: 1 })) }
      : { build_id: 'a', rows: [row(file, 10)] }));
  } });
  const retainedByCaller = await repository.loadFamily('b');
  await repository.loadFamily('a');
  await repository.loadFamily('c');
  await repository.loadFamily('a');
  await repository.loadFamily('d');
  assert.deepEqual([...repository.resolvedFamilies.keys()], ['family:c', 'family:a', 'family:d']);
  assert.ok(!repository.promises.has('family:b'));
  assert.equal(retainedByCaller.rows[0].values.chi, 10);
  assert.throws(() => { retainedByCaller.rows[0].values.chi = 20; }, TypeError);
  await repository.loadFamily('b');
  assert.equal(calls.get('b.json'), 2);
  assert.equal(calls.get('a.json'), 1);
  assert.equal(repository.resolvedFamilies.size, 3);
  assert.equal([...repository.promises.keys()].filter(key => key.startsWith('family:')).length, 3);
});

test('Family LRU never evicts an in-flight request and preserves concurrent request deduplication', async () => {
  let releaseSlow, startedSlow;
  const waiting = new Promise(resolve => { releaseSlow = resolve; });
  const started = new Promise(resolve => { startedSlow = resolve; });
  const calls = new Map();
  const repository = new RnaDataRepository({ manifestUrl: 'https://example.org/manifest.json', maxCachedFamilies: 1, fetchImpl: async url => {
    const file = url.split('/').at(-1);
    calls.set(file, (calls.get(file) || 0) + 1);
    if (file === 'manifest.json') return new Response(JSON.stringify({ schema_version: 'rna-explorer-1', molecule_type: 'RNA', build_id: 'a', families: ['slow', 'a', 'b'].map(id => ({ id, path: `${id}.json` })) }));
    if (file === 'slow.json') { startedSlow(); await waiting; }
    return new Response(JSON.stringify({ build_id: 'a', rows: [row(file, 10)] }));
  } });
  const slow = repository.loadFamily('slow');
  await started;
  await repository.loadFamily('a');
  await repository.loadFamily('b');
  assert.equal(repository.loadFamily('slow'), slow);
  assert.equal(repository.resolvedFamilies.size, 1);
  assert.ok(repository.promises.has('family:slow'));
  releaseSlow();
  await slow;
  assert.equal(calls.get('slow.json'), 1);
  assert.deepEqual([...repository.resolvedFamilies.keys()], ['family:slow']);
  assert.equal([...repository.promises.keys()].filter(key => key.startsWith('family:')).length, 1);
});

test('RNA pair policy defaults to exact and never filters residue or step observations', () => {
  const pair = (id, extra = {}) => ({ id, pdb_id: '1ABC', residue1_id: 'r1', residue2_id: 'r2', family: 'cWW', pair_label: 'G-U',
    near: false, alternative: false, stem_eligible: true, values: { opening: 1 }, ...extra });
  const rows = [pair('exact'), pair('near-flag', { near: true }), pair('near-name', { family: 'ncWW', near: undefined }),
    pair('alternative-flag', { alternative: true }), pair('alternative-name', { family: 'cWWa', alternative: undefined }),
    row('unpaired-residue', 3), { id: 'step', pdb_id: '1ABC', pair1_id: 'p1', pair2_id: 'p2', residue_ids: ['r1', 'r2', 'r3', 'r4'], values: { shift: 2 } }];
  const selected = selectRows(rows);
  assert.deepEqual(selected.rows.map(item => item.id), ['exact', 'unpaired-residue', 'step']);
  assert.equal(selected.spec.pairPolicy, 'exact');
  assert.equal(selected.coverage.pairPolicyExcluded, 4);
  assert.deepEqual(selectRows(rows, {}, { pairPolicy: 'exact' }).indices, selected.indices);
  assert.equal(selectRows(rows, {}, { pairPolicy: 'all' }).rows.length, rows.length);
  assert.deepEqual(selectRows(rows, {}, { pairPolicy: 'near' }).rows.map(item => item.id), ['near-flag', 'near-name', 'unpaired-residue', 'step']);
});

test('RNA pair family and stem filters preserve noncanonical classes and apply only to pairs', () => {
  const pair = (id, family, stem_eligible = false) => ({ id, pdb_id: '1ABC', residue1_id: 'r1', residue2_id: 'r2', family,
    near: family.startsWith('n'), alternative: false, stem_eligible, values: { opening: 1 } });
  const rows = [pair('watson', 'cWW', true), pair('wobble', 'cWW', true), pair('near', 'ncWW'), pair('hoogsteen', 'tWH'), row('residue', 3)];
  assert.deepEqual(selectRows(rows, {}, { interactionFamilies: ['tWH'] }).rows.map(item => item.id), ['hoogsteen', 'residue']);
  assert.deepEqual(selectRows(rows, {}, { pairPolicy: 'all', interactionFamilies: ['cWW'] }).rows.map(item => item.id), ['watson', 'wobble', 'residue']);
  assert.deepEqual(selectRows(rows, {}, { pairPolicy: 'near', interactionFamilies: ['ncWW'] }).rows.map(item => item.id), ['near', 'residue']);
  assert.deepEqual(selectRows(rows, {}, { pairPolicy: 'all', stemOnly: true }).rows.map(item => item.id), ['watson', 'wobble', 'residue']);
});

test('Interaction grouping keeps exact, near and alternative pair observations distinct', () => {
  const base = { pdb_id: '1ABC', residue1_id: 'r1', residue2_id: 'r2', family: 'cWW', pair_label: 'G-U', values: { opening: 1 } };
  const rows = selectRows([{ ...base, id: 'exact' }, { ...base, id: 'near', near: true }, { ...base, id: 'alternate', alternative: true }], {}, { pairPolicy: 'all' }).rows;
  const result = distribution(rows, { id: 'opening' }, { groupBy: 'interaction', sigma: 0 });
  assert.deepEqual(result.series.map(series => series.key), ['cWW', 'ncWW', 'cWW (alternative)']);
  assert.deepEqual(result.series.map(series => series.rowIds), [['exact'], ['near'], ['alternate']]);
  assert.equal(result.coverage.memberships, 3);
});

test('Pucker selection follows all recorded members of survey and pair observations', () => {
  const rows = [
    { id: 'local', pdb_id: '1ABC', pucker_class: "C3'-endo" },
    { id: 'survey', pdb_id: '1ABC', pucker_classes: ["C3'-endo"] },
    { id: 'pair', pdb_id: '1ABC', pucker_classes: ["C3'-endo", "C3'-endo"] },
    { id: 'mixed', pdb_id: '1ABC', pucker_classes: ["C3'-endo", "C2'-endo"] },
    { id: 'missing', pdb_id: '1ABC', pucker_classes: [null] },
  ];
  assert.deepEqual(selectRows(rows, {}, { puckerStates: ["C3'-endo"] }).rows.map(row => row.id), ['local', 'survey', 'pair']);
  assert.equal(selectRows(rows).rows.length, 5);
});
