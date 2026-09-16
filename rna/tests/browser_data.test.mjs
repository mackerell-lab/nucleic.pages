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
