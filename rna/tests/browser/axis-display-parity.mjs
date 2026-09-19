/** Axis defaults and seam controls preserve actual deposited RNA measurements. */
import assert from 'node:assert/strict';
import { mkdir, readFile, writeFile } from 'node:fs/promises';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { waitReady, downloadCsv } from './helpers.mjs';
const workspace = path.resolve(process.env.RNA_WORKSPACE || process.cwd());
const output = path.resolve(process.env.RNA_BROWSER_OUTPUT || path.join(workspace, 'data/pure_rna/browser_validation/axis-display-parity'));
await mkdir(output, { recursive: true });
const { chromium } = await import(pathToFileURL(process.env.PLAYWRIGHT_MODULE || path.join(workspace, 'data/pure_rna/browser-tools/node_modules/playwright-core/index.mjs')).href);
const browser = await chromium.launch({ headless: true });
const page = await browser.newPage({ viewport: { width: 1440, height: 1000 }, acceptDownloads: true });
const report = { startedAt: new Date().toISOString(), checks: [], pageErrors: [], failedRequests: [] };
page.on('pageerror', error => report.pageErrors.push(error.message));
page.on('requestfailed', request => report.failedRequests.push(request.url()));
const record = (name, evidence) => report.checks.push({ name, passed: true, ...evidence });
try {
  if (process.env.RNA_AXIS_LEGACY_ANALYSIS) {
    const body = await readFile(process.env.RNA_AXIS_LEGACY_ANALYSIS, 'utf8');
    await page.route('**/rna/core/analysis.js', route => route.fulfill({ body, contentType: 'text/javascript' }));
    report.legacyAnalysis = process.env.RNA_AXIS_LEGACY_ANALYSIS;
  }
  await page.goto(process.env.RNA_URL || 'http://127.0.0.1:8767/nucleic.pages/rna/', { waitUntil: 'domcontentloaded', timeout: 120000 }); await waitReady(page);
  const singleton = await page.evaluate(async () => {
    const app = window.rnaExplorer, rows = (await app.repository.loadFamily('backbone')).rows, groups = new Map();
    for (const row of rows) {
      if (!Number.isFinite(row.values.e_z) || !['available', 'ok'].includes(row.statuses.e_z)) continue;
      const key = `${row.pdb_id}|${row.comp_id}`;
      if (!groups.has(key)) groups.set(key, []); groups.get(key).push(row);
    }
    const row = [...groups.values()].find(rows => rows.length === 1 && Number.isFinite(rows[0].values.chi) && Math.abs(rows[0].values.e_z) < 180)?.[0];
    if (!row) throw new Error('Source lacks a singleton linear/circular subset');
    await app.setSelection({ components: 'all', methods: [], resolutionMax: null, contexts: [row.comp_id], search: row.pdb_id });
    return { id: row.id, pdb: row.pdb_id, base: row.comp_id, e_z: row.values.e_z, chi: row.values.chi };
  }); await waitReady(page);
  await page.selectOption('#parameterSelect', 'e_z'); await waitReady(page);
  await page.selectOption('#family2Select', 'backbone'); await waitReady(page);
  await page.selectOption('#parameter2Select', 'chi'); await waitReady(page);
  const single = await page.evaluate(() => {
    const app = window.rnaExplorer, d = app.snapshots.distribution.result, j = app.snapshots.joint.result;
    return { buildId: app.manifest.build_id, range: d.range, raw: d.series.flatMap(s => s.values),
      xRange: j.xRange, points: j.points.map(p => [p.left_id, p.x, p.y]), finite: d.coverage.finiteRows, plotted: d.coverage.plottedRows };
  });
  assert.deepEqual(single.range, [-180, 180], 'A singleton should retain the advisory linear range');
  assert.deepEqual(single.xRange, [-180, 180]); assert.deepEqual(single.raw, [singleton.e_z]);
  assert.deepEqual(single.points, [[singleton.id, singleton.e_z, singleton.chi]]); assert.equal(single.finite, single.plotted);
  record('Actual singleton retains advisory 1D and joint range', { singleton, ...single });
  await page.evaluate(() => window.rnaExplorer.setSelection({ search: 'NO_MATCHING_RNA_ENTRY_000' })); await waitReady(page);
  const empty = await page.evaluate(() => ({ range: window.rnaExplorer.snapshots.distribution.result.range, count: window.rnaExplorer.snapshots.distribution.result.coverage.plottedRows, points: window.rnaExplorer.snapshots.joint.result.points.length }));
  assert.deepEqual(empty.range, [-180, 180]); assert.equal(empty.count, 0); assert.equal(empty.points, 0);
  record('Empty subset clears old observations without collapsing default axis', empty);

  await page.click('#resetFilters'); await waitReady(page);
  await page.selectOption('#parameterSelect', 'chi'); await waitReady(page);
  await page.selectOption('#family2Select', 'ribose_2oh'); await waitReady(page);
  await page.selectOption('#parameter2Select', 'c2_o2_length'); await waitReady(page);
  await page.click('#baseGeometryLoad'); await waitReady(page);
  await page.selectOption('#baseGeometryTermSelect', 'c_n1_c2_n3_c4'); await waitReady(page);
  const baseline = new Map();
  for (const mode of ['wrap_360', 'signed_180', 'auto']) {
    await page.click(`#circularModeGroup button[data-value="${mode}"]`); await waitReady(page);
    const evidence = await page.evaluate(async () => {
      const app = window.rnaExplorer;
      const left = new Map((await app.repository.loadFamily('backbone')).rows.map(row => [row.id, row]));
      const right = new Map((await app.repository.loadFamily('ribose_2oh')).rows.map(row => [row.id, row]));
      const surveySource = await app.repository.loadSurveyScalars(app.state.survey.termId);
      const scalars = new Map((surveySource.rows ?? surveySource).map(row => [row.id, row]));
      const d = app.snapshots.distribution.result, j = app.snapshots.joint.result, s = app.snapshots.survey.result;
      const familyRaw = d.series.every(series => series.values.every((value, i) => Object.is(value, left.get(series.rowIds[i]).values.chi)));
      const surveyRaw = s.series.every(series => series.values.every((value, i) => Object.is(value, scalars.get(series.rowIds[i]).value)));
      const jointRaw = j.points.every(point => Object.is(point.x, left.get(point.left_id).values.chi) && Object.is(point.y, right.get(point.right_id).values.c2_o2_length));
      const rendered = (id, result, kind) => {
        const plot = document.getElementById(id), xRange = kind === 'joint' ? result.xRange : result.range;
        return { x: [...plot._fullLayout.xaxis.range], result: xRange,
          centersVisible: (kind === 'joint' ? result.x : result.series.flatMap(series => series.x)).every(x => x >= plot._fullLayout.xaxis.range[0] && x <= plot._fullLayout.xaxis.range[1]) };
      };
      return { mode: app.state.display.circularMode, term: app.state.survey.termId, familyRaw, surveyRaw, jointRaw,
        distribution: { range: d.range, finite: d.coverage.finiteRows, plotted: d.coverage.plottedRows, cut: d.displayCut },
        survey: { range: s.range, finite: s.coverage.finiteRows, plotted: s.coverage.plottedRows, cut: s.displayCut },
        joint: { xRange: j.xRange, yRange: j.yRange, finite: j.coverage.finitePoints, plotted: j.coverage.plottedPoints },
        rendered: [rendered('plot', d), rendered('baseGeometryPlot', s), rendered('jointPlot', j, 'joint')] };
    });
    assert(evidence.familyRaw && evidence.surveyRaw && evidence.jointRaw, 'Plotted numbers changed from independent source assets');
    for (const panel of [evidence.distribution, evidence.survey, evidence.joint]) { assert(panel.finite > 0); assert.equal(panel.finite, panel.plotted); }
    for (const panel of evidence.rendered) assert(panel.centersVisible);
    for (const range of [evidence.distribution.range, evidence.survey.range, evidence.joint.xRange]) {
      assert.equal(range[1] - range[0], 360);
      if (mode !== 'auto') assert.deepEqual(range, mode === 'signed_180' ? [-180, 180] : [0, 360]);
    }
    assert(evidence.joint.yRange[0] <= 1.2 && evidence.joint.yRange[1] >= 1.6, 'Advisory bond length range was shrunk');
    const exports = [];
    for (const [name, selector] of [['distribution', '#filteredCsvDownload'], ['survey', '#surveyCsvDownload'], ['joint', '#jointCsvDownload']]) {
      const csv = await downloadCsv(page, selector, path.join(output, `${name}-${mode}.csv`));
      const rawRows = csv.rows.map(({ snapshot_id, ...row }) => row);
      if (!baseline.has(name)) baseline.set(name, { headers: csv.headers, rows: rawRows });
      else assert.deepEqual({ headers: csv.headers, rows: rawRows }, baseline.get(name), `${name} raw CSV changed with circular display mode`);
      exports.push({ panel: name, rows: csv.rows.length });
    }
    record('Actual 1D Survey and cross-family joint axes preserve source and CSV', { ...evidence, exports });
  }
  // Independent deterministic oracles execute the integrated browser analysis module.
  const oracles = await page.evaluate(async () => {
    const { distribution, histogram2D } = await import('./core/analysis.js');
    const circular = { id: 'angle', period: 360 }, linear = { id: 'distance', display_range_default: [-10, 10] };
    const rows = Array.from({ length: 8 }, (_, i) => ({ id: `fixture-${i}`, values: { angle: i % 2 ? 357.5 : 2.5 } }));
    const auto = distribution(rows, circular, { bins: 72, circularMode: 'auto', sigma: 0 });
    const small = distribution(rows.slice(0, 7), circular, { bins: 72, circularMode: 'auto', sigma: 0 });
    const expanded = distribution([{ id: 'low', values: { distance: -4 } }, { id: 'high', values: { distance: 20 } }], linear, { sigma: 0 });
    const j = histogram2D([{ x: 2.5, y: -4 }, { x: 357.5, y: 20 }], circular, linear, { sigma: 0 });
    return { auto: auto.range, lowCount: small.range, expanded: expanded.range, raw: expanded.series.flatMap(s => s.values), joint: j.yRange, points: j.points.map(p => [p.x, p.y]) };
  });
  assert.deepEqual(oracles.auto, [180, 540]); assert.deepEqual(oracles.lowCount, [0, 360]);
  assert.deepEqual(oracles.expanded, [-10, 21.2]); assert.deepEqual(oracles.raw, [-4, 20]); assert.deepEqual(oracles.joint, [-10, 21.2]);
  assert.deepEqual(oracles.points, [[2.5, -4], [357.5, 20]]);
  record('Independent browser-module seam and outside-default range oracles', oracles);
  assert.deepEqual(report.pageErrors, []); assert.deepEqual(report.failedRequests, []);
  report.passed = true; console.log(`PASS ${report.checks.length} axis display checks`);
} catch (error) { report.passed = false; report.failure = error.stack; throw error; }
finally { report.finishedAt = new Date().toISOString(); await writeFile(path.join(output, 'rna-axis-display-parity.json'), JSON.stringify(report, null, 2) + '\n'); await browser.close(); }
