/** Actual summary presentation and Survey semantic colors against independent DNA oracles. */
import assert from 'node:assert/strict';
import { mkdir, readFile, writeFile } from 'node:fs/promises';
import path from 'node:path';
import vm from 'node:vm';
import { createHash } from 'node:crypto';
import { pathToFileURL } from 'node:url';
import { downloadCsv, waitReady } from './helpers.mjs';

const workspace = path.resolve(process.env.RNA_WORKSPACE || process.cwd());
const output = process.env.RNA_BROWSER_OUTPUT || path.join(workspace, 'data/pure_rna/browser_validation/statistic-display-parity');
await mkdir(output, { recursive: true });
const dna = await readFile(new URL('../../../js/pure-dna.js', import.meta.url), 'utf8');
const extract = name => { const start = dna.indexOf(`function ${name}(`), end = dna.indexOf('\nfunction ', start + 1); assert(start >= 0); return dna.slice(start, end); };
const oracle = vm.createContext({});
vm.runInContext([extract('wrapCircular'), extract('circularDisplayValue'), dna.match(/const BASE_GEOMETRY_BIN_META = \{[\s\S]*?\n\};/)[0], 'globalThis.binColors = BASE_GEOMETRY_BIN_META;'].join('\n'), oracle);
const expectedColors = JSON.parse(JSON.stringify(oracle.binColors));
const format = value => Number.isFinite(value) ? value.toLocaleString('en-US', { maximumFractionDigits: 3 }) : 'Undefined';
const expectedAngle = (value, mode, period = 360) => oracle.circularDisplayValue(value, mode, period);
const { chromium } = await import(pathToFileURL(process.env.PLAYWRIGHT_MODULE || path.join(workspace, 'data/pure_rna/browser-tools/node_modules/playwright-core/index.mjs')).href);
const browser = await chromium.launch({ headless: true });
const page = await browser.newPage({ viewport: { width: 1440, height: 1000 }, acceptDownloads: true });
const report = { started: new Date().toISOString(), dnaSha256: createHash('sha256').update(dna).digest('hex'), checks: [], errors: [] };
page.on('pageerror', e => report.errors.push(e.message));
page.on('console', m => { if (m.type() === 'error') report.errors.push(m.text()); });
const click = async (group, value) => { await page.click(`#${group} button[data-value="${value}"]`); await waitReady(page); };
const record = (name, evidence) => report.checks.push({ name, ...evidence });
const withoutSnapshot = rows => rows.map(({ snapshot_id, ...row }) => row);
async function openingEvidence(name) {
  const e = await page.evaluate(() => {
    const a = window.rnaExplorer, result = a.snapshots.survey.result;
    return { keys: result.series.map(s => s.key), traces: document.querySelector('#baseGeometryPlot').data.map(t => ({ key: t.name.split(' (n=')[0], color: t.line.color, fillcolor: t.fillcolor })),
      scientific: result.series.map(s => ({ key: s.key, rows: s.rowIds, values: s.values })), n: result.coverage.plottedRows };
  });
  const order = ['small', 'middle', 'large'].filter(key => e.keys.includes(key));
  assert.deepEqual(e.traces.map(t => t.key), order);
  for (const trace of e.traces) { assert.equal(trace.color, expectedColors[trace.key].color); assert.equal(trace.fillcolor, trace.color + '18'); }
  record(name, { keys: e.keys, traces: e.traces, n: e.n }); return e;
}
try {
  await page.goto(process.env.RNA_URL || 'http://127.0.0.1:8767/nucleic.pages/rna/'); await waitReady(page);
  report.release = await page.evaluate(() => ({ build: window.rnaExplorer.manifest.build_id, partial: window.rnaExplorer.manifest.partial })); assert.equal(report.release.partial, false);
  let csvReference;
  for (const mode of ['wrap_360', 'signed_180', 'auto']) {
    await click('circularModeGroup', mode);
    const e = await page.evaluate(async () => {
      const a = window.rnaExplorer, snapshot = a.snapshots.distribution, r = snapshot.result;
      const source = new Map((await a.repository.loadFamily('backbone')).rows.map(row => [row.id, row]));
      const cards = [...document.querySelectorAll('#seriesSummary > .card')].map(card => ({ title: card.querySelector('h3').textContent,
        metrics: Object.fromEntries([...card.querySelectorAll('.metric')].map(m => [m.querySelector('.metric-label').textContent, m.querySelector('.metric-value').textContent])) }));
      return { mode: a.state.display.circularMode, cut: r.displayCut, unit: r.parameter.unit, period: r.parameter.period,
        series: r.series.map((s, i) => {
          const raw = s.rowIds.map(id => source.get(id).values.chi);
          let sin = 0, cos = 0; for (const x of raw) { sin += Math.sin(x * Math.PI / 180); cos += Math.cos(x * Math.PI / 180); }
          const mean = (Math.atan2(sin, cos) * 180 / Math.PI + 360) % 360;
          return { key: s.key, statistics: s.statistics, independentMean: mean, exactRaw: raw.every((v, j) => Object.is(v, s.values[j])), card: cards[i] };
        }) };
    });
    for (const series of e.series) {
      assert(series.exactRaw); assert(Math.abs(series.independentMean - series.statistics.mean) < 1e-9);
      assert.equal(series.card.metrics.Mean, `${format(expectedAngle(series.independentMean, mode, e.period))} ${e.unit}`);
      assert.equal(series.card.metrics['Smoothed peak'], `${format(expectedAngle(series.statistics.peak, mode, e.period))} ${e.unit}`);
      assert.equal(series.card.metrics['Circular std. deviation'], `${format(series.statistics.std)} ${e.unit}`);
      assert.equal(series.card.metrics['Resultant length'], format(series.statistics.resultant));
    }
    const csv = await downloadCsv(page, '#filteredCsvDownload', path.join(output, `main-${mode}.csv`));
    if (csvReference) assert.deepEqual(withoutSnapshot(csv.rows), csvReference); else csvReference = withoutSnapshot(csv.rows);
    record(`Actual ${mode} summary cards match source mean and DNA angular display`, { ...e, csvRows: csv.rows.length });
  }
  await page.click('#baseGeometryLoad'); await waitReady(page);
  await page.selectOption('#baseGeometryTermSelect', 'c_n1_c2_n3_c4'); await waitReady(page);
  await page.selectOption('#surveyOpeningSelect', 'bins'); await waitReady(page);
  const baseline = await openingEvidence('Actual Survey opening groups use canonical colors/order');
  assert(baseline.n > 0);
  const before = await downloadCsv(page, '#surveyCsvDownload', path.join(output, 'survey-before-restyle.csv'));
  await page.evaluate(() => { window.presentationSnapshot = window.rnaExplorer.snapshots.survey; });
  await click('traceStyleGroup', 'line');
  assert(await page.evaluate(() => window.presentationSnapshot.result.series === window.rnaExplorer.snapshots.survey.result.series));
  const styled = await openingEvidence('Actual Survey trace-style redraw retains colors and scientific rows');
  assert.deepEqual(styled.scientific, baseline.scientific);
  const after = await downloadCsv(page, '#surveyCsvDownload', path.join(output, 'survey-after-restyle.csv'));
  assert.deepEqual(withoutSnapshot(after.rows), withoutSnapshot(before.rows));
  // Choose an actual release entry whose selected incidences lose one or more bins.
  const subset = await page.evaluate(() => {
    const byEntry = new Map();
    for (const series of window.rnaExplorer.snapshots.survey.result.series) for (const row of series.rows) {
      if (!byEntry.has(row.pdb_id)) byEntry.set(row.pdb_id, new Set()); byEntry.get(row.pdb_id).add(series.key);
    }
    return [...byEntry].find(([, bins]) => bins.has('large') && bins.size < 3)?.[0] ?? [...byEntry].find(([, bins]) => bins.size < 3)?.[0];
  });
  assert(subset, 'Real release provides a partial-bin subset');
  await page.evaluate(pdb => window.rnaExplorer.setSelection({ search: pdb }), subset); await waitReady(page);
  const sparse = await openingEvidence('Actual PDB subset retains colors when opening bins disappear');
  assert(sparse.keys.length < baseline.keys.length); assert(sparse.n > 0); report.realSparsePdb = subset;

  // Explicit synthetic seam and rare subset examples use the real DOM module.
  const fixture = await page.evaluate(async () => {
    const { summaryCards, distributionTraces } = await import('./views/panels.js');
    const { csv, createPlotSnapshot } = await import('./core/export.js');
    const host = document.createElement('div'); document.body.append(host);
    const metric = label => [...host.querySelectorAll('.metric')].find(m => m.querySelector('.metric-label').textContent === label)?.querySelector('.metric-value').textContent;
    const cards = [];
    for (const mode of ['auto', 'wrap_360', 'signed_180']) for (const mean of [5, 180, 270, 540, null]) {
      const result = { parameter: { id: 'fixture_angle', period: 360, unit: '°' }, displaySpec: { circularMode: mode }, displayCut: 185,
        series: [{ key: 'fixture', values: [5], statistics: { n: 1, mean, peak: mean, std: 7, resultant: 0.5 } }] };
      const original = JSON.stringify(result); summaryCards(host, result);
      cards.push({ mode, mean, displayedMean: metric('Mean'), peak: metric('Smoothed peak'), unchanged: JSON.stringify(result) === original });
    }
    const make = keys => ({ kind: 'distribution', parameter: { id: 'fixture_angle', period: 360 }, displaySpec: { groupBy: 'opening_bin' },
      series: keys.map((key, i) => ({ key, values: [i], rowIds: [key], rows: [{ id: key }], x: [10], y: [1] })) });
    const subsets = [];
    for (const keys of [['large', 'small', 'middle'], ['large', 'small'], ['large'], ['middle']]) {
      const snapshot = createPlotSnapshot({ result: make(keys) }), prior = csv(snapshot), serial = JSON.stringify(snapshot);
      const traces = distributionTraces(snapshot.result, { groupBy: 'base' });
      const plot = document.createElement('div'); host.replaceChildren(plot);
      await window.Plotly.newPlot(plot, traces, {}, { displayModeBar: false });
      subsets.push({ inputKeys: keys, keys: plot.data.map(t => t.name.split(' (n=')[0]), colors: plot.data.map(t => t.line.color),
        csvUnchanged: csv(snapshot) === prior, snapshotUnchanged: JSON.stringify(snapshot) === serial });
      window.Plotly.purge(plot);
    }
    host.remove(); return { cards, subsets };
  });
  for (const card of fixture.cards) { const value = expectedAngle(card.mean, card.mode); const expected = value === null ? 'Undefined' : `${format(value)} °`; assert.equal(card.displayedMean, expected); assert.equal(card.peak, expected); assert(card.unchanged); }
  for (const subset of fixture.subsets) { const order = ['small', 'middle', 'large'].filter(key => subset.inputKeys.includes(key)); assert.deepEqual(subset.keys, order); assert.deepEqual(subset.colors, order.map(key => expectedColors[key].color)); assert(subset.csvUnchanged && subset.snapshotUnchanged); }
  record('Explicit synthetic DOM seam/half-period/undefined and missing-middle examples', fixture);
  assert.deepEqual(report.errors, []); report.passed = true; console.log(`PASS ${report.checks.length} statistic display browser checks`);
} catch (error) { report.passed = false; report.failure = error.stack; throw error; }
finally { report.finished = new Date().toISOString(); await writeFile(path.join(output, 'report.json'), JSON.stringify(report, null, 2)); await browser.close(); }
