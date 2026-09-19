/** Actual RNA axes use DNA canonical labels without modifying raw observations. */
import assert from 'node:assert/strict';
import { mkdir, readFile, writeFile } from 'node:fs/promises';
import { createHash } from 'node:crypto';
import vm from 'node:vm';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { downloadCsv, waitReady } from './helpers.mjs';

const workspace = path.resolve(process.env.RNA_WORKSPACE || process.cwd());
const output = process.env.RNA_BROWSER_OUTPUT || path.join(workspace, 'data/pure_rna/browser_validation/circular-tick-parity');
await mkdir(output, { recursive: true });
// Execute the existing DNA tick functions as an independent display oracle.
const dna = await readFile(new URL('../../../js/pure-dna.js', import.meta.url), 'utf8');
const extract = name => {
  const start = dna.indexOf(`function ${name}(`); assert(start >= 0);
  const end = dna.indexOf('\nfunction ', start + 1); return dna.slice(start, end < 0 ? undefined : end);
};
const oracle = vm.createContext({});
vm.runInContext(['wrapCircular', 'formatCircularTickLabel', 'buildCircularTickSpec'].map(extract).join('\n'), oracle);
const ticks = (period, cut, compact, mode) => JSON.parse(JSON.stringify(oracle.buildCircularTickSpec(period, cut, compact, mode)));
const { chromium } = await import(pathToFileURL(process.env.PLAYWRIGHT_MODULE || path.join(workspace, 'data/pure_rna/browser-tools/node_modules/playwright-core/index.mjs')).href);
const browser = await chromium.launch({ headless: true });
const page = await browser.newPage({ viewport: { width: 1440, height: 1000 }, acceptDownloads: true });
const report = { started: new Date().toISOString(), dnaSourceSha256: createHash('sha256').update(dna).digest('hex'), checks: [], errors: [] };
page.on('pageerror', e => report.errors.push(e.message));
page.on('console', m => { if (m.type() === 'error') report.errors.push(m.text()); });
page.on('requestfailed', r => report.errors.push(r.failure()?.errorText));
const click = async (group, value) => { await page.click(`#${group} button[data-value="${value}"]`); await waitReady(page); };
const select = async (id, value) => { await page.selectOption(`#${id}`, value); await waitReady(page); };
const sameNumbers = (actual, expected, name) => { assert.equal(actual.length, expected.length, name); actual.forEach((value, i) => assert(Math.abs(value - expected[i]) < 1e-8, `${name}: ${value} versus ${expected[i]}`)); };
async function axes(name, circularJoint = true) {
  const e = await page.evaluate(() => {
    const a = window.rnaExplorer, d = a.snapshots.distribution.result, s = a.snapshots.survey.result, j = a.snapshots.joint.result;
    const axis = (plot, key, range, period, compact = false) => ({ key, expectedRange: range, period, compact,
      actualRange: [...plot._fullLayout[key].range], tickvals: plot.layout[key].tickvals,
      ticktext: plot.layout[key].ticktext, domTicks: [...plot.querySelectorAll(key === 'xaxis' ? '.xtick text' : '.ytick text')].map(node => node.textContent),
      traceCoordinates: plot.data.map(trace => trace[key === 'xaxis' ? 'x' : 'y']).flat().every(value => value >= range[0] && value <= range[1]) });
    const mini = document.querySelector('#familyOverview [data-parameter="chi"] .rna-mini-plot');
    const centers = mini.data[0].x, width = centers[1] - centers[0];
    return { mode: a.state.display.circularMode, axes: [
      { panel: 'main', ...axis(document.querySelector('#plot'), 'xaxis', d.range, d.parameter.period) },
      { panel: 'mini', ...axis(mini, 'xaxis', [centers[0] - width / 2, centers.at(-1) + width / 2], 360, true) },
      { panel: 'survey', ...axis(document.querySelector('#baseGeometryPlot'), 'xaxis', s.range, s.parameter.period) },
      { panel: 'joint-x', ...axis(document.querySelector('#jointPlot'), 'xaxis', j.xRange, j.xParameter.period) },
      { panel: 'joint-y', ...axis(document.querySelector('#jointPlot'), 'yaxis', j.yRange, j.yParameter.period) }],
      counts: { main: d.coverage.plottedRows, survey: s.coverage.plottedRows, joint: j.points.length } };
  });
  for (const axis of e.axes) {
    sameNumbers(axis.actualRange, axis.expectedRange, `${name} ${axis.panel} complete range`);
    assert(axis.traceCoordinates, `${name} ${axis.panel} centers inside range`);
    if (!axis.period) continue;
    const expected = ticks(axis.period, axis.expectedRange[0], axis.compact, e.mode);
    sameNumbers(axis.tickvals, expected.tickvals.map(tick => tick + axis.expectedRange[0]), `${name} ${axis.panel} positions`);
    assert.deepEqual(axis.ticktext, expected.ticktext, `${name} ${axis.panel} canonical labels`);
    assert.deepEqual(axis.domTicks, expected.ticktext, `${name} ${axis.panel} actual SVG text`);
    assert.equal(axis.actualRange[1] - axis.actualRange[0], axis.period);
  }
  if (circularJoint) assert(e.axes.filter(a => a.panel.startsWith('joint')).every(a => a.period === 360));
  report.checks.push({ name, ...e }); return e;
}
const csvBaselines = new Map();
async function exports(mode) {
  const evidence = [];
  for (const [panel, button] of [['main', '#filteredCsvDownload'], ['survey', '#surveyCsvDownload'], ['joint', '#jointCsvDownload']]) {
    const data = await downloadCsv(page, button, path.join(output, `${panel}-${mode}.csv`));
    const normalized = data.rows.map(({ snapshot_id, ...row }) => row);
    if (csvBaselines.has(panel)) assert.deepEqual(normalized, csvBaselines.get(panel), `${panel} exact CSV identity and values across modes`);
    else csvBaselines.set(panel, normalized);
    evidence.push({ panel, rows: data.rows.length });
  }
  report.checks.push({ name: `Raw CSV invariant ${mode}`, exports: evidence });
}
async function svg(name, selector, expectedX, expectedY) {
  const plot = page.locator(selector); await plot.scrollIntoViewIfNeeded(); await plot.hover({ position: { x: 20, y: 20 } });
  const downloadEvent = page.waitForEvent('download'); await plot.locator('.modebar-btn[data-title*="Download plot"]').click();
  const download = await downloadEvent; assert.equal(await download.failure(), null);
  const file = path.join(output, name + '.svg'); await download.saveAs(file);
  const xml = await readFile(file, 'utf8');
  const actual = await page.evaluate(xml => {
    const doc = new DOMParser().parseFromString(xml, 'image/svg+xml');
    return { x: [...doc.querySelectorAll('.xtick text')].map(n => n.textContent), y: [...doc.querySelectorAll('.ytick text')].map(n => n.textContent), errors: doc.querySelectorAll('parsererror').length };
  }, xml);
  assert.equal(actual.errors, 0); assert.deepEqual(actual.x, expectedX); if (expectedY) assert.deepEqual(actual.y, expectedY);
  report.checks.push({ name: `Standalone SVG canonical ticks ${name}`, file, ...actual });
}
try {
  await page.goto(process.env.RNA_URL || 'http://127.0.0.1:8767/nucleic.pages/rna/'); await waitReady(page);
  report.release = await page.evaluate(() => ({ build: window.rnaExplorer.manifest.build_id, partial: window.rnaExplorer.manifest.partial }));
  assert.equal(report.release.partial, false);
  await select('family2Select', 'backbone'); await select('parameter2Select', 'alpha');
  await page.click('#baseGeometryLoad'); await waitReady(page); await select('baseGeometryTermSelect', 'c_n1_c2_n3_c4');
  for (const mode of ['wrap_360', 'signed_180', 'auto']) {
    await click('circularModeGroup', mode); const e = await axes(`Complete circular axes ${mode}`); await exports(mode);
    await page.evaluate(() => { const a = window.rnaExplorer; window.axisPrevious = { main: a.snapshots.distribution.result, joint: a.snapshots.joint.result, survey: a.snapshots.survey.result }; });
    await click('traceStyleGroup', mode === 'signed_180' ? 'filled' : 'line');
    await click('jointPaletteGroup', mode === 'signed_180' ? 'warm' : 'viridis');
    assert(await page.evaluate(() => { const a = window.rnaExplorer, p = window.axisPrevious; return p.main.series === a.snapshots.distribution.result.series && p.joint.points === a.snapshots.joint.result.points && p.joint.z === a.snapshots.joint.result.z && p.survey.series === a.snapshots.survey.result.series; }), 'Visual restyles preserve scientific arrays');
    await axes(`Restyled complete axes ${mode}`);
    await svg(`main-${mode}`, '#plot', e.axes.find(a => a.panel === 'main').ticktext);
    await svg(`joint-${mode}`, '#jointPlot', e.axes.find(a => a.panel === 'joint-x').ticktext, e.axes.find(a => a.panel === 'joint-y').ticktext);
  }
  await select('parameterSelect', 'e_z');
  await axes('Linear main and joint-x retain full advisory range', false);
  const linear = await page.evaluate(() => ({ main: window.rnaExplorer.snapshots.distribution.result.range, joint: window.rnaExplorer.snapshots.joint.result.xRange }));
  assert(linear.main[0] <= -180 && linear.main[1] >= 180); assert(linear.joint[0] <= -180 && linear.joint[1] >= 180);
  assert.deepEqual(report.errors, []); report.passed = true;
  console.log(`PASS ${report.checks.length} circular tick browser checks`);
} catch (error) { report.passed = false; report.failure = error.stack; throw error; }
finally { report.finished = new Date().toISOString(); await writeFile(path.join(output, 'report.json'), JSON.stringify(report, null, 2)); await browser.close(); }
