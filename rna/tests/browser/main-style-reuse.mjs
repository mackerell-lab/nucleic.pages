/** Actual trace-style controls preserve observations, exports and unrelated panels. */
import assert from 'node:assert/strict';
import { mkdir, writeFile } from 'node:fs/promises';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { downloadCsv, waitReady } from './helpers.mjs';
const workspace = path.resolve(process.env.RNA_WORKSPACE || process.cwd());
const output = process.env.RNA_BROWSER_OUTPUT || path.join(workspace, 'data/pure_rna/browser_validation/main-style-reuse');
await mkdir(output, { recursive: true });
const { chromium } = await import(pathToFileURL(process.env.PLAYWRIGHT_MODULE || path.join(workspace, 'data/pure_rna/browser-tools/node_modules/playwright-core/index.mjs')).href);
const browser = await chromium.launch({ headless: true });
const page = await browser.newPage({ viewport: { width: 1440, height: 1000 }, acceptDownloads: true });
const report = { started: new Date().toISOString(), checks: [], errors: [] };
let expectedFailure = false;
page.on('pageerror', error => report.errors.push(error.message));
page.on('console', message => { if (message.type() === 'error' && !(expectedFailure && message.text().includes('Injected trace failure'))) report.errors.push(message.text()); });
page.on('requestfailed', request => report.errors.push(request.failure()?.errorText));
const click = async (group, value) => { await page.click(`#${group} button[data-value="${value}"]`); await waitReady(page); };
const remember = () => page.evaluate(() => {
  const a = window.rnaExplorer;
  window.traceBefore = { ...a.snapshots, cards: [...document.querySelectorAll('#familyOverview [data-parameter]')],
    loads: window.traceLoads, requests: performance.getEntriesByType('resource').length };
});
const clean = rows => rows.map(({ snapshot_id, ...row }) => row);
try {
  await page.goto(process.env.RNA_URL || 'http://127.0.0.1:8767/nucleic.pages/rna/'); await waitReady(page);
  report.build = await page.evaluate(() => window.rnaExplorer.manifest.build_id);
  assert.equal(report.build, 'full_packed_family_20260919');
  await page.evaluate(() => {
    const a = window.rnaExplorer; window.traceLoads = 0;
    for (const name of ['loadFamily', 'loadRelations', 'loadMetadata', 'loadSurveyScalars', 'loadSurveyCoordinates']) {
      const original = a.repository[name]?.bind(a.repository); if (!original) continue;
      a.repository[name] = (...args) => { window.traceLoads++; return original(...args); };
    }
  });
  await remember(); await click('traceStyleGroup', 'line');
  const lazy = await page.evaluate(() => ({ survey: window.rnaExplorer.state.survey.loaded,
    coordinates: window.rnaExplorer.state.survey.coordinatesLoaded, loads: window.traceLoads,
    reused: window.rnaExplorer.snapshots.distribution.result.series === window.traceBefore.distribution.result.series }));
  assert.deepEqual(lazy, { survey: false, coordinates: false, loads: 0, reused: true });
  report.checks.push({ name: 'Unopened Survey and coordinates remain lazy', ...lazy });
  await page.click('#baseGeometryLoad'); await waitReady(page);
  await page.selectOption('#family2Select', 'backbone'); await waitReady(page);
  await page.selectOption('#parameter2Select', 'alpha'); await waitReady(page);
  const baseline = {};
  for (const [key, selector] of [['distribution', '#filteredCsvDownload'], ['survey', '#surveyCsvDownload'], ['joint', '#jointCsvDownload']]) {
    baseline[key] = await downloadCsv(page, selector, path.join(output, `before-${key}.csv`));
  }
  for (const style of ['filled', 'line', 'filled']) {
    await remember(); await click('traceStyleGroup', style);
    const evidence = await page.evaluate(() => {
      const a = window.rnaExplorer, b = window.traceBefore;
      return { style: a.state.display.traceStyle, loads: window.traceLoads - b.loads,
        requests: performance.getEntriesByType('resource').length - b.requests,
        mainReuse: b.distribution.result.series === a.snapshots.distribution.result.series,
        surveyReuse: b.survey.result.series === a.snapshots.survey.result.series,
        jointReuse: b.joint.result.points === a.snapshots.joint.result.points && b.joint.result.z === a.snapshots.joint.result.z,
        fresh: ['distribution', 'survey', 'joint'].every(k => a.snapshots[k].snapshot_id !== b[k].snapshot_id && a.snapshots[k].display_spec.traceStyle === a.state.display.traceStyle && a.snapshots[k].result.displaySpec.traceStyle === a.state.display.traceStyle),
        overviewStable: b.cards.every((node, i) => node === document.querySelectorAll('#familyOverview [data-parameter]')[i]),
        cards: b.cards.length,
        fills: ['plot', 'baseGeometryPlot'].map(id => document.getElementById(id).data.map(t => t.fill)),
        snapshotIds: Object.fromEntries(['distribution', 'survey', 'joint'].map(k => [k, a.snapshots[k].snapshot_id])),
        jointKeyCurrent: a.completedJointKey === a.jointAnalysisKey(a.state) };
    });
    for (const key of ['mainReuse', 'surveyReuse', 'jointReuse', 'fresh', 'overviewStable', 'jointKeyCurrent']) assert(evidence[key], key);
    assert(evidence.cards > 0); assert.equal(evidence.loads, 0); assert.equal(evidence.requests, 0);
    assert(evidence.fills.flat().every(fill => fill === (style === 'filled' ? 'tozeroy' : 'none')));
    for (const [key, selector] of [['distribution', '#filteredCsvDownload'], ['survey', '#surveyCsvDownload'], ['joint', '#jointCsvDownload']]) {
      const result = await downloadCsv(page, selector, path.join(output, `${report.checks.length}-${style}-${key}.csv`));
      assert.deepEqual(clean(result.rows), clean(baseline[key].rows));
      assert(result.rows.every(row => row.snapshot_id === evidence.snapshotIds[key]));
    }
    report.checks.push({ name: `Trace style ${style}`, ...evidence });
  }
  await remember(); await click('jointPaletteGroup', 'ocean');
  const palette = await page.evaluate(() => ({ reused: window.traceBefore.joint.result === window.rnaExplorer.snapshots.joint.result,
    loads: window.traceLoads - window.traceBefore.loads }));
  assert.deepEqual(palette, { reused: true, loads: 0 }); report.checks.push({ name: 'Palette still reuses joint after trace style', ...palette });
  // Joint-only controls must update the completed full-state marker for subsequent trace changes.
  await remember(); await click('traceStyleGroup', 'line');
  assert(await page.evaluate(() => window.traceBefore.distribution.result.series === window.rnaExplorer.snapshots.distribution.result.series));
  report.checks.push({ name: 'Trace style still reuses after palette' });
  await remember(); await click('smoothingSigmaGroup', '0.8');
  assert(await page.evaluate(() => window.traceBefore.distribution.result.series !== window.rnaExplorer.snapshots.distribution.result.series));
  report.checks.push({ name: 'Numerical smoothing change recomputes' });
  // Fail the second plot so the already painted first panel cannot publish stale exports.
  await remember(); expectedFailure = true;
  await page.evaluate(() => {
    const a = window.rnaExplorer, original = a.plot.bind(a); window.traceOriginalPlot = original; let count = 0;
    a.plot = async (...args) => { if (++count === 2) throw Error('Injected trace failure'); return original(...args); };
  });
  await page.click('#traceStyleGroup button[data-value="filled"]');
  await page.waitForFunction(() => document.querySelector('#appStatus').dataset.state === 'error');
  assert(await page.evaluate(() => {
    const a = window.rnaExplorer;
    return !a.fullRenderComplete && a.snapshots.distribution === window.traceBefore.distribution
      && a.snapshots.survey === window.traceBefore.survey && document.querySelector('#filteredCsvDownload').disabled;
  }));
  await page.evaluate(() => { window.rnaExplorer.plot = window.traceOriginalPlot; });
  await page.click('#resetFilters'); await waitReady(page); expectedFailure = false;
  assert(await page.evaluate(() => window.rnaExplorer.fullRenderComplete && !document.querySelector('#filteredCsvDownload').disabled));
  report.checks.push({ name: 'Failed partial trace update blocks exports and Reset repairs all panels' });
  // Force two real button clicks to overlap at Plotly completion.
  await page.evaluate(() => {
    const a = window.rnaExplorer, original = a.plot.bind(a); let first = true;
    a.plot = async (...args) => {
      if (first) { first = false; window.traceEntered = true; await new Promise(resolve => { window.traceFinish = resolve; }); }
      return original(...args);
    };
  });
  await page.click('#traceStyleGroup button[data-value="line"]'); await page.waitForFunction(() => window.traceEntered);
  await page.click('#traceStyleGroup button[data-value="filled"]');
  await page.evaluate(() => window.traceFinish()); await waitReady(page);
  assert(await page.evaluate(() => {
    const a = window.rnaExplorer; return a.state.display.traceStyle === 'filled' && a.snapshots.distribution.display_spec.traceStyle === 'filled'
      && document.querySelector('#plot').data.every(t => t.fill === 'tozeroy') && a.fullRenderComplete;
  }));
  report.checks.push({ name: 'Overlapping style clicks finish with current data and latest style' });
  assert.deepEqual(report.errors, []); report.passed = true;
  console.log(`PASS ${report.checks.length} main style reuse browser checks`);
} catch (error) { report.failure = error.stack; report.passed = false; throw error; }
finally { report.finished = new Date().toISOString(); await writeFile(path.join(output, 'report.json'), JSON.stringify(report, null, 2) + '\n'); await browser.close(); }
