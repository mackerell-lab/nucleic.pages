/** Real timer-delivered input interrupts obsolete scientific rendering between stages. */
import assert from 'node:assert/strict';
import { mkdir, readFile, writeFile } from 'node:fs/promises';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { createHash } from 'node:crypto';
import { waitReady, downloadCsv } from './helpers.mjs';

const workspace = path.resolve(process.env.RNA_WORKSPACE || process.cwd());
const output = process.env.RNA_BROWSER_OUTPUT || path.join(workspace, 'data/pure_rna/browser_validation/cooperative-render');
const phase = process.env.RNA_COOPERATIVE_PHASE || 'after', source = process.env.RNA_HELD_SOURCE;
const inputMode = process.env.RNA_COOPERATIVE_INPUT || 'timer';
const base = process.env.RNA_URL || 'http://127.0.0.1:8767/nucleic.pages/rna/';
await mkdir(output, { recursive: true });
const { chromium } = await import(pathToFileURL(process.env.PLAYWRIGHT_MODULE || path.join(workspace, 'data/pure_rna/browser-tools/node_modules/playwright-core/index.mjs')).href);
const browser = await chromium.launch({ headless: true });
const page = await browser.newPage({ viewport: { width: 1440, height: 1000 }, acceptDownloads: true });
const report = { phase, inputMode, started: new Date().toISOString(), checks: [], errors: [], sourceSha256: {}, limitation: 'One matched shared-host experiment. Input delivery and obsolete completed stages measure interactivity, not total rendering speed. A single synchronous stage can still block.' };
for (const file of ['app/PureRnaExplorer.js', 'app/NucleicAcidExplorer.js']) report.sourceSha256[file] = createHash('sha256').update(await readFile(source ? path.join(source, file) : new URL('../../' + file, import.meta.url))).digest('hex');
page.on('pageerror', error => report.errors.push(error.message));
page.on('console', message => { if (message.type() === 'error' && !message.text().includes('Injected cooperative render failure')) report.errors.push(message.text()); });
if (source) await page.route(base + '**/*.js', async route => {
  const relative = new URL(route.request().url()).pathname.slice(new URL(base).pathname.length);
  await route.fulfill({ status: 200, contentType: 'application/javascript', body: await readFile(path.join(source, relative)) });
});
try {
  await page.goto(base); await waitReady(page);
  await page.selectOption('#family2Select', 'backbone'); await waitReady(page);
  await page.selectOption('#parameter2Select', 'alpha'); await waitReady(page);
  await page.click('#baseGeometryLoad'); await waitReady(page);
  if (inputMode === 'trusted') {
    const button = page.locator('#contextGroup button[data-value="U"]'); await button.scrollIntoViewIfNeeded();
    const bounds = await button.boundingBox(), cdp = await page.context().newCDPSession(page);
    await page.exposeFunction('cooperativeInputRequested', async () => {
      const point = { x: bounds.x + bounds.width / 2, y: bounds.y + bounds.height / 2, button: 'left', clickCount: 1 };
      await cdp.send('Input.dispatchMouseEvent', { ...point, type: 'mousePressed' });
      await cdp.send('Input.dispatchMouseEvent', { ...point, type: 'mouseReleased' });
    });
  }
  const run = await page.evaluate(async inputMode => {
    const a = window.rnaExplorer, started = performance.now(), events = [], heartbeats = [], longtasks = [];
    const observer = new PerformanceObserver(list => { for (const entry of list.getEntries()) longtasks.push({ start: entry.startTime - started, duration: entry.duration }); });
    observer.observe({ type: 'longtask' });
    const heartbeat = setInterval(() => heartbeats.push(performance.now() - started), 25);
    const broadRevision = a.revision + 1, originalPlot = a.plot.bind(a);
    let queued = false, delivered = false, queuedAt, deliveredAt, pendingButtons, broadMiniCount = 0, trusted = false;
    let resolveInput; const inputFinished = new Promise(resolve => { resolveInput = resolve; });
    const deliveredInput = event => {
      if (!queued || delivered || !event.target.closest('#contextGroup button[data-value="U"]')) return;
      delivered = true; trusted = event.isTrusted; deliveredAt = performance.now() - started;
      events.push({ stage: 'scheduled-context-click', at: deliveredAt, completedBroadMinis: broadMiniCount, trusted });
    };
    const received = async event => {
      if (!delivered || !event.target.closest('#contextGroup button[data-value="U"]')) return;
      // Bubble listener runs after the real target control captured its revision.
      pendingButtons = ['filteredCsvDownload', 'plotProvenanceDownload', 'jointCsvDownload', 'surveyCsvDownload'].map(id => ({ id, disabled: document.getElementById(id).disabled }));
      while (document.querySelector('#appStatus').dataset.state === 'loading') await new Promise(resolve => setTimeout(resolve, 10));
      resolveInput();
    };
    document.addEventListener('click', deliveredInput, true); document.addEventListener('click', received);
    a.plot = async function (node, ...args) {
      const revision = this.revision, mini = node.classList.contains('rna-mini-plot');
      if (mini && revision === broadRevision && !queued) {
        queued = true; queuedAt = performance.now() - started;
        // A browser task, not a directly awaited application call: baseline
        // microtask-only rendering cannot observe this until it yields to tasks.
        if (inputMode === 'trusted') window.cooperativeInputRequested();
        else setTimeout(() => document.querySelector('#contextGroup button[data-value="U"]').click(), 0);
      }
      await originalPlot(node, ...args);
      if (mini && revision === broadRevision) broadMiniCount++;
      events.push({ stage: mini ? 'mini' : node.id, revision, at: performance.now() - started, delivered });
    };
    try {
      await a.setSelection({ components: 'all', methods: [], resolutionMax: null, resolution: 'any', contexts: [], functions: [], subtypes: [], structures: [], puckerStates: [], search: '' });
      await inputFinished;
      await new Promise(resolve => setTimeout(resolve, 40));
    } finally { a.plot = originalPlot; clearInterval(heartbeat); observer.disconnect(); document.removeEventListener('click', deliveredInput, true); document.removeEventListener('click', received); }
    const d = a.snapshots.distribution, j = a.snapshots.joint, s = a.snapshots.survey;
    return { build: a.manifest.build_id, broadRevision, finalRevision: a.revision, familyParameters: a.parameters().length,
      queuedAt, deliveredAt, inputTaskDelayMs: deliveredAt - queuedAt, broadMiniCount, events, heartbeats, longtasks, pendingButtons, trusted,
      final: { state: a.state, status: document.querySelector('#appStatus').dataset.state, mainSelection: d.selection_spec,
        mainRows: d.result.coverage.selectedRows, mainFinite: d.result.coverage.plottedRows, mainBases: [...new Set(d.result.series.flatMap(s => s.rows.map(r => r.comp_id)))],
        jointSelection: j.selection_spec, jointPoints: j.result.points.length, jointBases: [...new Set(j.result.points.map(p => p.left.comp_id))],
        surveySelection: s.selection_spec, surveyFinite: s.result.coverage.plottedRows } };
  }, inputMode);
  report.run = run; await waitReady(page);
  assert.deepEqual(run.final.state.selection.contexts, ['U']); assert.deepEqual(run.final.mainSelection.contexts, ['U']);
  assert.deepEqual(run.final.mainBases, ['U']); assert.deepEqual(run.final.jointSelection.contexts, ['U']); assert.deepEqual(run.final.jointBases, ['U']);
  if (phase === 'after' || inputMode === 'timer') assert(run.pendingButtons.every(button => button.disabled), 'New scientific input disables stale exports immediately');
  assert(run.finalRevision > run.broadRevision); assert(run.final.mainFinite > 0); assert(run.final.jointPoints > 0);
  if (phase === 'after') assert(run.broadMiniCount < run.familyParameters, 'Input cancels obsolete overview before all parameters finish');
  if (inputMode === 'trusted') assert(run.trusted, 'CDP delivers a trusted browser input event');
  report.checks.push({ name: 'Scheduled input reaches latest selection with correctly pending exports', passed: true });
  const expected = await page.evaluate(async () => {
    const a = window.rnaExplorer, table = await a.repository.loadFamily('backbone');
    const map = new Map(table.rows.map(row => [row.id, row]));
    const d = a.snapshots.distribution, j = a.snapshots.joint;
    return { main: d.result.series.flatMap(s => s.rowIds.map(id => ({ id, x: map.get(id).values.chi, base: map.get(id).comp_id }))),
      joint: j.result.points.map(p => ({ id: p.left_id, yid: p.right_id, x: map.get(p.left_id).values.chi, y: map.get(p.right_id).values.alpha, base: map.get(p.left_id).comp_id })),
      mainId: d.snapshot_id, jointId: j.snapshot_id };
  });
  const main = await downloadCsv(page, '#filteredCsvDownload', path.join(output, `${phase}-${inputMode}-main.csv`));
  assert.equal(main.rows.length, expected.main.length);
  main.rows.forEach((row, i) => { const e = expected.main[i]; assert.equal(e.base, 'U'); assert.equal(row.id, e.id); assert.equal(Number(row.value), e.x); assert.equal(row.snapshot_id, expected.mainId); });
  const joint = await downloadCsv(page, '#jointCsvDownload', path.join(output, `${phase}-${inputMode}-joint.csv`));
  assert.equal(joint.rows.length, expected.joint.length);
  joint.rows.forEach((row, i) => { const e = expected.joint[i]; assert.equal(e.base, 'U'); assert.equal(row.x_id, e.id); assert.equal(row.y_id, e.yid); assert.equal(Number(row.x_value), e.x); assert.equal(Number(row.y_value), e.y); assert.equal(row.snapshot_id, expected.jointId); });
  report.checks.push({ name: 'Latest completed CSV retains exact raw identities and values', mainRows: main.rows.length, jointRows: joint.rows.length });
  // Inject failure through the actual completed plotting boundary, then use Reset.
  await page.evaluate(() => {
    const a = window.rnaExplorer, original = a.plot.bind(a); let armed = true;
    window.cooperativeOriginalPlot = original;
    a.plot = async (...args) => { if (armed && args[0].id === 'plot') { armed = false; throw Error('Injected cooperative render failure'); } return original(...args); };
  });
  await page.click('#contextGroup button[data-value="A"]');
  await page.waitForFunction(() => document.querySelector('#appStatus').dataset.state === 'error');
  assert(await page.locator('#filteredCsvDownload').isDisabled());
  await page.click('#resetFilters'); await waitReady(page);
  const reset = await page.evaluate(() => ({ selection: window.rnaExplorer.state.selection, complete: window.rnaExplorer.fullRenderComplete,
    plotted: window.rnaExplorer.snapshots.distribution.result.coverage.plottedRows, exportDisabled: document.querySelector('#filteredCsvDownload').disabled }));
  assert(reset.complete && reset.plotted > 0 && !reset.exportDisabled); assert.deepEqual(reset.selection.contexts, []);
  report.checks.push({ name: 'Injected render failure repairs through actual Reset button', ...reset });
  assert.deepEqual(report.errors, []); report.passed = true;
  console.log(JSON.stringify({ phase, checks: report.checks.length, inputTaskDelayMs: run.inputTaskDelayMs, broadMiniCount: run.broadMiniCount, totalParameters: run.familyParameters, heartbeatCount: run.heartbeats.length, longtasks: run.longtasks }));
} catch (error) { report.passed = false; report.failure = error.stack; throw error; }
finally { report.finished = new Date().toISOString(); await writeFile(path.join(output, `report-${phase}${inputMode === 'timer' ? '' : '-trusted'}.json`), JSON.stringify(report, null, 2)); await browser.close(); }
