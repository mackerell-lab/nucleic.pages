/** Queued actual input interrupts ownership transfer of a real full-release result. */
import assert from 'node:assert/strict';
import { mkdir, readFile, writeFile } from 'node:fs/promises';
import { createHash } from 'node:crypto';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { downloadCsv, waitReady } from './helpers.mjs';

const workspace = path.resolve(process.env.RNA_WORKSPACE || process.cwd());
const output = process.env.RNA_BROWSER_OUTPUT || path.join(workspace, 'data/pure_rna/browser_validation/cooperative-snapshot');
const phase = process.env.RNA_SNAPSHOT_PHASE || 'after', heldApp = process.env.RNA_SNAPSHOT_BASELINE_APP;
const base = process.env.RNA_URL || 'http://127.0.0.1:8767/nucleic.pages/rna/';
await mkdir(output, { recursive: true });
const { chromium } = await import(pathToFileURL(process.env.PLAYWRIGHT_MODULE || path.join(workspace, 'data/pure_rna/browser-tools/node_modules/playwright-core/index.mjs')).href);
const browser = await chromium.launch({ headless: true });
const page = await browser.newPage({ viewport: { width: 1440, height: 1000 }, acceptDownloads: true });
const report = { phase, started: new Date().toISOString(), sourceSha256: {}, checks: [], errors: [], limitation: 'One matched browser task-delivery experiment on a real result graph. No fabricated graph or substituted repository. A checkpoint time budget does not bound indivisible work or guarantee universal latency.' };
for (const file of ['app/PureRnaExplorer.js', 'app/NucleicAcidExplorer.js', 'core/export.js', 'core/repository.js', 'core/cooperative-freeze.js']) {
  const bytes = await readFile(file === 'app/PureRnaExplorer.js' && heldApp ? heldApp : new URL('../../' + file, import.meta.url));
  report.sourceSha256[file] = createHash('sha256').update(bytes).digest('hex');
}
page.on('pageerror', error => report.errors.push(error.message));
page.on('console', message => { if (message.type() === 'error' && !message.text().includes('Injected snapshot ownership failure')) report.errors.push(message.text()); });
if (heldApp) await page.route(base + 'app/PureRnaExplorer.js', async route => route.fulfill({ status: 200, contentType: 'application/javascript', body: await readFile(heldApp) }));
try {
  await page.goto(base); await waitReady(page);
  await page.selectOption('#family2Select', 'backbone'); await waitReady(page);
  await page.selectOption('#parameter2Select', 'alpha'); await waitReady(page);
  const run = await page.evaluate(async () => {
    const a = window.rnaExplorer, originalSnapshot = a.snapshot.bind(a), originalSnapshots = a.snapshots;
    const started = performance.now(), events = [], publications = [], beats = [], longtasks = [];
    const obsoleteRevision = a.revision + 1;
    let queued = false, delivered = false, queuedAt, deliveredAt, obsoleteOutcome, obsoleteId, pending;
    const observer = new PerformanceObserver(list => { for (const e of list.getEntries()) longtasks.push({ start: e.startTime - started, duration: e.duration }); });
    observer.observe({ type: 'longtask' }); const heartbeat = setInterval(() => beats.push(performance.now() - started), 25);
    let resolveLatest; const latestFinished = new Promise(resolve => { resolveLatest = resolve; });
    a.snapshots = new Proxy(originalSnapshots, { set(target, key, value) {
      publications.push({ key, at: performance.now() - started, delivered, id: value?.snapshot_id, contexts: value?.selection_spec?.contexts, selectedRows: value?.result?.coverage?.selectedRows });
      target[key] = value; return true;
    } });
    a.snapshot = function (options) {
      const old = options.revision === obsoleteRevision && options.result.kind === 'distribution';
      if (old && !queued) {
        queued = true; queuedAt = performance.now() - started;
        events.push({ stage: 'obsolete-snapshot-start', at: queuedAt, rows: options.result.coverage.selectedRows });
        setTimeout(async () => {
          delivered = true; deliveredAt = performance.now() - started;
          events.push({ stage: 'queued-context-input', at: deliveredAt, obsoleteOutcome });
          document.querySelector('#contextGroup button[data-value="U"]').click();
          pending = ['filteredCsvDownload', 'plotProvenanceDownload', 'jointCsvDownload'].map(id => ({ id, disabled: document.getElementById(id).disabled }));
          while (document.querySelector('#appStatus').dataset.state === 'loading') await new Promise(resolve => setTimeout(resolve, 10));
          resolveLatest();
        }, 0);
      }
      const completed = snapshot => {
        if (old) { obsoleteOutcome = 'completed'; obsoleteId = snapshot?.snapshot_id; events.push({ stage: 'obsolete-snapshot-completed', at: performance.now() - started, delivered }); }
        return snapshot;
      };
      const failed = error => { if (old) { obsoleteOutcome = error.name; events.push({ stage: 'obsolete-snapshot-aborted', name: error.name, at: performance.now() - started, delivered }); } throw error; };
      try { const value = originalSnapshot(options); return value?.then ? value.then(completed, failed) : completed(value); }
      catch (error) { return failed(error); }
    };
    try {
      await a.setSelection({ components: 'all', methods: [], resolutionMax: null, resolution: 'any', contexts: [], functions: [], subtypes: [], structures: [], puckerStates: [], search: '' });
      await latestFinished; await a.commitQueue; await new Promise(resolve => setTimeout(resolve, 40));
    } finally { a.snapshot = originalSnapshot; a.snapshots = originalSnapshots; clearInterval(heartbeat); observer.disconnect(); }
    return { build: a.manifest.build_id, obsoleteRevision, finalRevision: a.revision, inputDelayMs: deliveredAt - queuedAt, obsoleteOutcome,
      obsoleteId, pending, publications, events, beats, longtasks, obsoletePublished: Boolean(obsoleteId && publications.some(p => p.id === obsoleteId)),
      final: { state: a.state, complete: a.fullRenderComplete, main: { id: a.snapshots.distribution.snapshot_id, selection: a.snapshots.distribution.selection_spec, count: a.snapshots.distribution.result.coverage.plottedRows },
        joint: { id: a.snapshots.joint.snapshot_id, selection: a.snapshots.joint.selection_spec, count: a.snapshots.joint.result.points.length },
        frozen: Object.isFrozen(a.snapshots.distribution) && Object.isFrozen(a.snapshots.distribution.result) && Object.isFrozen(a.snapshots.joint.result.points) } };
  });
  report.run = run; await waitReady(page);
  assert(run.events.some(e => e.rows === 242078)); assert(run.final.complete && run.final.frozen);
  assert.deepEqual(run.final.main.selection.contexts, ['U']); assert.deepEqual(run.final.joint.selection.contexts, ['U']);
  assert(run.pending.every(button => button.disabled)); assert(run.finalRevision > run.obsoleteRevision);
  assert(!run.publications.some(p => p.delivered && p.selectedRows === 242078), 'Superseded broad result never publishes after latest input');
  if (phase === 'after') { assert.equal(run.obsoleteOutcome, 'AbortError'); assert.equal(run.obsoletePublished, false); }
  record('Actual input supersedes broad snapshot construction without late publication', { inputDelayMs: run.inputDelayMs, obsoleteOutcome: run.obsoleteOutcome, obsoletePublished: run.obsoletePublished });
  const source = await page.evaluate(async () => {
    const a = window.rnaExplorer, table = await a.repository.loadFamily('backbone');
    const valid = (row, parameter) => {
      const status = row.statuses?.[parameter], code = typeof status === 'object' ? status.code : status;
      return (!code || ['ok', 'available', 'valid', 'computed'].includes(code)) && Number.isFinite(row.values?.[parameter]);
    };
    const rows = table.rows.filter(row => row.comp_id === 'U');
    return { main: rows.filter(row => valid(row, 'chi')).map(row => ({ id: row.id, value: row.values.chi })),
      joint: rows.filter(row => valid(row, 'chi') && valid(row, 'alpha')).map(row => ({ id: row.id, x: row.values.chi, y: row.values.alpha })) };
  });
  const main = await downloadCsv(page, '#filteredCsvDownload', path.join(output, `${phase}-main.csv`));
  const expectedMain = new Map(source.main.map(row => [row.id, row.value]));
  assert.equal(main.rows.length, source.main.length); assert.equal(new Set(main.rows.map(row => row.id)).size, source.main.length);
  main.rows.forEach(row => { assert(expectedMain.has(row.id)); assert.equal(Number(row.value), expectedMain.get(row.id)); assert.equal(row.snapshot_id, run.final.main.id); });
  const joint = await downloadCsv(page, '#jointCsvDownload', path.join(output, `${phase}-joint.csv`));
  const expectedJoint = new Map(source.joint.map(row => [row.id, row]));
  assert.equal(joint.rows.length, source.joint.length); assert.equal(new Set(joint.rows.map(row => row.x_id)).size, source.joint.length);
  joint.rows.forEach(row => { const e = expectedJoint.get(row.x_id); assert(e); assert.equal(row.y_id, e.id); assert.equal(Number(row.x_value), e.x); assert.equal(Number(row.y_value), e.y); assert.equal(row.snapshot_id, run.final.joint.id); });
  record('Exact latest CSV against independent complete source subsets', { mainRows: main.rows.length, jointRows: joint.rows.length });
  await page.evaluate(() => {
    const a = window.rnaExplorer, original = a.snapshot.bind(a); let armed = true;
    a.snapshot = function (options) { if (armed) { armed = false; throw Error('Injected snapshot ownership failure'); } return original(options); };
  });
  await page.click('#contextGroup button[data-value="A"]'); await page.waitForFunction(() => document.querySelector('#appStatus').dataset.state === 'error');
  assert(await page.locator('#filteredCsvDownload').isDisabled());
  await page.click('#resetFilters'); await waitReady(page);
  const reset = await page.evaluate(() => ({ complete: window.rnaExplorer.fullRenderComplete, contexts: window.rnaExplorer.state.selection.contexts,
    finite: window.rnaExplorer.snapshots.distribution.result.coverage.plottedRows, frozen: Object.isFrozen(window.rnaExplorer.snapshots.distribution) }));
  assert(reset.complete && reset.finite > 0 && reset.frozen); assert.deepEqual(reset.contexts, []); assert(!(await page.locator('#filteredCsvDownload').isDisabled()));
  record('Snapshot failure recovers through actual Reset control', reset);
  assert.deepEqual(report.errors, []); report.passed = true;
  console.log(JSON.stringify({ phase, checks: report.checks.length, inputDelayMs: run.inputDelayMs, outcome: run.obsoleteOutcome, obsoletePublished: run.obsoletePublished, mainRows: main.rows.length, jointRows: joint.rows.length }));
} catch (error) { report.passed = false; report.failure = error.stack; throw error; }
finally { report.finished = new Date().toISOString(); await writeFile(path.join(output, `report-${phase}.json`), JSON.stringify(report, null, 2)); await browser.close(); }
function record(name, evidence) { report.checks.push({ name, ...evidence }); }
