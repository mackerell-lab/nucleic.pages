/** Real secondary selectors must keep the completed main/Survey marker current. */
import assert from 'node:assert/strict';
import { mkdir, readFile, writeFile } from 'node:fs/promises';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { createHash } from 'node:crypto';
import { waitReady } from './helpers.mjs';
const workspace = process.env.RNA_WORKSPACE || '/home/zhaomt/cmap/test15';
const output = process.env.RNA_BROWSER_OUTPUT || path.join(workspace, 'data/pure_rna/browser_validation/secondary-selection-trace-reuse');
const phase = process.env.RNA_TIMING_PHASE || 'after';
await mkdir(output, { recursive: true });
const { chromium } = await import(pathToFileURL(process.env.PLAYWRIGHT_MODULE || path.join(workspace, 'data/pure_rna/browser-tools/node_modules/playwright-core/index.mjs')).href);
const held = path.join(output, 'baseline-source');
const source = await readFile(phase === 'before' ? path.join(held, 'app/PureRnaExplorer.js') : new URL('../../app/PureRnaExplorer.js', import.meta.url));
const browser = await chromium.launch({ headless: true }); const page = await browser.newPage({ viewport: { width: 1440, height: 1000 } });
const report = { started: new Date().toISOString(), phase, sourceSha256: createHash('sha256').update(source).digest('hex'), measurements: [], errors: [], limitation: 'Three sequential samples on shared host. Timer surrounds real Playwright click and ready-state wait, not browser paint completion; no general speed claim.' };
page.on('pageerror', error => report.errors.push(error.message));
page.on('console', message => { if (message.type() === 'error') report.errors.push(message.text()); });
page.on('requestfailed', request => report.errors.push(request.failure()?.errorText));
const digest = () => page.evaluate(async () => {
  const a = window.rnaExplorer, series = result => result.series.map(s => ({ key: s.key, rowIds: s.rowIds, values: s.values, weights: s.weights, x: s.x, y: s.y, counts: s.counts, statistics: s.statistics }));
  const j = a.snapshots.joint?.result;
  const bytes = new TextEncoder().encode(JSON.stringify({ main: series(a.snapshots.distribution.result), survey: series(a.snapshots.survey.result), joint: j && { points: j.points.map(p => [p.left_id, p.right_id, p.x, p.y]), x: j.x, y: j.y, z: j.z, statistics: j.statistics } }));
  return [...new Uint8Array(await crypto.subtle.digest('SHA-256', bytes))].map(x => x.toString(16).padStart(2, '0')).join('');
});
try {
  const base = process.env.RNA_URL || 'http://127.0.0.1:8767/nucleic.pages/rna/';
  await page.route(base + '**/*.js', async route => {
    const relative = new URL(route.request().url()).pathname.slice(new URL(base).pathname.length);
    const body = relative === 'app/PureRnaExplorer.js' ? source : await readFile(path.join(held, relative));
    report.modules ??= {}; report.modules[relative] = createHash('sha256').update(body).digest('hex');
    await route.fulfill({ status: 200, contentType: 'application/javascript', body });
  });
  await page.goto(base); await waitReady(page);
  await page.evaluate(() => window.rnaExplorer.setSelection({ components: 'all', methods: [], resolutionMax: null, resolution: 'any', contexts: [], functions: [], subtypes: [], structures: [], puckerStates: [], search: '' }));
  await waitReady(page); await page.click('#baseGeometryLoad'); await waitReady(page);
  report.setup = await page.evaluate(() => ({ build: window.rnaExplorer.manifest.build_id, state: window.rnaExplorer.state, coverage: window.rnaExplorer.snapshots.distribution.result.coverage }));
  for (const [index, parameter] of ['alpha', 'beta', 'delta'].entries()) {
    await page.evaluate(() => { window.preSecondaryMain = window.rnaExplorer.snapshots.distribution.result; });
    const selectorStart = performance.now();
    if (index === 0) { await page.selectOption('#family2Select', 'backbone'); await waitReady(page); }
    await page.selectOption('#parameter2Select', parameter); await waitReady(page);
    const selectorMs = performance.now() - selectorStart;
    const selectorReusedMain = await page.evaluate(() => window.preSecondaryMain === window.rnaExplorer.snapshots.distribution.result);
    assert.equal(selectorReusedMain, phase === 'after');
    const beforeDigest = await digest();
    await page.evaluate(() => {
      const a = window.rnaExplorer;
      window.secondaryBefore = { revision: a.revision, main: a.snapshots.distribution.result, survey: a.snapshots.survey.result, joint: a.snapshots.joint.result,
        overview: [...document.querySelectorAll('#familyOverview [data-parameter]')], loads: 0 };
      if (!window.secondaryInstrumented) {
        for (const name of ['loadFamily', 'loadRelations', 'loadSurveyScalars']) {
          const original = a.repository[name].bind(a.repository);
          a.repository[name] = (...args) => { window.secondaryBefore.loads++; return original(...args); };
        }
        window.secondaryInstrumented = true;
      }
    });
    const style = index % 2 ? 'filled' : 'line'; const start = performance.now();
    await page.click(`#traceStyleGroup button[data-value="${style}"]`); await waitReady(page);
    const ms = performance.now() - start;
    const evidence = await page.evaluate(() => {
      const a = window.rnaExplorer, b = window.secondaryBefore;
      return { reused: a.snapshots.distribution.result.series === b.main.series && a.snapshots.survey.result.series === b.survey.series && a.snapshots.joint.result.points === b.joint.points && a.snapshots.joint.result.z === b.joint.z,
        overview: b.overview.every((node, i) => node === document.querySelectorAll('#familyOverview [data-parameter]')[i]),
        loads: b.loads, revisionDelta: a.revision - b.revision, parameter: a.snapshots.joint.result.yParameter.id,
        points: a.snapshots.joint.result.points.length, style: a.state.display.traceStyle,
        metadata: Object.values(a.snapshots).every(s => s.display_spec.traceStyle === a.state.display.traceStyle), keyCurrent: a.completedTraceKey === a.traceAnalysisKey(a.state) };
    });
    const afterDigest = await digest(); assert.equal(beforeDigest, afterDigest); assert.equal(evidence.parameter, parameter); assert.equal(evidence.revisionDelta, 1); assert(evidence.metadata); assert(evidence.keyCurrent);
    if (phase === 'after') { assert(evidence.reused); assert(evidence.overview); assert.equal(evidence.loads, 0); }
    else { assert(evidence.reused); assert.equal(evidence.loads, 0); }
    report.measurements.push({ index, parameter, selectorMs, selectorReusedMain, ms, beforeDigest, afterDigest, ...evidence });
    console.log(JSON.stringify({ parameter, selectorMs, selectorReusedMain, ms, ...evidence }));
  }
  assert.deepEqual(report.errors, []); report.passed = true;
} catch (error) { report.failure = error.stack; report.passed = false; throw error; }
finally { report.finished = new Date().toISOString(); await writeFile(path.join(output, `${phase}.json`), JSON.stringify(report, null, 2)); await browser.close(); }
