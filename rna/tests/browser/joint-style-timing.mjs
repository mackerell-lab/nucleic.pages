/** Matched real-browser timings; retained source directory pins the before run. */
import assert from 'node:assert/strict';
import { mkdir, readFile, writeFile } from 'node:fs/promises';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { createHash } from 'node:crypto';
import { waitReady } from './helpers.mjs';

const workspace = path.resolve(process.env.RNA_WORKSPACE || process.cwd());
const output = process.env.RNA_BROWSER_OUTPUT || path.join(workspace, 'data/pure_rna/browser_validation/joint-style-reuse');
const phase = process.env.RNA_TIMING_PHASE || 'after';
const source = process.env.RNA_HELD_SOURCE;
const base = process.env.RNA_URL || 'http://127.0.0.1:8767/nucleic.pages/rna/';
await mkdir(output, { recursive: true });
const { chromium } = await import(pathToFileURL(process.env.PLAYWRIGHT_MODULE || path.join(workspace, 'data/pure_rna/browser-tools/node_modules/playwright-core/index.mjs')).href);
const browser = await chromium.launch({ headless: true });
const page = await browser.newPage({ viewport: { width: 1440, height: 1000 } });
const report = { phase, started: new Date().toISOString(), errors: [], measurements: [], sourceSha256: {}, limitation: 'Three warm measurements on a shared host. Matched browser, viewport, build, selection, parameters, histogram settings and palette sequence; not a universal speed claim.' };
for (const name of ['app/PureRnaExplorer.js', 'core/analysis.js', 'core/selection.js', 'core/joints.js', 'core/export.js']) {
  report.sourceSha256[name] = createHash('sha256').update(await readFile(source ? path.join(source, name) : new URL('../../' + name, import.meta.url))).digest('hex');
}
page.on('pageerror', error => report.errors.push(error.message));
if (source) await page.route(base + '**/*.js', async route => {
  const relative = new URL(route.request().url()).pathname.slice(new URL(base).pathname.length);
  await route.fulfill({ status: 200, contentType: 'application/javascript', body: await readFile(path.join(source, relative)) });
});
try {
  await page.goto(base); await waitReady(page);
  await page.evaluate(async () => {
    const a = window.rnaExplorer;
    await a.setSelection({ components: 'all', methods: [], resolutionMax: null, resolution: 'any', contexts: [], functions: [], subtypes: [], structures: [], puckerStates: [], search: '' });
    a.state.family2Id = 'backbone'; a.state.parameter2Id = 'alpha'; a.updateSelectors();
    await a.requestJointOnly();
  });
  await waitReady(page);
  report.setup = await page.evaluate(() => {
    const a = window.rnaExplorer;
    window.styleTimingLoads = 0;
    for (const name of ['loadFamily', 'loadRelations', 'loadMetadata']) {
      const original = a.repository[name].bind(a.repository);
      a.repository[name] = (...args) => { window.styleTimingLoads++; return original(...args); };
    }
    return { build: a.manifest.build_id, state: a.state, coverage: a.snapshots.distribution.result.coverage, points: a.snapshots.joint.result.points.length, grid: [a.snapshots.joint.result.x.length, a.snapshots.joint.result.y.length] };
  });
  assert.equal(report.setup.build, 'full_packed_family_20260919');
  for (const [index, palette] of ['warm', 'viridis', 'cividis', 'ocean'].entries()) {
    const measurement = await page.evaluate(async ({ palette, index }) => {
      const a = window.rnaExplorer, previous = a.snapshots.joint.result;
      const loads = window.styleTimingLoads, start = performance.now();
      a.state.joint.palette = palette; await a.requestJointOnly();
      return { index, palette, ms: performance.now() - start, loads: window.styleTimingLoads - loads, reused: previous === a.snapshots.joint.result, points: a.snapshots.joint.result.points.length };
    }, { palette, index });
    await waitReady(page);
    if (phase === 'after') { assert(measurement.reused); assert.equal(measurement.loads, 0); }
    report.measurements.push({ ...measurement, warmup: index === 0 });
    console.log(JSON.stringify(measurement));
    await writeFile(path.join(output, `timing-${phase}.json`), JSON.stringify(report, null, 2));
  }
  report.numericDigest = await page.evaluate(async () => {
    const r = window.rnaExplorer.snapshots.joint.result;
    const bytes = new TextEncoder().encode(JSON.stringify({ points: r.points.map(p => [p.left_id, p.right_id, p.x, p.y]), x: r.x, y: r.y, z: r.z, statistics: r.statistics }));
    return [...new Uint8Array(await crypto.subtle.digest('SHA-256', bytes))].map(x => x.toString(16).padStart(2, '0')).join('');
  });
  assert.deepEqual(report.errors, []); report.passed = true;
} catch (error) { report.failure = error.stack; report.passed = false; throw error; }
finally { report.finished = new Date().toISOString(); await writeFile(path.join(output, `timing-${phase}.json`), JSON.stringify(report, null, 2)); await browser.close(); }
