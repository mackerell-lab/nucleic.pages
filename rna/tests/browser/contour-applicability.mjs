/** Contour controls apply only to actual contour traces and retain preferences. */
import assert from 'node:assert/strict';
import { mkdir, writeFile } from 'node:fs/promises';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { waitReady } from './helpers.mjs';
const workspace = path.resolve(process.env.RNA_WORKSPACE || process.cwd());
const output = path.resolve(process.env.RNA_BROWSER_OUTPUT || path.join(workspace, 'data/pure_rna/browser_validation/contour-applicability'));
await mkdir(output, { recursive: true });
const { chromium } = await import(pathToFileURL(process.env.PLAYWRIGHT_MODULE || path.join(workspace, 'data/pure_rna/browser-tools/node_modules/playwright-core/index.mjs')).href);
const browser = await chromium.launch({ headless: true });
const page = await browser.newPage({ viewport: { width: 1440, height: 1000 } });
const report = { startedAt: new Date().toISOString(), checks: [], pageErrors: [], failedRequests: [] };
page.on('pageerror', error => report.pageErrors.push(error.message));
page.on('requestfailed', request => report.failedRequests.push(request.url()));
try {
  await page.goto(process.env.RNA_URL || 'http://127.0.0.1:8767/nucleic.pages/rna/', { waitUntil: 'domcontentloaded', timeout: 120000 });
  await waitReady(page);
  await page.selectOption('#family2Select', 'backbone'); await waitReady(page);
  await page.selectOption('#parameter2Select', 'delta'); await waitReady(page);
  await page.evaluate(() => {
    const app = window.rnaExplorer;
    window.contourBaseline = { distribution: app.snapshots.distribution, points: app.snapshots.joint.result.points.map(p => [p.left_id, p.right_id, p.x, p.y]) };
  });
  async function check(type, labels, count) {
    const evidence = await page.evaluate(() => {
      const app = window.rnaExplorer, groups = ['jointContourLabelsGroup', 'jointContourWidthGroup'].map(id => document.getElementById(id));
      const contours = document.getElementById('jointPlot').data.filter(trace => trace.type === 'contour');
      return { type: app.state.joint.type, labels: app.state.joint.labels, count: app.state.joint.contourCount,
        groups: groups.map(group => ({ hidden: group.parentElement.hidden, disabled: [...group.querySelectorAll('button')].every(button => button.disabled), active: group.querySelector('[aria-pressed="true"]')?.dataset.value })),
        contours: contours.map(trace => ({ labels: trace.contours.showlabels, count: trace.ncontours })),
        distributionUnchanged: app.snapshots.distribution === window.contourBaseline.distribution,
        pointsUnchanged: JSON.stringify(app.snapshots.joint.result.points.map(p => [p.left_id, p.right_id, p.x, p.y])) === JSON.stringify(window.contourBaseline.points),
        points: app.snapshots.joint.result.points.length };
    });
    assert.equal(evidence.type, type); assert.equal(evidence.labels, labels); assert.equal(evidence.count, count);
    for (const group of evidence.groups) { assert.equal(group.hidden, type === 'heatmap'); assert.equal(group.disabled, type === 'heatmap'); }
    assert.equal(evidence.groups[0].active, labels ? 'on' : 'off'); assert.equal(evidence.groups[1].active, String(count));
    assert.equal(evidence.contours.length, type === 'heatmap' ? 0 : 1);
    for (const contour of evidence.contours) { assert.equal(contour.labels, labels); assert.equal(contour.count, count); }
    assert(evidence.distributionUnchanged); assert(evidence.pointsUnchanged); assert(evidence.points > 0);
    report.checks.push({ passed: true, ...evidence });
  }
  await check('heatmap', false, 12);
  const inactive = await page.evaluate(() => {
    const app = window.rnaExplorer, revision = app.revision, snapshot = app.snapshots.joint;
    document.querySelector('#jointContourLabelsGroup button[data-value="on"]').click();
    document.querySelector('#jointContourWidthGroup button[data-value="24"]').click();
    return { revisionUnchanged: revision === app.revision, snapshotUnchanged: snapshot === app.snapshots.joint, labels: app.state.joint.labels, count: app.state.joint.contourCount };
  });
  assert.deepEqual(inactive, { revisionUnchanged: true, snapshotUnchanged: true, labels: false, count: 12 });
  report.checks.push({ name: 'Inactive disabled buttons do not render or mutate preferences', passed: true, ...inactive });
  await page.click('#jointPlotTypeGroup button[data-value="contour"]'); await waitReady(page);
  await page.click('#jointContourLabelsGroup button[data-value="on"]'); await waitReady(page);
  await page.click('#jointContourWidthGroup button[data-value="24"]'); await waitReady(page);
  await check('contour', true, 24);
  for (const type of ['heatmap', 'contour', 'filled_contour', 'heatmap_contour', 'heatmap', 'filled_contour']) {
    await page.click(`#jointPlotTypeGroup button[data-value="${type}"]`); await waitReady(page);
    await check(type, true, 24);
  }
  assert.deepEqual(report.pageErrors, []); assert.deepEqual(report.failedRequests, []);
  report.passed = true;
  console.log(`PASS ${report.checks.length} contour applicability checks`);
} catch (error) { report.passed = false; report.failure = error.stack; throw error; }
finally { report.finishedAt = new Date().toISOString(); await writeFile(path.join(output, 'rna-contour-applicability.json'), JSON.stringify(report, null, 2) + '\n'); await browser.close(); }
