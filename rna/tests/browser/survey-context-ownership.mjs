/** A stale Survey context cannot filter the next term after a failed load. */
import assert from 'node:assert/strict';
import { mkdir, readFile, writeFile } from 'node:fs/promises';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { waitReady } from './helpers.mjs';
const workspace = process.env.RNA_WORKSPACE || '/home/zhaomt/cmap/test15';
const output = process.env.RNA_BROWSER_OUTPUT || path.join(workspace, 'data/pure_rna/browser_validation/survey-context-ownership');
await mkdir(output, { recursive: true });
const { chromium } = await import(pathToFileURL(process.env.PLAYWRIGHT_MODULE || path.join(workspace, 'data/pure_rna/browser-tools/node_modules/playwright-core/index.mjs')).href);
const browser = await chromium.launch({ headless: true }); const page = await browser.newPage();
const report = { started: new Date().toISOString(), checks: [], pageErrors: [] };
page.on('pageerror', error => report.pageErrors.push(error.message));
const state = () => page.evaluate(() => {
  const a = window.rnaExplorer; return { revision: a.revision, term: a.state.survey.termId, contexts: a.state.survey.contexts,
    status: document.querySelector('#appStatus').dataset.state, snapshot: a.snapshots.survey.snapshot_id,
    plotted: a.snapshots.survey.result.coverage.finiteRows, exportDisabled: document.querySelector('#surveyCsvDownload').disabled };
});
try {
  const base = process.env.RNA_URL || 'http://127.0.0.1:8767/nucleic.pages/rna/';
  if (process.env.RNA_HELD_APP) await page.route(base + 'app/PureRnaExplorer.js', async route => route.fulfill({ status: 200, contentType: 'application/javascript', body: await readFile(process.env.RNA_HELD_APP) }));
  await page.goto(base); await waitReady(page); await page.click('#baseGeometryLoad'); await waitReady(page);
  assert.equal((await state()).term, 'c_n1_c2_n3');
  assert(await page.locator('#baseGeometryContextGroup button[data-value="C"]').count());
  const url = await page.evaluate(() => {
    const a = window.rnaExplorer; return new URL(a.manifest.survey.scalars.terms.g_n1_c2_n3.path, a.repository.releaseUrl).href;
  });
  await page.route(url, route => route.fulfill({ status: 503, body: 'Injected Survey scalar load failure' }), { times: 1 });
  await page.selectOption('#baseGeometryTermSelect', 'g_n1_c2_n3');
  await page.waitForFunction(() => document.querySelector('#appStatus').dataset.state === 'error');
  report.before = await state();
  await page.evaluate(() => { window.oldSurveyContext = document.querySelector('#baseGeometryContextGroup button[data-value="C"]'); });
  await page.click('#baseGeometryContextGroup button[data-value="C"]');
  if (process.env.RNA_EXPECT_RED) {
    await waitReady(page); report.after = await state();
    assert.deepEqual(report.after.contexts, ['C']); assert.equal(report.after.plotted, 0);
    report.reproduced = true; report.passed = true;
  } else {
    report.after = await state(); assert.deepEqual(report.after, report.before);
    assert.equal(await page.locator('#baseGeometryContextGroup button[data-value="C"]').getAttribute('aria-pressed'), 'false');
    assert.equal(await page.locator('#baseGeometryContextGroup button[data-all="true"]').getAttribute('aria-pressed'), 'true');
    report.checks.push('Old C context cannot mutate failed G term or reenable stale exports');
    await page.click('#baseGeometryContextGroup button[data-all="true"]'); await waitReady(page);
    const recovered = await state(); assert.equal(recovered.term, 'g_n1_c2_n3'); assert.deepEqual(recovered.contexts, []); assert(recovered.plotted > 0); assert.equal(recovered.exportDisabled, false);
    report.checks.push('Old All contexts safely retries G term with nonzero observations');
    const beforeDetached = await state();
    await page.evaluate(() => { const old = window.oldSurveyContext; old.disabled = false; old.click(); });
    await waitReady(page);
    assert.deepEqual(await state(), beforeDetached);
    assert.equal(await page.locator('#baseGeometryContextGroup button[data-value="G"]').isDisabled(), false);
    report.checks.push('Detached old handler cannot disable the replacement context control');
    await page.click('#baseGeometryContextGroup button[data-value="G"]'); await waitReady(page);
    assert.deepEqual((await state()).contexts, ['G']); report.checks.push('Current G context remains usable');
    await page.click('#traceStyleGroup button[data-value="line"]'); await waitReady(page);
    await page.click('#baseGeometryContextGroup button[data-value="G"]'); await waitReady(page);
    assert.deepEqual((await state()).contexts, []); report.checks.push('Same-term reused context remains usable after trace display revision');
    // Group changes clear term ID before the new data arrives. The old specific
    // control must not populate contexts in that transitional state either.
    const target = await page.evaluate(() => {
      const a = window.rnaExplorer, term = a.surveyTerms().find(t => t.group !== 'base_internal_angles');
      return { group: term.group, id: term.id, url: new URL(a.manifest.survey.scalars.terms[term.id].path, a.repository.releaseUrl).href };
    });
    await page.route(target.url, route => route.fulfill({ status: 503, body: 'Injected Survey group load failure' }), { times: 1 });
    await page.selectOption('#surveyGroupSelect', target.group); await page.waitForFunction(() => document.querySelector('#appStatus').dataset.state === 'error');
    const beforeGroup = await state(); await page.click('#baseGeometryContextGroup button[data-value="G"]'); assert.deepEqual(await state(), beforeGroup);
    await page.click('#baseGeometryContextGroup button[data-all="true"]'); await waitReady(page);
    assert.equal((await state()).term, target.id); assert.deepEqual((await state()).contexts, []);
    report.checks.push('Group-change failure rejects old context and All safely recovers resolved term');
    assert.deepEqual(report.pageErrors, []); report.passed = true;
  }
} catch (error) { report.failure = error.stack; report.passed = false; throw error; }
finally { report.finished = new Date().toISOString(); await writeFile(path.join(output, 'report.json'), JSON.stringify(report, null, 2)); console.log(JSON.stringify(report, null, 2)); await browser.close(); }
