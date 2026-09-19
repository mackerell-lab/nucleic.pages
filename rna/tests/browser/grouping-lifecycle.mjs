/** Pair-only grouping cannot leak through old controls into a residue family. */
import assert from 'node:assert/strict';
import { mkdir, readFile, writeFile } from 'node:fs/promises';
import { createHash } from 'node:crypto';
import { pathToFileURL } from 'node:url';
import { waitReady } from './helpers.mjs';

const { chromium } = await import(pathToFileURL(process.env.PLAYWRIGHT_MODULE).href);
const browser = await chromium.launch({ headless: true });
const report = { startedAt: new Date().toISOString(), checks: [], pageErrors: [] };
const sourceHashes = async () => Object.fromEntries(await Promise.all(['app/PureRnaExplorer.js', 'views/panels.js'].map(async file =>
  [file, createHash('sha256').update(await readFile(new URL(`../../${file}`, import.meta.url))).digest('hex')])));
report.sourceBefore = await sourceHashes();
try {
  const page = await browser.newPage();
  page.on('pageerror', error => report.pageErrors.push(error.message));
  await page.goto(process.env.RNA_URL || 'http://127.0.0.1:8767/nucleic.pages/rna/'); await waitReady(page);
  await page.selectOption('#familySelect', 'base_pair'); await waitReady(page);
  const targetUrl = await page.evaluate(() => {
    const app = window.rnaExplorer;
    window.oldGrouping = document.querySelector('#groupingGroup');
    return new URL(app.families.find(family => family.id === 'ribose_2oh').path, app.repository.releaseUrl).href;
  });
  await page.route(targetUrl, route => route.fulfill({ status: 503, body: 'Deliberate grouping regression failure' }), { times: 1 });
  await page.selectOption('#familySelect', 'ribose_2oh');
  await page.waitForFunction(() => document.querySelector('#appStatus').dataset.state === 'error');
  const before = await page.evaluate(() => ({ revision: window.rnaExplorer.revision, grouping: window.rnaExplorer.state.display.groupBy }));
  await page.click('#groupingGroup button[data-value="interactionFamily"]');
  assert.deepEqual(await page.evaluate(() => ({ revision: window.rnaExplorer.revision, grouping: window.rnaExplorer.state.display.groupBy })), before);
  assert.equal(before.grouping, 'base');
  assert.equal(await page.locator('#groupingGroup button[data-value="base"]').getAttribute('aria-pressed'), 'true');
  assert.equal(await page.locator('#groupingGroup button[data-value="interactionFamily"]').getAttribute('aria-pressed'), 'false');
  assert(await page.locator('#filteredCsvDownload').isDisabled());
  report.checks.push('Old connected pair-only grouping is rejected and its appearance rolled back');

  await page.click('#groupingGroup button[data-value="method"]'); await waitReady(page);
  assert(await page.evaluate(() => {
    const app = window.rnaExplorer;
    return app.state.familyId === 'ribose_2oh' && app.state.display.groupBy === 'method'
      && app.snapshots.distribution.display_spec.groupBy === 'method'
      && app.snapshots.distribution.result.coverage.plottedRows > 0;
  }));
  assert(await page.locator('#filteredCsvDownload').isEnabled());
  assert.equal(await page.locator('#groupingGroup button[data-value="interactionFamily"]').count(), 0);
  report.checks.push('Supported old Method grouping retries the residue family and restores exports');

  const detached = await page.evaluate(() => {
    const app = window.rnaExplorer, current = document.querySelector('#groupingGroup');
    const before = { revision: app.revision, grouping: app.state.display.groupBy, snapshot: app.snapshots.distribution.snapshot_id };
    if (window.oldGrouping.isConnected) throw Error('Saved grouping was not detached');
    window.oldGrouping.querySelector('[data-value="interactionFamily"]').click();
    window.oldGrouping.querySelector('[data-value="base"]').click();
    return { before, after: { revision: app.revision, grouping: app.state.display.groupBy, snapshot: app.snapshots.distribution.snapshot_id }, sameDom: current === document.querySelector('#groupingGroup') };
  });
  assert.deepEqual(detached.after, detached.before); assert(detached.sameDom);
  report.checks.push('Detached old specific and portable grouping buttons cannot change current data');

  await page.selectOption('#familySelect', 'base_pair'); await waitReady(page);
  await page.evaluate(() => { window.reusedGrouping = document.querySelector('#groupingGroup'); });
  await page.click('#traceStyleGroup button[data-value="line"]'); await waitReady(page);
  assert(await page.evaluate(() => window.reusedGrouping === document.querySelector('#groupingGroup')));
  await page.click('#groupingGroup button[data-value="interactionFamily"]'); await waitReady(page);
  assert(await page.evaluate(() => window.rnaExplorer.state.display.groupBy === 'interactionFamily'
    && window.rnaExplorer.snapshots.distribution.display_spec.groupBy === 'interactionFamily'));
  await page.click('#resetFilters'); await waitReady(page);
  assert.equal(await page.evaluate(() => window.rnaExplorer.state.display.groupBy), 'base');
  assert.equal(await page.locator('#groupingGroup button[data-value="base"]').getAttribute('aria-pressed'), 'true');
  assert.deepEqual(report.pageErrors, []);
  report.checks.push('Current reused pair grouping remains usable and Reset restores base grouping');
  report.sourceAfter = await sourceHashes();
  assert.deepEqual(report.sourceAfter, report.sourceBefore, 'Source changed during grouping regression');
  report.passed = true;
} catch (error) { report.passed = false; report.failure = error.stack; throw error; }
finally {
  report.finishedAt = new Date().toISOString();
  if (process.env.RNA_BROWSER_OUTPUT) {
    await mkdir(process.env.RNA_BROWSER_OUTPUT, { recursive: true });
    await writeFile(`${process.env.RNA_BROWSER_OUTPUT}/report.json`, JSON.stringify(report, null, 2) + '\n');
  }
  console.log(JSON.stringify(report, null, 2)); await browser.close();
}
