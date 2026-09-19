/** Old family overview cards cannot corrupt a failed or pending selection. */
import assert from 'node:assert/strict';
import { mkdir, writeFile } from 'node:fs/promises';
import { pathToFileURL } from 'node:url';
import { waitReady } from './helpers.mjs';

const { chromium } = await import(pathToFileURL(process.env.PLAYWRIGHT_MODULE).href);
const browser = await chromium.launch({ headless: true });
const report = { startedAt: new Date().toISOString(), checks: [], pageErrors: [] };
try {
  const page = await browser.newPage();
  page.on('pageerror', error => report.pageErrors.push(error.message));
  await page.goto(process.env.RNA_URL || 'http://127.0.0.1:8767/nucleic.pages/rna/');
  await waitReady(page);
  const familyUrl = await page.evaluate(() => {
    const app = window.rnaExplorer;
    return new URL(app.families.find(family => family.id === 'ribose_2oh').path, app.repository.releaseUrl).href;
  });
  await page.route(familyUrl, route => route.fulfill({ status: 503, body: 'Temporary regression-test failure' }), { times: 1 });
  await page.selectOption('#familySelect', 'ribose_2oh');
  await page.waitForFunction(() => document.querySelector('#appStatus').dataset.state === 'error');
  const before = await page.evaluate(() => {
    const app = window.rnaExplorer;
    return { revision: app.revision, family: app.state.familyId, parameter: app.state.parameterId,
      snapshot: app.snapshots.distribution.snapshot_id };
  });
  assert.equal(before.family, 'ribose_2oh');
  assert(await page.locator('#filteredCsvDownload').isDisabled());
  assert(await page.locator('#plotProvenanceDownload').isDisabled());
  await page.click('#familyOverview [data-parameter="beta"]');
  const after = await page.evaluate(() => {
    const app = window.rnaExplorer;
    return { revision: app.revision, family: app.state.familyId, parameter: app.state.parameterId,
      snapshot: app.snapshots.distribution.snapshot_id };
  });
  report.staleClick = { before, after };
  assert.deepEqual(after, before, 'Stale overview card mutated the current family selection');
  assert.deepEqual(report.pageErrors, []);
  assert(await page.locator('#filteredCsvDownload').isDisabled());
  report.checks.push('Failed-family stale overview click leaves selection and exports unchanged');

  await page.click('#resetFilters');
  await waitReady(page);
  report.recovered = await page.evaluate(() => {
    const app = window.rnaExplorer;
    return { family: app.state.familyId, parameter: app.state.parameterId,
      plotted: app.snapshots.distribution.result.parameter.id, buildId: app.manifest.build_id };
  });
  assert.equal(report.recovered.family, 'ribose_2oh');
  assert.equal(report.recovered.plotted, report.recovered.parameter);
  assert(await page.locator('#filteredCsvDownload').isEnabled());
  assert(await page.locator('#plotProvenanceDownload').isEnabled());
  report.checks.push('Reset retries the failed family and restores matching plot and exports');

  await page.selectOption('#familySelect', 'backbone');
  await waitReady(page);
  await page.evaluate(() => {
    window.staleOverviewCurrentCard = document.querySelector('#familyOverview [data-parameter="beta"]');
    window.staleOverviewReusedCard = document.querySelector('#familyOverview [data-parameter="alpha"]');
  });
  await page.click('#familyOverview [data-parameter="beta"]');
  await waitReady(page);
  assert(await page.evaluate(() => {
    const app = window.rnaExplorer;
    return app.state.parameterId === 'beta' && app.snapshots.distribution.result.parameter.id === 'beta'
      && window.staleOverviewCurrentCard === document.querySelector('#familyOverview [data-parameter="beta"]')
      && window.staleOverviewCurrentCard.getAttribute('aria-pressed') === 'true';
  }));
  assert.deepEqual(report.pageErrors, []);
  report.checks.push('Current-family overview click still selects beta and reuses its card');
  await page.click('#familyOverview [data-parameter="alpha"]');
  await waitReady(page);
  assert(await page.evaluate(() => {
    const app = window.rnaExplorer;
    return app.state.parameterId === 'alpha' && app.snapshots.distribution.result.parameter.id === 'alpha'
      && window.staleOverviewReusedCard === document.querySelector('#familyOverview [data-parameter="alpha"]')
      && window.staleOverviewReusedCard.getAttribute('aria-pressed') === 'true';
  }));
  assert.deepEqual(report.pageErrors, []);
  report.checks.push('Reused same-family card remains live after another completed revision');
  report.passed = true;
} catch (error) {
  report.passed = false; report.failure = error.stack; throw error;
} finally {
  report.finishedAt = new Date().toISOString();
  if (process.env.RNA_BROWSER_OUTPUT) {
    await mkdir(process.env.RNA_BROWSER_OUTPUT, { recursive: true });
    await writeFile(`${process.env.RNA_BROWSER_OUTPUT}/result.json`, JSON.stringify(report, null, 2) + '\n');
  }
  console.log(JSON.stringify(report, null, 2));
  await browser.close();
}
