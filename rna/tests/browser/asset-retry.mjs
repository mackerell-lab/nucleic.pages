import assert from 'node:assert/strict';
import { pathToFileURL } from 'node:url';
import { waitReady } from './helpers.mjs';
const { chromium } = await import(pathToFileURL(process.env.PLAYWRIGHT_MODULE).href);
const browser = await chromium.launch({ headless: true });
try {
  const page = await browser.newPage();
  const errors = [];
  page.on('pageerror', error => errors.push(error.message));
  await page.goto('http://127.0.0.1:8767/nucleic.pages/rna/'); await waitReady(page);
  const familyUrl = await page.evaluate(() => {
    const app = window.rnaExplorer;
    return new URL(app.families.find(family => family.id === 'ribose_2oh').path, app.repository.releaseUrl).href;
  });
  await page.route(familyUrl, route => route.fulfill({ status: 503, body: 'Temporary test failure' }), { times: 1 });
  await page.selectOption('#familySelect', 'ribose_2oh');
  await page.waitForFunction(() => document.querySelector('#appStatus').dataset.state === 'error');
  assert(await page.locator('#filteredCsvDownload').isDisabled());
  assert(await page.locator('#plotProvenanceDownload').isDisabled());
  await page.click('#resetFilters'); await waitReady(page);
  assert(await page.evaluate(() => {
    const app = window.rnaExplorer;
    return app.parameters('ribose_2oh').some(parameter => parameter.id === app.snapshots.distribution.result.parameter.id);
  }));
  assert(await page.locator('#filteredCsvDownload').isEnabled());
  console.log('PASS Family fetch failure disables stale exports and Reset retries');

  await page.route('**/survey/scalars/**', route => route.fulfill({ status: 503, body: 'Temporary test failure' }), { times: 1 });
  await page.click('#baseGeometryLoad');
  await page.waitForFunction(() => document.querySelector('#appStatus').dataset.state === 'error');
  assert(await page.locator('#surveyCsvDownload').isDisabled());
  assert(await page.locator('#baseGeometryLoad').isEnabled());
  await page.click('#baseGeometryLoad'); await waitReady(page);
  assert(await page.locator('#surveyCsvDownload').isEnabled());
  assert(await page.evaluate(() => Boolean(window.rnaExplorer.snapshots.survey)));
  assert.deepEqual(errors, []);
  console.log('PASS Survey fetch failure can be retried with the load button');
} finally { await browser.close(); }
