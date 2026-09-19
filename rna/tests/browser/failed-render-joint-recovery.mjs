/** A joint control must repair a failed full refresh before declaring ready. */
import assert from 'node:assert/strict';
import {mkdir, writeFile} from 'node:fs/promises';
import {pathToFileURL} from 'node:url';
import {waitReady} from './helpers.mjs';

const {chromium} = await import(pathToFileURL(process.env.PLAYWRIGHT_MODULE).href);
const browser = await chromium.launch({headless: true});
try {
  const page = await browser.newPage();
  const errors = [];
  page.on('pageerror', error => errors.push(error.message));
  await page.goto(process.env.RNA_URL || 'http://127.0.0.1:8767/nucleic.pages/rna/');
  await waitReady(page);
  const familyUrl = await page.evaluate(() => {
    const app = window.rnaExplorer;
    return new URL(app.families.find(family => family.id === 'ribose_2oh').path, app.repository.releaseUrl).href;
  });
  await page.route(familyUrl, route => route.fulfill({status: 503, body: 'Temporary regression-test failure'}), {times: 1});
  await page.selectOption('#familySelect', 'ribose_2oh');
  await page.waitForFunction(() => document.querySelector('#appStatus').dataset.state === 'error');
  assert(await page.locator('#filteredCsvDownload').isDisabled());
  assert.equal(await page.evaluate(() => window.rnaExplorer.snapshots.distribution.result.parameter.id), 'chi');

  await page.click('#jointPaletteGroup button[data-value="viridis"]');
  await waitReady(page);
  const recovered = await page.evaluate(() => {
    const app = window.rnaExplorer;
    return {selectedParameter: app.state.parameterId, plottedParameter: app.snapshots.distribution.result.parameter.id,
      snapshotId: app.snapshots.distribution.snapshot_id, family: app.state.familyId, buildId: app.manifest.build_id};
  });
  assert.equal(recovered.plottedParameter, recovered.selectedParameter,
    'Joint click declared ready while the main plot retained the failed family selection');
  assert(await page.locator('#filteredCsvDownload').isEnabled(), 'Full recovery left the primary CSV disabled');
  assert(await page.locator('#plotProvenanceDownload').isEnabled());

  await page.click('#jointPaletteGroup button[data-value="ocean"]');
  await waitReady(page);
  const optimizedSnapshotId = await page.evaluate(() => window.rnaExplorer.snapshots.distribution.snapshot_id);
  assert.equal(optimizedSnapshotId, recovered.snapshotId, 'Healthy joint controls unnecessarily redrew the completed main plot');
  assert.deepEqual(errors, []);
  const result = {passed: true, recovered, healthyJointPreservedSnapshot: true, pageErrors: errors};
  if (process.env.RNA_BROWSER_OUTPUT) {
    await mkdir(process.env.RNA_BROWSER_OUTPUT, {recursive: true});
    await writeFile(`${process.env.RNA_BROWSER_OUTPUT}/result.json`, JSON.stringify(result, null, 2) + '\n');
  }
  console.log(JSON.stringify(result, null, 2));
} finally {await browser.close();}
