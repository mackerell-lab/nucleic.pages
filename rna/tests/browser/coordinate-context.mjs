/** Exercise coordinate context controls with real release data and UI changes. */
import assert from 'node:assert/strict';
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
  await page.click('#coordinatesLoad');
  await waitReady(page);
  const selected = await page.locator('#coordinateContextSelect option').evaluateAll(options => options.find(option => option.value !== 'all')?.value);
  assert(selected, 'The default coordinate frame needs an observed context for this regression');
  await page.selectOption('#coordinateContextSelect', selected);
  await waitReady(page);
  assert(await page.evaluate(() => window.rnaExplorer.coordinateSummary.length > 0));

  const emptyPopulation = await page.evaluate(async () => {
    const {selectRows} = await import('./core/selection.js');
    const app = window.rnaExplorer;
    for (const method of ['em', 'other', 'nmr']) for (const profile of ['relaxed', 'conservative']) {
      const selection = {...app.state.selection, methods: [method], components: profile, resolutionMax: 1.5};
      if (selectRows([], app.metadata, selection).entryIds.length === 0) return {method, profile};
    }
    return null;
  });
  assert(emptyPopulation, 'Real release needs a representable empty entry population for this regression');
  await page.click(`#methodGroup [data-value="${emptyPopulation.method}"]`);
  await waitReady(page);
  await page.click('#methodGroup [data-value="xray"]');
  await waitReady(page);
  await page.click(`#cleanlinessGroup [data-value="${emptyPopulation.profile}"]`);
  await waitReady(page);
  await page.click('#resolutionGroup [data-value="1.5"]');
  await waitReady(page);
  assert.equal(await page.evaluate(() => window.rnaExplorer.coordinateSummary.length), 0,
    'The empty entry population must have no coordinate observations');
  assert.equal(await page.locator('#coordinateContextSelect').inputValue(), selected,
    'Filtered-out context was visually changed to All while its filter remained active');
  assert.equal(await page.evaluate(() => window.rnaExplorer.state.survey.coordinateContext), selected);
  assert.match(await page.locator('#coordinateContextSelect option:checked').textContent(), /no observations/i);

  // An actual control change to All must now clear the hidden restriction.
  await page.selectOption('#coordinateContextSelect', 'all');
  await waitReady(page);
  assert.equal(await page.evaluate(() => window.rnaExplorer.state.survey.coordinateContext), 'all');
  assert.equal(await page.locator('#coordinateContextSelect').inputValue(), 'all');
  await page.click('#resolutionGroup [data-value="3"]');
  await waitReady(page);
  await page.click('#cleanlinessGroup [data-value="relaxed"]');
  await waitReady(page);
  await page.click('#methodGroup [data-value="xray"]');
  await waitReady(page);
  await page.click(`#methodGroup [data-value="${emptyPopulation.method}"]`);
  await waitReady(page);
  const restored = await page.evaluate(() => ({buildId: window.rnaExplorer.manifest.build_id,
    context: window.rnaExplorer.state.survey.coordinateContext, atomRows: window.rnaExplorer.coordinateSummary.length}));
  assert(restored.atomRows > 0, 'Restoring X-ray did not restore coordinate observations');
  assert.equal(restored.context, 'all');
  assert.equal(await page.locator('#coordinateContextSelect').inputValue(), 'all');
  assert.deepEqual(errors, []);
  console.log(JSON.stringify({passed: true, selectedContext: selected, emptyPopulation, restored, pageErrors: errors}, null, 2));
} finally {await browser.close();}
