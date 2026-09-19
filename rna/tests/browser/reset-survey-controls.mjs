import assert from 'node:assert/strict';
import { pathToFileURL } from 'node:url';
import { waitReady } from './helpers.mjs';
const { chromium } = await import(pathToFileURL(process.env.PLAYWRIGHT_MODULE).href);
const browser = await chromium.launch({ headless: true });
try {
  const page = await browser.newPage();
  await page.goto('http://127.0.0.1:8767/nucleic.pages/rna/');
  await waitReady(page);
  await page.click('#baseGeometryLoad'); await waitReady(page);
  await page.selectOption('#surveyOpeningSelect', 'bins'); await waitReady(page);
  await page.click('#baseGeometryMinObsGroup [data-value="100"]'); await waitReady(page);
  await page.click('#resetFilters'); await waitReady(page);
  assert.equal(await page.locator('#surveyOpeningSelect').inputValue(), 'all');
  assert.equal(await page.locator('#baseGeometryMinObsGroup [aria-pressed="true"]').getAttribute('data-value'), '20');
  await page.click('#baseGeometryMinObsGroup [data-value="100"]'); await waitReady(page);
  assert.equal(await page.evaluate(() => window.rnaExplorer.state.survey.minimum), 100);
  console.log('PASS Survey reset synchronizes controls and subsequent clicks');
} finally { await browser.close(); }
