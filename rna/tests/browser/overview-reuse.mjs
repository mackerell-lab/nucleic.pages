import assert from 'node:assert/strict';
import { pathToFileURL } from 'node:url';
import { waitReady } from './helpers.mjs';
const { chromium } = await import(pathToFileURL(process.env.PLAYWRIGHT_MODULE).href);
const browser = await chromium.launch({ headless: true });
try {
  const page = await browser.newPage();
  await page.goto('http://127.0.0.1:8767/nucleic.pages/rna/'); await waitReady(page);
  await page.selectOption('#familySelect', 'backbone'); await waitReady(page);
  const ids = await page.locator('#familyOverview [data-parameter]').evaluateAll(nodes => nodes.map(node => node.dataset.parameter));
  assert(ids.length > 2);
  await page.evaluate(() => { window.savedOverview = [...document.querySelector('#familyOverview').children]; });
  await page.click(`#familyOverview [data-parameter="${ids[1]}"]`); await waitReady(page);
  assert(await page.evaluate(() => window.savedOverview.every((node, i) => document.querySelector('#familyOverview').children[i] === node)));
  assert.equal(await page.locator('#familyOverview [aria-pressed="true"]').getAttribute('data-parameter'), ids[1]);
  await page.selectOption('#parameterSelect', ids[2]); await waitReady(page);
  assert(await page.evaluate(() => window.savedOverview.every((node, i) => document.querySelector('#familyOverview').children[i] === node)));
  assert.equal(await page.locator('#familyOverview [aria-pressed="true"]').getAttribute('data-parameter'), ids[2]);
  await page.click('#contextGroup [data-value="U"]'); await waitReady(page);
  assert(await page.evaluate(() => window.savedOverview.every(node => !node.isConnected)));
  assert.equal(await page.locator('#familyOverview [data-parameter]').count(), ids.length);
  console.log('PASS Parameter clicks reuse overview; population changes rebuild it');
} finally { await browser.close(); }
