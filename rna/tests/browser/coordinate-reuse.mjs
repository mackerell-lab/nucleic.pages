import assert from 'node:assert/strict';
import { pathToFileURL } from 'node:url';
import { waitReady } from './helpers.mjs';
const { chromium } = await import(pathToFileURL(process.env.PLAYWRIGHT_MODULE).href);
const browser = await chromium.launch({ headless: true });
try {
  const page = await browser.newPage();
  let requests = 0;
  page.on('request', request => { if (request.url().includes('/survey/coordinates/')) requests++; });
  await page.goto('http://127.0.0.1:8767/nucleic.pages/rna/'); await waitReady(page);
  await page.click('#coordinatesLoad'); await waitReady(page);
  assert(requests > 0);
  const initial = requests;
  const summary = await page.locator('#baseGeometryCoordBody').textContent();
  await page.click('#traceStyleGroup [data-value="line"]'); await waitReady(page);
  assert.equal(requests, initial, 'Display-only change reloaded coordinate partitions');
  assert.equal(await page.locator('#baseGeometryCoordBody').textContent(), summary);
  const opening = await page.locator('#coordinateOpeningSelect option').evaluateAll(nodes => nodes.map(node => node.value).find(value => value !== 'all'));
  assert(opening);
  await page.selectOption('#coordinateOpeningSelect', opening); await waitReady(page);
  assert(requests > initial, 'Coordinate population change failed to reload');
  assert.equal(await page.locator('#coordinateOpeningSelect').inputValue(), opening);
  console.log('PASS Display changes retain coordinates; opening changes reload partitions');
} finally { await browser.close(); }
