/** Read-only DNA smoke and byte-preservation audit. No baseline files are rewritten. */
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { readFile, mkdir, writeFile } from 'node:fs/promises';
import path from 'node:path';
import { pathToFileURL } from 'node:url';

const workspace = path.resolve(process.env.RNA_WORKSPACE || process.cwd());
const output = path.resolve(process.env.RNA_BROWSER_OUTPUT || path.join(workspace, 'data/pure_rna/browser_validation'));
const baseline = JSON.parse(await readFile(path.join(workspace, 'data/pure_rna/dna_baseline.json'), 'utf8'));
await mkdir(output, { recursive: true });
const checks = [];
for (const [file, expected] of Object.entries(baseline.files)) {
  const actual = createHash('sha256').update(await readFile(path.join(workspace, 'nucleic.pages', file))).digest('hex');
  checks.push({ file, expected, actual, unchanged: actual === expected });
}
const report = { checkedAt: new Date().toISOString(), baselineCreatedAt: baseline.created_at, checks };
await writeFile(path.join(output, 'dna-preservation.json'), JSON.stringify(report, null, 2) + '\n');
assert(checks.every(check => check.unchanged), 'Protected DNA bytes changed; see dna-preservation.json');

if (process.argv.includes('--hashes-only')) {
  console.log(`Protected DNA files unchanged: ${checks.length}`);
  process.exit(0);
}
const playwrightPath = process.env.PLAYWRIGHT_MODULE || path.join(workspace, 'data/pure_rna/browser-tools/node_modules/playwright/index.mjs');
const { chromium } = await import(pathToFileURL(playwrightPath).href);
const browser = await chromium.launch({ headless: true });
try {
  const page = await browser.newPage({ viewport: { width: 1440, height: 1000 } });
  const errors = [], failedRequests = [];
  page.on('pageerror', error => errors.push(error.message));
  page.on('requestfailed', request => failedRequests.push({ url: request.url(), error: request.failure()?.errorText }));
  const url = process.env.DNA_URL || 'http://127.0.0.1:8767/nucleic.pages/';
  await page.goto(url, { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => document.querySelector('#parameterSelect')?.options.length > 0, null, { timeout: 60000 });
  await page.waitForFunction(() => document.querySelector('#filteredObservationCount')?.textContent.trim() !== '-', null, { timeout: 60000 });
  await page.screenshot({ path: path.join(output, 'dna-baseline.png'), fullPage: true });
  const initial = await page.evaluate(() => ({ title: document.title, family: document.querySelector('#familySelect').value, parameter: document.querySelector('#parameterSelect').value, count: document.querySelector('#filteredObservationCount').textContent, plotTraces: document.querySelector('#plot').data?.length || 0 }));
  assert.equal(initial.title, 'Pure DNA Explorer');
  assert(initial.plotTraces > 0, 'DNA distribution did not draw');
  await page.click('#universeToggle');
  assert(await page.locator('#universeDrawer').isVisible());
  await page.click('#filteredToggle');
  assert(await page.locator('#filteredDrawer').isVisible());
  const parameters = await page.locator('#parameterSelect option').evaluateAll(options => options.map(option => option.value));
  if (parameters.length > 1) {
    await page.selectOption('#parameterSelect', parameters.find(value => value !== initial.parameter));
    await page.waitForFunction(parameter => document.querySelector('#parameterSelect').value === parameter, parameters.find(value => value !== initial.parameter));
  }
  await page.click('#baseGeometryLoad');
  await page.waitForFunction(() => document.querySelector('#baseGeometryTermSelect')?.options.length > 0, null, { timeout: 60000 });
  const survey = await page.evaluate(() => ({ terms: document.querySelector('#baseGeometryTermSelect').options.length, scalarRows: document.querySelector('#baseGeometryScalarRows').textContent, coordinateRows: document.querySelector('#baseGeometryCoordBody').rows.length }));
  await writeFile(path.join(output, 'dna-browser-baseline.json'), JSON.stringify({ checkedAt: new Date().toISOString(), url, initial, survey, errors, failedRequests, note: 'Existing DNA behavior is recorded without repair; hash checks protect the supplied working tree.' }, null, 2) + '\n');
  console.log(JSON.stringify({ protectedFiles: checks.length, initial, survey, errors, failedRequests }));
} finally { await browser.close(); }
