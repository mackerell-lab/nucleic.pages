/** Entity-scoped NAKB annotations must be searchable in the PDB drawer. */
import assert from 'node:assert/strict';
import { mkdir, writeFile } from 'node:fs/promises';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { waitReady } from './helpers.mjs';

const workspace = path.resolve(process.env.RNA_WORKSPACE || process.cwd());
const output = path.join(workspace, 'data/pure_rna/browser_validation/entity-annotation-search');
await mkdir(output, { recursive: true });
const { chromium } = await import(pathToFileURL(process.env.PLAYWRIGHT_MODULE || path.join(workspace, 'data/pure_rna/browser-tools/node_modules/playwright/index.mjs')).href);
const browser = await chromium.launch({ headless: true });
const page = await browser.newPage();
const report = { startedAt: new Date().toISOString(), checks: [], errors: [] };
page.on('pageerror', error => report.errors.push(error.message));
page.on('console', message => { if (message.type() === 'error') report.errors.push(message.text()); });
page.on('requestfailed', request => report.errors.push(request.failure()?.errorText));
try {
  await page.goto(process.env.RNA_URL || 'http://127.0.0.1:8767/nucleic.pages/rna/', { waitUntil: 'domcontentloaded', timeout: 120000 });
  await waitReady(page);
  await page.click('#universeToggle');
  const probe = await page.evaluate(() => {
    const app = window.rnaExplorer;
    for (const entity of app.metadata.entities ?? []) {
      const pdbId = String(entity.pdb_id ?? entity.entry_id ?? '').toUpperCase();
      const entry = app.entries.find(item => String(item.pdb_id ?? item.entry_id ?? '').toUpperCase() === pdbId);
      if (!entry || !pdbId) continue;
      const entryText = JSON.stringify(entry).toLowerCase();
      const values = [entity.functions, entity.function_tags, entity.structures, entity.structural_tags, entity.subtypes, entity.rna_types, entity.annotation_tags]
        .flatMap(value => Array.isArray(value) ? value : value == null ? [] : [value])
        .filter(value => typeof value === 'string' && value.length > 2);
      const value = values.find(item => !entryText.includes(item.toLowerCase()));
      if (value) return { pdbId, value };
    }
    return null;
  });
  assert(probe, 'No entity-only annotation exists in the full release');
  await page.fill('#universeSearch', probe.value);
  const rows = await page.locator('#universeTableBody tr').evaluateAll(items => items.map(row => ({ text: row.textContent, pdb: row.querySelector('a')?.textContent })));
  assert(rows.some(row => row.pdb === probe.pdbId), `Search missed ${probe.pdbId} for entity annotation ${probe.value}`);
  assert(rows.some(row => row.text.toLowerCase().includes(probe.value.toLowerCase())), 'Matching annotation was not rendered in the table');
  assert.deepEqual(report.errors, []);
  report.checks.push({ name: 'Entity annotation search and rendering', ...probe });
  report.passed = true;
} catch (error) {
  report.passed = false;
  report.failure = { message: error.message, stack: error.stack };
  throw error;
} finally {
  report.finishedAt = new Date().toISOString();
  await writeFile(path.join(output, 'report.json'), JSON.stringify(report, null, 2) + '\n');
  await browser.close();
}
console.log(JSON.stringify(report));
