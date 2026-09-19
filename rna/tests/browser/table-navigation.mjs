/** Real table navigation against archive metadata and raw finite observations. */
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { mkdir, readFile, writeFile } from 'node:fs/promises';
import path from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';
import { gunzipSync } from 'node:zlib';
import { waitReady, numericText } from './helpers.mjs';

const repository = fileURLToPath(new URL('../../../', import.meta.url));
const workspace = path.resolve(process.env.RNA_WORKSPACE || path.join(repository, '..'));
const output = path.resolve(process.env.RNA_BROWSER_OUTPUT || path.join(workspace, 'data/pure_rna/browser_validation/table-navigation'));
await mkdir(output, { recursive: true });
const pointerPath = path.join(repository, 'assets/pure_rna/manifest.json');
const pointer = JSON.parse(await readFile(pointerPath, 'utf8'));
const manifestPath = path.resolve(path.dirname(pointerPath), pointer.manifest);
const manifest = JSON.parse(await readFile(manifestPath, 'utf8'));
const metadataBytes = await readFile(path.resolve(path.dirname(manifestPath), manifest.metadata.path));
assert.equal(createHash('sha256').update(metadataBytes).digest('hex'), manifest.metadata.sha256);
const metadata = JSON.parse(gunzipSync(metadataBytes));
const id = entry => entry.pdb_id.toUpperCase();
const allIds = metadata.entries.map(id);
assert.equal(new Set(allIds).size, allIds.length);
// Deliberately independent of selectRows(), renderTable(), and page state.
const eligible = maximum => metadata.entries.filter(entry => entry.methods.includes('X-RAY DIFFRACTION')
  && entry.profiles.relaxed === true && Number.isFinite(entry.resolution) && entry.resolution <= maximum).map(id);
const { chromium } = await import(pathToFileURL(process.env.PLAYWRIGHT_MODULE || path.join(workspace, 'data/pure_rna/browser-tools/node_modules/playwright-core/index.mjs')).href);
const browser = await chromium.launch({ headless: true });
const page = await browser.newPage({ viewport: { width: 1440, height: 1000 } });
const report = { startedAt: new Date().toISOString(), buildId: manifest.build_id, metadataSha256: manifest.metadata.sha256, checks: [], errors: [] };
page.on('pageerror', error => report.errors.push(error.message));
page.on('console', message => { if (message.type() === 'error') report.errors.push(message.text()); });
page.on('requestfailed', request => report.errors.push(`${request.url()}: ${request.failure()?.errorText}`));
page.on('response', response => { if (response.status() >= 400) report.errors.push(`${response.status()} ${response.url()}`); });
const record = (name, evidence = {}) => report.checks.push({ name, passed: true, ...evidence });
async function table(name, expected, pageIndex, contributing = null) {
  const evidence = await page.evaluate(name => {
    const q = suffix => document.querySelector(`#${name}${suffix}`);
    return { label: q('PageLabel').textContent, previousDisabled: q('Prev').disabled, nextDisabled: q('Next').disabled,
      emptyText: q('TableBody').textContent, rows: [...q('TableBody').rows].filter(row => row.cells[0].querySelector('a')).map(row => ({ id: row.cells[0].querySelector('a').textContent,
        href: row.cells[0].querySelector('a').href, contributes: row.cells[6].textContent })) };
  }, name);
  const pages = Math.max(1, Math.ceil(expected.length / 100));
  if (!expected.length) assert.equal(evidence.emptyText, 'No entries match this selection.');
  assert.deepEqual(evidence.rows.map(row => row.id), expected.slice(pageIndex * 100, (pageIndex + 1) * 100));
  assert.equal(evidence.label, `Page ${pageIndex + 1} / ${pages} · ${expected.length.toLocaleString('en-US')} entries`);
  assert.equal(evidence.previousDisabled, pageIndex === 0);
  assert.equal(evidence.nextDisabled, pageIndex === pages - 1);
  for (const row of evidence.rows) {
    assert.equal(row.href, `https://www.rcsb.org/structure/${row.id}`);
    if (contributing) assert.equal(row.contributes, contributing.has(row.id) ? 'Yes' : 'No', `${row.id} finite contribution`);
  }
  return evidence.rows.map(row => row.id);
}
async function traverse(name, expected, contributing = null) {
  while (!await page.locator(`#${name}Prev`).isDisabled()) await page.click(`#${name}Prev`);
  const seen = [], pageCount = Math.max(1, Math.ceil(expected.length / 100));
  for (let index = 0; index < pageCount; index++) {
    seen.push(...await table(name, expected, index, contributing));
    if (index < pageCount - 1) await page.click(`#${name}Next`);
  }
  assert.deepEqual(seen, expected); assert.equal(new Set(seen).size, seen.length);
  // Exercise reverse navigation from the last page as well.
  if (pageCount > 1) { await page.click(`#${name}Prev`); await table(name, expected, pageCount - 2, contributing); }
  record(`${name}: every forward page and reverse boundary`, { pages: pageCount, entries: seen.length });
}
try {
  await page.goto(process.env.RNA_URL || 'http://127.0.0.1:8767/nucleic.pages/rna/', { waitUntil: 'domcontentloaded', timeout: 120000 });
  await waitReady(page);
  assert.equal(await page.evaluate(() => window.rnaExplorer.manifest.build_id), manifest.build_id);
  assert.equal(await page.locator('#universeDrawer').isVisible(), false);
  assert.equal(await page.locator('#filteredDrawer').isVisible(), false);
  await page.click('#universeToggle'); await page.click('#filteredToggle');
  assert.equal(await page.getAttribute('#universeToggle', 'aria-expanded'), 'true');
  assert.equal(await page.getAttribute('#filteredToggle', 'aria-expanded'), 'true');
  await traverse('universe', allIds);
  // Only the source loader is shared. Expected contribution bypasses production
  // selection, distribution, snapshot, and contributing-set implementations.
  const rawRows = await page.evaluate(async () => {
    const family = await window.rnaExplorer.repository.loadFamily('backbone');
    return (family.rows ?? family).map(row => ({ pdb: row.pdb_id.toUpperCase(), base: row.comp_id, chi: row.values.chi }));
  });
  const contributors = (entries, base = null) => {
    const allowed = new Set(entries);
    return new Set(rawRows.filter(row => allowed.has(row.pdb) && (!base || row.base === base) && Number.isFinite(row.chi)).map(row => row.pdb));
  };
  const defaultEntries = eligible(3);
  await traverse('filtered', defaultEntries, contributors(defaultEntries));
  // An exact ID chosen from the end must reset a later Universe page and trim.
  const probeId = allIds.at(-1);
  await page.fill('#universeSearch', `  ${probeId.toLowerCase()}  `);
  await table('universe', [probeId], 0);
  record('Whitespace and lowercase exact-ID search resets the page', { probeId });
  await page.fill('#universeSearch', 'RNA_TABLE_NO_SUCH_ENTRY_20260919');
  await table('universe', [], 0);
  record('Empty search result disables both pager buttons');
  await page.fill('#universeSearch', '   ');
  await table('universe', allIds, 0);
  await page.click('#universeNext'); await table('universe', allIds, 1);
  await page.click('#universeToggle'); assert.equal(await page.locator('#universeDrawer').isVisible(), false);
  assert.equal(await page.getAttribute('#universeToggle', 'aria-expanded'), 'false');
  await page.click('#universeToggle'); await table('universe', allIds, 1);
  record('Clearing search restores all entries and drawer reopening retains page');
  while (!await page.locator('#filteredNext').isDisabled()) await page.click('#filteredNext');
  await page.click('#resolutionGroup button[data-value="1.5"]'); await waitReady(page);
  const highResolution = eligible(1.5);
  assert(highResolution.length > 0 && highResolution.length < 100, 'Fixture must shrink to one page');
  await table('filtered', highResolution, 0, contributors(highResolution));
  await table('universe', allIds, 1);
  record('Entry filter shrink clamps old last page without changing Universe', { entries: highResolution.length });
  await page.click('#resetFilters'); await waitReady(page);
  await table('filtered', defaultEntries, 0, contributors(defaultEntries));
  await page.click('#contextGroup button[data-value="U"]'); await waitReady(page);
  const uracilContributors = contributors(defaultEntries, 'U');
  assert(uracilContributors.size > 0 && uracilContributors.size < defaultEntries.length, 'Fixture must distinguish eligible and contributing entries');
  await traverse('filtered', defaultEntries, uracilContributors);
  assert.equal(numericText(await page.textContent('#filteredPdbCount')), defaultEntries.length);
  assert.equal(numericText(await page.textContent('#contributingPdbCount')), uracilContributors.size);
  record('Uracil selection keeps eligible entries and marks noncontributors No', { eligible: defaultEntries.length, contributing: uracilContributors.size });
  await page.click('#resetFilters'); await waitReady(page);
  await traverse('filtered', defaultEntries, contributors(defaultEntries));
  assert.equal(numericText(await page.textContent('#contributingPdbCount')), contributors(defaultEntries).size);
  record('Reset restores default eligible and finite contribution populations');
  // Real touch interactions in a separate mobile context; no programmatic state
  // mutation is used to navigate or apply selection anywhere in this script.
  const mobile = await browser.newContext({ viewport: { width: 390, height: 844 }, isMobile: true, hasTouch: true });
  const phone = await mobile.newPage();
  phone.on('pageerror', error => report.errors.push(`mobile: ${error.message}`));
  phone.on('console', message => { if (message.type() === 'error') report.errors.push(`mobile: ${message.text()}`); });
  phone.on('requestfailed', request => report.errors.push(`mobile: ${request.failure()?.errorText}`));
  phone.on('response', response => { if (response.status() >= 400) report.errors.push(`mobile: ${response.status()} ${response.url()}`); });
  await phone.goto(process.env.RNA_URL || 'http://127.0.0.1:8767/nucleic.pages/rna/', { waitUntil: 'domcontentloaded', timeout: 120000 }); await waitReady(phone);
  await phone.tap('#universeToggle'); await phone.tap('#universeNext');
  assert.deepEqual(await phone.locator('#universeTableBody tr').evaluateAll(rows => rows.map(row => row.cells[0].querySelector('a').textContent)), allIds.slice(100, 200));
  await phone.fill('#universeSearch', probeId.toLowerCase());
  assert.equal(await phone.locator('#universeTableBody tr').count(), 1);
  assert.equal(await phone.locator('#universePrev').isDisabled(), true); assert.equal(await phone.locator('#universeNext').isDisabled(), true);
  await phone.tap('#filteredToggle'); await phone.tap('#filteredNext');
  assert.deepEqual(await phone.locator('#filteredTableBody tr').evaluateAll(rows => rows.map(row => row.cells[0].querySelector('a').textContent)), defaultEntries.slice(100, 200));
  await phone.tap('#filteredPrev');
  assert.equal(await phone.locator('#filteredPrev').isDisabled(), true);
  await phone.screenshot({ path: path.join(output, 'mobile.png') });
  await mobile.close();
  record('Mobile touch drawers, both pagers and exact search operate');
  assert.deepEqual(report.errors, []); report.passed = true;
  console.log(`PASS ${report.checks.length} table navigation browser checks`);
} catch (error) { report.passed = false; report.failure = error.stack; throw error; }
finally { report.finishedAt = new Date().toISOString(); await writeFile(path.join(output, 'report.json'), JSON.stringify(report, null, 2) + '\n'); await browser.close(); }
