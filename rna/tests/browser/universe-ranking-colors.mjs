/** Full-release inventory, ranking lifecycle counts and stable semantic curve colors. */
import assert from 'node:assert/strict';
import { mkdir, writeFile } from 'node:fs/promises';
import path from 'node:path';
import { createHash } from 'node:crypto';
import { pathToFileURL } from 'node:url';
import { downloadCsv, waitReady } from './helpers.mjs';

const workspace = path.resolve(process.env.RNA_WORKSPACE || process.cwd());
const output = process.env.RNA_BROWSER_OUTPUT || path.join(workspace, 'data/pure_rna/browser_validation/universe-ranking-colors');
const base = process.env.RNA_URL || 'http://127.0.0.1:8767/nucleic.pages/rna/';
await mkdir(output, { recursive: true });
const { chromium } = await import(pathToFileURL(process.env.PLAYWRIGHT_MODULE || path.join(workspace, 'data/pure_rna/browser-tools/node_modules/playwright-core/index.mjs')).href);
const browser = await chromium.launch({ headless: true });
const page = await browser.newPage({ viewport: { width: 1440, height: 1000 }, acceptDownloads: true });
const report = { started: new Date().toISOString(), checks: [], errors: [] };
const sourceReads = []; report.servedSourceSha256 = {};
page.on('response', response => {
  if (response.url().startsWith(base) && response.url().endsWith('.js')) sourceReads.push(response.body().then(bytes => {
    report.servedSourceSha256[response.url().slice(base.length)] = createHash('sha256').update(bytes).digest('hex');
  }));
});
page.on('pageerror', e => report.errors.push(e.message));
page.on('console', m => { if (m.type() === 'error' && !m.text().includes('Injected ranking count failure')) report.errors.push(m.text()); });
await page.route(base + 'main.js', async route => {
  const response = await route.fetch(); let body = await response.text();
  body = body.replace('const root =', `window.inventoryFamilyCalls=[];
const originalFamily=RnaDataRepository.prototype.loadFamily;
RnaDataRepository.prototype.loadFamily=function(id,...args){ window.inventoryFamilyCalls.push(id); return originalFamily.call(this,id,...args); };
const root =`);
  await route.fulfill({ response, body });
});
const click = async (group, value) => { await page.click(`#${group} button[data-value="${value}"]`); await waitReady(page); };
const select = async (id, value) => { await page.selectOption(`#${id}`, value); await waitReady(page); };
const number = text => Number(text.replaceAll(',', ''));
const record = (name, value = {}) => report.checks.push({ name, ...value });
async function inventory() { return page.locator('#overviewCards').textContent(); }
async function rankingText() { return page.evaluate(() => Object.fromEntries(['baseGeometrySurveyTerms', 'baseGeometryRankRows', 'baseGeometrySufficientRankRows'].map(id => [id, document.getElementById(id)?.textContent]))); }
async function colors() { return page.evaluate(() => window.rnaExplorer.snapshots.distribution.result.series.map((series, i) => ({ key: series.key, color: document.querySelector('#plot').data[i].line.color }))); }
const omitSnapshot = rows => rows.map(({ snapshot_id, ...row }) => row);
try {
  await page.goto(base); await waitReady(page);
  report.release = await page.evaluate(() => ({ build: window.rnaExplorer.manifest.build_id, partial: window.rnaExplorer.manifest.partial })); assert.equal(report.release.partial, false);
  const e = await page.evaluate(() => {
    const a = window.rnaExplorer, metadata = a.metadata, entries = metadata.entries, ids = new Set(entries.map(row => row.pdb_id));
    const entities = metadata.entities.filter(row => row.type === 'polymer' && row.polymer_type === 'polyribonucleotide' && ids.has(row.pdb_id));
    const axes = ['functions', 'subtypes', 'structures'].map(axis => {
      const tagSets = new Map(), known = new Set(), knownPdb = new Set();
      for (const entity of entities) {
        const identity = JSON.stringify([entity.pdb_id, String(entity.entity_id)]);
        const tags = new Set(entity[axis] || []);
        if (tags.size) { known.add(identity); knownPdb.add(entity.pdb_id); }
        for (const tag of tags) { if (!tagSets.has(tag)) tagSets.set(tag, { entities: new Set(), pdbs: new Set() }); const set = tagSets.get(tag); set.entities.add(identity); set.pdbs.add(entity.pdb_id); }
      }
      const details = document.querySelector(`#universeAnnotationInventory [data-annotation-axis="${axis}"]`);
      return { axis, expected: [...tagSets].map(([tag, set]) => ({ tag, entities: set.entities.size, entries: set.pdbs.size })),
        unknownEntities: entities.length - known.size, noAnnotatedPdb: entries.length - knownPdb.size,
        text: details.querySelector('p').textContent, rows: [...details.querySelectorAll('tbody tr')].map(row => [...row.children].map(cell => cell.textContent)) };
    });
    return { entries: entries.length, entities: entities.length, uniqueEntries: ids.size,
      expectedProfiles: ['conservative', 'relaxed', 'mw100'].map(key => entries.filter(row => row.profiles[key] === true).length),
      actualProfiles: [...document.querySelectorAll('#universeProfileCounts .metric-value')].map(node => node.textContent),
      families: a.families.map(family => ({ label: family.label ?? family.name ?? family.id.replaceAll('_', ' '), count: family.row_count })),
      familyRows: [...document.querySelectorAll('#universeFamilyInventory tbody tr')].map(row => [...row.children].map(cell => cell.textContent)),
      axes, familyCalls: window.inventoryFamilyCalls, familyHelp: document.querySelector('#universeFamilyInventory .meta').textContent,
      profileHelp: document.querySelector('#universeProfileCounts .meta').textContent };
  });
  assert.equal(e.entries, e.uniqueEntries); assert.deepEqual(e.actualProfiles.map(number), e.expectedProfiles);
  assert.deepEqual(e.familyRows, e.families.map(f => [f.label, f.count.toLocaleString('en-US')]));
  for (const axis of e.axes) {
    assert.equal(axis.rows.length, axis.expected.length);
    for (const tag of axis.expected) { const row = axis.rows.find(row => row[0].endsWith(`[${tag.tag}]`)); assert(row); assert.equal(number(row[1]), tag.entities); assert.equal(number(row[2]), tag.entries); }
    assert(axis.text.includes(`${axis.unknownEntities.toLocaleString('en-US')} RNA entities`));
    assert(axis.text.includes(`${axis.noAnnotatedPdb.toLocaleString('en-US')} PDB entries`));
  }
  assert.deepEqual([...new Set(e.familyCalls)], ['backbone'], 'Inventory must not eagerly load other scientific families');
  assert.match(e.familyHelp, /not unique observations/); assert.match(e.profileHelp, /must not be added/);
  record('Inventory equals independent raw metadata and manifest counts without eager family loads', e);
  const fixedInventory = await inventory();
  const expectedBaseColors = { A: '#174a7e', C: '#8c3b2a', G: '#146c43', U: '#8659a1' };
  for (const item of await colors()) assert.equal(item.color, expectedBaseColors[item.key]);
  const allBases = await downloadCsv(page, '#filteredCsvDownload', path.join(output, 'bases-all.csv'));
  await click('contextGroup', 'U');
  assert.deepEqual(await colors(), [{ key: 'U', color: expectedBaseColors.U }]);
  const uracil = await downloadCsv(page, '#filteredCsvDownload', path.join(output, 'bases-U.csv'));
  assert.deepEqual(omitSnapshot(uracil.rows), omitSnapshot(allBases.rows.filter(row => row.group === 'U')));
  assert.equal(await inventory(), fixedInventory); record('Uracil subset keeps purple and exact source CSV population; universe remains fixed', { allRows: allBases.rows.length, uracilRows: uracil.rows.length });
  await click('contextGroup', 'U'); await click('groupingGroup', 'method'); await click('methodGroup', 'nmr');
  const methodColors = Object.fromEntries((await colors()).map(item => [item.key, item.color]));
  assert.equal(methodColors.xray, '#174a7e'); assert.equal(methodColors.nmr, '#8c3b2a');
  const mixedMethods = await downloadCsv(page, '#filteredCsvDownload', path.join(output, 'methods-mixed.csv'));
  await click('methodGroup', 'xray'); assert.deepEqual(await colors(), [{ key: 'nmr', color: '#8c3b2a' }]);
  const nmrOnly = await downloadCsv(page, '#filteredCsvDownload', path.join(output, 'methods-nmr.csv'));
  assert.deepEqual(omitSnapshot(nmrOnly.rows), omitSnapshot(mixedMethods.rows.filter(row => row.group === 'nmr')));
  assert.equal(await inventory(), fixedInventory); record('NMR subset retains red and exact CSV population', { mixed: mixedMethods.rows.length, nmr: nmrOnly.rows.length });
  await page.click('#resetFilters'); await waitReady(page);

  // Native disclosures must work from keyboard and on a narrow viewport.
  const summary = page.locator('#universeFamilyInventory summary'); await summary.focus(); await page.keyboard.press('Enter');
  assert(await page.locator('#universeFamilyInventory details').evaluate(node => node.open));
  await page.keyboard.press('Enter'); assert(!(await page.locator('#universeFamilyInventory details').evaluate(node => node.open)));
  await page.setViewportSize({ width: 390, height: 844 });
  // Plotly schedules responsive resizing asynchronously after viewport changes.
  await page.waitForFunction(() => Math.abs(document.querySelector('#plot')._fullLayout.width - document.querySelector('#plot').clientWidth) < 4);
  for (const selector of ['#universeFamilyInventory details', ...['functions', 'subtypes', 'structures'].map(axis => `#universeAnnotationInventory [data-annotation-axis="${axis}"]`)]) {
    await page.locator(selector + ' summary').click(); assert(await page.locator(selector).evaluate(node => node.open));
  }
  assert(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth + 2), 'Mobile inventory avoids page overflow');
  await page.locator('#overviewCards').screenshot({ path: path.join(output, 'mobile-inventory.png') });
  assert.equal(await inventory(), fixedInventory); record('Inventory disclosures work with keyboard and mobile viewport');
  await page.setViewportSize({ width: 1440, height: 1000 });

  await page.click('#baseGeometryLoad'); await waitReady(page); await select('surveyGroupSelect', 'major_groove_distances');
  const totalDefinitions = await page.evaluate(() => window.rnaExplorer.surveyTerms().length);
  const initialRank = await rankingText(); assert.equal(number(initialRank.baseGeometrySurveyTerms), totalDefinitions);
  assert.equal(initialRank.baseGeometryRankRows, 'Not computed'); assert.equal(initialRank.baseGeometrySufficientRankRows, 'Not computed');
  record('Survey separates definitions from uncomputed ranking counts', initialRank);
  await page.evaluate(() => {
    const a = window.rnaExplorer, original = a.repository.loadSurveyScalars.bind(a.repository); let armed = true;
    a.repository.loadSurveyScalars = async function (term, ...args) {
      if (armed && a.state.survey.ranking && term !== a.lastSurveyTerm) {
        armed = false; await new Promise((resolve, reject) => { window.rejectRankingCount = reject; });
      }
      return original(term, ...args);
    };
  });
  await page.click('#surveyRankingLoad'); await page.waitForFunction(() => typeof window.rejectRankingCount === 'function');
  const pending = await rankingText(); assert.equal(pending.baseGeometryRankRows, 'Computing…'); assert.equal(pending.baseGeometrySufficientRankRows, 'Computing…');
  await page.evaluate(() => window.rejectRankingCount(Error('Injected ranking count failure')));
  await page.waitForFunction(() => document.querySelector('#appStatus').dataset.state === 'error');
  const failed = await rankingText(); assert.equal(failed.baseGeometryRankRows, 'Unavailable'); assert.equal(failed.baseGeometrySufficientRankRows, 'Unavailable');
  record('Real pending ranking and injected load failure never display stale numeric counts', { pending, failed });
  await page.click('#resetFilters'); await waitReady(page);
  const resetCounts = await rankingText(); assert.equal(resetCounts.baseGeometryRankRows, 'Not computed'); assert.equal(resetCounts.baseGeometrySufficientRankRows, 'Not computed');
  await select('surveyGroupSelect', 'major_groove_distances'); await page.click('#surveyRankingLoad'); await waitReady(page);
  const independent = await page.evaluate(async () => {
    const a = window.rnaExplorer, eligible = new Set(a.filteredEntries.map(row => row.pdb_id ?? row.accession));
    const pairs = new Map((await a.repository.loadFamily('base_pair')).rows.filter(row => eligible.has(row.pdb_id) && !row.near && !row.alternative && !/^n[ct]/i.test(row.family) && !/^[n]?[ct][WHS]{2}a$/i.test(row.family)).map(row => [row.id, row]));
    const expected = [];
    for (const term of a.surveyTerms().filter(term => term.group === 'major_groove_distances')) {
      const table = await a.repository.loadSurveyScalars(term.id), contexts = new Map();
      for (const row of table.rows ?? table) {
        if (!eligible.has(row.pdb_id) || !pairs.has(row.pair_id) || !Number.isFinite(row.value) || !['available', 'ok', 'computed', 'valid'].includes(row.status)) continue;
        const opening = pairs.get(row.pair_id).values.opening;
        const bin = opening >= -16 && opening < -8 ? 0 : opening >= -8 && opening < 2 ? 1 : opening >= 2 && opening <= 10 ? 2 : -1;
        if (bin < 0) continue;
        const context = row.context ?? row.sequence_context ?? row.base;
        if (!contexts.has(context)) contexts.set(context, [0, 0, 0]); contexts.get(context)[bin]++;
      }
      for (const [context, counts] of contexts) expected.push({ term: term.id, context, counts });
    }
    return expected;
  });
  for (const minimum of [1, 20]) {
    await click('baseGeometryMinObsGroup', String(minimum));
    const text = await rankingText(), ranks = await page.evaluate(() => window.rnaExplorer.surveyRanks.map(rank => ({ term: rank.term.id, context: rank.context, counts: rank.counts, sufficient: rank.sufficient })));
    assert.equal(ranks.length, independent.length); assert.equal(number(text.baseGeometryRankRows), independent.length);
    const sufficient = independent.filter(rank => rank.counts.every(n => n >= minimum)).length;
    assert.equal(number(text.baseGeometrySufficientRankRows), sufficient);
    for (const expected of independent) { const actual = ranks.find(rank => rank.term === expected.term && rank.context === expected.context); assert(actual); assert.deepEqual(actual.counts, expected.counts); assert.equal(actual.sufficient, expected.counts.every(n => n >= minimum)); }
    record(`Actual minimum ${minimum} ranking population and sufficient counts match raw incidences`, { text, ranks });
  }
  await page.evaluate(() => window.rnaExplorer.setSelection({ search: 'NO_RNA_MATCH_INVENTORY_000' })); await waitReady(page);
  const empty = await rankingText(); assert.equal(number(empty.baseGeometryRankRows), 0); assert.equal(number(empty.baseGeometrySufficientRankRows), 0);
  assert.equal(await page.locator('#baseGeometryRankingBody tr[data-term]').count(), 0);
  assert.equal(await inventory(), fixedInventory); record('Empty ranking has zero real rows despite placeholder and fixed universe', empty);
  await page.click('#resetFilters'); await waitReady(page); assert.equal((await rankingText()).baseGeometryRankRows, 'Not computed');
  record('Final Reset clears obsolete counts and ranking state');
  assert.deepEqual(report.errors, []); report.passed = true; console.log(`PASS ${report.checks.length} universe/ranking/color browser checks`);
} catch (error) { report.passed = false; report.failure = error.stack; throw error; }
finally { await Promise.all(sourceReads); report.finished = new Date().toISOString(); await writeFile(path.join(output, 'report.json'), JSON.stringify(report, null, 2)); await browser.close(); }
