/** Actual Survey means, units, row ownership, and stale-group recovery. */
import assert from 'node:assert/strict';
import { mkdir, readFile, writeFile } from 'node:fs/promises';
import path from 'node:path';
import vm from 'node:vm';
import { pathToFileURL } from 'node:url';
import { downloadCsv, waitReady } from './helpers.mjs';

const workspace = path.resolve(process.env.RNA_WORKSPACE || process.cwd());
const output = path.resolve(process.env.RNA_BROWSER_OUTPUT || path.join(workspace, 'data/pure_rna/browser_validation/survey-ranking-display'));
await mkdir(output, { recursive: true });
const dna = await readFile(new URL('../../../js/pure-dna.js', import.meta.url), 'utf8');
const dnaFunction = name => { const start = dna.indexOf(`function ${name}(`); return dna.slice(start, dna.indexOf('\nfunction ', start + 1)); };
const oracle = vm.createContext({ state: {} });
vm.runInContext(['wrapCircular', 'circularDisplayValue', 'displayBaseGeometryValue', 'unitLabel'].map(dnaFunction).join('\n'), oracle);
const { chromium } = await import(pathToFileURL(process.env.PLAYWRIGHT_MODULE || path.join(workspace, 'data/pure_rna/browser-tools/node_modules/playwright-core/index.mjs')).href);
const browser = await chromium.launch({ headless: true }); const page = await browser.newPage();
const report = { startedAt: new Date().toISOString(), checks: [], errors: [], expectedErrors: [] };
page.on('pageerror', error => report.errors.push(error.message));
page.on('console', message => { if (message.type() === 'error') (message.text().includes('INJECTED ranking group load failure') ? report.expectedErrors : report.errors).push(message.text()); });
async function select(selector, value) { await page.selectOption(selector, value); await waitReady(page); }
async function check(name, mode, cacheUnchanged = false) {
  const evidence = await page.evaluate(() => {
    const app = window.rnaExplorer;
    return { group: app.state.survey.group, selected: app.state.survey.termId, contexts: app.state.survey.contexts,
      cacheUnchanged: !window.rankingDisplayCache || window.rankingDisplayCache.every(([key, cache, entries]) => app.rankingCache.get(key) === cache && entries.every(([term, ranks]) => cache.get(term) === ranks)),
      rawUnchanged: !window.rankingDisplayRaw || JSON.stringify(app.surveyRanks) === window.rankingDisplayRaw,
      ranks: app.surveyRanks.map(rank => {
        const row = [...document.querySelectorAll('#baseGeometryRankingBody tr[data-term]')].find(row => row.dataset.term === rank.term.id && row.dataset.context === rank.context);
        return { rank, cells: [...row.cells].map(cell => cell.textContent), unit: row.cells[0].querySelector('.meta')?.textContent.trim(), active: row.classList.contains('active-row'), current: row.getAttribute('aria-current'), pressed: row.querySelector('button').getAttribute('aria-pressed') };
      }) };
  });
  assert(evidence.ranks.length > 0);
  oracle.state.circularMode = mode;
  let negativeMeans = 0, currentRows = 0;
  for (const item of evidence.ranks) {
    assert.equal(item.cells.length, 9);
    const term = item.rank.term;
    assert.equal(item.unit, oracle.unitLabel(term.unit));
    for (let index = 0; index < 3; index++) {
      const expected = oracle.displayBaseGeometryValue(item.rank.means[index], { isCircular: Boolean(term.period), period: term.period }, 3);
      assert.equal(item.cells[3 + index], expected); if (expected.startsWith('-') && expected !== '-') negativeMeans++;
    }
    assert.equal(item.cells[6], Number.isFinite(item.rank.difference) ? item.rank.difference.toFixed(4) : '-');
    const active = term.id === evidence.selected && evidence.contexts.length === 1 && evidence.contexts[0] === item.rank.context;
    assert.equal(item.active, active); assert.equal(item.current, String(active)); assert.equal(item.pressed, String(active)); if (active) currentRows++;
  }
  if (cacheUnchanged) { assert(evidence.cacheUnchanged); assert(evidence.rawUnchanged); }
  report.checks.push({ name, passed: true, rows: evidence.ranks.length, negativeMeans, currentRows, cacheUnchanged: evidence.cacheUnchanged, rawUnchanged: evidence.rawUnchanged });
  return evidence;
}
try {
  await page.goto(process.env.RNA_URL || 'http://127.0.0.1:8767/nucleic.pages/rna/'); await waitReady(page);
  await page.click('#baseGeometryLoad'); await waitReady(page); await select('#surveyGroupSelect', 'major_groove_distances');
  await page.click('#surveyRankingLoad'); await waitReady(page);
  report.buildId = await page.evaluate(() => window.rnaExplorer.manifest.build_id);
  await check('Linear distances show declared angstrom units', 'wrap_360');
  const chosen = await page.evaluate(() => {
    const row = document.querySelector('#baseGeometryRankingBody tr[data-term]'); window.staleRankingRow = row; window.staleRankingButton = row.querySelector('button');
    const app = window.rnaExplorer, original = app.repository.loadSurveyScalars.bind(app.repository);
    app.repository.loadSurveyScalars = async (...args) => { app.repository.loadSurveyScalars = original; throw Error('INJECTED ranking group load failure'); };
    return { term: row.dataset.term, context: row.dataset.context, group: app.state.survey.group };
  });
  await page.selectOption('#surveyGroupSelect', 'base_local_dihedrals'); await page.waitForFunction(() => document.querySelector('#appStatus').dataset.state === 'error');
  await page.evaluate(() => window.staleRankingButton.click()); await waitReady(page);
  const recovered = await page.evaluate(() => ({ group: window.rnaExplorer.state.survey.group, visibleGroup: document.querySelector('#surveyGroupSelect').value, term: window.rnaExplorer.snapshots.survey.provenance.survey_term, contexts: window.rnaExplorer.state.survey.contexts, count: window.rnaExplorer.snapshots.survey.result.series.reduce((n, series) => n + series.values.length, 0) }));
  assert.equal(recovered.group, chosen.group); assert.equal(recovered.visibleGroup, chosen.group); assert.equal(recovered.term, chosen.term); assert.deepEqual(recovered.contexts, [chosen.context]); assert(recovered.count > 0);
  const recoveredCsv = await downloadCsv(page, '#surveyCsvDownload', path.join(output, 'recovered.csv'));
  assert.equal(recoveredCsv.rows.length, recovered.count); assert(recoveredCsv.rows.every(row => row.parameter === chosen.term && row.context === chosen.context));
  await check('Stale ranking click restores exact term group and selected row', 'wrap_360');
  report.checks.push({ name: 'Recovered CSV selects intended nonempty term and context', passed: true, ...recovered });
  assert(await page.evaluate(() => { const app = window.rnaExplorer, revision = app.revision; window.staleRankingButton.click(); return revision === app.revision; }));
  report.checks.push({ name: 'Detached old ranking button does not mutate state', passed: true });
  await select('#surveyGroupSelect', 'base_local_dihedrals');
  await page.evaluate(() => {
    const app = window.rnaExplorer;
    window.rankingDisplayCache = [...app.rankingCache].map(([key, cache]) => [key, cache, [...cache]]);
    window.rankingDisplayRaw = JSON.stringify(app.surveyRanks);
  });
  await check('Circular canonical means match independent DNA formatter', 'wrap_360', true);
  const before = await downloadCsv(page, '#surveyCsvDownload', path.join(output, 'canonical.csv'));
  for (const mode of ['signed_180', 'auto', 'wrap_360']) {
    await page.click(`#circularModeGroup button[data-value="${mode}"]`); await waitReady(page);
    await check(`Circular ${mode} preserves cached raw statistics`, mode, true);
    if (mode === 'signed_180') assert(report.checks.at(-1).negativeMeans > 0);
  }
  const after = await downloadCsv(page, '#surveyCsvDownload', path.join(output, 'wrapped.csv'));
  const raw = rows => rows.map(({ snapshot_id, ...row }) => row);
  assert.deepEqual(raw(after.rows), raw(before.rows)); report.checks.push({ name: 'Display mode changes preserve raw Survey CSV', passed: true, rows: after.rows.length });
  assert.deepEqual(report.errors, []); assert.equal(report.expectedErrors.length, 1); report.passed = true;
  console.log(`PASS ${report.checks.length} Survey ranking display checks`);
} catch (error) { report.passed = false; report.failure = error.stack; throw error; }
finally { report.finishedAt = new Date().toISOString(); await writeFile(path.join(output, 'report.json'), JSON.stringify(report, null, 2) + '\n'); await browser.close(); }
