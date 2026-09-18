/** Check context ranking against raw real-release scalars and opening values. */
import assert from 'node:assert/strict';
import { mkdir, writeFile } from 'node:fs/promises';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { downloadCsv, waitReady } from './helpers.mjs';

const workspace = path.resolve(process.env.RNA_WORKSPACE || process.cwd());
const output = path.join(workspace, 'data/pure_rna/browser_validation/survey-context-ranking');
await mkdir(output, { recursive: true });
const { chromium } = await import(pathToFileURL(process.env.PLAYWRIGHT_MODULE || path.join(workspace, 'data/pure_rna/browser-tools/node_modules/playwright/index.mjs')));
const browser = await chromium.launch({ headless: true });
const page = await browser.newPage({ acceptDownloads: true });
const report = { startedAt: new Date().toISOString(), checks: [], errors: [] };
page.on('pageerror', error => report.errors.push(error.message));
page.on('console', message => { if (message.type() === 'error') report.errors.push(message.text()); });
page.on('requestfailed', request => report.errors.push(request.failure()?.errorText));
try {
  await page.goto(process.env.RNA_URL || 'http://127.0.0.1:8767/nucleic.pages/rna/'); await waitReady(page);
  await page.click('#baseGeometryLoad'); await waitReady(page);
  await page.selectOption('#surveyGroupSelect', 'major_groove_distances'); await waitReady(page);
  await page.click('#surveyRankingLoad'); await waitReady(page);
  const evidence = await page.evaluate(async () => {
    const app = window.rnaExplorer;
    if (app.manifest.partial) throw new Error('Full RNA release required');
    const eligible = new Set(app.filteredEntries.map(entry => entry.pdb_id ?? entry.accession));
    const pairs = new Map((await app.repository.loadFamily('base_pair')).rows.filter(row => eligible.has(row.pdb_id)
      && !row.near && !row.alternative && !/^n[ct]/i.test(row.family) && !/^[n]?[ct][WHS]{2}a$/i.test(row.family)).map(row => [row.id, row]));
    const expected = [];
    for (const term of app.surveyTerms().filter(term => term.group === 'major_groove_distances')) {
      const table = await app.repository.loadSurveyScalars(term.id), source = Array.isArray(table) ? table : table.rows;
      const contexts = new Map();
      for (const row of source) {
        if (!eligible.has(row.pdb_id) || !pairs.has(row.pair_id) || !Number.isFinite(row.value) || !['available', 'ok', 'computed', 'valid'].includes(row.status)) continue;
        const opening = pairs.get(row.pair_id).values.opening;
        if (!Number.isFinite(opening)) continue;
        const bin = opening >= -16 && opening < -8 ? 0 : opening >= -8 && opening < 2 ? 1 : opening >= 2 && opening <= 10 ? 2 : -1;
        if (bin < 0) continue;
        const context = row.context ?? row.sequence_context ?? row.base;
        if (!contexts.has(context)) contexts.set(context, [[], [], []]);
        contexts.get(context)[bin].push(row.value);
      }
      for (const [context, groups] of contexts) expected.push({ term: term.id, context, counts: groups.map(values => values.length),
        means: groups.map(values => values.length ? values.reduce((a, b) => a + b, 0) / values.length : null) });
    }
    return { build: app.manifest.build_id, expected, actual: app.surveyRanks.map(rank => ({ ...rank, term: rank.term.id })) };
  });
  assert(evidence.expected.length >= 2, 'Need multiple real term/context comparisons');
  assert.equal(evidence.actual.length, evidence.expected.length);
  for (const expected of evidence.expected) {
    const actual = evidence.actual.find(rank => rank.term === expected.term && rank.context === expected.context);
    assert(actual); assert.deepEqual(actual.counts, expected.counts);
    expected.means.forEach((mean, index) => mean === null ? assert.equal(actual.means[index], null) : assert(Math.abs(actual.means[index] - mean) < 1e-10));
  }
  report.checks.push({ name: 'Term/context counts and means match raw scalars', ...evidence });
  for (const minimum of [5, 100]) {
    await page.click(`#baseGeometryMinObsGroup button[data-value="${minimum}"]`); await waitReady(page);
    const ranks = await page.evaluate(() => window.rnaExplorer.surveyRanks);
    let insufficientSeen = false;
    for (const rank of ranks) {
      assert.equal(rank.sufficient, rank.counts.every(n => n >= minimum));
      if (!rank.sufficient) insufficientSeen = true; else assert(!insufficientSeen, 'Sufficient row follows insufficient row');
    }
    report.checks.push({ name: `Minimum ${minimum} coverage ordering`, rows: ranks.length });
  }
  const selected = await page.evaluate(() => {
    const rank = window.rnaExplorer.surveyRanks.find(rank => rank.counts.every(n => n > 0));
    return { term: rank.term.id, context: rank.context, n: rank.counts.reduce((a, b) => a + b, 0) };
  });
  await page.locator(`#baseGeometryRankingBody tr[data-term="${selected.term}"][data-context="${selected.context}"] button`).click(); await waitReady(page);
  const exported = await downloadCsv(page, '#surveyCsvDownload', path.join(output, 'ranked-context.csv'));
  assert.equal(exported.rows.length, selected.n);
  assert(exported.rows.every(row => row.context === selected.context && row.parameter === selected.term && row.pair_id));
  const snapshot = await page.evaluate(() => ({ spec: window.rnaExplorer.snapshots.survey.selection_spec, provenance: window.rnaExplorer.snapshots.survey.provenance }));
  assert.deepEqual(snapshot.spec.contexts, [selected.context]);
  assert.deepEqual(snapshot.provenance.survey_contexts, [selected.context]);
  report.checks.push({ name: 'Ranking click selects exact context and CSV population', ...selected });
  await page.click(`#baseGeometryContextGroup button[data-value="${selected.context}"]`); await waitReady(page);
  const allContexts = await page.evaluate(() => window.rnaExplorer.snapshots.survey.selection_spec.contexts);
  assert.deepEqual(allContexts, []);
  await page.click('#resetFilters'); await waitReady(page);
  const reset = await page.evaluate(() => ({ ranking: window.rnaExplorer.state.survey.ranking,
    rows: [...document.querySelectorAll('#baseGeometryRankingBody tr[data-term]')].length,
    notice: document.querySelector('#baseGeometryRankingBody')?.textContent.trim() }));
  assert.equal(reset.ranking, false);
  assert.equal(reset.rows, 0, 'Reset retained stale term-ranking rows');
  assert.match(reset.notice, /Compute term ranking/);
  await page.setViewportSize({ width: 390, height: 844 });
  await page.screenshot({ path: path.join(output, 'mobile-survey.png'), fullPage: true });
  assert.deepEqual(report.errors, []); report.passed = true;
} catch (error) { report.passed = false; report.failure = error.stack; throw error; }
finally { report.finishedAt = new Date().toISOString(); await writeFile(path.join(output, 'report.json'), JSON.stringify(report, null, 2)); await browser.close(); }
console.log(JSON.stringify({ passed: report.passed, checks: report.checks.map(check => check.name), errors: report.errors }));
