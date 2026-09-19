/** Actual browser decode of a scalar-only candidate; never activates it. */
import assert from 'node:assert/strict';
import { pathToFileURL } from 'node:url';
import { writeFile } from 'node:fs/promises';
const { chromium } = await import(pathToFileURL(process.env.PLAYWRIGHT_MODULE).href);
const browser = await chromium.launch({ headless: true });
try {
  const page = await browser.newPage();
  const errors = [];
  page.on('pageerror', error => errors.push(error.message));
  await page.goto('http://127.0.0.1:8767/');
  const evidence = await page.evaluate(async () => {
    const { RnaDataRepository } = await import('/nucleic.pages/rna/core/repository.js');
    const { decodeSurveyRows } = await import('/nucleic.pages/rna/core/survey-codec.js');
    const url = new URL('/data/pure_rna/shared_survey_candidate_20260918/candidate.json', location.href).href;
    let columnRequests = 0;
    const repository = new RnaDataRepository({ manifestUrl: url, fetchImpl: async (input, options) => {
      if (String(input).includes('/survey/columns/')) columnRequests++;
      const response = await fetch(input, options);
      if (String(input) !== url) return response;
      const candidate = await response.json();
      // Adapt only the standalone candidate index for this isolated repository.
      // The production release pointer and production explorer are not involved.
      return new Response(JSON.stringify({ ...candidate, molecule_type: 'RNA', schema_version: 'rna-explorer-1' }));
    } });
    const started = performance.now();
    const term = 'g_n1_c2_n3';
    const full = await repository.loadSurveyScalars(term);
    const original = await repository.readJson(new URL(`/nucleic.pages/assets/pure_rna/releases/full_columnar_interaction_20260918/survey/scalars/${term}.json.gz`, location.href).href);
    const expected = decodeSurveyRows(original);
    const equal = JSON.stringify(full.rows) === JSON.stringify(expected);
    const projected = await repository.loadSurveyScalars(term, { fields: ['id', 'value', 'status'] });
    const projectionEqual = JSON.stringify(projected.rows) === JSON.stringify(decodeSurveyRows(original, ['id', 'value', 'status']));
    return { term, rows: full.rows.length, equal, projectionEqual, columnRequests,
      rawColumnsRetained: Object.hasOwn(full, 'columns') || Object.hasOwn(projected, 'columns'),
      elapsedMs: performance.now() - started };
  });
  assert.equal(evidence.rows, 70800);
  assert.equal(evidence.equal, true);
  assert.equal(evidence.projectionEqual, true);
  assert.equal(evidence.rawColumnsRetained, false);
  assert.deepEqual(errors, []);
  const report = { completedAt: new Date().toISOString(), ...evidence, errors };
  if (process.env.RNA_SHARED_REPORT) await writeFile(process.env.RNA_SHARED_REPORT, JSON.stringify(report, null, 2));
  console.log(JSON.stringify(report, null, 2));
} finally { await browser.close(); }
