/** Compare every projected candidate term in Chromium, retaining bounded bundles. */
import assert from 'node:assert/strict';
import { pathToFileURL } from 'node:url';
import { writeFile } from 'node:fs/promises';
const { chromium } = await import(pathToFileURL(process.env.PLAYWRIGHT_MODULE).href);
const browser = await chromium.launch({ headless: true });
try {
  const page = await browser.newPage();
  const errors = [];
  page.on('pageerror', error => errors.push(error.message));
  page.on('console', message => { if (message.text().startsWith('Verified browser')) console.log(message.text()); });
  const latencyMs = Number(process.env.RNA_REQUEST_DELAY_MS ?? 0);
  assert(Number.isFinite(latencyMs) && latencyMs >= 0 && latencyMs <= 1000);
  if (latencyMs) await page.route('**/data/pure_rna/bundled_survey_candidate_20260918/**', async route => {
    await new Promise(resolve => setTimeout(resolve, latencyMs)); await route.continue();
  });
  await page.goto('http://127.0.0.1:8767/');
  const evidence = await page.evaluate(async () => {
    const { RnaDataRepository } = await import('/nucleic.pages/rna/core/repository.js');
    const { decodeSurveyRows } = await import('/nucleic.pages/rna/core/survey-codec.js');
    const url = new URL('/data/pure_rna/bundled_survey_candidate_20260918/candidate.json', location.href).href;
    const bundleRequests = [], scalarRequests = [];
    const repository = new RnaDataRepository({ manifestUrl: url, fetchImpl: async (input, options) => {
      if (String(input).includes('/survey/bundles/')) bundleRequests.push(String(input));
      if (String(input).startsWith(new URL('.', url).href) && String(input).includes('/survey/scalars/')) scalarRequests.push(String(input));
      const response = await fetch(input, options);
      if (String(input) !== url) return response;
      const candidate = await response.json();
      return new Response(JSON.stringify({ ...candidate, molecule_type: 'RNA', schema_version: 'rna-explorer-1' }));
    } });
    const manifest = await repository.loadManifest();
    const sourceRoot = new URL('/nucleic.pages/assets/pure_rna/releases/full_columnar_interaction_20260918/', location.href);
    const fields = ['id', 'value', 'status'], terms = Object.keys(manifest.survey.scalars.terms);
    let count = 0, rows = 0, peakSerializedCacheBytes = 0, candidateLoadMs = 0;
    for (const term of terms) {
      const start = performance.now();
      const table = await repository.loadSurveyScalars(term, { fields });
      candidateLoadMs += performance.now() - start;
      const original = await repository.readJson(new URL(`survey/scalars/${term}.json.gz`, sourceRoot).href);
      if (JSON.stringify(table.rows) !== JSON.stringify(decodeSurveyRows(original, fields))) throw new Error(`Projected rows differ: ${term}`);
      if (Object.hasOwn(table, 'columns')) throw new Error('Retained raw columns');
      if (repository.bundleCacheBytes > repository.maxBundleCacheBytes) throw new Error('Bundle cache exceeded policy');
      peakSerializedCacheBytes = Math.max(peakSerializedCacheBytes, repository.bundleCacheBytes);
      rows += table.rows.length;
      repository.releaseSurvey('scalars', term);
      if (++count % 10 === 0) console.log(`Verified browser ${count}/${terms.length}`);
    }
    if (repository.bundleRequests.size) throw new Error('Completed bundle requests were retained');
    if ([...repository.promises.keys()].some(key => key.startsWith('survey:scalars:'))) throw new Error('Completed scalar tables were retained');
    return { terms: terms.length, rows, fields, exactProjectedEquality: true, candidateLoadMs,
      scalarRequests: scalarRequests.length, bundleRequests: bundleRequests.length,
      uniqueBundles: new Set(bundleRequests).size, peakSerializedCacheBytes,
      cacheBudgetBytes: repository.maxBundleCacheBytes, heapUsedBytes: performance.memory?.usedJSHeapSize ?? null };
  });
  assert.equal(evidence.terms, 98);
  assert.equal(evidence.rows, 6014017);
  assert.equal(evidence.scalarRequests, 98);
  assert(evidence.bundleRequests < 98, 'Bundling did not reduce additional requests below one per term');
  assert.deepEqual(errors, []);
  const report = { completedAt: new Date().toISOString(), latencyMs, ...evidence, errors,
    limitation: 'Projected-row correctness and request/cache probe, not a full Explorer UI acceptance or JavaScript heap bound.' };
  if (process.env.RNA_BUNDLED_REPORT) await writeFile(process.env.RNA_BUNDLED_REPORT, JSON.stringify(report, null, 2));
  console.log(JSON.stringify(report, null, 2));
} finally { await browser.close(); }
