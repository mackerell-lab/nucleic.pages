/** All-field Chromium acceptance for an authenticated, unactivated scalar candidate. */
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { mkdir, writeFile } from 'node:fs/promises';
import { configurePackedSurveyCandidate } from './packed-survey-candidate-routing.mjs';
import { downloadCsv, waitReady } from './helpers.mjs';

const workspace = process.env.RNA_WORKSPACE || '/home/zhaomt/cmap/test15';
const base = process.env.RNA_BROWSER_ORIGIN || 'http://127.0.0.1:8767';
const sourceRelative = 'nucleic.pages/assets/pure_rna/releases/full_packed_family_20260919/manifest.json';
const candidateUrl = new URL(process.env.RNA_SURVEY_CANDIDATE_URL
  || '/data/pure_rna/packed_survey_candidate_20260919/candidate.json', base).href;
const output = process.env.RNA_BROWSER_OUTPUT || path.join(workspace, 'data/pure_rna/browser_validation/packed-survey-candidate-20260919');
await mkdir(output, { recursive: true });
assert(process.env.PLAYWRIGHT_MODULE, 'PLAYWRIGHT_MODULE is required');
const { chromium } = await import(pathToFileURL(process.env.PLAYWRIGHT_MODULE).href);
const browser = await chromium.launch({ headless: true });
const report = { startedAt: new Date().toISOString(), errors: [] };
const moduleResponses = [];
try {
  const page = await browser.newPage({ acceptDownloads: true });
  page.on('pageerror', error => report.errors.push(error.message));
  page.on('response', response => {
    if (/\/rna\/core\/[^/]+\.js$/.test(response.url())) moduleResponses.push(response.body().then(bytes => ({
      url: response.url(), sha256: createHash('sha256').update(bytes).digest('hex'),
    })));
  });
  page.on('console', message => { if (message.text().startsWith('Verified browser Survey')) console.log(message.text()); });
  const routingOptions = {
    sourceManifestUrl: new URL(`/${sourceRelative}`, base).href,
    sourceManifestPath: path.join(workspace, sourceRelative),
  };
  const indexResponse = await page.request.get(candidateUrl);
  assert(indexResponse.ok());
  const candidateIndex = await indexResponse.json(); await indexResponse.dispose();
  const firstTerm = Object.keys(candidateIndex.survey.scalars.terms)[0];
  const mutations = {
    'source-hash': index => { index.source_manifest.sha256 = '0'.repeat(64); },
    'build-id': index => { index.build_id = 'wrong-build'; },
    'term-count': index => { index.survey.scalars.terms[firstTerm].row_count++; },
    'path-traversal': index => { index.survey.scalars.terms[firstTerm].path = '../escape.json.gz'; },
    'bundle-registry': index => { delete index.survey.bundles[Object.keys(index.survey.bundles)[0]]; },
    'resource-hash': index => { index.survey.scalars.terms[firstTerm].sha256 = '0'.repeat(64); },
  };
  report.preflightRejections = [];
  for (const [name, mutate] of Object.entries(mutations)) {
    const index = structuredClone(candidateIndex); mutate(index);
    let installedRoutes = 0, failure;
    const probe = { request: { get: href => href === candidateUrl ? Promise.resolve({
      ok: () => true, status: () => 200, url: () => href,
      body: async () => Buffer.from(JSON.stringify(index)), dispose: async () => {},
    }) : page.request.get(href) }, route: async () => { installedRoutes++; } };
    try { await configurePackedSurveyCandidate(probe, candidateUrl, routingOptions); } catch (error) { failure = error.message; }
    assert(failure, `Malformed candidate index accepted: ${name}`);
    assert.equal(installedRoutes, 0, `Malformed candidate installed routes: ${name}`);
    report.preflightRejections.push({ name, failure, installedRoutes });
  }
  const routing = await configurePackedSurveyCandidate(page, candidateUrl, routingOptions);
  report.routing = routing;
  assert.equal(routing.termCount, 98);
  await page.goto(base);
  report.repository = await page.evaluate(async ({ routing }) => {
    const { RnaDataRepository } = await import('/nucleic.pages/rna/core/repository.js');
    const scalarRequests = [], bundleRequests = [];
    const candidateRoot = new URL('.', routing.candidateUrl).href;
    const repository = new RnaDataRepository({ manifestUrl: routing.sourceManifestUrl, fetchImpl: (input, options) => {
      const url = String(input);
      if (url.startsWith(candidateRoot) && url.includes('/survey/scalars/')) scalarRequests.push(url);
      if (url.includes('/survey/bundles/')) bundleRequests.push(url);
      return fetch(input, options);
    } });
    const source = new RnaDataRepository({ manifestUrl: routing.originalManifestUrl, fetchImpl: (input, options) => {
      const url = new URL(input);
      if (url.pathname.includes('/survey/bundles/')) url.searchParams.set('rna-survey-source', '1');
      return fetch(url.href, options);
    } });
    const manifest = await repository.loadManifest(), original = await source.loadManifest();
    if (Object.values(original.survey.scalars.terms).some(item => item.encoding !== 'rna-survey-bundled-columns-1' || /^https?:/.test(item.path))) {
      throw new Error('Independent source Survey was substituted');
    }
    if (manifest.build_id !== routing.sourceBuildId || original.build_id !== routing.sourceBuildId) throw new Error('Browser source build mismatch');
    function equal(a, b, location) {
      if (Object.is(a, b)) return;
      if (!a || !b || typeof a !== 'object' || typeof b !== 'object' || Array.isArray(a) !== Array.isArray(b)) throw new Error(`Survey value differs: ${location}`);
      const ak = Object.keys(a).sort(), bk = Object.keys(b).sort();
      if (ak.length !== bk.length || ak.some((key, index) => key !== bk[index])) throw new Error(`Survey fields differ: ${location}`);
      for (const key of ak) equal(a[key], b[key], `${location}.${key}`);
    }
    const terms = []; let peakSerializedCacheBytes = 0, totalRows = 0;
    await repository.loadSurveyScalars();
    if (scalarRequests.length || bundleRequests.length) throw new Error('Survey inventory eagerly loaded resources');
    const fields = ['id', 'value', 'status'];
    for (const [id, descriptor] of Object.entries(manifest.survey.scalars.terms)) {
      const before = { scalars: scalarRequests.length, bundles: bundleRequests.length };
      const start = performance.now(), candidate = await repository.loadSurveyScalars(id);
      const loadMs = performance.now() - start, expected = await source.loadSurveyScalars(id);
      equal(candidate.rows, expected.rows, id);
      if (Object.hasOwn(candidate, 'columns') || Object.hasOwn(candidate, 'missing')) throw new Error('Decoded Survey retained transport columns');
      const allowedBundles = new Set(routing.termBundles[id].map(reference => new URL(`survey/bundles/${reference}.json.gz`, routing.sourceManifestUrl).href));
      if (scalarRequests.length - before.scalars !== 1 || scalarRequests.at(-1) !== descriptor.path) throw new Error('Survey load fetched unrelated scalar payloads');
      if (bundleRequests.slice(before.bundles).some(url => !allowedBundles.has(url))) throw new Error('Survey load fetched unrelated bundles');
      if (!terms.length && bundleRequests.length !== allowedBundles.size) throw new Error('Cold Survey omitted or repeated bundles');
      const requestCount = scalarRequests.length - before.scalars + bundleRequests.length - before.bundles;
      if (requestCount > 3) throw new Error('Survey term exceeds three cold requests');
      // All terms exercise downstream ranking projection against independent rows.
      const projected = await repository.loadSurveyScalars(id, { fields });
      equal(projected.rows, expected.rows.map(row => Object.fromEntries(fields.filter(key => Object.hasOwn(row, key)).map(key => [key, row[key]]))), `projected:${id}`);
      const warmBefore = scalarRequests.length + bundleRequests.length;
      if (await repository.loadSurveyScalars(id, { fields }) !== projected || scalarRequests.length + bundleRequests.length !== warmBefore) throw new Error('Warm projected Survey load repeated work');
      if (repository.bundleCacheBytes > repository.maxBundleCacheBytes || repository.bundleRequests.size) throw new Error('Survey bundle cache policy failed');
      peakSerializedCacheBytes = Math.max(peakSerializedCacheBytes, repository.bundleCacheBytes);
      totalRows += candidate.rows.length;
      terms.push({ id, rows: candidate.rows.length, loadMs, fullLoadRequests: requestCount,
        projectionAdditionalRequests: scalarRequests.length - before.scalars - 1 });
      repository.releaseSurvey('scalars', id); source.releaseSurvey('scalars', id);
      if ([...repository.promises.keys(), ...source.promises.keys()].some(key => key.startsWith(`survey:scalars:${id}`))) throw new Error('Released term retains full or projected tables');
      console.log(`Verified browser Survey ${terms.length}/98: ${id}`);
    }
    const first = Object.keys(manifest.survey.scalars.terms)[0];
    const expected = await source.loadSurveyScalars(first), retries = [];
    for (const kind of ['bundle', 'payload']) {
      let failedUrl = null, attempts = 0, injected = false;
      const retry = new RnaDataRepository({ manifestUrl: routing.sourceManifestUrl, fetchImpl: (input, options) => {
        const url = String(input);
        if (url === failedUrl) attempts++;
        const match = kind === 'bundle' ? url.includes('/survey/bundles/') : url === manifest.survey.scalars.terms[first].path;
        if (!injected && match) {
          injected = true; failedUrl = url; attempts = 1;
          return Promise.resolve(new Response('Intentional packed Survey retry', { status: 503 }));
        }
        return fetch(input, options);
      } });
      let failure;
      try { await retry.loadSurveyScalars(first); } catch (error) { failure = error.message; }
      if (!injected || !failure?.includes('503') || retry.promises.has(`survey:scalars:${first}`)) throw new Error('Survey failure poisoned retry cache');
      // A healthy sibling bundle can still be running after Promise.all rejects.
      // Await those exact requests before testing completed-request eviction.
      const pendingSiblings = retry.bundleRequests.size;
      await Promise.allSettled([...retry.bundleRequests.values()]);
      if (retry.bundleRequests.size) throw new Error('Completed sibling bundle requests retained');
      const recovered = await retry.loadSurveyScalars(first);
      equal(recovered.rows, expected.rows, `retry:${kind}`);
      if (attempts !== 2 || retry.bundleRequests.size) throw new Error('Survey retry attempts differ');
      retries.push({ kind, failure, attempts, pendingSiblings, exactAllFieldEquality: true });
    }
    // Projection cannot hide corrupt unselected values or mismatched transports.
    const rejections = [];
    for (const corruption of ['null-mask', 'value-bytes', 'encoding', 'build', 'derived-id', 'bundle-content']) {
      const bad = new RnaDataRepository({ manifestUrl: routing.sourceManifestUrl, fetchImpl: async (input, options) => {
        const url = String(input);
        if (url !== manifest.survey.scalars.terms[first].path && !(corruption === 'bundle-content' && url.includes('/survey/bundles/'))) return fetch(input, options);
        const reader = new RnaDataRepository({ manifestUrl: routing.originalManifestUrl });
        const payload = await reader.readJson(url);
        if (url.includes('/survey/bundles/')) payload.columns[Object.keys(payload.columns)[0]][0] = 'intentional-corruption';
        else if (corruption === 'null-mask') payload.columns.value.nulls = [payload.row_count];
        else if (corruption === 'value-bytes') payload.columns.value.values.data = '!';
        else if (corruption === 'encoding') payload.encoding = 'rna-survey-bundled-columns-1';
        else if (corruption === 'build') payload.build_id = 'wrong-build';
        else if (corruption === 'derived-id') payload.columns.id.column = 'id';
        return new Response(JSON.stringify(payload));
      } });
      let failure;
      try { await bad.loadSurveyScalars(first, { fields: ['id', 'status'] }); } catch (error) { failure = error.message; }
      if (!failure || [...bad.promises.keys()].some(key => key.startsWith(`survey:scalars:${first}`))) throw new Error(`Malformed Survey accepted or retained: ${corruption}`);
      rejections.push({ corruption, failure });
    }
    return { terms, totalRows, exactAllFieldEquality: true, exactAllTermProjection: fields,
      scalarRequests: scalarRequests.length, bundleRequests: bundleRequests.length,
      uniqueBundleRequests: new Set(bundleRequests).size, peakSerializedCacheBytes, cacheBudgetBytes: repository.maxBundleCacheBytes,
      retries, rejections, heapUsedBytes: performance.memory?.usedJSHeapSize ?? null };
  }, { routing });
  assert.equal(report.repository.terms.length, 98);
  assert.equal(report.repository.totalRows, 6014017);
  // Open the ordinary Explorer, substituting only its pinned release manifest.
  await page.route(input => input.href === new URL('/nucleic.pages/assets/pure_rna/manifest.json', base).href,
    route => route.fulfill({ status: 200, contentType: 'application/json', body: JSON.stringify({
      build_id: routing.sourceBuildId, manifest: routing.sourceManifestUrl,
    }) }));
  await page.goto(new URL('/nucleic.pages/rna/', base).href); await waitReady(page);
  await page.click('#baseGeometryLoad'); await waitReady(page);
  const probes = await page.evaluate(() => {
    const terms = window.rnaExplorer.surveyTerms();
    const choose = test => terms.find(test);
    return [choose(term => term.id.startsWith('u_') && term.group === 'base_internal_angles'),
      choose(term => term.group === 'base_local_dihedrals'), choose(term => term.group === 'major_groove_distances')]
      .map(term => { if (!term) throw new Error('RNA Survey UI probe unavailable'); return { id: term.id, group: term.group }; });
  });
  report.ui = [];
  for (const term of probes) {
    await page.selectOption('#surveyGroupSelect', term.group); await waitReady(page);
    await page.selectOption('#baseGeometryTermSelect', term.id); await waitReady(page);
    const evidence = await page.evaluate(async ({ routing, term }) => {
      const { RnaDataRepository } = await import('/nucleic.pages/rna/core/repository.js');
      const source = new RnaDataRepository({ manifestUrl: routing.originalManifestUrl, fetchImpl: (input, options) => {
        const url = new URL(input); if (url.pathname.includes('/survey/bundles/')) url.searchParams.set('rna-survey-source', '1');
        return fetch(url.href, options);
      } });
      const raw = await source.loadSurveyScalars(term.id), rows = new Map(raw.rows.map(row => [row.id, row]));
      const app = window.rnaExplorer, snapshot = app.snapshots.survey;
      if (app.manifest.survey.scalars.terms[term.id].encoding !== 'rna-survey-float64-1') throw new Error('UI failed to pin packed Survey');
      if (snapshot.result.parameter.id !== term.id) throw new Error('UI retained wrong Survey term');
      const expected = snapshot.result.series.flatMap(series => series.values.map((value, index) => {
        const row = rows.get(series.rowIds[index]);
        if (!row || !Object.is(value, row.value)) throw new Error('UI measurement differs from independent original');
        return { id: row.id, value, group: series.key, weight: series.weights[index] };
      }));
      if (!expected.length) throw new Error('UI Survey probe empty');
      return { expected, snapshot: snapshot.snapshot_id, build: snapshot.build_id };
    }, { routing, term });
    const csv = await downloadCsv(page, '#surveyCsvDownload', path.join(output, `${term.id}.csv`));
    assert.equal(csv.rows.length, evidence.expected.length);
    for (let index = 0; index < csv.rows.length; index++) {
      const actual = csv.rows[index], expected = evidence.expected[index];
      assert.equal(actual.id, expected.id); assert.equal(Number(actual.value), expected.value);
      assert.equal(actual.group, expected.group); assert.equal(Number(actual.weight), expected.weight);
      assert.equal(actual.parameter, term.id); assert.equal(actual.snapshot_id, evidence.snapshot); assert.equal(actual.build_id, evidence.build);
    }
    report.ui.push({ ...term, csvRows: csv.rows.length, exactIndependentMeasurements: true, snapshot: evidence.snapshot });
  }
  await page.screenshot({ path: path.join(output, 'packed-survey-ui.png'), fullPage: true });
  assert.deepEqual(report.errors, []);
  report.passed = true;
} catch (error) { report.passed = false; report.failure = error.stack; throw error; }
finally {
  report.moduleResponses = await Promise.all(moduleResponses);
  report.finishedAt = new Date().toISOString();
  report.limitation = 'Full scalar all-field and projected comparisons plus three selected-term UI/CSV probes. No activation, full Explorer gate, memory ceiling, or matched performance claim. Preflight traffic excluded from browser request counts.';
  await writeFile(path.join(output, 'report.json'), JSON.stringify(report, null, 2));
  await browser.close();
}
console.log(JSON.stringify({ passed: report.passed, terms: report.repository.terms.length, rows: report.repository.totalRows, ui: report.ui, errors: report.errors }));
