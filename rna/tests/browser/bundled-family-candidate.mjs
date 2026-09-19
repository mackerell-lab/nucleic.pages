/** Browser repository acceptance for an unactivated, family-only candidate. */
import assert from 'node:assert/strict';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { writeFile } from 'node:fs/promises';
import { configureFamilyCandidate } from './family-candidate-routing.mjs';

const workspace = process.env.RNA_WORKSPACE || '/home/zhaomt/cmap/test15';
const base = process.env.RNA_BROWSER_ORIGIN || 'http://127.0.0.1:8767';
const sourceRelative = 'nucleic.pages/assets/pure_rna/releases/full_bundled_survey_20260919/manifest.json';
const candidateUrl = new URL(process.env.RNA_FAMILY_CANDIDATE_URL
  || '/data/pure_rna/bundled_family_candidate_20260919/candidate.json', base).href;
assert(process.env.PLAYWRIGHT_MODULE, 'PLAYWRIGHT_MODULE is required');
const { chromium } = await import(pathToFileURL(process.env.PLAYWRIGHT_MODULE).href);
const browser = await chromium.launch({ headless: true });
try {
  const page = await browser.newPage();
  const errors = [];
  page.on('pageerror', error => errors.push(error.message));
  page.on('console', message => { if (message.text().startsWith('Verified browser family')) console.log(message.text()); });
  const routing = await configureFamilyCandidate(page, candidateUrl, {
    sourceManifestUrl: new URL(`/${sourceRelative}`, base).href,
    sourceManifestPath: path.join(workspace, sourceRelative),
  });
  assert.equal(routing.familyCount, 14);
  await page.goto(base);
  const evidence = await page.evaluate(async ({ routing }) => {
    const { RnaDataRepository } = await import('/nucleic.pages/rna/core/repository.js');
    const familyRequests = [], bundleRequests = [];
    const candidateRoot = new URL('.', routing.candidateUrl).href;
    const trackedFetch = async (input, options) => {
      const url = String(input);
      if (url.startsWith(candidateRoot) && /\/families\/[^/]+\.json\.gz$/.test(url)) familyRequests.push(url);
      if (url.includes('/families/bundles/')) bundleRequests.push(url);
      return fetch(input, options);
    };
    // One completed family per repository bounds the independent comparison.
    const repository = new RnaDataRepository({ manifestUrl: routing.sourceManifestUrl, fetchImpl: trackedFetch, maxCachedFamilies: 1 });
    const source = new RnaDataRepository({ manifestUrl: routing.originalManifestUrl, maxCachedFamilies: 1 });
    const manifest = await repository.loadManifest(), original = await source.loadManifest();
    if (original.family_bundles) throw new Error('Source comparison unexpectedly uses bundled families');
    if (original.families.some(item => /^https?:/.test(item.path))) throw new Error('Source family paths were substituted');
    if (manifest.build_id !== routing.sourceBuildId || original.build_id !== routing.sourceBuildId) throw new Error('Browser source build mismatch');
    function equal(a, b, location) {
      if (Object.is(a, b)) return;
      if (!a || !b || typeof a !== 'object' || typeof b !== 'object' || Array.isArray(a) !== Array.isArray(b)) throw new Error(`Family value differs: ${location}`);
      const ak = Object.keys(a).sort(), bk = Object.keys(b).sort();
      if (ak.length !== bk.length || ak.some((key, index) => key !== bk[index])) throw new Error(`Family fields differ: ${location}`);
      for (const key of ak) equal(a[key], b[key], `${location}.${key}`);
    }
    const families = []; let peakSerializedCacheBytes = 0, totalRows = 0;
    for (const descriptor of manifest.families) {
      const before = { families: familyRequests.length, bundles: bundleRequests.length };
      const start = performance.now(), candidate = await repository.loadFamily(descriptor.id);
      const loadMs = performance.now() - start, expected = await source.loadFamily(descriptor.id);
      if (candidate.rows.length !== expected.rows.length) throw new Error(`Family row count differs: ${descriptor.id}`);
      for (let index = 0; index < candidate.rows.length; index++) equal(candidate.rows[index], expected.rows[index], `${descriptor.id}[${index}]`);
      if (Object.hasOwn(candidate, 'columns') || Object.hasOwn(candidate, 'missing')) throw new Error('Decoded family retained transport columns');
      if (repository.bundleCacheBytes > repository.maxBundleCacheBytes || repository.resolvedFamilies.size > 1) throw new Error('Candidate repository exceeded cache policy');
      if (repository.bundleRequests.size) throw new Error('Completed bundle request retained');
      peakSerializedCacheBytes = Math.max(peakSerializedCacheBytes, repository.bundleCacheBytes);
      totalRows += candidate.rows.length;
      const resourceRequests = familyRequests.length - before.families + bundleRequests.length - before.bundles;
      if (resourceRequests > 3) throw new Error(`Family exceeds three cold resource requests: ${descriptor.id}`);
      families.push({ id: descriptor.id, rows: candidate.rows.length, loadMs, resourceRequests });
      console.log(`Verified browser family ${families.length}/${manifest.families.length}: ${descriptor.id}`);
    }
    const last = manifest.families.at(-1).id;
    const beforeWarm = { families: familyRequests.length, bundles: bundleRequests.length };
    const previous = await repository.loadFamily(last), warmStart = performance.now();
    const warm = await repository.loadFamily(last), warmMs = performance.now() - warmStart;
    if (warm !== previous || familyRequests.length !== beforeWarm.families || bundleRequests.length !== beforeWarm.bundles) throw new Error('Warm family load repeated requests or decoded rows');

    // Inject one real repository network failure before any bundle can be cached.
    let failureInjected = false, failedBundleUrl = null, failedBundleAttempts = 0;
    const retry = new RnaDataRepository({ manifestUrl: routing.sourceManifestUrl, maxCachedFamilies: 1, fetchImpl: async (input, options) => {
      const url = String(input);
      if (url === failedBundleUrl) failedBundleAttempts++;
      if (!failureInjected && url.includes('/families/bundles/')) {
        failureInjected = true; failedBundleUrl = url; failedBundleAttempts = 1;
        return new Response('Intentional family bundle retry probe', { status: 503 });
      }
      return fetch(input, options);
    } });
    const first = manifest.families[0].id;
    let failedMessage = null;
    try { await retry.loadFamily(first); } catch (error) { failedMessage = error.message; }
    if (!failureInjected || !failedMessage?.includes('503')) throw new Error('Bundle failure injection did not reject family load');
    if (retry.promises.has(`family:${first}`) || retry.bundleRequests.size) throw new Error('Failed bundle/family promise poisoned retry cache');
    const recovered = await retry.loadFamily(first), expected = await source.loadFamily(first);
    for (let index = 0; index < expected.rows.length; index++) equal(recovered.rows[index], expected.rows[index], `retry:${first}[${index}]`);
    if (recovered.rows.length !== expected.rows.length || failedBundleAttempts !== 2 || retry.bundleRequests.size) throw new Error('Family bundle retry failed');
    return { families, totalRows, exactAllFieldEquality: true, independentSourceManifest: routing.originalManifestUrl,
      coldFamilyRequests: beforeWarm.families, coldBundleRequests: beforeWarm.bundles,
      uniqueBundleRequests: new Set(bundleRequests).size, peakSerializedCacheBytes,
      cacheBudgetBytes: repository.maxBundleCacheBytes, completedFamilyCacheLimit: repository.maxCachedFamilies,
      warmMs, warmAdditionalRequests: 0, retry: { failureInjected, failedMessage, failedBundleAttempts, exactAllFieldEquality: true },
      heapUsedBytes: performance.memory?.usedJSHeapSize ?? null };
  }, { routing });
  assert.equal(evidence.families.length, 14);
  assert.equal(evidence.coldFamilyRequests, 14);
  assert(evidence.coldBundleRequests > 0);
  assert.deepEqual(errors, []);
  const report = { completedAt: new Date().toISOString(), routing, ...evidence, errors,
    limitation: 'Family-only candidate repository acceptance. No activation, full Explorer UI acceptance, or JavaScript heap bound is claimed. Preflight traffic is excluded from browser request counts.' };
  if (process.env.RNA_FAMILY_REPORT) await writeFile(process.env.RNA_FAMILY_REPORT, JSON.stringify(report, null, 2));
  console.log(JSON.stringify(report, null, 2));
} finally { await browser.close(); }
