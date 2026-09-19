/** Stream every candidate coordinate row against an independent source reader. */
import assert from 'node:assert/strict';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { writeFile } from 'node:fs/promises';
import { configureCoordinateCandidate } from './coordinate-candidate-routing.mjs';

const workspace = process.env.RNA_WORKSPACE || '/home/zhaomt/cmap/test15';
const base = process.env.RNA_BROWSER_ORIGIN || 'http://127.0.0.1:8767';
const sourceRelative = 'nucleic.pages/assets/pure_rna/releases/full_bundled_family_20260919/manifest.json';
const candidateUrl = new URL(process.env.RNA_COORDINATE_CANDIDATE_URL
  || '/data/pure_rna/packed_coordinate_candidate_20260919/candidate.json', base).href;
assert(process.env.PLAYWRIGHT_MODULE, 'PLAYWRIGHT_MODULE is required');
const { chromium } = await import(pathToFileURL(process.env.PLAYWRIGHT_MODULE).href);
const browser = await chromium.launch({ headless: true });
try {
  const page = await browser.newPage(), errors = [];
  page.on('pageerror', error => errors.push(error.message));
  page.on('console', message => { if (message.text().startsWith('Verified browser coordinates')) console.log(message.text()); });
  const routing = await configureCoordinateCandidate(page, candidateUrl, {
    sourceManifestUrl: new URL(`/${sourceRelative}`, base).href,
    sourceManifestPath: path.join(workspace, sourceRelative),
  });
  assert.equal(routing.groupCount, 5);
  assert.equal(routing.verifiedResources, 340);
  await page.goto(base);
  const evidence = await page.evaluate(async ({ routing }) => {
    const { RnaDataRepository } = await import('/nucleic.pages/rna/core/repository.js');
    const requests = [], root = new URL('.', routing.candidateUrl).href;
    const isCoordinate = url => String(url).startsWith(root) && String(url).includes('/survey/coordinates/');
    const repository = new RnaDataRepository({ manifestUrl: routing.sourceManifestUrl, fetchImpl: (input, options) => {
      if (isCoordinate(input)) requests.push(String(input));
      return fetch(input, options);
    } });
    const source = new RnaDataRepository({ manifestUrl: routing.originalManifestUrl });
    const manifest = await repository.loadManifest(), original = await source.loadManifest();
    if (manifest.build_id !== routing.sourceBuildId || original.build_id !== routing.sourceBuildId) throw new Error('Coordinate browser source build mismatch');
    if (Object.values(original.survey.coordinates.groups).some(group => group.partitions.some(part => /^https?:/.test(part.path)))) throw new Error('Independent source coordinates were substituted');
    function equal(a, b, location) {
      if (Object.is(a, b)) return;
      if (!a || !b || typeof a !== 'object' || typeof b !== 'object' || Array.isArray(a) !== Array.isArray(b)) throw new Error(`Coordinate value differs: ${location}`);
      const ak = Object.keys(a).sort(), bk = Object.keys(b).sort();
      if (ak.length !== bk.length || ak.some((key, index) => key !== bk[index])) throw new Error(`Coordinate fields differ: ${location}`);
      for (const key of ak) equal(a[key], b[key], `${location}.${key}`);
    }
    function noRetention(repo) {
      if (repo.coordinateRequests.size) throw new Error('Completed coordinate requests retained');
      if ([...repo.promises.keys()].some(key => key.startsWith('survey:coordinates:'))) throw new Error('Streaming coordinate rows retained');
      if (repo.bundleCacheBytes !== 0) throw new Error('Coordinate stream unexpectedly populated bundle cache');
    }
    const groups = []; let totalRows = 0, totalPartitions = 0, numericCoordinates = 0;
    async function compareGroup(group, options = {}) {
      const actual = repository.iterateSurveyCoordinates(group, options), expected = source.iterateSurveyCoordinates(group, options);
      let rows = 0, partitions = 0;
      try {
        while (true) {
          const [a, b] = await Promise.all([actual.next(), expected.next()]);
          if (a.done !== b.done) throw new Error(`Coordinate partition count differs: ${group}`);
          if (a.done) break;
          if (a.value.rows.length !== b.value.rows.length || a.value.rows.length > 10000) throw new Error('Coordinate row count differs or exceeds limit');
          if (Object.hasOwn(a.value, 'columns') || Object.hasOwn(a.value, 'missing')) throw new Error('Coordinate stream retained packed transport');
          equal(a.value.source_row_count, b.value.source_row_count, 'source row count');
          equal(a.value.selected_row_count, b.value.selected_row_count, 'selected row count');
          for (let index = 0; index < a.value.rows.length; index++) {
            const row = a.value.rows[index];
            equal(row, b.value.rows[index], `${group}/${partitions}/${index}`);
            numericCoordinates += ['x', 'y', 'z'].filter(key => typeof row[key] === 'number').length;
          }
          rows += a.value.rows.length; partitions++;
          noRetention(repository); noRetention(source);
          if (!options.entryIds && (totalPartitions + partitions) % 25 === 0) console.log(`Verified browser coordinates ${totalPartitions + partitions}/340`);
        }
      } finally { await actual.return(); await expected.return(); }
      return { group, rows, partitions };
    }
    for (const group of Object.keys(manifest.survey.coordinates.groups)) {
      const result = await compareGroup(group);
      const declared = manifest.survey.coordinates.groups[group];
      if (result.rows !== declared.row_count || result.partitions !== declared.partitions.length) throw new Error('Coordinate group totals differ');
      groups.push(result); totalRows += result.rows; totalPartitions += result.partitions;
    }
    if (requests.length !== totalPartitions || new Set(requests).size !== totalPartitions) throw new Error('Coordinate full sweep repeated or omitted partition requests');
    const allCoordinateRequests = requests.length, allNumericCoordinates = numericCoordinates;
    const firstGroup = Object.keys(manifest.survey.coordinates.groups)[0];
    const beforeEmpty = requests.length;
    const empty = await compareGroup(firstGroup, { entryIds: [] });
    if (empty.rows || empty.partitions || requests.length !== beforeEmpty) throw new Error('Explicit empty selection loaded coordinates');
    const occurrences = new Map();
    for (const part of manifest.survey.coordinates.groups[firstGroup].partitions) {
      for (const id of part.entry_ids) occurrences.set(id, (occurrences.get(id) || 0) + 1);
    }
    const [entry, eligible] = [...occurrences].sort((a, b) => a[1] - b[1] || a[0].localeCompare(b[0]))[0];
    const beforeTiny = requests.length;
    const tiny = await compareGroup(firstGroup, { entryIds: [entry.toLowerCase()] });
    if (!tiny.rows || tiny.partitions !== eligible || requests.length - beforeTiny !== eligible) throw new Error('Tiny entry selection did not prune coordinate partitions');
    const expectedTinyUrls = manifest.survey.coordinates.groups[firstGroup].partitions.filter(part => part.entry_ids.includes(entry)).map(part => part.path).sort();
    equal(requests.slice(beforeTiny).sort(), expectedTinyUrls, 'tiny eligible request URLs');

    let failureInjected = false, failureAttempts = 0;
    const retry = new RnaDataRepository({ manifestUrl: routing.sourceManifestUrl, fetchImpl: (input, options) => {
      if (isCoordinate(input)) {
        failureAttempts++;
        if (!failureInjected) { failureInjected = true; return Promise.resolve(new Response('Intentional coordinate retry probe', { status: 503 })); }
      }
      return fetch(input, options);
    } });
    let failureMessage = null;
    try { await retry.iterateSurveyCoordinates(firstGroup).next(); } catch (error) { failureMessage = error.message; }
    if (!failureInjected || !failureMessage?.includes('503')) throw new Error('Coordinate failure injection did not reject');
    noRetention(retry);
    const recovered = retry.iterateSurveyCoordinates(firstGroup), baseline = source.iterateSurveyCoordinates(firstGroup);
    try {
      const [a, b] = await Promise.all([recovered.next(), baseline.next()]);
      equal(a.value.rows, b.value.rows, 'retried coordinate rows');
    } finally { await recovered.return(); await baseline.return(); }
    if (failureAttempts !== 2) throw new Error('Coordinate retry did not refetch exactly once');
    noRetention(retry);

    const preAborted = new AbortController(); preAborted.abort(new DOMException('Pre-aborted test', 'AbortError'));
    const beforeAborted = requests.length;
    let preAbortName = null;
    try { await repository.iterateSurveyCoordinates(firstGroup, { signal: preAborted.signal }).next(); } catch (error) { preAbortName = error.name; }
    if (preAbortName !== 'AbortError' || requests.length !== beforeAborted) throw new Error('Pre-aborted stream fetched coordinates');

    let started, sawSignal = false, cancelAttempts = 0;
    const entered = new Promise(resolve => { started = resolve; }), controller = new AbortController();
    const cancelled = new RnaDataRepository({ manifestUrl: routing.sourceManifestUrl, fetchImpl: (input, options) => {
      if (!isCoordinate(input)) return fetch(input, options);
      cancelAttempts++; sawSignal ||= options?.signal === controller.signal;
      const response = fetch(input, options);
      started(); return response;
    } });
    const pending = cancelled.iterateSurveyCoordinates(firstGroup, { signal: controller.signal }).next();
    const observed = pending.then(() => null, error => error.name);
    await entered; controller.abort(new DOMException('In-flight test', 'AbortError'));
    if (await observed !== 'AbortError' || !sawSignal || cancelAttempts !== 1) throw new Error('In-flight coordinate cancellation failed');
    noRetention(cancelled);
    const afterCancel = cancelled.iterateSurveyCoordinates(firstGroup, { entryIds: [entry] });
    let afterCancelRows = 0;
    for await (const chunk of afterCancel) afterCancelRows += chunk.rows.length;
    if (afterCancelRows !== tiny.rows || cancelAttempts !== 1 + eligible) throw new Error('Cancellation poisoned subsequent coordinate loads');
    noRetention(cancelled);
    return { groups, totalRows, totalPartitions, numericCoordinates: allNumericCoordinates,
      exactAllFieldEquality: true, objectIsCoordinateComparison: true, allCoordinateRequests,
      noRetainedCoordinateRequestsOrTables: true, emptySelectionAdditionalRequests: 0,
      tinySelection: { entry, eligiblePartitions: eligible, ...tiny, exactAllFieldEquality: true },
      retry: { failureMessage, failureAttempts, exactAllFieldEquality: true },
      cancellation: { preAbortedRequests: 0, inFlightSignalForwarded: sawSignal, cancelledAttempts: 1, recoveryRows: afterCancelRows },
      heapUsedBytes: performance.memory?.usedJSHeapSize ?? null };
  }, { routing });
  assert.equal(evidence.totalPartitions, 340);
  assert.equal(evidence.groups.length, 5);
  assert.equal(evidence.allCoordinateRequests, 340);
  assert(evidence.numericCoordinates > 0);
  assert.deepEqual(errors, []);
  const report = { completedAt: new Date().toISOString(), routing, ...evidence, errors,
    limitation: 'Coordinate-only streamed repository acceptance. No activation or full Explorer UI acceptance. Preflight requests excluded; heap reading is not a heap bound.' };
  if (process.env.RNA_COORDINATE_REPORT) await writeFile(process.env.RNA_COORDINATE_REPORT, JSON.stringify(report, null, 2));
  console.log(JSON.stringify(report, null, 2));
} finally { await browser.close(); }
