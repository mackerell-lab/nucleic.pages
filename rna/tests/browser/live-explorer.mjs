/** Real-release end-to-end validation, with optional authenticated candidate routing. */
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { mkdir, writeFile } from 'node:fs/promises';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { downloadCsv, waitReady, numericText } from './helpers.mjs';
import { checkPairControls, checkPuckerSurvey, checkBroadResidueScope, checkPalettes } from './rna-specific-controls.mjs';
import { configureSurveyCandidate } from './candidate-routing.mjs';
import { configureReleaseCandidate } from './release-routing.mjs';
import { configurePackedSurveyCandidate } from './packed-survey-candidate-routing.mjs';

const workspace = path.resolve(process.env.RNA_WORKSPACE || process.cwd());
const output = path.resolve(process.env.RNA_BROWSER_OUTPUT || path.join(workspace, 'data/pure_rna/browser_validation'));
await mkdir(output, { recursive: true });
const { chromium } = await import(pathToFileURL(process.env.PLAYWRIGHT_MODULE || path.join(workspace, 'data/pure_rna/browser-tools/node_modules/playwright/index.mjs')).href);
const report = { startedAt: new Date().toISOString(), url: process.env.RNA_URL || 'http://127.0.0.1:8767/nucleic.pages/rna/', checks: [], pageErrors: [], consoleErrors: [], failedRequests: [], responses: [] };
const browser = await chromium.launch({ headless: true });
const page = await browser.newPage({ viewport: { width: 1440, height: 1000 }, acceptDownloads: true });
const packedModuleResponses = [];
await page.addInitScript(() => {
  performance.setResourceTimingBufferSize(10000);
  globalThis.rnaBrowserLongTasks = { count: 0, totalMs: 0, maximumMs: 0 };
  if (PerformanceObserver.supportedEntryTypes.includes('longtask')) {
    new PerformanceObserver(list => { for (const task of list.getEntries()) {
      const report = globalThis.rnaBrowserLongTasks; report.count++; report.totalMs += task.duration; report.maximumMs = Math.max(report.maximumMs, task.duration);
    } }).observe({ type: 'longtask', buffered: true });
  }
});
page.on('pageerror', error => { report.pageErrors.push(error.message); console.error(`PAGE ERROR ${error.message}`); });
page.on('console', message => { if (message.type() === 'error') { report.consoleErrors.push(message.text()); console.error(`CONSOLE ERROR ${message.text()}`); } });
page.on('requestfailed', request => report.failedRequests.push({ url: request.url(), error: request.failure()?.errorText }));
page.on('response', response => report.responses.push({ url: response.url(), status: response.status() }));
page.on('response', response => {
  if (process.env.RNA_PACKED_SURVEY_CANDIDATE_URL && /\/rna\/core\/[^/]+\.js$/.test(response.url())) {
    packedModuleResponses.push(response.body().then(bytes => ({ url: response.url(),
      sha256: createHash('sha256').update(bytes).digest('hex') })));
  }
});
const record = (name, evidence) => { report.checks.push({ name, passed: true, elapsedMs: Date.now() - Date.parse(report.startedAt), evidence }); console.log(`PASS ${name}`); };

async function snapshot(kind = 'distribution') {
  return page.evaluate(kind => {
    const source = window.rnaExplorer.snapshots[kind];
    if (!source) return null;
    const result = source.result;
    return { snapshot_id: source.snapshot_id, build_id: source.build_id, selection_spec: source.selection_spec,
      result: { ...result,
        series: result.series?.map(series => ({ ...series, rows: series.rows?.map(row => ({ id: row.id, comp_id: row.comp_id,
          model_id: row.model_id, insertion_code: row.insertion_code, altloc: row.altloc,
          context: row.context ?? row.context_id ?? row.sequence_context ?? row.pair_label ?? row.step_label ?? row.comp_id ?? '' })) })),
        points: result.points?.map(point => ({ x: point.x, y: point.y, left_id: point.left_id, right_id: point.right_id })) } };
  }, kind);
}
async function verifyDistribution(name, button = '#filteredCsvDownload', kind = 'distribution') {
  const current = await snapshot(kind);
  assert(current?.result?.series, `${kind} has no completed snapshot`);
  const exported = await downloadCsv(page, button, path.join(output, `${name}.csv`));
  const expected = current.result.series.flatMap(series => series.values.map((value, index) => ({
    id: series.rowIds[index], group: series.key, value, weight: series.weights[index], context: series.rows[index]?.context ?? '',
    model_id: series.rows[index]?.model_id, insertion_code: series.rows[index]?.insertion_code, altloc: series.rows[index]?.altloc,
  })));
  assert.equal(exported.rows.length, expected.length, 'CSV membership count differs from plotted snapshot');
  for (let i = 0; i < expected.length; i++) {
    const actual = exported.rows[i], wanted = expected[i];
    assert.equal(actual.snapshot_id, current.snapshot_id);
    assert.equal(actual.build_id, current.build_id);
    assert.equal(actual.id, wanted.id);
    assert.equal(actual.group, wanted.group);
    assert.equal(actual.context, wanted.context, 'CSV lost the displayed sequence context');
    for (const field of ['model_id', 'insertion_code', 'altloc']) assert.equal(actual[field], wanted[field] == null ? '' : String(wanted[field]), `CSV lost ${field}`);
    assert.equal(Number(actual.value), wanted.value, 'Raw measurement lost numeric precision');
    assert.equal(Number(actual.weight), wanted.weight);
    assert.equal(actual.parameter, current.result.parameter.id);
  }
  if (kind === 'distribution') {
    const uiCount = numericText(await page.locator('#filteredObservationCount').textContent());
    assert.equal(uiCount, new Set(expected.map(row => row.id)).size);
  }
  record(name, { buildId: current.build_id, parameter: current.result.parameter.id, csvRows: expected.length, uniqueObservations: new Set(expected.map(row => row.id)).size,
    rowsWithInsertionCodes: expected.filter(row => row.insertion_code).length, rowsWithAltlocs: expected.filter(row => row.altloc).length });
  return current;
}

try {
  if ([process.env.RNA_RELEASE_CANDIDATE_URL, process.env.RNA_SURVEY_CANDIDATE_URL,
    process.env.RNA_PACKED_SURVEY_CANDIDATE_URL].filter(Boolean).length > 1) throw new Error('Choose only one release, scalar, or packed scalar candidate mode');
  if (process.env.RNA_RELEASE_CANDIDATE_URL) report.releaseCandidate = await configureReleaseCandidate(page, process.env.RNA_RELEASE_CANDIDATE_URL);
  if (process.env.RNA_SURVEY_CANDIDATE_URL) report.surveyCandidate = await configureSurveyCandidate(page, process.env.RNA_SURVEY_CANDIDATE_URL);
  if (process.env.RNA_PACKED_SURVEY_CANDIDATE_URL) {
    const sourceRelative = 'nucleic.pages/assets/pure_rna/releases/full_packed_family_20260919/manifest.json';
    report.packedSurveyCandidate = await configurePackedSurveyCandidate(page, process.env.RNA_PACKED_SURVEY_CANDIDATE_URL, {
      sourceManifestUrl: new URL(`/${sourceRelative}`, report.url).href,
      sourceManifestPath: path.join(workspace, sourceRelative),
    });
  }
  const navigationStarted = Date.now();
  await page.goto(report.url, { waitUntil: 'domcontentloaded', timeout: 120000 });
  await waitReady(page);
  report.initialReadyMs = Date.now() - navigationStarted;
  assert.equal(await page.title(), 'Pure RNA Explorer');
  const release = await page.evaluate(() => window.rnaExplorer.repository.loadManifest());
  if (report.releaseCandidate) {
    assert.equal(release.build_id, report.releaseCandidate.buildId, 'Explorer did not select the staged release');
    assert.equal(await page.evaluate(() => window.rnaExplorer.repository.releaseUrl), report.releaseCandidate.candidateUrl);
  }
  if (report.packedSurveyCandidate) {
    const routing = report.packedSurveyCandidate, candidateRoot = new URL('.', routing.candidateUrl).href;
    assert.equal(release.build_id, routing.sourceBuildId, 'Explorer did not select the packed scalar source build');
    assert.equal(await page.evaluate(() => window.rnaExplorer.repository.releaseUrl), routing.sourceManifestUrl);
    assert.equal(Object.keys(release.survey.scalars.terms).length, routing.termCount);
    assert(Object.values(release.survey.scalars.terms).every(term => term.encoding === 'rna-survey-float64-1'
      && term.path.startsWith(candidateRoot)), 'Explorer did not select all packed scalar descriptors');
    assert(release.families.every(family => !/^https?:/.test(family.path)), 'Scalar-only candidate changed family paths');
  }
  report.release = { buildId: release.build_id, partial: release.partial === true, counts: release.counts, capabilities: release.capabilities };
  assert(!release.partial || process.env.RNA_ALLOW_PARTIAL === '1', 'Partial releases require RNA_ALLOW_PARTIAL=1 and establish integration evidence only');
  await page.waitForFunction(() => document.querySelector('#plot')?.data?.length > 0);
  const initial = await verifyDistribution('initial-distribution');
  assert(initial.result.coverage.plottedRows > 0, 'Default real-data selection is empty');
  assert.equal(initial.result.parameter.id, 'chi');
  const summaryLabels = await page.locator('#seriesSummary .metric-label').allTextContents();
  assert(summaryLabels.includes('Rows') && summaryLabels.includes('PDBs'), 'Series summary lacks DNA-compatible row and PDB counts');
  assert(!report.responses.some(response => /assets\/pure_dna\//.test(response.url)), 'RNA fetched DNA datasets');
  const initiallyFetched = [...report.responses];
  const contexts = await page.locator('#contextGroup button').allTextContents();
  assert(contexts.some(value => value.trim() === 'U'), 'Uracil context missing');
  assert(!contexts.some(value => value.trim() === 'T'), 'Thymine context exposed for canonical RNA');
  await page.screenshot({ path: path.join(output, 'rna-initial-desktop.png'), fullPage: true });
  record('RNA identity and default controls', { contexts, initialRequests: initiallyFetched.length });
  if (!release.partial) await checkBroadResidueScope(page, record);

  await page.click('#universeToggle');
  assert(await page.locator('#universeDrawer').isVisible());
  assert(await page.locator('#universeTableBody tr').count() > 0);
  const entityAnnotation = await page.evaluate(() => {
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
  assert(entityAnnotation, 'Full RNA metadata has no entity-only annotation probe');
  await page.fill('#universeSearch', entityAnnotation.value);
  const searchedPdbs = await page.locator('#universeTableBody a').allTextContents();
  assert(searchedPdbs.includes(entityAnnotation.pdbId), 'Entity-scoped annotation search missed its PDB entry');
  await page.fill('#universeSearch', 'THIS_ID_DOES_NOT_EXIST');
  assert.equal(await page.locator('#universeTableBody a').count(), 0);
  await page.fill('#universeSearch', '');
  await page.click('#filteredToggle');
  assert(await page.locator('#filteredDrawer').isVisible());
  record('Entry table drawers and entity annotation search', entityAnnotation);

  // Controls capture a fresh render revision. Two unresolved renders must never
  // commit the first selection after the second one completes.
  await page.evaluate(async () => {
    const explorer = window.rnaExplorer;
    await Promise.all([explorer.setSelection({ contexts: ['A'] }), explorer.setSelection({ contexts: ['U'] })]);
  });
  await waitReady(page);
  const uracil = await verifyDistribution('rapid-change-uracil');
  assert.deepEqual(uracil.selection_spec.contexts, ['U']);
  assert(uracil.result.series.every(series => series.rows.every(row => (row.comp_id || row.context) === 'U')));
  assert(uracil.result.coverage.plottedRows > 0, 'Real release has no finite uracil chi values');

  await page.click('#resetFilters');
  await waitReady(page);
  const resetState = await page.evaluate(() => ({ selection: window.rnaExplorer.state.selection, display: window.rnaExplorer.state.display, joint: window.rnaExplorer.state.joint, family2: window.rnaExplorer.state.family2Id, survey: window.rnaExplorer.state.survey }));
  assert.deepEqual(resetState.selection.contexts, [], 'Reset retained sequence context');
  assert.equal(resetState.display.groupBy, 'base', 'Reset retained display grouping');
  assert.equal(resetState.joint.mode, 'identity', 'Reset retained joint mode');
  assert.equal(resetState.family2, '', 'Reset retained secondary family');
  assert.equal(resetState.survey.group, 'all', 'Reset retained Survey group');
  assert.equal(resetState.survey.opening, 'all', 'Reset retained Survey opening');
  record('Reset RNA filters restores defaults', resetState);
  await waitReady(page);
  await page.selectOption('#family2Select', 'backbone');
  await page.waitForFunction(() => Array.from(document.querySelector('#parameter2Select').options).some(option => option.value === 'delta'));
  await page.selectOption('#parameter2Select', 'delta');
  await waitReady(page);
  await page.waitForFunction(() => window.rnaExplorer.snapshots.joint?.result?.points?.length > 0);
  const joint = await snapshot('joint');
  const jointCsv = await downloadCsv(page, '#jointCsvDownload', path.join(output, 'joint-chi-delta.csv'));
  assert.equal(jointCsv.rows.length, joint.result.points.length);
  for (let i = 0; i < jointCsv.rows.length; i++) {
    const row = jointCsv.rows[i], point = joint.result.points[i];
    assert.equal(row.snapshot_id, joint.snapshot_id);
    assert.equal(row.x_id, row.y_id, 'Same-residue joint joined different observation identities');
    assert.equal(Number(row.x_value), point.x);
    assert.equal(Number(row.y_value), point.y);
  }
  record('Same-residue 2D plot and exact export', { points: jointCsv.rows.length, x: joint.result.xParameter.id, y: joint.result.yParameter.id });
  await checkPalettes(page, record);

  await page.selectOption('#familySelect', 'base_pair');
  await waitReady(page);
  await page.selectOption('#parameterSelect', 'opening');
  await waitReady(page);
  await page.locator('#jointJoinModeGroup button').filter({ hasText: 'Pair' }).click();
  await waitReady(page);
  await page.selectOption('#family2Select', 'backbone');
  await waitReady(page);
  await page.selectOption('#parameter2Select', 'chi');
  await waitReady(page);
  await page.waitForFunction(() => window.rnaExplorer.snapshots.joint?.result?.points?.length > 0, null, { timeout: 120000 });
  const endpointJoint = await snapshot('joint');
  assert(endpointJoint.result.points.length > 0, 'Real supported pair-to-residue view is empty');
  const endpointCsv = await downloadCsv(page, '#jointCsvDownload', path.join(output, 'joint-opening-chi-endpoints.csv'));
  assert.equal(endpointCsv.rows.length, endpointJoint.result.points.length);
  const incidences = new Set();
  for (const row of endpointCsv.rows) {
    assert(row.pair_id && row.residue_id && row.endpoint_role, 'Endpoint CSV lost explicit relationship identity');
    assert.notEqual(row.x_id, row.y_id);
    const key = JSON.stringify([row.pair_id, row.residue_id, row.endpoint_role]);
    assert(!incidences.has(key), 'Duplicate endpoint incidence'); incidences.add(key);
  }
  record('Pair-to-residue relation join', { incidences: endpointCsv.rows.length, pairs: new Set(endpointCsv.rows.map(row => row.pair_id)).size });
  await checkPairControls(page, record);

  // Exercise each released family through the ordinary controls. An intentionally
  // unavailable family is not invented to make this sweep pass.
  await page.selectOption('#family2Select', '');
  await waitReady(page);
  const families = await page.locator('#familySelect option').evaluateAll(options => options.map(option => ({ id: option.value, label: option.textContent })));
  for (const family of families) {
    await page.selectOption('#familySelect', family.id);
    await waitReady(page);
    await verifyDistribution(`family-${family.id}`);
    if (['step', 'helical', 'step_position', 'same_strand', 'helix_radius'].includes(family.id)) {
      const contexts = await page.locator('#contextGroup button').allTextContents();
      assert(contexts.length && contexts.every(value => value.trim() !== 'Unknown'), `${family.id} lost its step sequence contexts`);
      record(`Step contexts: ${family.id}`, { contexts });
    }
  }
  await page.selectOption('#familySelect', 'backbone');
  await waitReady(page);
  await page.selectOption('#parameterSelect', 'chi');
  await waitReady(page);
  record('All released family controls and exports', { families });

  const manifest = release;
  const paths = descriptor => typeof descriptor === 'string' ? [descriptor] : descriptor?.path ? [descriptor.path]
    : Object.values(descriptor?.terms || descriptor?.groups || descriptor?.partitions || {}).flatMap(paths);
  const coordinatesPaths = paths(manifest.survey?.coordinates), scalarPaths = paths(manifest.survey?.scalars);
  const fetched = (responses, paths) => responses.some(response => paths.some(path => response.url.endsWith(path)));
  assert(scalarPaths.length, 'Real release has no scalar survey');
  assert(coordinatesPaths.length, 'Real release has no coordinate survey');
  assert(!fetched(initiallyFetched, [...scalarPaths, ...coordinatesPaths]), 'Survey loaded before user request');
  await page.click('#baseGeometryLoad');
  await page.waitForFunction(() => window.rnaExplorer.snapshots.survey?.result?.series?.length > 0, null, { timeout: 120000 });
  await verifyDistribution('survey', '#surveyCsvDownload', 'survey');
  await checkPuckerSurvey(page, record);
  assert(!fetched(report.responses, coordinatesPaths), 'Scalar request eagerly loaded coordinates');
  await page.selectOption('#surveyOpeningSelect', 'bins');
  await waitReady(page);
  const conditioned = await verifyDistribution('survey-opening-conditioned', '#surveyCsvDownload', 'survey');
  assert(conditioned.result.series.every(series => ['small', 'middle', 'large'].includes(series.key)));
  await page.click('#surveyRankingLoad');
  await waitReady(page);
  record('Opening-conditioned survey and term ranking', { incidences: conditioned.result.coverage.plottedRows, rankedTerms: await page.locator('#baseGeometryRankingBody tr').count() });
  assert(!fetched(report.responses, coordinatesPaths), 'Term ranking eagerly loaded coordinate assets');
  await page.click('#coordinatesLoad');
  await page.waitForFunction(() => document.querySelector('#baseGeometryCoordBody')?.rows.length > 0, null, { timeout: 120000 });
  assert(fetched(report.responses, coordinatesPaths), 'Coordinate request did not fetch real asset');
  const coordinateHeaders = await page.locator('#coordinateBody thead th').allTextContents();
  assert.deepEqual(coordinateHeaders, ['Context', 'Atom', 'Observations', 'Residues', 'Pairs', 'PDB entries', 'Mean x (Å)', 'Mean y (Å)', 'Mean z (Å)', 'RMS spread (Å)']);
  const coordinateBinNote = await page.locator('#coordinateBinNote').textContent();
  assert.match(coordinateBinNote, /Opening bins:/, 'Coordinate opening boundaries are not explained');
  const coordinateRows = await page.locator('#baseGeometryCoordBody tr').count();
  if (coordinateRows) assert((await page.locator('#baseGeometryCoordBody tr').first().locator('td').count()) === 10, 'Coordinate table columns shifted');
  record('Scalar and coordinate surveys load independently', { scalarPartitions: scalarPaths.length, coordinatePartitions: coordinatesPaths.length, coordinateRows, coordinateHeaders, coordinateBinNote });
  await page.screenshot({ path: path.join(output, 'rna-survey-desktop.png'), fullPage: true });

  await page.setViewportSize({ width: 390, height: 844 });
  await page.screenshot({ path: path.join(output, 'rna-mobile.png'), fullPage: true });
  const dimensions = await page.evaluate(() => ({ viewport: innerWidth, body: document.body.scrollWidth, root: document.documentElement.scrollWidth }));
  assert(dimensions.root <= dimensions.viewport + 2, `Mobile document overflows horizontally: ${JSON.stringify(dimensions)}`);
  record('Mobile layout', dimensions);
  report.performance = await page.evaluate(() => ({
    longTasks: globalThis.rnaBrowserLongTasks,
    heap: performance.memory ? { used: performance.memory.usedJSHeapSize, total: performance.memory.totalJSHeapSize, limit: performance.memory.jsHeapSizeLimit } : null,
    resources: performance.getEntriesByType('resource').map(entry => ({ url: entry.name, durationMs: entry.duration, transferBytes: entry.transferSize, encodedBytes: entry.encodedBodySize, decodedBytes: entry.decodedBodySize })),
    note: 'One Chromium session on this shared host; includes the complete interaction sweep, all cached families, ranking, exports, and screenshots. Not a controlled hardware benchmark.',
  }));
  assert.deepEqual(report.pageErrors, [], 'Uncaught browser exceptions');
  assert.deepEqual(report.consoleErrors, [], 'Browser console errors');
  assert.deepEqual(report.failedRequests, [], 'Failed network requests');
  assert(report.responses.every(response => response.status < 400), 'HTTP error responses');
  if (report.releaseCandidate) assert(!report.responses.some(response => /\/assets\/pure_rna\/releases\//.test(response.url)), 'Staged release test fetched published release resources');
  if (report.packedSurveyCandidate) {
    const routing = report.packedSurveyCandidate, candidateRoot = new URL('.', routing.candidateUrl).href;
    const scalarRequests = report.responses.filter(response => response.url.startsWith(candidateRoot)
      && response.url.includes('/survey/scalars/'));
    assert(scalarRequests.length > 0, 'Packed scalar UI sweep never fetched candidate payloads');
    assert(!report.responses.some(response => /\/assets\/pure_rna\/releases\/[^/]+\/survey\/scalars\//.test(response.url)),
      'Packed scalar UI sweep fell back to original scalar payloads');
    assert(report.responses.some(response => response.url.startsWith(new URL('.', routing.sourceManifestUrl).href)
      && response.url.includes('/survey/coordinates/')), 'Mixed candidate sweep did not load original coordinates');
    report.packedSurveyIntegration = { scalarRequests: scalarRequests.length,
      uniqueScalarRequests: new Set(scalarRequests.map(response => response.url)).size,
      noOriginalScalarFallback: true, originalFamiliesAndCoordinates: true,
      limitation: 'Mixed scalar candidate and original installed release integration; not an immutable release activation gate.' };
  }
  report.passed = true;
} catch (error) {
  report.passed = false; report.failure = { message: error.message, stack: error.stack };
  await page.screenshot({ path: path.join(output, 'rna-failure.png'), fullPage: true }).catch(() => {});
  throw error;
} finally {
  if (process.env.RNA_PACKED_SURVEY_CANDIDATE_URL) report.packedModuleResponses = await Promise.all(packedModuleResponses);
  report.finishedAt = new Date().toISOString();
  await writeFile(path.join(output, 'rna-live-browser.json'), JSON.stringify(report, null, 2) + '\n');
  await browser.close();
}
