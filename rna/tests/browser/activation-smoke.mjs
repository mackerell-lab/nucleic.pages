/** Unmocked smoke of the installed RNA pointer and retained previous release. */
import assert from 'node:assert/strict';
import { mkdir, writeFile } from 'node:fs/promises';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { downloadCsv, waitReady } from './helpers.mjs';

const expectedBuild = process.env.RNA_EXPECT_BUILD_ID;
const previousBuild = process.env.RNA_PREVIOUS_BUILD_ID;
const safeBuild = value => typeof value === 'string' && /^[A-Za-z0-9][A-Za-z0-9_-]*$/.test(value);
assert(safeBuild(expectedBuild), 'RNA_EXPECT_BUILD_ID must identify the newly installed release');
assert(!previousBuild || (safeBuild(previousBuild) && previousBuild !== expectedBuild), 'Previous build must be a distinct safe build ID');
assert(!process.env.RNA_RELEASE_CANDIDATE_URL && !process.env.RNA_SURVEY_CANDIDATE_URL, 'Activation smoke requires the real installed pointer');
const workspace = path.resolve(process.env.RNA_WORKSPACE || process.cwd());
const output = path.resolve(process.env.RNA_BROWSER_OUTPUT || path.join(workspace, 'data/pure_rna/browser_validation/activation-smoke'));
await mkdir(output, { recursive: true });
const { chromium } = await import(pathToFileURL(process.env.PLAYWRIGHT_MODULE || path.join(workspace, 'data/pure_rna/browser-tools/node_modules/playwright/index.mjs')).href);
const report = { startedAt: new Date().toISOString(), url: process.env.RNA_URL || 'http://127.0.0.1:8767/nucleic.pages/rna/',
  expectedBuildId: expectedBuild, previousBuildId: previousBuild ?? null, checks: [], responses: [], pageErrors: [], consoleErrors: [], failedRequests: [] };
const browser = await chromium.launch({ headless: true });
const page = await browser.newPage({ viewport: { width: 1440, height: 1000 }, acceptDownloads: true });
page.on('pageerror', error => report.pageErrors.push(error.message));
page.on('console', message => { if (message.type() === 'error') report.consoleErrors.push(message.text()); });
page.on('requestfailed', request => report.failedRequests.push({ url: request.url(), error: request.failure()?.errorText }));
page.on('response', response => report.responses.push({ url: response.url(), status: response.status() }));
const record = (name, evidence) => { report.checks.push({ name, passed: true, evidence }); console.log(`PASS ${name}`); };

try {
  await page.goto(report.url, { waitUntil: 'domcontentloaded', timeout: 120000 });
  await waitReady(page);
  const loaded = await page.evaluate(() => {
    const app = window.rnaExplorer;
    return { buildId: app.manifest.build_id, partial: app.manifest.partial, counts: app.manifest.counts,
      releaseUrl: app.repository.releaseUrl, snapshotBuild: app.snapshots.distribution?.build_id,
      plottedRows: app.snapshots.distribution?.result.coverage.plottedRows };
  });
  const installedBase = new URL(`../assets/pure_rna/releases/${expectedBuild}/`, report.url).href;
  const pointerUrl = new URL('../assets/pure_rna/manifest.json', report.url).href;
  assert.equal(loaded.buildId, expectedBuild);
  assert.equal(loaded.snapshotBuild, expectedBuild);
  assert.equal(loaded.partial, false);
  assert.equal(loaded.releaseUrl, new URL('manifest.json', installedBase).href);
  assert(loaded.counts.entries > 0 && loaded.plottedRows > 0);
  record('Unmocked installed pointer and initial distribution', loaded);

  await page.selectOption('#familySelect', 'base_pair');
  await waitReady(page);
  const family = await page.evaluate(() => ({ selected: window.rnaExplorer.state.familyId,
    snapshotBuild: window.rnaExplorer.snapshots.distribution.build_id,
    rows: window.rnaExplorer.snapshots.distribution.result.coverage.plottedRows }));
  assert.equal(family.selected, 'base_pair');
  assert.equal(family.snapshotBuild, expectedBuild);
  assert(family.rows > 0);
  assert(report.responses.some(item => item.url.startsWith(installedBase) && /\/families\/base_pair\./.test(item.url)), 'Additional family did not load from installed release');
  record('Additional family loads through ordinary controls', family);

  await page.click('#baseGeometryLoad');
  await waitReady(page);
  await page.waitForFunction(() => window.rnaExplorer.snapshots.survey?.result?.series?.length > 0, null, { timeout: 120000 });
  const survey = await page.evaluate(() => {
    const snapshot = window.rnaExplorer.snapshots.survey;
    return { buildId: snapshot.build_id, snapshotId: snapshot.snapshot_id, term: snapshot.result.parameter.id,
      rows: snapshot.result.series.flatMap(series => series.values.map((value, index) => ({ id: series.rowIds[index], value, group: series.key }))) };
  });
  const csv = await downloadCsv(page, '#surveyCsvDownload', path.join(output, 'activation-survey.csv'));
  assert.equal(survey.buildId, expectedBuild);
  assert(survey.rows.length > 0);
  assert.equal(csv.rows.length, survey.rows.length);
  for (let index = 0; index < csv.rows.length; index++) {
    const row = csv.rows[index], expected = survey.rows[index];
    assert.equal(row.build_id, expectedBuild);
    assert.equal(row.snapshot_id, survey.snapshotId);
    assert.equal(row.id, expected.id);
    assert.equal(row.group, expected.group);
    assert.equal(Number(row.value), expected.value);
  }
  record('Installed Survey values and exact CSV', { buildId: survey.buildId, term: survey.term, rows: survey.rows.length });
  assert(!report.responses.some(item => /\/survey\/coordinates\//.test(item.url)), 'Coordinates loaded before request');

  await page.click('#coordinatesLoad');
  await waitReady(page);
  await page.waitForFunction(() => window.rnaExplorer.coordinateSummary?.length > 0, null, { timeout: 120000 });
  const coordinates = await page.evaluate(() => ({ group: window.rnaExplorer.state.survey.coordinateGroup,
    atoms: window.rnaExplorer.coordinateSummary.length, tableRows: document.querySelector('#baseGeometryCoordBody').rows.length }));
  const coordinateRequests = report.responses.filter(item => item.url.startsWith(installedBase) && /\/survey\/coordinates\//.test(item.url));
  assert(coordinateRequests.length > 0, 'Coordinates did not load installed partitions');
  assert.equal(coordinates.tableRows, coordinates.atoms);
  record('Installed coordinate partitions and populated table', { ...coordinates, partitions: coordinateRequests.length });

  const installedRequests = report.responses.filter(item => item.url.includes('/assets/pure_rna/'));
  assert(installedRequests.length > 0);
  assert(installedRequests.every(item => item.url === pointerUrl || item.url.startsWith(installedBase)), 'Current explorer mixed release directories');
  if (previousBuild) {
    const requestStart = report.responses.length;
    const previous = await page.evaluate(async ({ buildId, term }) => {
      const { RnaDataRepository } = await import(new URL('./core/repository.js', location.href).href);
      const manifestUrl = new URL(`../${buildId}/manifest.json`, window.rnaExplorer.repository.releaseUrl).href;
      // A new repository has no cached term; loading it proves old URLs survive.
      const repository = new RnaDataRepository({ manifestUrl });
      const manifest = await repository.loadManifest();
      const table = await repository.loadSurveyScalars(term, { fields: ['id', 'value', 'status'] });
      return { buildId: manifest.build_id, tableBuildId: table.build_id, releaseUrl: repository.releaseUrl,
        rows: table.rows.length, expectedRows: manifest.survey.scalars.terms[term].row_count,
        currentExplorerBuildId: window.rnaExplorer.manifest.build_id };
    }, { buildId: previousBuild, term: survey.term });
    const previousBase = new URL(`../${previousBuild}/`, installedBase).href;
    assert.equal(previous.buildId, previousBuild);
    assert.equal(previous.tableBuildId, previousBuild);
    assert.equal(previous.rows, previous.expectedRows);
    assert(previous.rows > 0);
    assert.equal(previous.currentExplorerBuildId, expectedBuild);
    const previousRequests = report.responses.slice(requestStart).filter(item => item.url.includes('/assets/pure_rna/'));
    assert(previousRequests.some(item => item.url.includes('/survey/scalars/')));
    assert(previousRequests.every(item => item.url.startsWith(previousBase)), 'Previous repository mixed release directories');
    record('Previous immutable release still supports uncached lazy loading', previous);
  }
  assert(!report.responses.some(item => new URL(item.url).pathname.startsWith('/data/')), 'Smoke fetched staging resources');
  assert(!report.responses.some(item => item.url.includes('/assets/pure_dna/')), 'Smoke fetched DNA assets');
  assert(report.responses.every(item => item.status < 400), 'HTTP errors during activation smoke');
  assert.deepEqual(report.pageErrors, []);
  assert.deepEqual(report.consoleErrors, []);
  assert.deepEqual(report.failedRequests, []);
  report.passed = true;
} catch (error) {
  report.passed = false;
  report.failure = { message: error.message, stack: error.stack };
  throw error;
} finally {
  report.finishedAt = new Date().toISOString();
  await writeFile(path.join(output, 'rna-activation-smoke.json'), JSON.stringify(report, null, 2) + '\n');
  await browser.close();
}
