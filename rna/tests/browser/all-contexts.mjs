/** Context controls restore the exact full population, without resetting other filters. */
import assert from 'node:assert/strict';
import { mkdir, writeFile } from 'node:fs/promises';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { downloadCsv, waitReady } from './helpers.mjs';
const workspace = path.resolve(process.env.RNA_WORKSPACE || process.cwd());
const output = path.join(workspace, 'data/pure_rna/browser_validation/all-contexts');
await mkdir(output, { recursive: true });
const { chromium } = await import(pathToFileURL(process.env.PLAYWRIGHT_MODULE || path.join(workspace, 'data/pure_rna/browser-tools/node_modules/playwright/index.mjs')));
const browser = await chromium.launch({ headless: true });
const page = await browser.newPage({ acceptDownloads: true });
const report = { startedAt: new Date().toISOString(), checks: [], errors: [] };
page.on('pageerror', error => report.errors.push(error.message));
page.on('console', message => { if (message.type() === 'error') report.errors.push(message.text()); });
page.on('requestfailed', request => report.errors.push(request.failure()?.errorText));
async function click(selector) {
  await page.evaluate(() => { globalThis.contextRender = null; });
  await page.click(selector);
  await page.evaluate(async () => {
    if (!globalThis.contextRender) throw Error('Control did not request a render');
    await globalThis.contextRender;
  });
  await waitReady(page);
}
async function evidence(group, key) {
  return page.evaluate(({ group, key }) => {
    const app = window.rnaExplorer, snapshot = app.snapshots[key];
    const state = key === 'joint' ? app.state.joint.residueContexts : key === 'survey' ? app.state.survey.contexts : app.state.selection.contexts;
    const ids = key === 'joint' ? snapshot.result.points.map(p => JSON.stringify([p.pair_id, p.residue_id, p.endpoint_role, p.x, p.y]))
      : snapshot.result.series.flatMap(s => s.rowIds.map((id, i) => JSON.stringify([s.key, id, s.values[i]])));
    return { state, ids: ids.sort(), all: document.querySelector(`#${group} [data-all]`).getAttribute('aria-pressed'),
      active: Array.from(document.querySelectorAll(`#${group} [data-value][aria-pressed="true"]`), node => node.dataset.value),
      primaryId: app.snapshots.distribution.snapshot_id, selection: app.state.selection };
  }, { group, key });
}
try {
  await page.goto(process.env.RNA_URL || 'http://127.0.0.1:8767/nucleic.pages/rna/'); await waitReady(page);
  await page.evaluate(() => {
    const app = window.rnaExplorer;
    for (const name of ['requestRender', 'requestJointOnly']) {
      const original = app[name].bind(app);
      app[name] = (...args) => (globalThis.contextRender = original(...args));
    }
  });
  const baseline = await evidence('contextGroup', 'distribution');
  assert.equal(baseline.all, 'true'); assert(baseline.ids.length > 0);
  await click('#contextGroup [data-value="U"]');
  let selected = await evidence('contextGroup', 'distribution');
  assert.deepEqual(selected.state, ['U']); assert.equal(selected.all, 'false');
  assert(selected.ids.length > 0 && selected.ids.length < baseline.ids.length);
  await click('#contextGroup [data-value="A"]');
  selected = await evidence('contextGroup', 'distribution');
  assert.deepEqual(selected.state, ['U', 'A']); assert.deepEqual(selected.active, ['A', 'U']);
  await click('#contextGroup [data-all]');
  assert.deepEqual((await evidence('contextGroup', 'distribution')).ids, baseline.ids);
  await click('#contextGroup [data-value="U"]'); await click('#contextGroup [data-value="U"]');
  const restored = await evidence('contextGroup', 'distribution');
  assert.deepEqual(restored.ids, baseline.ids); assert.equal(restored.all, 'true');
  assert.deepEqual(restored.selection, baseline.selection);
  const exported = await downloadCsv(page, '#filteredCsvDownload', path.join(output, 'all-contexts.csv'));
  assert.deepEqual(exported.rows.map(row => JSON.stringify([row.group, row.id, Number(row.value)])).sort(), baseline.ids);
  report.checks.push({ name: 'Main multi-select, last deselection, All and exact CSV', rows: baseline.ids.length });
  await page.evaluate(async () => {
    const app = window.rnaExplorer;
    Object.assign(app.state, { familyId: 'base_pair', parameterId: 'opening', family2Id: 'backbone', parameter2Id: 'chi' });
    app.state.joint.mode = 'relation'; app.updateSelectors(); app.renderJointControls(); await app.requestRender();
  });
  await waitReady(page);
  const joint = await evidence('jointResidueContextGroup', 'joint'); assert(joint.ids.length > 0);
  await click('#jointResidueContextGroup [data-value="U"]');
  assert.deepEqual((await evidence('jointResidueContextGroup', 'joint')).state, ['U']);
  await click('#jointResidueContextGroup [data-all]');
  const allJoint = await evidence('jointResidueContextGroup', 'joint');
  assert.deepEqual(allJoint.ids, joint.ids); assert.equal(allJoint.primaryId, joint.primaryId); assert.equal(allJoint.all, 'true');
  report.checks.push({ name: 'Joint All restores endpoints and preserves 1D snapshot', rows: joint.ids.length });
  await page.evaluate(async () => {
    const app = window.rnaExplorer;
    app.state.familyId = 'backbone'; app.state.parameterId = 'chi'; app.state.family2Id = '';
    app.state.survey.loaded = true; app.$('baseGeometryBody').hidden = false;
    app.updateSelectors(); await app.requestRender();
  });
  await waitReady(page);
  const survey = await evidence('baseGeometryContextGroup', 'survey'); assert(survey.ids.length > 0);
  const choice = await page.locator('#baseGeometryContextGroup [data-value]').first().getAttribute('data-value');
  await click(`#baseGeometryContextGroup [data-value="${choice}"]`);
  assert.deepEqual((await evidence('baseGeometryContextGroup', 'survey')).state, [choice]);
  await click('#baseGeometryContextGroup [data-all]');
  const allSurvey = await evidence('baseGeometryContextGroup', 'survey');
  assert.deepEqual(allSurvey.ids, survey.ids); assert.equal(allSurvey.all, 'true');
  report.checks.push({ name: 'Survey All restores scalar population', rows: survey.ids.length });
  assert.deepEqual(report.errors, []); report.passed = true;
} catch (error) { report.passed = false; report.failure = error.stack; throw error; }
finally { report.finishedAt = new Date().toISOString(); await writeFile(path.join(output, 'report.json'), JSON.stringify(report, null, 2)); await browser.close(); }
console.log(JSON.stringify(report));
