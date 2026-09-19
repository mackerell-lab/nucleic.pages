/** Full-release endpoint checks compare UI subsets to raw deposited relationships. */
import assert from 'node:assert/strict';
import { mkdir, writeFile } from 'node:fs/promises';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { downloadCsv, waitReady } from './helpers.mjs';

const workspace = path.resolve(process.env.RNA_WORKSPACE || process.cwd());
const output = path.join(workspace, 'data/pure_rna/browser_validation/joint-residue-filters');
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
  await page.selectOption('#familySelect', 'base_pair'); await waitReady(page);
  await page.selectOption('#parameterSelect', 'opening'); await waitReady(page);
  await page.click('#jointJoinModeGroup button[data-value="relation"]'); await waitReady(page);
  await page.selectOption('#family2Select', 'backbone'); await waitReady(page);
  await page.selectOption('#parameter2Select', 'chi'); await waitReady(page);
  const baseline = await page.evaluate(async () => {
    const app = window.rnaExplorer;
    if (app.manifest.partial) throw new Error('Full RNA release required');
    const pairs = new Map((await app.repository.loadFamily('base_pair')).rows.map(row => [row.id, row]));
    const residues = new Map((await app.repository.loadFamily('backbone')).rows.map(row => [row.id, row]));
    const relationTable = await app.repository.loadRelations('observations');
    const links = Array.isArray(relationTable) ? relationTable : relationTable.rows;
    const identity = link => JSON.stringify([link.pair_id, link.residue_id, link.endpoint_role ?? (link.side === 1 ? 'first' : 'second')]);
    const rawRelations = new Set(links.filter(link => link.kind === 'pair_residue').map(identity));
    const points = app.snapshots.joint.result.points;
    if (!points.length) throw new Error('Empty initial endpoint plot');
    for (const point of points) {
      if (!rawRelations.has(identity(point))) throw new Error('Invented endpoint relationship');
      if (point.x !== pairs.get(point.pair_id).values.opening || point.y !== residues.get(point.residue_id).values.chi) throw new Error('Raw value mismatch');
    }
    globalThis.endpointBaseline = points.map(point => ({ pair: point.pair_id, residue: point.residue_id, side: point.endpoint_role,
      base: residues.get(point.residue_id).comp_id, pucker: residues.get(point.residue_id).pucker_class, x: point.x, y: point.y }));
    globalThis.primaryBaseline = app.snapshots.distribution.result.series.flatMap(series => series.rowIds).sort().join('\n');
    return { build: app.manifest.build_id, points: points.length };
  });
  report.checks.push({ name: 'Raw relation and scalar identities', ...baseline });
  await page.click('#jointResidueContextGroup button[data-value="U"]'); await waitReady(page);
  const pucker = await page.evaluate(() => globalThis.endpointBaseline.find(point => point.base === 'U' && point.pucker)?.pucker);
  assert(pucker, 'No real uracil pucker observations');
  await page.selectOption('#jointResiduePuckerGroup', pucker); await waitReady(page);
  for (const [choice, side] of [['both', null], ['nt1', 'first'], ['nt2', 'second']]) {
    await page.click(`#jointResidueSideGroup button[data-value="${choice}"]`); await waitReady(page);
    const evidence = await page.evaluate(({ side, pucker }) => {
      const app = window.rnaExplorer, snapshot = app.snapshots.joint;
      const expected = globalThis.endpointBaseline.filter(point => point.base === 'U' && point.pucker === pucker && (!side || point.side === side));
      const actual = snapshot.result.points.map(point => ({ pair: point.pair_id, residue: point.residue_id, side: point.endpoint_role, base: point.right.comp_id, pucker: point.right.pucker_class, x: point.x, y: point.y }));
      return { expected, actual, oneDimensionalUnchanged: globalThis.primaryBaseline === app.snapshots.distribution.result.series.flatMap(series => series.rowIds).sort().join('\n'),
        axis: snapshot.provenance.axis_selections.y, snapshotId: snapshot.snapshot_id };
    }, { side, pucker });
    assert.deepEqual(evidence.actual, evidence.expected); assert(evidence.oneDimensionalUnchanged);
    assert.deepEqual(evidence.axis.contexts, ['U']); assert.deepEqual(evidence.axis.puckerStates, [pucker]);
    const exported = await downloadCsv(page, '#jointCsvDownload', path.join(output, `${choice}.csv`));
    assert.equal(exported.rows.length, evidence.expected.length);
    exported.rows.forEach((row, index) => {
      const expected = evidence.expected[index];
      assert.equal(row.pair_id, expected.pair); assert.equal(row.residue_id, expected.residue);
      assert.equal(Number(row.x_value), expected.x); assert.equal(Number(row.y_value), expected.y);
      assert.equal(row.snapshot_id, evidence.snapshotId);
    });
    report.checks.push({ name: `Uracil pucker ${choice} subset and CSV`, points: evidence.actual.length, pucker });
  }
  assert.equal(await page.locator('#family2Select option[value="base_pair"]').count(), 0);
  await page.selectOption('#family2Select', ''); await waitReady(page);
  assert(await page.locator('#jointCsvDownload').isDisabled());
  assert((await page.locator('#jointPlot').textContent()).includes('Select a second parameter'));
  assert.equal(await page.evaluate(() => window.rnaExplorer.snapshots.joint), null);
  report.checks.push({ name: 'Incompatible relation unavailable and None disables stale export' });
  assert.deepEqual(report.errors, []); report.passed = true;
} catch (error) { report.passed = false; report.failure = error.stack; throw error; }
finally { report.finishedAt = new Date().toISOString(); await writeFile(path.join(output, 'report.json'), JSON.stringify(report, null, 2)); await browser.close(); }
console.log(JSON.stringify(report, null, 2));
