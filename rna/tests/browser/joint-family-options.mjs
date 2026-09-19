/** Compatible families follow scientific levels across real selector clicks. */
import assert from 'node:assert/strict';
import { mkdir, writeFile } from 'node:fs/promises';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { waitReady } from './helpers.mjs';
const workspace = path.resolve(process.env.RNA_WORKSPACE || process.cwd());
const output = path.resolve(process.env.RNA_BROWSER_OUTPUT || path.join(workspace, 'data/pure_rna/browser_validation/joint-family-options'));
await mkdir(output, { recursive: true });
const { chromium } = await import(pathToFileURL(process.env.PLAYWRIGHT_MODULE || path.join(workspace, 'data/pure_rna/browser-tools/node_modules/playwright-core/index.mjs')).href);
const browser = await chromium.launch({ headless: true });
const page = await browser.newPage();
const report = { startedAt: new Date().toISOString(), checks: [], errors: [] };
page.on('pageerror', error => report.errors.push(error.message));
page.on('console', message => { if (message.type() === 'error') report.errors.push(message.text()); });
async function select(id, value) { await page.selectOption(`#${id}`, value); await waitReady(page); }
async function mode(value) { await page.click(`#jointJoinModeGroup button[data-value="${value}"]`); await waitReady(page); }
async function check(name, expected, plotted = false) {
  const evidence = await page.evaluate(() => {
    const app = window.rnaExplorer, snapshot = app.snapshots.joint, points = snapshot?.result.points ?? [];
    const x = snapshot?.result.xParameter, y = snapshot?.result.yParameter;
    return { familyId: app.state.familyId, family2Id: app.state.family2Id, parameter2Id: app.state.parameter2Id,
      families: [...document.querySelector('#family2Select').options].map(option => option.value).filter(Boolean),
      familyDisabled: document.querySelector('#family2Select').disabled, parameterDisabled: document.querySelector('#parameter2Select').disabled,
      exportDisabled: document.querySelector('#jointCsvDownload').disabled, snapshot: Boolean(snapshot), points: points.length,
      rawValuesMatch: points.every(point => point.x === (point.left.values?.[x.id] ?? point.left[x.id]) && point.y === (point.right.values?.[y.id] ?? point.right[y.id])),
      correctIdentities: points.every(point => app.state.joint.mode === 'identity' ? point.left_id === point.right_id : Boolean(point.pair_id && point.residue_id && point.endpoint_role)),
      x: x?.id, y: y?.id, note: document.querySelector('#jointNote').textContent };
  });
  for (const [key, value] of Object.entries(expected)) assert.deepEqual(evidence[key], value, `${name}: ${key}`);
  assert.equal(evidence.snapshot, plotted, `${name}: snapshot`); assert.equal(evidence.exportDisabled, !plotted, `${name}: export`);
  if (plotted) { assert(evidence.points > 0); assert(evidence.rawValuesMatch); assert(evidence.correctIdentities); }
  report.checks.push({ name, passed: true, ...evidence });
}
try {
  await page.goto(process.env.RNA_URL || 'http://127.0.0.1:8767/nucleic.pages/rna/'); await waitReady(page);
  const residues = ['backbone', 'pseudo_torsion', 'sugar_torsion', 'pucker', 'glycosidic_sugar_angles', 'glycosidic_base_angles', 'ribose_2oh'];
  await check('Residue identity excludes pair and step families', { families: residues, family2Id: '', parameterDisabled: true });
  await select('family2Select', 'backbone'); await select('parameter2Select', 'delta');
  await check('Same-family identity keeps finite raw values', { family2Id: 'backbone', parameter2Id: 'delta' }, true);
  await page.evaluate(() => { window.compatibilityOverviewCard = document.querySelector('#familyOverview [data-parameter="beta"]'); });
  await page.click('#familyOverview [data-parameter="beta"]'); await waitReady(page);
  await check('Overview card refresh retains compatible secondary parameter', { family2Id: 'backbone', parameter2Id: 'delta', x: 'beta' }, true);
  assert(await page.evaluate(() => window.compatibilityOverviewCard === document.querySelector('#familyOverview [data-parameter="beta"]')));
  await select('familySelect', 'pseudo_torsion');
  await check('Primary same-level change retains valid secondary parameter', { family2Id: 'backbone', parameter2Id: 'delta' }, true);
  await mode('relation');
  await check('Mode change clears stale identity axis and exports', { families: ['base_pair', 'pair_quality'], family2Id: '', parameterDisabled: true });
  await select('family2Select', 'base_pair'); await select('parameter2Select', 'opening');
  await check('Residue to pair relation uses explicit endpoints', { family2Id: 'base_pair', y: 'opening' }, true);
  await select('familySelect', 'base_pair'); await select('parameterSelect', 'opening');
  await check('Primary level change clears incompatible relation', { families: residues, family2Id: '' });
  await select('family2Select', 'backbone'); await select('parameter2Select', 'chi');
  await check('Pair to residue relation uses explicit endpoints', { family2Id: 'backbone', x: 'opening', y: 'chi' }, true);
  await select('family2Select', '');
  await check('Explicit None clears plot and export without errors', { family2Id: '', parameter2Id: '', parameterDisabled: true });
  await select('familySelect', 'step');
  await check('Step relation has no compatible family', { families: [], familyDisabled: true, parameterDisabled: true });
  assert.match(report.checks.at(-1).note, /one pair parameter and one residue parameter/);
  await mode('identity');
  await check('Identity restores compatible step families', { families: ['step', 'helical', 'step_position', 'same_strand', 'helix_radius'], familyDisabled: false });
  await select('family2Select', 'helical');
  await check('Cross-family step identity works', { family2Id: 'helical' }, true);
  await mode('relation'); await mode('identity');
  await check('Returning to identity requires explicit new secondary choice', { family2Id: '', parameter2Id: '' });
  assert.deepEqual(report.errors, []); report.passed = true;
  console.log(`PASS ${report.checks.length} joint compatibility browser checks`);
} catch (error) { report.passed = false; report.failure = error.stack; throw error; }
finally { report.finishedAt = new Date().toISOString(); await writeFile(path.join(output, 'report.json'), JSON.stringify(report, null, 2) + '\n'); await browser.close(); }
