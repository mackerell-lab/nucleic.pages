/** Scientific help survives dynamic controls without changing the plotted population. */
import assert from 'node:assert/strict';
import { mkdir, writeFile } from 'node:fs/promises';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { waitReady } from './helpers.mjs';

const workspace = path.resolve(process.env.RNA_WORKSPACE || process.cwd());
const output = path.resolve(process.env.RNA_BROWSER_OUTPUT || path.join(workspace, 'data/pure_rna/browser_validation/control-help-content'));
await mkdir(output, { recursive: true });
const { chromium } = await import(pathToFileURL(process.env.PLAYWRIGHT_MODULE || path.join(workspace, 'data/pure_rna/browser-tools/node_modules/playwright-core/index.mjs')).href);
const browser = await chromium.launch({ headless: true });
const report = { startedAt: new Date().toISOString(), checks: [], pageErrors: [] };

async function verifyHelp(page, id, fragments, mobile, phase) {
  const group = page.locator(`#${id}`), disclosure = group.locator('..').locator('details.rna-control-help');
  assert.equal(await disclosure.count(), 1, `${id} needs exactly one native help disclosure`);
  const summary = disclosure.locator('summary'), note = disclosure.locator('p');
  const text = await note.textContent();
  for (const fragment of fragments) assert(text.includes(fragment), `${id} is missing scientific meaning: ${fragment}`);
  assert.equal(await disclosure.evaluate(node => node.open), false);
  await page.evaluate(() => {
    const app = window.rnaExplorer;
    window.rnaHelpContentBaseline = { state: JSON.stringify(app.state), revision: app.revision,
      distribution: app.snapshots.distribution, joint: app.snapshots.joint, plotData: document.querySelector('#plot').data };
  });
  if (mobile) await summary.tap();
  else { await summary.focus(); await page.keyboard.press('Enter'); }
  assert(await disclosure.evaluate(node => node.open));
  assert(await note.isVisible());
  assert((await summary.getAttribute('aria-label'))?.endsWith(' help'));
  const evidence = await page.evaluate(id => {
    const app = window.rnaExplorer, before = window.rnaHelpContentBaseline;
    const cluster = document.getElementById(id).parentElement;
    const rect = cluster.querySelector('.rna-control-help').getBoundingClientRect();
    return { stateUnchanged: before.state === JSON.stringify(app.state), revisionUnchanged: before.revision === app.revision,
      distributionUnchanged: before.distribution === app.snapshots.distribution, jointUnchanged: before.joint === app.snapshots.joint,
      plotUnchanged: before.plotData === document.querySelector('#plot').data,
      withinViewport: rect.left >= 0 && rect.right <= innerWidth,
      renderedText: cluster.querySelector('.rna-control-help p').textContent };
  }, id);
  assert(evidence.stateUnchanged && evidence.revisionUnchanged && evidence.distributionUnchanged && evidence.jointUnchanged && evidence.plotUnchanged);
  assert(evidence.withinViewport, `${id} help exceeds viewport`);
  if (mobile) await summary.tap();
  else await page.keyboard.press('Space');
  assert.equal(await disclosure.evaluate(node => node.open), false);
  report.checks.push({ id, mobile, phase, passed: true, ...evidence });
}

try {
  for (const mobile of [false, true]) {
    const context = await browser.newContext({ viewport: mobile ? { width: 390, height: 844 } : { width: 1440, height: 1000 }, hasTouch: mobile, isMobile: mobile });
    const page = await context.newPage();
    page.on('pageerror', error => report.pageErrors.push(error.message));
    await page.goto(process.env.RNA_URL || 'http://127.0.0.1:8767/nucleic.pages/rna/', { waitUntil: 'domcontentloaded', timeout: 120000 });
    await waitReady(page);
    report.build = await page.evaluate(() => window.rnaExplorer.manifest.build_id);
    const controls = [
      ['circularModeGroup', ['declared circular', 'Bond angles', 'remain linear', 'Raw measurements are unchanged']],
      ['displayScaleGroup', ['Probability sums to 1', 'bin width in 1D or bin area in 2D', 'integral is 1']],
      ['binDetailGroup', ['64 linear or 72 circular', '36 bins per axis', 'Fine uses 72']],
      ['jointJoinModeGroup', ['explicit IDs', 'recorded endpoint links', 'not independent samples']],
      ['jointColorScaleGroup', ['log₁₀', '10⁻⁸', 'including zero']],
      ['groupingGroup', ['Each curve is normalized separately', 'can overlap', 'double-count']],
      ['contextGroup', ['base, pair, or step contexts', 'does not establish observation identity']],
    ];
    for (const [id, fragments] of controls) await verifyHelp(page, id, fragments, mobile, 'initial');
    await page.evaluate(() => { window.rnaHelpOldControls = { context: document.getElementById('contextGroup'), grouping: document.getElementById('groupingGroup') }; });
    await page.selectOption('#familySelect', 'ribose_2oh');
    await waitReady(page);
    assert(await page.evaluate(() => window.rnaHelpOldControls.context !== document.getElementById('contextGroup')
      && window.rnaHelpOldControls.grouping !== document.getElementById('groupingGroup')), 'Family change did not exercise rebuilt controls');
    for (const [id, fragments] of controls.filter(([id]) => ['contextGroup', 'groupingGroup'].includes(id))) {
      await verifyHelp(page, id, fragments, mobile, 'family-change');
    }
    await context.close();
  }
  assert.deepEqual(report.pageErrors, []);
  report.passed = true;
  console.log(`PASS ${report.checks.length} scientific help checks across desktop and touch`);
} catch (error) { report.passed = false; report.failure = error.stack; throw error; }
finally {
  report.finishedAt = new Date().toISOString();
  await writeFile(path.join(output, 'rna-control-help-content.json'), JSON.stringify(report, null, 2) + '\n');
  await browser.close();
}
