/** Static selector disclosures preserve labels, selections, and plotted snapshots. */
import assert from 'node:assert/strict';
import { mkdir, writeFile } from 'node:fs/promises';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { waitReady } from './helpers.mjs';

const workspace = path.resolve(process.env.RNA_WORKSPACE || process.cwd());
const output = path.resolve(process.env.RNA_BROWSER_OUTPUT || path.join(workspace, 'data/pure_rna/browser_validation/static-control-help'));
await mkdir(output, { recursive: true });
const { chromium } = await import(pathToFileURL(process.env.PLAYWRIGHT_MODULE || path.join(workspace, 'data/pure_rna/browser-tools/node_modules/playwright-core/index.mjs')).href);
const browser = await chromium.launch({ headless: true });
const report = { startedAt: new Date().toISOString(), checks: [], pageErrors: [], consoleErrors: [] };
const controls = [
  ['familySelect', 'Family', ['residue, pair, or step', 'finite values']],
  ['parameterSelect', 'Parameter', ['missing or unsupported', 'not independent RNA molecules']],
  ['family2Select', 'Family 2', ['Explicit identities', 'endpoint links']],
  ['parameter2Select', 'Parameter 2', ['Pearson r and r²', 'mixed circular/linear', 'statistical significance']],
  ['surveyOpeningSelect', 'Opening conditioning', ['declared pair-opening boundaries', 'explicit pair endpoint', 'not RNA conformation thresholds']],
];

async function capture(page) {
  await page.evaluate(() => {
    const app = window.rnaExplorer;
    window.staticHelpBaseline = { state: JSON.stringify(app.state), revision: app.revision, snapshots: { ...app.snapshots },
      plots: ['plot', 'jointPlot', 'baseGeometryPlot'].map(id => document.getElementById(id).data),
      selects: [...document.querySelectorAll('[data-rna-control-help] select')].map(node => ({ node, value: node.value, label: node.labels[0] })) };
  });
}

async function unchanged(page) {
  const evidence = await page.evaluate(() => {
    const app = window.rnaExplorer, baseline = window.staticHelpBaseline;
    return { stateUnchanged: baseline.state === JSON.stringify(app.state), revisionUnchanged: baseline.revision === app.revision,
      snapshotsUnchanged: Object.keys(baseline.snapshots).every(key => baseline.snapshots[key] === app.snapshots[key]),
      plotsUnchanged: ['plot', 'jointPlot', 'baseGeometryPlot'].every((id, i) => document.getElementById(id).data === baseline.plots[i]),
      selectorsPreserved: baseline.selects.every(({ node, value, label }) => node.isConnected && node.value === value && node.labels[0] === label),
      noOverflow: document.documentElement.scrollWidth <= innerWidth };
  });
  assert(Object.values(evidence).every(Boolean), JSON.stringify(evidence));
  return evidence;
}

try {
  for (const mobile of [false, true]) {
    const context = await browser.newContext({ viewport: mobile ? { width: 390, height: 844 } : { width: 1440, height: 1000 }, hasTouch: mobile, isMobile: mobile });
    const page = await context.newPage();
    page.on('pageerror', error => report.pageErrors.push(error.message));
    page.on('console', message => { if (message.type() === 'error') report.consoleErrors.push(message.text()); });
    await page.goto(process.env.RNA_URL || 'http://127.0.0.1:8767/nucleic.pages/rna/', { waitUntil: 'domcontentloaded', timeout: 120000 });
    await waitReady(page);
    report.build = await page.evaluate(() => window.rnaExplorer.manifest.build_id);
    await page.locator('#baseGeometryLoad').click(); await waitReady(page);
    assert(await page.evaluate(() => Boolean(window.rnaExplorer.snapshots.distribution && window.rnaExplorer.snapshots.survey)));
    const cdp = await context.newCDPSession(page);
    const accessibility = await cdp.send('Accessibility.getFullAXTree');
    for (const [id, title, fragments] of controls) {
      const cluster = page.locator(`[data-rna-control-help="${id}"]`), disclosure = cluster.locator('details.rna-control-help');
      assert.equal(await disclosure.count(), 1);
      assert.equal(await cluster.locator('label details').count(), 0, 'Help must not activate the select label');
      const summary = disclosure.locator('summary'), note = disclosure.locator('p');
      assert.equal(await summary.getAttribute('aria-label'), `${title} help`);
      assert(accessibility.nodes.some(node => node.name?.value === `${title} help` && ['DisclosureTriangle', 'button'].includes(node.role?.value)));
      for (const fragment of fragments) assert((await note.textContent()).includes(fragment));
      await capture(page);
      if (mobile) await summary.tap();
      else { await summary.focus(); await page.keyboard.press('Enter'); }
      assert(await disclosure.evaluate(node => node.open)); assert(await note.isVisible());
      assert(await disclosure.evaluate(node => { const box = node.getBoundingClientRect(); return box.left >= 0 && box.right <= innerWidth; }));
      const evidence = await unchanged(page);
      if (mobile) await summary.tap();
      else {
        await page.keyboard.press('Space');
        assert.equal(await disclosure.evaluate(node => node.open), false);
        await page.keyboard.press('Tab'); assert(await summary.evaluate(node => document.activeElement !== node));
        await page.keyboard.press('Shift+Tab'); assert(await summary.evaluate(node => document.activeElement === node));
        await summary.click(); assert(await disclosure.evaluate(node => node.open));
        await summary.click();
      }
      assert.equal(await disclosure.evaluate(node => node.open), false); await unchanged(page);
      report.checks.push({ id, mobile, passed: true, ...evidence });
    }
    await capture(page);
    await page.evaluate(async () => {
      const { mountStaticControlHelp } = await import('./views/panels.js');
      mountStaticControlHelp(document.getElementById('rnaExplorer'));
      mountStaticControlHelp(document.getElementById('rnaExplorer'));
    });
    assert.equal(await page.locator('[data-rna-control-help] details.rna-control-help').count(), 5);
    report.checks.push({ name: 'Idempotent help mounting', mobile, passed: true, ...await unchanged(page) });
    // Real selection and reset actions prove existing listeners still work.
    await page.selectOption('#familySelect', 'ribose_2oh'); await waitReady(page);
    assert.equal(await page.evaluate(() => window.rnaExplorer.state.familyId), 'ribose_2oh');
    await page.selectOption('#family2Select', 'ribose_2oh'); await waitReady(page);
    assert(await page.evaluate(() => Boolean(window.rnaExplorer.snapshots.joint)));
    await page.locator('#resetFilters').click(); await waitReady(page);
    assert.equal(await page.locator('[data-rna-control-help] details.rna-control-help').count(), 5);
    assert(await page.evaluate(() => window.rnaExplorer.state.familyId === 'ribose_2oh' && window.rnaExplorer.state.family2Id === '' && window.rnaExplorer.state.survey.opening === 'all'));
    report.checks.push({ name: 'Selectors and reset preserve help', mobile, passed: true });
    await page.screenshot({ path: path.join(output, mobile ? 'static-help-mobile.png' : 'static-help-desktop.png'), fullPage: false });
    await context.close();
  }
  assert.deepEqual(report.pageErrors, []); assert.deepEqual(report.consoleErrors, []);
  report.passed = true;
  console.log(`PASS ${report.checks.length} static selector help checks`);
} catch (error) { report.passed = false; report.failure = error.stack; throw error; }
finally {
  report.finishedAt = new Date().toISOString();
  await writeFile(path.join(output, 'rna-static-control-help.json'), JSON.stringify(report, null, 2) + '\n');
  await browser.close();
}
