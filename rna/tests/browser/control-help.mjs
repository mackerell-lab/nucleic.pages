/** Native help disclosures work with keyboard, pointer and touch without filtering. */
import assert from 'node:assert/strict';
import { mkdir, writeFile } from 'node:fs/promises';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { waitReady } from './helpers.mjs';
const workspace = path.resolve(process.env.RNA_WORKSPACE || process.cwd());
const output = path.resolve(process.env.RNA_BROWSER_OUTPUT || path.join(workspace, 'data/pure_rna/browser_validation/control-help'));
await mkdir(output, { recursive: true });
const { chromium } = await import(pathToFileURL(process.env.PLAYWRIGHT_MODULE || path.join(workspace, 'data/pure_rna/browser-tools/node_modules/playwright-core/index.mjs')).href);
const browser = await chromium.launch({ headless: true });
const report = { startedAt: new Date().toISOString(), checks: [], pageErrors: [] };
try {
  for (const mobile of [false, true]) {
    const context = await browser.newContext({ viewport: mobile ? { width: 390, height: 844 } : { width: 1440, height: 1000 }, hasTouch: mobile, isMobile: mobile });
    const page = await context.newPage(); page.on('pageerror', error => report.pageErrors.push(error.message));
    await page.goto(process.env.RNA_URL || 'http://127.0.0.1:8767/nucleic.pages/rna/', { waitUntil: 'domcontentloaded', timeout: 120000 }); await waitReady(page);
    const summary = page.locator('summary[aria-label="Cleanliness help"]');
    assert.equal(await summary.count(), 1, 'Help needs a named keyboard-focusable disclosure');
    const cdp = await context.newCDPSession(page);
    const accessibility = await cdp.send('Accessibility.getFullAXTree');
    assert(accessibility.nodes.some(node => node.name?.value === 'Cleanliness help' && ['DisclosureTriangle', 'button'].includes(node.role?.value)), 'Disclosure has no accessible name/role');
    const disclosure = page.locator('#cleanlinessGroup').locator('..').locator('details.rna-control-help');
    const note = disclosure.locator('p');
    await page.evaluate(() => { const app = window.rnaExplorer; window.helpBaseline = { revision: app.revision, snapshot: app.snapshots.distribution, selection: JSON.stringify(app.state.selection) }; });
    assert.equal(await disclosure.evaluate(node => node.open), false);
    if (mobile) {
      await summary.tap(); assert.equal(await disclosure.evaluate(node => node.open), true); assert(await note.isVisible());
      await summary.tap(); assert.equal(await disclosure.evaluate(node => node.open), false);
    } else {
      await summary.focus(); await page.keyboard.press('Enter'); assert.equal(await disclosure.evaluate(node => node.open), true); assert(await note.isVisible());
      assert(await summary.evaluate(node => document.activeElement === node));
      await page.keyboard.press('Space'); assert.equal(await disclosure.evaluate(node => node.open), false);
      assert(await summary.evaluate(node => document.activeElement === node));
      await page.keyboard.press('Space'); assert.equal(await disclosure.evaluate(node => node.open), true);
      await page.keyboard.press('Enter'); assert.equal(await disclosure.evaluate(node => node.open), false);
      await page.keyboard.press('Tab'); assert(await summary.evaluate(node => document.activeElement !== node), 'Disclosure traps keyboard focus');
      await page.keyboard.press('Shift+Tab'); assert(await summary.evaluate(node => document.activeElement === node), 'Disclosure absent from sequential keyboard order');
      await summary.click(); assert.equal(await disclosure.evaluate(node => node.open), true);
      await summary.click(); assert.equal(await disclosure.evaluate(node => node.open), false);
    }
    // All current explanations must be present as readable text, including select controls.
    const coverage = await page.evaluate(() => [...document.querySelectorAll('.filter-cluster')].filter(cluster => cluster.querySelector('.cluster-title[title]')).map(cluster => {
      const label = cluster.querySelector('.cluster-title'), details = cluster.querySelector('.rna-control-help');
      if (details) details.open = true;
      return { title: label.textContent, matched: details?.querySelector('p')?.textContent === label.title, named: details?.querySelector('summary')?.getAttribute('aria-label') === `${label.textContent} help` };
    }));
    assert(coverage.length > 5); assert(coverage.every(item => item.matched && item.named));
    const evidence = await page.evaluate(() => {
      const app = window.rnaExplorer, baseline = window.helpBaseline;
      return { viewport: innerWidth, documentWidth: document.documentElement.scrollWidth, helpCount: document.querySelectorAll('.rna-control-help').length,
        stateUnchanged: baseline.selection === JSON.stringify(app.state.selection), revisionUnchanged: baseline.revision === app.revision,
        snapshotUnchanged: baseline.snapshot === app.snapshots.distribution,
        helpFits: [...document.querySelectorAll('.rna-control-help')].every(node => { const r = node.getBoundingClientRect(); return !r.width || r.left >= 0 && r.right <= innerWidth; }) };
    });
    assert(evidence.stateUnchanged && evidence.revisionUnchanged && evidence.snapshotUnchanged); assert(evidence.helpFits); assert(evidence.documentWidth <= evidence.viewport);
    report.checks.push({ passed: true, mobile, coverage, ...evidence });
    await page.screenshot({ path: path.join(output, mobile ? 'help-mobile.png' : 'help-desktop.png'), fullPage: false });
    await context.close();
  }
  assert.deepEqual(report.pageErrors, []); report.passed = true; console.log('PASS keyboard, pointer and touch control help');
} catch (error) { report.passed = false; report.failure = error.stack; throw error; }
finally { report.finishedAt = new Date().toISOString(); await writeFile(path.join(output, 'rna-control-help.json'), JSON.stringify(report, null, 2) + '\n'); await browser.close(); }
