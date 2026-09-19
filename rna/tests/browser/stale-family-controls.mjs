/** Family-owned context/pucker controls must not install invisible filters. */
import assert from 'node:assert/strict';
import { mkdir, writeFile } from 'node:fs/promises';
import { pathToFileURL } from 'node:url';
import { waitReady } from './helpers.mjs';

const { chromium } = await import(pathToFileURL(process.env.PLAYWRIGHT_MODULE).href);
const browser = await chromium.launch({ headless: true });
const report = { startedAt: new Date().toISOString(), checks: [], pageErrors: [] };
const url = process.env.RNA_URL || 'http://127.0.0.1:8767/nucleic.pages/rna/';
async function failedSwitch(source, target) {
  const page = await browser.newPage();
  page.on('pageerror', error => report.pageErrors.push(error.message));
  await page.goto(url); await waitReady(page);
  if (source !== 'backbone') { await page.selectOption('#familySelect', source); await waitReady(page); }
  const targetUrl = await page.evaluate(target => {
    const app = window.rnaExplorer;
    return new URL(app.families.find(family => family.id === target).path, app.repository.releaseUrl).href;
  }, target);
  await page.route(targetUrl, route => route.fulfill({ status: 503, body: 'Deliberate family-control regression failure' }), { times: 1 });
  await page.selectOption('#familySelect', target);
  await page.waitForFunction(() => document.querySelector('#appStatus').dataset.state === 'error');
  assert(await page.locator('#filteredCsvDownload').isDisabled());
  return page;
}
async function selection(page) {
  return page.evaluate(() => {
    const app = window.rnaExplorer;
    return { family: app.state.familyId, contexts: app.state.selection.contexts,
      puckers: app.state.selection.puckerStates, revision: app.revision };
  });
}
async function recovered(page, family) {
  await waitReady(page);
  assert.equal((await selection(page)).family, family);
  assert(await page.locator('#filteredCsvDownload').isEnabled());
  assert(await page.evaluate(() => {
    const app = window.rnaExplorer;
    return app.state.parameterId === app.snapshots.distribution.result.parameter.id
      && Number(document.querySelector('#filteredObservationCount').textContent.replaceAll(',', '')) > 0;
  }));
}
try {
  for (const repair of ['trace', 'reset', 'all']) {
    const page = await failedSwitch('base_pair', 'ribose_2oh');
    const before = await selection(page);
    await page.evaluate(() => {
      window.detachedContextButtons = [document.querySelector('#contextGroup button[data-value="A-A"]'), document.querySelector('#contextGroup button[data-all]')];
    });
    await page.click('#contextGroup button[data-value="A-A"]');
    assert.deepEqual(await selection(page), before, 'Old pair context changed the residue-family selection');
    assert.equal(await page.locator('#contextGroup button[data-all]').getAttribute('aria-pressed'), 'true');
    assert.equal(await page.locator('#contextGroup button[data-value="A-A"]').getAttribute('aria-pressed'), 'false');
    assert(await page.locator('#filteredCsvDownload').isDisabled());
    if (repair === 'trace') await page.click('#traceStyleGroup button[data-value="line"]');
    if (repair === 'reset') await page.click('#resetFilters');
    if (repair === 'all') await page.click('#contextGroup button[data-all]');
    await recovered(page, 'ribose_2oh');
    assert.deepEqual((await selection(page)).contexts, []);
    assert.equal(await page.locator('#contextGroup button[data-all]').getAttribute('aria-pressed'), 'true');
    const beforeDetached = await selection(page);
    assert(await page.evaluate(() => {
      const current = document.querySelector('#contextGroup');
      for (const button of window.detachedContextButtons) {
        if (button.isConnected) return false;
        button.click();
      }
      return current === document.querySelector('#contextGroup');
    }), 'Detached context callback replaced current controls');
    assert.deepEqual(await selection(page), beforeDetached, 'Detached context callback mutated selection or revision');
    if (repair === 'trace') {
      await page.evaluate(() => { window.savedContextButton = document.querySelector('#contextGroup button[data-value="U"]'); });
      await page.click('#traceStyleGroup button[data-value="filled"]'); await waitReady(page);
      assert(await page.evaluate(() => window.savedContextButton === document.querySelector('#contextGroup button[data-value="U"]')));
      await page.click('#contextGroup button[data-value="U"]'); await recovered(page, 'ribose_2oh');
      assert.deepEqual((await selection(page)).contexts, ['U']);
      assert.equal(await page.locator('#contextGroup button[data-value="U"]').getAttribute('aria-pressed'), 'true');
    }
    report.checks.push(`Stale context ignored; ${repair} restores current data${repair === 'trace' ? '; reused current-family context remains active' : ''}`);
    await page.close();
  }
  for (const repair of ['reset', 'all']) {
    const page = await failedSwitch('backbone', 'base_pair');
    const before = await selection(page);
    await page.evaluate(() => { window.detachedPuckerControl = document.querySelector('#puckerGroup'); });
    await page.selectOption('#puckerGroup', "C1'-endo");
    assert.deepEqual(await selection(page), before, 'Old residue pucker changed the pair-family selection');
    assert.equal(await page.locator('#puckerGroup').inputValue(), 'all');
    if (repair === 'reset') await page.click('#resetFilters');
    else await page.selectOption('#puckerGroup', 'all');
    await recovered(page, 'base_pair');
    assert.deepEqual((await selection(page)).puckers, []);
    assert.equal(await page.locator('#puckerGroup').count(), 0);
    const beforeDetached = await selection(page);
    assert(await page.evaluate(() => {
      const input = window.detachedPuckerControl;
      if (input.isConnected) return false;
      for (const value of ["C1'-endo", 'all']) { input.value = value; input.dispatchEvent(new Event('change')); }
      return !document.querySelector('#puckerGroup');
    }), 'Detached pucker callback recreated an inapplicable control');
    assert.deepEqual(await selection(page), beforeDetached, 'Detached pucker callback mutated selection or revision');
    await page.selectOption('#familySelect', 'backbone'); await waitReady(page);
    await page.evaluate(() => { window.savedPuckerControl = document.querySelector('#puckerGroup'); });
    await page.click('#traceStyleGroup button[data-value="line"]'); await waitReady(page);
    assert(await page.evaluate(() => window.savedPuckerControl === document.querySelector('#puckerGroup')));
    await page.selectOption('#puckerGroup', "C3'-endo"); await recovered(page, 'backbone');
    assert.deepEqual((await selection(page)).puckers, ["C3'-endo"]);
    report.checks.push(`Stale pucker ignored; ${repair} restores pair data; reused current-family pucker remains active`);
    await page.close();
  }
  assert.deepEqual(report.pageErrors, []); report.passed = true;
} catch (error) { report.passed = false; report.failure = error.stack; throw error; }
finally {
  report.finishedAt = new Date().toISOString();
  if (process.env.RNA_BROWSER_OUTPUT) {
    await mkdir(process.env.RNA_BROWSER_OUTPUT, { recursive: true });
    await writeFile(`${process.env.RNA_BROWSER_OUTPUT}/report.json`, JSON.stringify(report, null, 2) + '\n');
  }
  console.log(JSON.stringify(report, null, 2)); await browser.close();
}
