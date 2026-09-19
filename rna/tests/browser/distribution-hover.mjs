/** Real distribution, miniature and Survey hover labels preserve numerical data. */
import assert from 'node:assert/strict';
import { mkdir, writeFile } from 'node:fs/promises';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { waitReady } from './helpers.mjs';

const workspace = path.resolve(process.env.RNA_WORKSPACE || process.cwd());
const output = path.resolve(process.env.RNA_BROWSER_OUTPUT || path.join(workspace, 'data/pure_rna/browser_validation/distribution-hover'));
await mkdir(output, { recursive: true });
const { chromium } = await import(pathToFileURL(process.env.PLAYWRIGHT_MODULE || path.join(workspace, 'data/pure_rna/browser-tools/node_modules/playwright-core/index.mjs')).href);
const browser = await chromium.launch({ headless: true });
const page = await browser.newPage({ viewport: { width: 1440, height: 1000 } });
const report = { startedAt: new Date().toISOString(), checks: [], pageErrors: [] };
page.on('pageerror', error => report.pageErrors.push(error.message));
async function probe(kind, selector) {
  const evidence = await page.evaluate(async ({ kind, selector }) => {
    const app = window.rnaExplorer, result = app.snapshots[kind].result, plot = document.querySelector(selector);
    const trace = plot.data[0];
    let peak = 0;
    for (let index = 1; index < trace.y.length; index++) if (trace.y[index] > trace.y[peak]) peak = index;
    const { csv } = await import('./core/export.js');
    const before = csv(app.snapshots[kind]);
    Plotly.Fx.hover(plot, [{ curveNumber: 0, pointNumber: peak }]);
    await new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(resolve)));
    const hoverText = plot.querySelector('.hoverlayer')?.textContent;
    Plotly.Fx.unhover(plot);
    const same = plot.data.every((item, index) => item.x.length === result.series[index].x.length
      && item.x.every((value, bin) => value === result.series[index].x[bin])
      && item.y.every((value, bin) => value === result.series[index].y[bin]));
    const canonicalCorrect = !result.parameter.period || plot.data.every(item => item.x.every((value, index) => item.customdata?.[index]?.[0] === value
      && item.customdata[index][1] === ((value % result.parameter.period) + result.parameter.period) % result.parameter.period));
    const miniPlots = kind === 'distribution' ? [...document.querySelectorAll('.rna-mini-plot')].filter(node => node.data?.[0]) : [];
    const miniTemplates = miniPlots.map(node => node.data[0].hovertemplate);
    let miniHoverText = null;
    if (miniPlots.length) {
      const mini = miniPlots[0], y = mini.data[0].y, peak = y.indexOf(Math.max(...y));
      Plotly.Fx.hover(mini, [{ curveNumber: 0, pointNumber: peak }]);
      await new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(resolve)));
      miniHoverText = mini.querySelector('.hoverlayer')?.textContent;
      Plotly.Fx.unhover(mini);
    }
    return { kind, build: app.manifest.build_id, parameter: result.parameter, normalization: result.displaySpec.normalization, circularMode: result.displaySpec.circularMode,
      hoverText, template: trace.hovertemplate, rawDataUnchanged: same, csvUnchanged: before === csv(app.snapshots[kind]), canonicalCorrect, miniTemplates,
      miniHoverText, peak: trace.y[peak], view: trace.x[peak], canonical: trace.customdata?.[peak]?.[1], points: result.coverage.plottedRows };
  }, { kind, selector });
  assert(evidence.rawDataUnchanged && evidence.csvUnchanged && evidence.canonicalCorrect);
  assert(evidence.points > 0);
  const intensity = evidence.normalization === 'density' ? 'Probability density (smoothed)' : 'Probability (smoothed)';
  assert(evidence.hoverText?.includes(intensity), `Rendered hover missing ${intensity}: ${evidence.hoverText}`);
  assert(evidence.hoverText.includes(evidence.parameter.label));
  if (evidence.parameter.unit) assert(evidence.hoverText.includes(`(${evidence.parameter.unit})`));
  if (evidence.parameter.period) {
    assert(evidence.hoverText.includes(`View ${evidence.view.toFixed(3)}`));
    assert(evidence.hoverText.includes(`Angle ${evidence.canonical.toFixed(3)}`));
  } else assert(!evidence.template.includes('Angle'));
  assert(evidence.miniTemplates.every(template => template.includes(intensity)), 'Miniature normalization differs from result');
  if (evidence.miniTemplates.length) assert(evidence.miniHoverText?.includes(intensity), 'Rendered miniature hover lost normalization');
  report.checks.push(evidence);
}
try {
  await page.goto(process.env.RNA_URL || 'http://127.0.0.1:8767/nucleic.pages/rna/', { waitUntil: 'domcontentloaded', timeout: 120000 });
  await waitReady(page);
  await page.evaluate(() => {
    const result = window.rnaExplorer.snapshots.distribution.result;
    window.rnaHoverOriginalObservations = JSON.stringify(result.series.map(series => ({ ids: series.rowIds, values: series.values })));
  });
  for (const normalization of ['probability', 'density']) {
    await page.click(`#displayScaleGroup button[data-value="${normalization}"]`); await waitReady(page);
    for (const mode of ['signed_180', 'wrap_360', 'auto']) {
      await page.click(`#circularModeGroup button[data-value="${mode}"]`); await waitReady(page);
      await probe('distribution', '#plot');
      assert(await page.evaluate(() => JSON.stringify(window.rnaExplorer.snapshots.distribution.result.series.map(series => ({ ids: series.rowIds, values: series.values }))) === window.rnaHoverOriginalObservations), 'Display control changed source identities or measurements');
    }
  }
  await page.selectOption('#parameterSelect', 'e_z'); await waitReady(page);
  await probe('distribution', '#plot');
  await page.selectOption('#familySelect', 'ribose_2oh'); await waitReady(page);
  await page.selectOption('#parameterSelect', 'c2_o2_length'); await waitReady(page);
  await probe('distribution', '#plot');
  await page.click('#baseGeometryLoad'); await waitReady(page);
  await probe('survey', '#baseGeometryPlot');
  const circularTerm = await page.evaluate(() => window.rnaExplorer.surveyTerms().find(term => term.period)?.id);
  assert(circularTerm, 'Release has no circular Survey term');
  await page.selectOption('#baseGeometryTermSelect', circularTerm); await waitReady(page);
  await page.click('#circularModeGroup button[data-value="signed_180"]'); await waitReady(page);
  await probe('survey', '#baseGeometryPlot');
  assert.deepEqual(report.pageErrors, []);
  report.passed = true;
  console.log(`PASS ${report.checks.length} real distribution and Survey hover cases`);
} catch (error) { report.error = error.stack; throw error; }
finally { report.finishedAt = new Date().toISOString(); await writeFile(path.join(output, 'report.json'), JSON.stringify(report, null, 2) + '\n'); await browser.close(); }
