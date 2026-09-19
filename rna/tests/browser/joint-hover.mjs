/** Actual control clicks and hover text preserve raw 2D probabilities in log view. */
import assert from 'node:assert/strict';
import { mkdir, writeFile } from 'node:fs/promises';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { waitReady } from './helpers.mjs';

const workspace = path.resolve(process.env.RNA_WORKSPACE || process.cwd());
const output = path.resolve(process.env.RNA_BROWSER_OUTPUT || path.join(workspace, 'data/pure_rna/browser_validation/joint-hover'));
await mkdir(output, { recursive: true });
const { chromium } = await import(pathToFileURL(process.env.PLAYWRIGHT_MODULE || path.join(workspace, 'data/pure_rna/browser-tools/node_modules/playwright/index.mjs')).href);
const browser = await chromium.launch({ headless: true });
const page = await browser.newPage({ viewport: { width: 1440, height: 1000 } });
const report = { startedAt: new Date().toISOString(), checks: [], pageErrors: [], failedRequests: [] };
page.on('pageerror', error => report.pageErrors.push(error.message));
page.on('requestfailed', request => report.failedRequests.push({ url: request.url(), error: request.failure()?.errorText }));
try {
  await page.goto(process.env.RNA_URL || 'http://127.0.0.1:8767/nucleic.pages/rna/', { waitUntil: 'domcontentloaded', timeout: 120000 });
  await waitReady(page);
  await page.selectOption('#family2Select', 'backbone'); await waitReady(page);
  await page.selectOption('#parameter2Select', 'delta'); await waitReady(page);
  await page.click('#circularModeGroup button[data-value="signed_180"]'); await waitReady(page);
  await page.evaluate(() => { window.rnaHoverBaseline = window.rnaExplorer.snapshots.joint.result.points.map(point => [point.left_id, point.right_id, point.x, point.y]); });
  for (const normalization of ['probability', 'density']) {
    await page.click(`#displayScaleGroup button[data-value="${normalization}"]`); await waitReady(page);
    for (const type of ['heatmap', 'contour', 'filled_contour', 'heatmap_contour']) {
      await page.click(`#jointPlotTypeGroup button[data-value="${type}"]`); await waitReady(page);
      for (const scale of ['linear', 'log']) {
        await page.click(`#jointColorScaleGroup button[data-value="${scale}"]`); await waitReady(page);
        const evidence = await page.evaluate(async () => {
          const app = window.rnaExplorer, result = app.snapshots.joint.result, plot = document.querySelector('#jointPlot');
          const pointsUnchanged = result.points.length === window.rnaHoverBaseline.length
            && result.points.every((point, index) => [point.left_id, point.right_id, point.x, point.y].every((value, field) => value === window.rnaHoverBaseline[index][field]));
          let maximum = { value: -1, x: 0, y: 0 }, allValuesEqual = true, zero = null, zeroBins = 0, belowFloorBins = 0;
          for (let y = 0; y < result.y.length; y++) for (let x = 0; x < result.x.length; x++) {
            const value = result.z[y][x];
            if (value > maximum.value) maximum = { value, x, y };
            if (value === 0) { zeroBins++; zero ??= { x, y }; }
            else if (value < 1e-8) belowFloorBins++;
            for (const trace of plot.data) {
              const hover = trace.customdata?.[y]?.[x];
              if (!hover || hover[0] !== result.x[x] || hover[2] !== result.y[y] || hover[4] !== value) allValuesEqual = false;
              const displayed = app.state.joint.colorScale === 'log' ? Math.log10(Math.max(value, 1e-8)) : value;
              if (trace.z[y][x] !== displayed) allValuesEqual = false;
              if (result.xParameter.period && hover?.[1] !== ((result.x[x] % result.xParameter.period) + result.xParameter.period) % result.xParameter.period) allValuesEqual = false;
              if (result.yParameter.period && hover?.[3] !== ((result.y[y] % result.yParameter.period) + result.yParameter.period) % result.yParameter.period) allValuesEqual = false;
            }
          }
          const log = app.state.joint.colorScale === 'log';
          const colorRangeCorrect = plot._fullData.every(trace => trace.zmin === (log ? -8 : 0)
            && trace.zmax === (log ? Math.log10(Math.max(maximum.value, 1e-8)) : maximum.value));
          let hoverText = null, zeroHoverText = null;
          if (app.state.joint.type === 'heatmap') {
            Plotly.Fx.hover(plot, [{ curveNumber: 0, xval: result.x[maximum.x], yval: result.y[maximum.y] }]);
            await new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(resolve)));
            hoverText = plot.querySelector('.hoverlayer')?.textContent;
            Plotly.Fx.unhover(plot);
            if (zero) {
              Plotly.Fx.hover(plot, [{ curveNumber: 0, xval: result.x[zero.x], yval: result.y[zero.y] }]);
              await new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(resolve)));
              zeroHoverText = plot.querySelector('.hoverlayer')?.textContent;
              Plotly.Fx.unhover(plot);
            }
          }
          return { buildId: app.manifest.build_id, normalization: app.state.display.normalization, type: app.state.joint.type,
            scale: app.state.joint.colorScale, points: result.points.length, pointsUnchanged, allValuesEqual,
            hoverTemplates: plot.data.map(trace => trace.hovertemplate), hoverText, zeroHoverText, zeroBins, belowFloorBins, colorRangeCorrect,
            peakIntensity: maximum.value, peakFormatted: maximum.value.toPrecision(4) };
        });
        assert(evidence.points > 0);
        assert(evidence.pointsUnchanged);
        assert(evidence.allValuesEqual, 'Hover fields disagree with untransformed histogram');
        assert(evidence.colorRangeCorrect, 'Rendered color range differs from the DNA display floor');
        assert(evidence.zeroBins > 0, 'Real histogram did not exercise zero bins');
        for (const template of evidence.hoverTemplates) {
          assert.match(template, /chi \(deg\)/);
          assert.match(template, /delta \(deg\)/);
          assert.match(template, /customdata\[4\]/);
        }
        if (type === 'heatmap') {
          assert.match(evidence.hoverText, normalization === 'density' ? /Probability density \(smoothed\)/ : /Probability \(smoothed\)/);
          assert(evidence.hoverText.includes(evidence.peakFormatted), `Rendered hover lost raw peak ${evidence.peakFormatted}: ${evidence.hoverText}`);
          assert.match(evidence.zeroHoverText, /Probability(?: density)? \(smoothed\): 0(?:\.0+)?$/,
            'Zero-bin hover must report actual zero, not the color floor');
        }
        report.checks.push({ passed: true, ...evidence });
      }
    }
  }
  assert.deepEqual(report.pageErrors, []);
  assert.deepEqual(report.failedRequests, []);
  report.passed = true;
  console.log(`PASS ${report.checks.length} raw joint hover control combinations`);
} catch (error) {
  report.passed = false; report.failure = { message: error.message, stack: error.stack }; throw error;
} finally {
  report.finishedAt = new Date().toISOString();
  await writeFile(path.join(output, 'rna-joint-hover.json'), JSON.stringify(report, null, 2) + '\n');
  await browser.close();
}
