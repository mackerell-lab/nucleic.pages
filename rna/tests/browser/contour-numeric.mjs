/** Real controls and Plotly-coerced numeric levels match independent DNA code. */
import assert from 'node:assert/strict';
import { mkdir, readFile, writeFile } from 'node:fs/promises';
import path from 'node:path';
import vm from 'node:vm';
import { pathToFileURL } from 'node:url';
import { downloadCsv, waitReady } from './helpers.mjs';

const workspace = path.resolve(process.env.RNA_WORKSPACE || process.cwd());
const output = path.resolve(process.env.RNA_BROWSER_OUTPUT || path.join(workspace, 'data/pure_rna/browser_validation/contour-numeric'));
await mkdir(output, { recursive: true });
const dna = await readFile(new URL('../../../js/pure-dna.js', import.meta.url), 'utf8');
const oracle = vm.createContext({ state: {} });
const dnaOptions = dna.match(/const JOINT_CONTOUR_WIDTH_OPTIONS = (\[[\s\S]*?\]);/)[1];
vm.runInContext(`const JOINT_CONTOUR_WIDTH_OPTIONS = ${dnaOptions};\n${dna.slice(dna.indexOf('function currentJointContourTargetLevels('), dna.indexOf('function nextNiceStepAtLeast('))}`, oracle);
vm.runInContext(dna.slice(dna.indexOf('function nextNiceStepAtLeast('), dna.indexOf('function buildJointPlotTraces(')), oracle);
const { chromium } = await import(pathToFileURL(process.env.PLAYWRIGHT_MODULE || path.join(workspace, 'data/pure_rna/browser-tools/node_modules/playwright-core/index.mjs')).href);
const browser = await chromium.launch({ headless: true });
const page = await browser.newPage({ viewport: { width: 1440, height: 1000 } });
const report = { startedAt: new Date().toISOString(), checks: [], errors: [] };
page.on('pageerror', error => report.errors.push(error.message));
page.on('console', message => { if (message.type() === 'error') report.errors.push(message.text()); });
const rawCsv = rows => rows.map(({ snapshot_id, ...row }) => row);
async function click(group, value) { await page.click(`#${group} button[data-value="${value}"]`); await waitReady(page); }
try {
  await page.goto(process.env.RNA_URL || 'http://127.0.0.1:8767/nucleic.pages/rna/', { waitUntil: 'domcontentloaded', timeout: 120000 }); await waitReady(page);
  await page.selectOption('#family2Select', 'backbone'); await waitReady(page);
  await page.selectOption('#parameter2Select', 'delta'); await waitReady(page);
  report.buildId = await page.evaluate(() => window.rnaExplorer.manifest.build_id);
  const baselineCsv = await downloadCsv(page, '#jointCsvDownload', path.join(output, 'before.csv'));
  await page.evaluate(() => {
    const app = window.rnaExplorer, plot = document.getElementById('jointPlot');
    window.numericContourBaseline = { distribution: app.snapshots.distribution, z: JSON.stringify(app.snapshots.joint.result.z), customdata: JSON.stringify(plot.data[0].customdata), points: JSON.stringify(app.snapshots.joint.result.points.map(point => [point.left_id, point.right_id, point.x, point.y])) };
  });
  for (const type of ['contour', 'filled_contour', 'heatmap_contour']) {
    await click('jointPlotTypeGroup', type);
    await click('jointContourLabelsGroup', 'on');
    for (const scale of ['linear', 'log']) {
      await click('jointColorScaleGroup', scale);
      const sizes = [];
      for (const [count, preset] of [[6, 'wide'], [12, 'standard'], [24, 'tight']]) {
        await click('jointContourWidthGroup', count);
        const evidence = await page.evaluate(() => {
          const app = window.rnaExplorer, plot = document.getElementById('jointPlot');
          const trace = plot.data.find(trace => trace.type === 'contour'), full = plot._fullData.find(trace => trace.type === 'contour');
          const baseline = window.numericContourBaseline;
          return { type: app.state.joint.type, scale: app.state.joint.colorScale, count: app.state.joint.contourCount,
            zmin: trace.zmin, zmax: trace.zmax, automatic: trace.autocontour, levels: trace.contours,
            fullAutomatic: full.autocontour, fullColoring: full.contours.coloring, fullLevels: { start: full.contours.start, end: full.contours.end, size: full.contours.size, showlabels: full.contours.showlabels },
            distributionUnchanged: app.snapshots.distribution === baseline.distribution,
            rawGridUnchanged: JSON.stringify(app.snapshots.joint.result.z) === baseline.z,
            rawPointsUnchanged: JSON.stringify(app.snapshots.joint.result.points.map(point => [point.left_id, point.right_id, point.x, point.y])) === baseline.points,
            hoverUnchanged: JSON.stringify(trace.customdata) === baseline.customdata,
            points: app.snapshots.joint.result.points.length,
            displayGridCorrect: trace.z.every((row, y) => row.every((value, x) => value === (app.state.joint.colorScale === 'log' ? Math.log10(Math.max(app.snapshots.joint.result.z[y][x], 1e-8)) : app.snapshots.joint.result.z[y][x]))) };
        });
        oracle.state = { jointContourWidth: preset, jointContourLabels: 'on' };
        const expected = JSON.parse(JSON.stringify(oracle.buildJointContourConfig(evidence.zmin, evidence.zmax)));
        assert.equal(evidence.automatic, false); assert.equal(evidence.fullAutomatic, false);
        assert.equal(evidence.levels.coloring, type === 'filled_contour' ? 'heatmap' : 'none');
        assert.equal(evidence.fullColoring, type === 'filled_contour' ? 'heatmap' : 'none');
        for (const key of ['start', 'end', 'size', 'showlabels']) assert.equal(evidence.levels[key], expected.contours[key], `${type}/${scale}/${preset}: input ${key}`);
        assert.deepEqual(evidence.fullLevels, expected.contours);
        for (const key of ['distributionUnchanged', 'rawGridUnchanged', 'rawPointsUnchanged', 'hoverUnchanged', 'displayGridCorrect']) assert(evidence[key], key);
        assert(evidence.points > 0); if (scale === 'log') assert.equal(evidence.zmin, -8);
        sizes.push(evidence.fullLevels.size); report.checks.push({ name: `${type}/${scale}/${preset}`, passed: true, ...evidence });
      }
      assert(sizes[0] > sizes[1] && sizes[1] > sizes[2]);
    }
  }
  await click('jointPlotTypeGroup', 'filled_contour');
  await click('jointContourWidthGroup', 12);
  for (const scale of ['linear', 'log']) {
    await click('jointColorScaleGroup', scale);
    const plot = page.locator('#jointPlot'); await plot.scrollIntoViewIfNeeded(); await plot.hover({ position: { x: 20, y: 20 } });
    const pending = page.waitForEvent('download'); await plot.locator('.modebar-btn[data-title*="Download plot"]').click();
    const download = await pending; assert.equal(await download.failure(), null); assert.match(download.suggestedFilename(), /\.svg$/);
    const artifact = path.join(output, `filled-${scale}.svg`); await download.saveAs(artifact);
    const svg = await readFile(artifact, 'utf8');
    const inspection = await page.evaluate(async svg => {
      const doc = new DOMParser().parseFromString(svg, 'image/svg+xml');
      const images = [...doc.querySelectorAll('image')];
      const standalone = new Image(); standalone.src = `data:image/svg+xml;base64,${btoa(unescape(encodeURIComponent(svg)))}`; await standalone.decode();
      const canvas = document.createElement('canvas'); canvas.width = standalone.naturalWidth; canvas.height = standalone.naturalHeight;
      const context = canvas.getContext('2d'); context.drawImage(standalone, 0, 0);
      const pixels = context.getImageData(0, 0, canvas.width, canvas.height).data;
      let colored = 0; for (let i = 0; i < pixels.length; i += 4) if (Math.max(...pixels.slice(i, i + 3)) - Math.min(...pixels.slice(i, i + 3)) > 30) colored++;
      return { parseErrors: doc.querySelectorAll('parsererror').length, paths: doc.querySelectorAll('path').length, embeddedRaster: images.some(image => (image.getAttribute('href') ?? image.getAttributeNS('http://www.w3.org/1999/xlink', 'href') ?? '').startsWith('data:image/png')), coloredPixels: colored };
    }, svg);
    assert.equal(inspection.parseErrors, 0); assert(inspection.paths > 0); assert(inspection.embeddedRaster); assert(inspection.coloredPixels > 100);
    const preview = await browser.newPage(); await preview.setContent('<img id="exported" style="max-width:100%">');
    await preview.locator('#exported').evaluate(async (image, svg) => { image.src = `data:image/svg+xml;base64,${btoa(unescape(encodeURIComponent(svg)))}`; await image.decode(); }, svg);
    await preview.locator('#exported').screenshot({ path: path.join(output, `filled-${scale}.png`) }); await preview.close();
    report.checks.push({ name: `Actual ${scale} filled contour SVG export`, passed: true, artifact, ...inspection });
  }
  const afterCsv = await downloadCsv(page, '#jointCsvDownload', path.join(output, 'after.csv'));
  assert.deepEqual(rawCsv(afterCsv.rows), rawCsv(baselineCsv.rows));
  report.checks.push({ name: 'Actual CSV download preserves all raw observations', passed: true, rows: afterCsv.rows.length });
  const edges = await page.evaluate(async () => {
    const { jointContourConfig } = await import('./core/contours.js');
    const cases = [
      { name: 'empty', min: 0, max: undefined, z: [] },
      { name: 'flat', min: 0, max: 0, z: [[0, 0], [0, 0]] },
      { name: 'flat log floor', min: -8, max: -8, z: [[-8, -8], [-8, -8]] },
      { name: 'tiny positive intensity', min: 0, max: 1e-20, z: [[0, 1e-20], [1e-21, 0]] },
    ];
    const records = [];
    for (const item of cases) {
      const div = document.createElement('div'); document.body.append(div);
      const config = jointContourConfig(item.min, item.max, { contourCount: 24 });
      const explicitSize = config.contours.size;
      await window.Plotly.newPlot(div, [{ type: 'contour', z: item.z, zmin: item.min, zmax: item.max, ...config }], { width: 400, height: 300 });
      const full = div._fullData[0];
      records.push({ name: item.name, automatic: config.autocontour, fallbackLevels: config.ncontours, explicitSize, visible: full.visible, fullSize: full.contours?.size });
      window.Plotly.purge(div); div.remove();
    }
    return records;
  });
  for (const edge of edges) {
    if (edge.name === 'tiny positive intensity') { assert.equal(edge.automatic, false); assert.equal(edge.fullSize, edge.explicitSize); assert(edge.fullSize > 0); }
    else { assert.equal(edge.automatic, true); assert.equal(edge.fallbackLevels, 14); if (edge.fullSize !== undefined) assert(Number.isFinite(edge.fullSize) && edge.fullSize > 0); }
    report.checks.push({ ...edge, name: `Plotly edge: ${edge.name}`, passed: true });
  }
  assert.deepEqual(report.errors, []); report.passed = true;
  console.log(`PASS ${report.checks.length} numeric contour checks`);
} catch (error) { report.passed = false; report.failure = error.stack; throw error; }
finally { report.finishedAt = new Date().toISOString(); await writeFile(path.join(output, 'report.json'), JSON.stringify(report, null, 2) + '\n'); await browser.close(); }
