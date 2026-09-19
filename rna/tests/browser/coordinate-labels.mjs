import assert from 'node:assert/strict';
import { mkdir, readFile, writeFile } from 'node:fs/promises';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { waitReady } from './helpers.mjs';

const workspace = process.env.RNA_WORKSPACE || '/home/zhaomt/cmap/test15';
const output = process.env.RNA_BROWSER_OUTPUT || path.join(workspace, 'data/pure_rna/browser_validation/coordinate-labels');
await mkdir(output, { recursive: true });
const { chromium } = await import(pathToFileURL(process.env.PLAYWRIGHT_MODULE || path.join(workspace, 'data/pure_rna/browser-tools/node_modules/playwright-core/index.mjs')).href);
const browser = await chromium.launch({ headless: true });
const report = { startedAt: new Date().toISOString(), checks: [], errors: [], failedRequests: [], exports: [] };
const context = await browser.newContext({ viewport: { width: 1440, height: 1000 }, acceptDownloads: true });
const page = await context.newPage();
let requests = 0;
page.on('request', request => { if (request.url().includes('/survey/coordinates/')) requests++; });
page.on('pageerror', error => report.errors.push(error.message));
page.on('requestfailed', request => report.failedRequests.push(request.url()));
const check = description => { report.checks.push(description); console.log(`PASS ${description}`); };
async function waitMode(mode) {
  await page.waitForFunction(mode => document.querySelector('#coordinatePlot').data?.[0]?.mode === mode && window.rnaExplorer.completedCoordinateLabels === (mode === 'markers' ? 'none' : 'all'), mode);
}
async function exportMarkers(name) {
  const plot = page.locator('#coordinatePlot'); await plot.scrollIntoViewIfNeeded();
  const screenshot = await plot.screenshot(), box = await plot.boundingBox();
  const markerPixels = await page.evaluate(async encoded => {
    const image = new Image(); image.src = `data:image/png;base64,${encoded}`; await image.decode();
    const canvas = document.createElement('canvas'); canvas.width = image.width; canvas.height = image.height;
    const context = canvas.getContext('2d'); context.drawImage(image, 0, 0);
    const pixels = context.getImageData(0, 0, image.width, image.height).data, candidates = [];
    for (let y = 0; y < image.height; y++) for (let x = 0; x < image.width; x++) {
      const index = 4 * (y * image.width + x);
      if (Math.abs(pixels[index] - 23) < 4 && Math.abs(pixels[index + 1] - 74) < 4 && Math.abs(pixels[index + 2] - 126) < 4
          && candidates.every(point => Math.hypot(point.x - x, point.y - y) > 12)) candidates.push({ x, y });
    }
    return candidates.slice(0, 40);
  }, screenshot.toString('base64'));
  let hover = '';
  for (const point of markerPixels) {
    await page.mouse.move(box.x + point.x, box.y + point.y);
    try { await page.waitForFunction(() => /Context:.*Atom:/.test(document.querySelector('#coordinatePlot')?.textContent ?? ''), null, { timeout: 350 }); } catch {}
    hover = await plot.textContent();
    if (hover.includes('Context:') && hover.includes('RMS spread:')) break;
  }
  report.hoverDebug = { name, markerPixels, box, hover };
  assert.match(hover, /Context:.*Atom:.*Observations:.*Residues:.*Pairs:.*PDB entries:.*RMS spread:/);
  report.checks.push(`${name} actual marker hover includes identity and counts`);
  await plot.hover({ position: { x: 20, y: 20 } });
  const pending = page.waitForEvent('download');
  await plot.locator('.modebar-btn[data-title*="Download plot"]').click();
  const download = await pending; assert.equal(await download.failure(), null);
  const file = path.join(output, `${name}.svg`); await download.saveAs(file);
  const svg = await readFile(file, 'utf8'); assert.match(svg, /data:image\/png;base64,/);
  const standalone = await page.context().newPage();
  await standalone.setContent('<body style="margin:0;background:white"><img id="export"></body>');
  await standalone.evaluate(async svg => { const image = document.querySelector('img'); image.src = `data:image/svg+xml;base64,${btoa(unescape(encodeURIComponent(svg)))}`; await image.decode(); }, svg);
  await standalone.locator('img').screenshot({ path: path.join(output, `${name}.png`) });
  await standalone.close(); report.exports.push(file); check(`${name} standalone SVG export`);
}
try {
  await page.goto(process.env.RNA_URL || 'http://127.0.0.1:8767/nucleic.pages/rna/', { waitUntil: 'domcontentloaded' }); await waitReady(page);
  await page.selectOption('#family2Select', 'backbone'); await waitReady(page);
  await page.click('#baseGeometryLoad'); await waitReady(page);
  await page.click('#coordinatesLoad'); await waitReady(page);
  report.release = await page.evaluate(() => window.rnaExplorer.manifest.build_id);
  assert(await page.evaluate(() => window.rnaExplorer.coordinateSummary.length > 0));
  const before = await page.evaluate(() => {
    const app = window.rnaExplorer, plot = document.querySelector('#coordinatePlot');
    window.coordinateLabelSnapshotRefs = { ...app.snapshots };
    return { summary: JSON.stringify(app.coordinateSummary), table: document.querySelector('#baseGeometryCoordBody').textContent,
      xyz: JSON.stringify([plot.data[0].x, plot.data[0].y, plot.data[0].z]), revision: app.revision };
  });
  const initialRequests = requests;
  await page.evaluate(() => Plotly.relayout('coordinatePlot', { 'scene.camera.eye': { x: 1.4, y: 1.8, z: 2.5 } }));
  for (const labels of ['none', 'all', 'none']) {
    await page.selectOption('#coordinateLabelsSelect', labels); await waitMode(labels === 'none' ? 'markers' : 'markers+text');
    const after = await page.evaluate(() => {
      const app = window.rnaExplorer, plot = document.querySelector('#coordinatePlot');
      return { summary: JSON.stringify(app.coordinateSummary), table: document.querySelector('#baseGeometryCoordBody').textContent,
        xyz: JSON.stringify([plot.data[0].x, plot.data[0].y, plot.data[0].z]), revision: app.revision,
        sameSnapshots: Object.entries(window.coordinateLabelSnapshotRefs).every(([key, value]) => app.snapshots[key] === value),
        camera: plot._fullLayout.scene.camera.eye };
    });
    assert.deepEqual({ summary: after.summary, table: after.table, xyz: after.xyz, revision: after.revision }, before);
    assert(after.sameSnapshots); assert.equal(requests, initialRequests);
    assert.deepEqual(after.camera, { x: 1.4, y: 1.8, z: 2.5 });
  }
  check('Labels toggle without changing coordinates, counts, camera, snapshots, revision, or partition requests');
  const identity = await page.evaluate(() => {
    const app = window.rnaExplorer, trace = document.querySelector('#coordinatePlot').data[0];
    return { summary: app.coordinateSummary[0], customdata: trace.customdata[0], hover: trace.hovertemplate,
      cells: [...document.querySelector('#baseGeometryCoordBody tr').cells].map(cell => cell.textContent) };
  });
  assert.equal(identity.cells.length, 10); assert.equal(identity.cells[0], identity.summary.context);
  assert.equal(identity.cells[1], identity.summary.atom_label);
  assert.deepEqual(identity.customdata, [identity.summary.context, identity.summary.atom_label, identity.summary.n, identity.summary.residues ?? 'Unavailable', identity.summary.pairs ?? 'Not applicable', identity.summary.entries, identity.summary.rms]);
  assert.match(identity.hover, /RMS spread/); assert.match(identity.hover, /PDB entries/);
  check('Table and hover expose context, exact atom identity, population counts, and RMS spread');
  for (const mobile of [false, true]) {
    await page.setViewportSize({ width: mobile ? 390 : 1440, height: 1000 });
    await page.waitForFunction(() => Math.abs(document.querySelector('#coordinatePlot')._fullLayout.width - document.querySelector('#coordinatePlot').clientWidth) < 4);
    await exportMarkers(`pair-markers-${mobile ? 'mobile' : 'desktop'}`);
  }
  // Pause an actual pending Plotly update to exercise a second label choice
  // during the coordinate commit, without injecting new scientific data.
  await page.evaluate(() => {
    const app = window.rnaExplorer, original = app.plot.bind(app);
    window.coordinatePlotHeld = false;
    app.plot = async (...args) => {
      if (args[0].id === 'coordinatePlot' && !window.coordinatePlotHeld) {
        window.coordinatePlotHeld = true;
        await new Promise(resolve => { window.releaseCoordinatePlot = resolve; });
      }
      return original(...args);
    };
  });
  await page.selectOption('#coordinateLabelsSelect', 'all');
  await page.waitForFunction(() => window.coordinatePlotHeld);
  await page.selectOption('#coordinateLabelsSelect', 'none');
  await page.evaluate(() => window.releaseCoordinatePlot()); await waitMode('markers');
  assert.equal(requests, initialRequests);
  check('Latest label choice wins during an awaited Plotly commit');
  await page.evaluate(() => {
    const app = window.rnaExplorer, original = app.renderJoint.bind(app);
    window.jointLabelHeld = false;
    app.renderJoint = async (...args) => {
      if (!window.jointLabelHeld) {
        window.jointLabelHeld = true;
        await new Promise(resolve => { window.releaseJointLabel = resolve; });
      }
      return original(...args);
    };
  });
  await page.click('#jointPaletteGroup [data-value="viridis"]');
  await page.waitForFunction(() => window.jointLabelHeld);
  await page.selectOption('#coordinateLabelsSelect', 'all');
  await page.evaluate(() => window.releaseJointLabel()); await waitReady(page); await waitMode('markers+text');
  assert.equal(requests, initialRequests);
  check('Label choice made during a joint-only update is rendered before ready');
  await page.evaluate(() => {
    window.rnaExplorer.commitQueue = new Promise(resolve => { window.releaseLabelQueue = resolve; });
  });
  await page.selectOption('#coordinateLabelsSelect', 'none');
  await page.click('#jointPaletteGroup [data-value="hotspots"]');
  await page.evaluate(() => window.releaseLabelQueue()); await waitReady(page); await waitMode('markers');
  assert.equal(requests, initialRequests);
  assert(await page.evaluate(() => ['distribution', 'survey'].every(key => window.rnaExplorer.snapshots[key] === window.coordinateLabelSnapshotRefs[key])));
  check('A newer joint revision drains queued labels while preserving main and Survey snapshots');
  const baseGroup = await page.locator('#coordinateGroupSelect option').evaluateAll(options => options.find(option => option.value.includes('rna_standard_base')).value);
  await page.selectOption('#coordinateGroupSelect', baseGroup); await waitReady(page); await waitMode('markers');
  assert(requests > initialRequests);
  assert(await page.evaluate(() => window.rnaExplorer.coordinateSummary.every(row => row.pairs === null)));
  for (const mobile of [false, true]) {
    await page.setViewportSize({ width: mobile ? 390 : 1440, height: 1000 });
    await page.waitForFunction(() => Math.abs(document.querySelector('#coordinatePlot')._fullLayout.width - document.querySelector('#coordinatePlot').clientWidth) < 4);
    await exportMarkers(`base-markers-${mobile ? 'mobile' : 'desktop'}`);
  }
  await page.click('#resetFilters'); await waitReady(page); await waitMode('markers+text');
  assert.equal(await page.locator('#coordinateLabelsSelect').inputValue(), 'all');
  check('Reset restores context and atom labels');
  assert.deepEqual(report.errors, []); assert.deepEqual(report.failedRequests, []); report.passed = true;
} catch (error) { report.error = error.stack; throw error; }
finally { report.finishedAt = new Date().toISOString(); await writeFile(path.join(output, 'report.json'), JSON.stringify(report, null, 2) + '\n'); await browser.close(); }
