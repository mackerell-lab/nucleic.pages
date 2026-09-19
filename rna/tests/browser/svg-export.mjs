/** Exercise actual Plotly SVG download buttons and inspect standalone artifacts. */
import assert from 'node:assert/strict';
import { mkdir, readFile, writeFile } from 'node:fs/promises';
import { createHash } from 'node:crypto';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { waitReady } from './helpers.mjs';

const workspace = path.resolve(process.env.RNA_WORKSPACE || process.cwd());
const output = path.resolve(process.env.RNA_BROWSER_OUTPUT || path.join(workspace, 'data/pure_rna/browser_validation/svg-export'));
await mkdir(output, { recursive: true });
const { chromium } = await import(pathToFileURL(process.env.PLAYWRIGHT_MODULE || path.join(workspace, 'data/pure_rna/browser-tools/node_modules/playwright-core/index.mjs')).href);
const report = { startedAt: new Date().toISOString(), url: process.env.RNA_URL || 'http://127.0.0.1:8767/nucleic.pages/rna/', checks: [], pageErrors: [], failedRequests: [] };
const browser = await chromium.launch({ headless: true });
const context = await browser.newContext({ viewport: { width: 1440, height: 1000 }, acceptDownloads: true });
const page = await context.newPage();
page.on('pageerror', error => report.pageErrors.push(error.message));
page.on('requestfailed', request => report.failedRequests.push({ url: request.url(), error: request.failure()?.errorText }));
async function exported(name, selector, is3D = false) {
  const before = await page.evaluate(selector => {
    const plot = document.querySelector(selector);
    const axes = ['x', 'y', 'z'];
    const frame = plot.layout.scene ? axes.map(axis => ({ axis, range: plot.layout.scene[`${axis}axis`].range,
      ratio: plot.layout.scene.aspectratio?.[axis], values: plot.data[0][axis] })) : null;
    return { data: JSON.stringify(plot.data), frame,
      x: plot._fullLayout.xaxis?.title?.text, y: plot._fullLayout.yaxis?.title?.text,
      scene: plot.layout.scene, atoms: Array.isArray(plot.data[0]?.text) ? plot.data[0].text : [] };
  }, selector);
  if (is3D) {
    const scales = before.frame.map(({ range, ratio, values }) => {
      assert(range?.every(Number.isFinite) && range[1] > range[0]);
      assert(values.every(value => value > range[0] && value < range[1]), 'Clipped coordinate domain');
      return ratio / (range[1] - range[0]);
    });
    assert(scales.every(value => Number.isFinite(value) && value > 0));
    assert(scales.every(value => Math.abs(value - scales[0]) < 1e-12), 'Unequal physical axis scale');
  }
  const plot = page.locator(selector);
  await plot.scrollIntoViewIfNeeded();
  await plot.hover({ position: { x: 20, y: 20 } });
  const button = plot.locator('.modebar-btn[data-title*="Download plot"]');
  await button.waitFor({ state: 'visible' });
  const title = await button.getAttribute('data-title');
  const pending = page.waitForEvent('download', { timeout: 120000 });
  await button.click();
  const download = await pending;
  assert.equal(await download.failure(), null);
  assert.match(download.suggestedFilename(), /\.svg$/i);
  const artifact = path.join(output, `${name}.svg`);
  await download.saveAs(artifact);
  const svg = await readFile(artifact, 'utf8');
  const inspection = await page.evaluate(svg => {
    const document = new DOMParser().parseFromString(svg, 'image/svg+xml');
    const root = document.documentElement;
    return { root: root.localName, namespace: root.namespaceURI, parseErrors: document.querySelectorAll('parsererror').length,
      width: root.getAttribute('width'), height: root.getAttribute('height'),
      text: [...document.querySelectorAll('text')].map(node => node.textContent),
      paths: document.querySelectorAll('path').length,
      images: [...document.querySelectorAll('image')].map(node => ({ href: node.getAttribute('href') ?? node.getAttributeNS('http://www.w3.org/1999/xlink', 'href'), width: node.getAttribute('width'), height: node.getAttribute('height') })),
      externalReferences: [...document.querySelectorAll('[href], [*|href]')].map(node => node.getAttribute('href') ?? node.getAttributeNS('http://www.w3.org/1999/xlink', 'href')).filter(value => value && !value.startsWith('data:') && !value.startsWith('#')) };
  }, svg);
  assert.equal(inspection.root, 'svg'); assert.equal(inspection.namespace, 'http://www.w3.org/2000/svg'); assert.equal(inspection.parseErrors, 0);
  assert(Number(inspection.width) > 0 && Number(inspection.height) > 0);
  assert.deepEqual(inspection.externalReferences, [], 'Export relies on nonembedded external references');
  if (is3D) {
    assert(inspection.images.some(image => image.href?.startsWith('data:image/png;base64,')), 'Plotly WebGL export lacks embedded PNG');
  } else {
    assert(inspection.paths > 0, 'Plot has no vector paths');
    assert(inspection.text.some(text => text.includes(before.x)), `Missing x-axis title ${before.x}`);
    assert(inspection.text.some(text => text.includes(before.y)), `Missing y-axis title ${before.y}`);
    assert(inspection.text.some(text => /[0-9]/.test(text)), 'Missing numerical tick text');
  }
  const standalone = await context.newPage();
  await standalone.setContent('<html><body style="margin:0;background:white"><img id="exported"></body></html>');
  const pixels = await standalone.evaluate(async svg => {
    const image = document.querySelector('#exported');
    image.src = `data:image/svg+xml;base64,${btoa(unescape(encodeURIComponent(svg)))}`;
    await image.decode();
    const canvas = document.createElement('canvas'); canvas.width = image.naturalWidth; canvas.height = image.naturalHeight;
    const context = canvas.getContext('2d'); context.drawImage(image, 0, 0);
    const data = context.getImageData(0, 0, canvas.width, canvas.height).data;
    let nonwhite = 0;
    for (let index = 0; index < data.length; index += 4) if (data[index + 3] > 0 && Math.min(data[index], data[index + 1], data[index + 2]) < 230) nonwhite++;
    return { width: canvas.width, height: canvas.height, nonwhite };
  }, svg);
  assert(pixels.nonwhite > 100, 'Standalone exported plot is blank');
  await standalone.locator('#exported').screenshot({ path: path.join(output, `${name}.png`) });
  await standalone.close();
  const unchanged = await page.evaluate(({ selector, data }) => JSON.stringify(document.querySelector(selector).data) === data, { selector, data: before.data });
  assert(unchanged, 'SVG export mutated plotted data');
  const evidence = { name, buttonTitle: title, suggestedFilename: download.suggestedFilename(), artifact, bytes: Buffer.byteLength(svg), sha256: createHash('sha256').update(svg).digest('hex'),
    is3D, expectedLabels: is3D ? { scene: before.scene, atoms: before.atoms } : { x: before.x, y: before.y },
    physicalFrame: before.frame, ...inspection, images: inspection.images.map(image => ({ ...image, href: image.href?.slice(0, 30), embeddedBytes: image.href?.length })), pixels, plottedDataUnchanged: unchanged };
  report.checks.push(evidence);
  console.log(`PASS ${name} actual SVG download`);
}
try {
  await page.goto(report.url, { waitUntil: 'domcontentloaded', timeout: 120000 }); await waitReady(page);
  report.release = await page.evaluate(() => ({ buildId: window.rnaExplorer.manifest.build_id, partial: window.rnaExplorer.manifest.partial, releaseUrl: window.rnaExplorer.repository.releaseUrl }));
  assert.equal(report.release.partial, false);
  if (process.env.RNA_EXPECT_BUILD_ID) assert.equal(report.release.buildId, process.env.RNA_EXPECT_BUILD_ID);
  await exported('distribution', '#plot');
  await page.selectOption('#family2Select', 'backbone'); await waitReady(page);
  await page.selectOption('#parameter2Select', 'delta'); await waitReady(page);
  await exported('joint', '#jointPlot');
  await page.click('#baseGeometryLoad'); await waitReady(page);
  await exported('survey', '#baseGeometryPlot');
  await page.click('#coordinatesLoad'); await waitReady(page);
  assert(await page.evaluate(() => window.rnaExplorer.coordinateSummary.length > 0));
  const groups = await page.locator('#coordinateGroupSelect option').evaluateAll(options => options.map(option => option.value));
  for (const group of [groups.find(group => group.includes('cytosine_standard_pair')), groups.find(group => group.includes('rna_standard_base'))]) {
    assert(group, 'Release needs pair and single-base coordinate groups');
    await page.selectOption('#coordinateGroupSelect', group); await waitReady(page);
    for (const mobile of [false, true]) {
      await page.setViewportSize({ width: mobile ? 390 : 1440, height: 1000 });
      await page.waitForFunction(() => Math.abs(document.querySelector('#coordinatePlot')._fullLayout.width - document.querySelector('#coordinatePlot').clientWidth) < 4);
      await exported(`coordinates-${group}-${mobile ? 'mobile' : 'desktop'}`, '#coordinatePlot', true);
    }
  }
  assert.deepEqual(report.pageErrors, []); assert.deepEqual(report.failedRequests, []);
  report.passed = true;
} catch (error) { report.error = error.stack; throw error; }
finally { report.finishedAt = new Date().toISOString(); await writeFile(path.join(output, 'report.json'), JSON.stringify(report, null, 2) + '\n'); await browser.close(); }
