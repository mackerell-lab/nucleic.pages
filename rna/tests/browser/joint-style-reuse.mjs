/** Real controls must reuse completed joint analysis without retaining stale filters. */
import assert from 'node:assert/strict';
import { mkdir, writeFile } from 'node:fs/promises';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { downloadCsv, waitReady } from './helpers.mjs';
import { jointColorscale } from '../../views/palettes.js';

const workspace = path.resolve(process.env.RNA_WORKSPACE || process.cwd());
const output = process.env.RNA_BROWSER_OUTPUT || path.join(workspace, 'data/pure_rna/browser_validation/joint-style-reuse');
await mkdir(output, { recursive: true });
const { chromium } = await import(pathToFileURL(process.env.PLAYWRIGHT_MODULE || path.join(workspace, 'data/pure_rna/browser-tools/node_modules/playwright-core/index.mjs')).href);
const browser = await chromium.launch({ headless: true });
const page = await browser.newPage({ viewport: { width: 1440, height: 1000 }, acceptDownloads: true });
const report = { started: new Date().toISOString(), checks: [], errors: [] };
page.on('pageerror', error => report.errors.push(error.message));
page.on('console', message => { if (message.type() === 'error') report.errors.push(message.text()); });
page.on('requestfailed', request => report.errors.push(request.failure()?.errorText));
const click = async (group, value) => { await page.click(`#${group} button[data-value="${value}"]`); await waitReady(page); };
const select = async (id, value) => { await page.selectOption(`#${id}`, value); await waitReady(page); };
async function remember() {
  await page.evaluate(() => {
    const a = window.rnaExplorer;
    window.styleBefore = { joint: a.snapshots.joint, main: a.snapshots.distribution, survey: a.snapshots.survey, loads: window.styleLoads };
  });
}
async function evidence(name, reuse, unchangedMain = true) {
  const e = await page.evaluate(() => {
    const a = window.rnaExplorer, previous = window.styleBefore, s = a.snapshots.joint;
    const traces = document.querySelector('#jointPlot').data;
    return { reused: s.result === previous.joint.result, newSnapshot: s !== previous.joint, loads: window.styleLoads - previous.loads,
      sameMain: previous.main === a.snapshots.distribution, sameSurvey: previous.survey === a.snapshots.survey,
      styleMatches: JSON.stringify(s.join_spec) === JSON.stringify(a.state.joint),
      points: s.result.points.length, types: traces.map(t => t.type), colorScale: a.state.joint.colorScale,
      renderedColorscale: traces[0].colorscale,
      displayValuesExact: traces.every(trace => trace.z.every((row, y) => row.every((value, x) => value === (a.state.joint.colorScale === 'log' ? Math.log10(Math.max(s.result.z[y][x], 1e-8)) : s.result.z[y][x])))),
      hoverValuesExact: traces.every(trace => trace.customdata.every((row, y) => row.every((value, x) => value[4] === s.result.z[y][x]))),
      palette: a.state.joint.palette, snapshotId: s.snapshot_id, contours: traces.find(t => t.type === 'contour')?.contours,
      valuesExact: s.result.points.every(p => p.x === (p.left.values?.[s.result.xParameter.id] ?? p.left[s.result.xParameter.id]) && p.y === (p.right.values?.[s.result.yParameter.id] ?? p.right[s.result.yParameter.id])) };
  });
  assert.equal(e.reused, reuse, name); assert(e.newSnapshot, `${name}: export metadata refreshed`); assert(e.styleMatches, `${name}: current join metadata`);
  assert(e.valuesExact, `${name}: raw scalar values retained`); assert(e.points > 0, `${name}: finite observations`);
  assert(e.displayValuesExact, `${name}: displayed intensity`); assert(e.hoverValuesExact, `${name}: raw hover intensity`);
  assert.deepEqual(e.renderedColorscale, jointColorscale(e.palette), `${name}: rendered palette`);
  if (unchangedMain) { assert(e.sameMain, name + ': 1D snapshot stable'); assert(e.sameSurvey, name + ': Survey snapshot stable'); }
  if (reuse) assert.equal(e.loads, 0, name + ': no repository calls');
  report.checks.push({ name, ...e }); return e;
}
async function mutate(name, action, unchangedMain = true) {
  await remember(); await action(); return evidence(name, false, unchangedMain);
}
const withoutSnapshot = rows => rows.map(({ snapshot_id, ...row }) => row);
try {
  await page.goto(process.env.RNA_URL || 'http://127.0.0.1:8767/nucleic.pages/rna/'); await waitReady(page);
  await page.click('#baseGeometryLoad'); await waitReady(page);
  await select('family2Select', 'backbone'); await select('parameter2Select', 'alpha');
  report.build = await page.evaluate(() => window.rnaExplorer.manifest.build_id);
  assert.equal(report.build, 'full_packed_family_20260919');
  const raw = await page.evaluate(async () => {
    const a = window.rnaExplorer;
    const table = await a.repository.loadFamily('backbone');
    const rows = new Map(table.rows.map(row => [row.id, row]));
    const expected = a.snapshots.joint.result.points.map(p => {
      if (p.left_id !== p.right_id || !rows.has(p.left_id)) throw Error('Invented identity');
      return { id: p.left_id, x: rows.get(p.left_id).values.chi, y: rows.get(p.right_id).values.alpha };
    });
    window.styleLoads = 0;
    for (const name of ['loadFamily', 'loadRelations', 'loadMetadata', 'loadSurveyScalars']) {
      const original = a.repository[name].bind(a.repository);
      a.repository[name] = (...args) => { window.styleLoads++; return original(...args); };
    }
    return expected;
  });
  const baseline = await downloadCsv(page, '#jointCsvDownload', path.join(output, 'baseline.csv'));
  assert.equal(baseline.rows.length, raw.length);
  baseline.rows.forEach((row, i) => { assert.equal(row.x_id, raw[i].id); assert.equal(row.y_id, raw[i].id); assert.equal(Number(row.x_value), raw[i].x); assert.equal(Number(row.y_value), raw[i].y); });
  // Every palette is exercised through its actual button, including returning to Hotspots.
  for (const palette of ['warm', 'viridis', 'cividis', 'ocean', 'forest', 'greys', 'hotspots']) {
    await remember(); await click('jointPaletteGroup', palette); const e = await evidence(`Palette ${palette}`, true);
    const exported = await downloadCsv(page, '#jointCsvDownload', path.join(output, `palette-${palette}.csv`));
    assert.deepEqual(withoutSnapshot(exported.rows), withoutSnapshot(baseline.rows));
    assert(exported.rows.every(row => row.snapshot_id === e.snapshotId));
  }
  for (const type of ['contour', 'filled_contour', 'heatmap_contour', 'heatmap', 'contour']) {
    await remember(); await click('jointPlotTypeGroup', type); const e = await evidence(`Plot type ${type}`, true);
    assert.deepEqual(e.types, type === 'heatmap_contour' ? ['heatmap', 'contour'] : [type === 'heatmap' ? 'heatmap' : 'contour']);
  }
  for (const count of ['6', '24', '12']) {
    await remember(); await click('jointContourWidthGroup', count); const e = await evidence(`Contour spacing ${count}`, true);
    assert(e.contours.size > 0);
  }
  for (const value of ['on', 'off']) {
    await remember(); await click('jointContourLabelsGroup', value); const e = await evidence(`Contour labels ${value}`, true);
    assert.equal(e.contours.showlabels, value === 'on');
  }
  for (const value of ['log', 'linear']) {
    await remember(); await click('jointColorScaleGroup', value); await evidence(`Color scale ${value}`, true);
  }
  await mutate('Secondary parameter invalidates completed analysis', () => select('parameter2Select', 'delta'), false);
  for (const [group, value] of [['circularModeGroup', 'signed_180'], ['smoothingSigmaGroup', '0.8'], ['displayScaleGroup', 'density'], ['binDetailGroup', 'standard']]) {
    await mutate(`Mathematical display ${group} invalidates`, () => click(group, value), false);
  }
  await mutate('Primary context invalidates', () => click('contextGroup', 'U'), false);
  await click('contextGroup', 'U');
  await select('familySelect', 'base_pair'); await select('parameterSelect', 'opening');
  await click('jointJoinModeGroup', 'relation'); await select('family2Select', 'backbone'); await select('parameter2Select', 'chi');
  for (const orientation of ['pair-to-residue', 'residue-to-pair']) {
    if (orientation === 'residue-to-pair') {
      await select('familySelect', 'backbone'); await select('parameterSelect', 'chi');
      await select('family2Select', 'base_pair'); await select('parameter2Select', 'opening');
    }
    // Independent source relationship membership rules out accidental identity joins.
    assert(await page.evaluate(async () => {
      const a = window.rnaExplorer, data = await a.repository.loadRelations('observations');
      const key = p => JSON.stringify([p.pair_id, p.residue_id, p.endpoint_role ?? (p.side === 1 ? 'first' : 'second')]);
      const links = new Set((Array.isArray(data) ? data : data.rows).filter(p => p.kind === 'pair_residue').map(key));
      return a.snapshots.joint.result.points.every(p => links.has(key(p)));
    }));
    await mutate(`${orientation} endpoint invalidates`, () => click('jointResidueSideGroup', 'nt1'));
    await mutate(`${orientation} residue context invalidates`, () => click('jointResidueContextGroup', 'U'));
    const pucker = await page.evaluate(() => window.rnaExplorer.snapshots.joint.result.points.map(p => p.left.pucker_class || p.right.pucker_class).find(Boolean));
    assert(pucker, `${orientation}: finite endpoint pucker`);
    await mutate(`${orientation} endpoint pucker invalidates`, () => select('jointResiduePuckerGroup', pucker));
    await remember(); await click('jointPaletteGroup', 'warm'); await evidence(`${orientation} restyle reuses filtered incidences`, true);
    const points = await page.evaluate(() => window.rnaExplorer.snapshots.joint.result.points.map(p => ({ pair: p.pair_id, residue: p.residue_id, role: p.endpoint_role, left: p.left_id, right: p.right_id, x: p.x, y: p.y })));
    const exported = await downloadCsv(page, '#jointCsvDownload', path.join(output, `${orientation}.csv`));
    assert.equal(exported.rows.length, points.length);
    exported.rows.forEach((row, i) => {
      const p = points[i]; assert.equal(row.pair_id, p.pair); assert.equal(row.residue_id, p.residue); assert.equal(row.endpoint_role, p.role);
      assert.equal(row.x_id, p.left); assert.equal(row.y_id, p.right); assert.equal(Number(row.x_value), p.x); assert.equal(Number(row.y_value), p.y);
    });
    await click('jointResidueSideGroup', 'both'); await click('jointResidueContextGroup', 'U'); await select('jointResiduePuckerGroup', 'all');
    await click('jointPaletteGroup', 'hotspots');
  }
  for (const endpoint of ['nt1', 'nt2']) {
    await click('jointResidueSideGroup', endpoint);
    await click('jointJoinModeGroup', 'identity');
    await select('family2Select', 'backbone'); await select('parameter2Select', 'alpha');
    const identity = await page.evaluate(() => {
      const a = window.rnaExplorer, r = a.snapshots.joint.result;
      return { endpoint: a.state.joint.endpoint, points: r.points.length, exact: r.points.every(p => p.left_id === p.right_id && p.x === p.left.values.chi && p.y === p.right.values.alpha) };
    });
    assert.equal(identity.endpoint, endpoint); assert(identity.points > 0); assert(identity.exact);
    await click('jointJoinModeGroup', 'relation'); await select('family2Select', 'base_pair'); await select('parameter2Select', 'opening');
    const relation = await page.evaluate(() => {
      const a = window.rnaExplorer; return { endpoint: a.state.joint.endpoint, points: a.snapshots.joint.result.points.length, roles: [...new Set(a.snapshots.joint.result.points.map(p => p.endpoint_role))] };
    });
    assert.equal(relation.endpoint, endpoint); assert(relation.points > 0); assert.deepEqual(relation.roles, [endpoint === 'nt1' ? 'first' : 'second']);
    report.checks.push({ name: `Identity ignores and preserves ${endpoint} relation preference`, identity, relation });
  }
  assert.deepEqual(report.errors, []); report.passed = true;
  console.log(`PASS ${report.checks.length} joint style reuse browser checks`);
} catch (error) { report.passed = false; report.failure = error.stack; throw error; }
finally { report.finished = new Date().toISOString(); await writeFile(path.join(output, 'report.json'), JSON.stringify(report, null, 2) + '\n'); await browser.close(); }
