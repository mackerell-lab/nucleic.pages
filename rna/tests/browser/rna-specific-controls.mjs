import assert from 'node:assert/strict';
import { waitReady } from './helpers.mjs';

// Audited protected DNA palette: js/pure-dna.js:1981. This fixture never imports
// or executes the DNA application's unconditional bootstrap.
const dnaHotspots = [[0, '#00205b'], [0.0526, '#003b8e'], [0.1053, '#0051a8'], [0.1579, '#0069b4'],
  [0.2105, '#0080b9'], [0.2632, '#0097bd'], [0.3158, '#00afb8'], [0.3684, '#00c6a7'],
  [0.4211, '#00dd8c'], [0.4737, '#2ff272'], [0.5263, '#7fff5a'], [0.5789, '#c7ff54'],
  [0.6316, '#ffe64c'], [0.6842, '#ffb53b'], [0.7368, '#ff7b2e'], [0.7895, '#ff4b2c'],
  [0.8421, '#f2252e'], [0.8947, '#d8002c'], [0.9474, '#b30027'], [1, '#7f001d']];

export async function checkPalettes(page, record) {
  const expected = { hotspots: dnaHotspots, warm: 'YlOrRd', viridis: 'Viridis', cividis: 'Cividis', ocean: 'YlGnBu', forest: 'Greens', greys: 'Greys' };
  const choices = await page.locator('#jointPaletteGroup button').evaluateAll(buttons => buttons.map(button => button.dataset.value));
  assert.deepEqual(choices, Object.keys(expected), 'RNA palette choices differ from DNA');
  const saved = await page.evaluate(() => {
    globalThis.rnaPaletteBaseline = window.rnaExplorer.snapshots.joint.result.points.map(point => [point.left_id, point.right_id, point.x, point.y]);
    return window.rnaExplorer.state.joint.palette;
  });
  const rendered = [];
  try {
    for (const [choice, colorscale] of Object.entries(expected)) {
      await page.click(`#jointPaletteGroup button[data-value="${choice}"]`); await waitReady(page);
      const evidence = await page.evaluate(() => {
        const points = window.rnaExplorer.snapshots.joint.result.points, baseline = globalThis.rnaPaletteBaseline;
        const mismatch = points.length !== baseline.length || points.some((point, index) => [point.left_id, point.right_id, point.x, point.y].some((value, key) => value !== baseline[index][key]));
        const plot = document.querySelector('#jointPlot');
        return { inputScale: plot.data[0].colorscale, renderedScale: plot._fullData[0].colorscale, mismatch, points: points.length };
      });
      assert.deepEqual(evidence.inputScale, colorscale, `${choice} does not use the DNA colorscale`);
      assert(Array.isArray(evidence.renderedScale) && evidence.renderedScale.length > 1, 'Plotly did not resolve the colorscale');
      assert.equal(evidence.mismatch, false, `${choice} changed raw joined observations`);
      rendered.push({ choice, colorscale: evidence.renderedScale, points: evidence.points });
    }
    assert.equal(new Set(rendered.map(item => JSON.stringify(item.colorscale))).size, 7, 'Different palette choices produced identical scales');
    record('Seven DNA-matching palettes preserve joint observations', { rendered });
  } finally {
    await page.evaluate(() => { delete globalThis.rnaPaletteBaseline; });
    await page.click(`#jointPaletteGroup button[data-value="${saved}"]`); await waitReady(page);
  }
}

export async function checkBroadResidueScope(page, record) {
  const started = Date.now();
  const evidence = await page.evaluate(async () => {
    const app = window.rnaExplorer, saved = structuredClone(app.state);
    app.state.familyId = 'backbone'; app.state.parameterId = 'chi';
    app.state.family2Id = 'backbone'; app.state.parameter2Id = 'delta'; app.state.joint.mode = 'identity';
    app.state.selection = { components: 'all', methods: [], resolutionMax: null, contexts: [], functions: [], subtypes: [], structures: [], includeEnds: true, puckerStates: [], pairPolicy: 'exact', interactionFamilies: [], stemOnly: false };
    app.updateSelectors();
    const renderingStarted = performance.now();
    await app.requestRender();
    const renderMs = performance.now() - renderingStarted;
    if (document.querySelector('#appStatus').dataset.state !== 'ready') throw new Error(document.querySelector('#appStatus').textContent);
    const rows = (await app.repository.loadFamily('backbone')).rows;
    const finite = (row, parameter) => Number.isFinite(row.values?.[parameter]) && (!row.statuses?.[parameter] || ['available', 'computed', 'ok', 'valid'].includes(row.statuses[parameter]));
    let expectedChi = 0, expectedJoint = 0;
    for (const row of rows) { const chi = finite(row, 'chi'); if (chi) expectedChi++; if (chi && finite(row, 'delta')) expectedJoint++; }
    const result = { loadedRows: rows.length, expectedChi, actualChi: app.snapshots.distribution.result.coverage.plottedRows,
      expectedJoint, actualJoint: app.snapshots.joint.result.points.length, renderMs,
      heapUsedBytes: performance.memory?.usedJSHeapSize ?? null, semantics: 'All released entries, all methods and component profiles; no CSV duplication for this broad probe.' };
    app.state = saved; app.updateSelectors(); await app.requestRender();
    return result;
  });
  await waitReady(page);
  assert.equal(evidence.actualChi, evidence.expectedChi, 'Broad chi view lost finite observations');
  assert.equal(evidence.actualJoint, evidence.expectedJoint, 'Broad same-residue joint lost finite identity matches');
  assert(evidence.actualJoint > 0);
  record('Broad all-method residue and joint scope', { ...evidence, probeAndRestoreMs: Date.now() - started });
}

export async function checkPairControls(page, record) {
  const saved = await page.evaluate(() => structuredClone(window.rnaExplorer.state));
  await page.evaluate(async () => {
    const app = window.rnaExplorer;
    app.state.family2Id = ''; app.updateSelectors();
    await app.setSelection({ components: 'all', methods: [], resolutionMax: null, resolution: 'any', resolutionKnown: false,
      functions: [], structures: [], subtypes: [], contexts: [], puckerStates: [], includeEnds: true,
      pairPolicy: 'all', interactionFamilies: [], stemOnly: false });
  });
  await waitReady(page);
  const expectedAndActual = async (policy, stemOnly = false, family = null) => page.evaluate(async ({ policy, stemOnly, family }) => {
    const app = window.rnaExplorer, table = await app.repository.loadFamily('base_pair');
    // Independent oracle over deposited flags and raw finite opening values.
    const rows = table.rows.filter(row => {
      const label = String(row.family || row.interaction_family || '');
      const near = row.near === true || /^n[ct]/i.test(label);
      const alternative = row.alternative === true || /^[n]?[ct][WHS]{2}a$/i.test(label);
      return Number.isFinite(row.values.opening) && (!row.statuses?.opening || ['available', 'computed', 'ok', 'valid'].includes(row.statuses.opening))
        && (policy === 'all' || policy === 'near' ? policy !== 'near' || near : !near && !alternative)
        && (!stemOnly || row.stem_eligible === true) && (!family || label === family);
    });
    return { expected: rows.map(row => row.id).sort(), actual: [...new Set(app.snapshots.distribution.result.series.flatMap(series => series.rowIds))].sort() };
  }, { policy, stemOnly, family });
  const counts = {};
  for (const policy of ['exact', 'all', 'near']) {
    await page.click(`#pairPolicyGroup button[data-value="${policy}"]`); await waitReady(page);
    const values = await expectedAndActual(policy); assert.deepEqual(values.actual, values.expected, `${policy} pair policy disagrees with raw flags`);
    counts[policy] = values.actual.length;
  }
  if (!await page.evaluate(() => window.rnaExplorer.manifest.partial)) assert(counts.near > 0, 'Full-release near-assignment branch has no positive observations');
  await page.click('#pairPolicyGroup button[data-value="all"]'); await waitReady(page);
  await page.click('#stemScopeGroup button[data-value="stem"]'); await waitReady(page);
  let values = await expectedAndActual('all', true); assert.deepEqual(values.actual, values.expected);
  counts.stem = values.actual.length;
  await page.click('#stemScopeGroup button[data-value="all"]'); await waitReady(page);
  const family = await page.locator('#interactionFamilyGroup option').evaluateAll(options => options.find(option => option.value !== 'all')?.value);
  assert(family, 'No recorded LW family choices');
  await page.selectOption('#interactionFamilyGroup', family); await waitReady(page);
  values = await expectedAndActual('all', false, family); assert.deepEqual(values.actual, values.expected);
  await page.selectOption('#interactionFamilyGroup', 'all'); await waitReady(page);
  await page.click('#groupingGroup button[data-value="interactionFamily"]'); await waitReady(page);
  const groups = await page.evaluate(() => window.rnaExplorer.snapshots.distribution.result.series.map(series => ({ key: series.key, n: series.rowIds.length })));
  assert(groups.length && groups.every(group => group.key !== 'Unknown'), 'Known interactions became Unknown groups');
  assert.equal(groups.reduce((sum, group) => sum + group.n, 0), counts.all);
  record('Exact, near, LW-family and stem controls match raw data', { counts, filteredFamily: family, groups });

  // Restore the endpoint join, then exclude terminal endpoints using the UI.
  await page.evaluate(async saved => { const app = window.rnaExplorer; app.state = saved; app.updateSelectors(); await app.requestRender(); }, saved);
  await waitReady(page);
  await page.click('#terminalGroup button[data-value="exclude"]'); await waitReady(page);
  const endpoints = await page.evaluate(() => window.rnaExplorer.snapshots.joint.result.points.map(point => ({
    pairTerminal: point.left.is_terminal_any, residueTerminal: point.right.is_terminal_any, pairId: point.pair_id, residueId: point.residue_id,
  })));
  assert(endpoints.length > 0, 'Terminal-exclusion test did not exercise any endpoint observations');
  assert(endpoints.every(point => point.pairTerminal !== true && point.residueTerminal !== true), 'Terminal pair or residue survived endpoint exclusion');
  record('Terminal exclusion applies to endpoint joins', { incidences: endpoints.length });
  await page.click('#terminalGroup button[data-value="include"]'); await waitReady(page);
}

export async function checkPuckerSurvey(page, record) {
  const choices = await page.locator('#puckerGroup option').evaluateAll(options => options.map(option => option.value).filter(value => value !== 'all'));
  assert(choices.length, 'Pucker classes are missing');
  const chosen = choices.find(value => value === "C3'-endo") || choices[0];
  await page.selectOption('#puckerGroup', chosen); await waitReady(page);
  const evidence = await page.evaluate(async chosen => {
    const app = window.rnaExplorer, snapshot = app.snapshots.survey;
    const table = await app.repository.loadSurveyScalars(snapshot.result.parameter.id);
    const raw = new Map((Array.isArray(table) ? table : table.rows).map(row => [row.id, row]));
    const plotted = snapshot.result.series.flatMap(series => series.rows);
    return { chosen, selection: snapshot.selection_spec.puckerStates, plottedRows: plotted.length,
      bad: plotted.filter(row => {
        const source = raw.get(row.id);
        return !source || !(source.pucker_classes || [source.pucker_class || source.pucker_state]).every(value => value === chosen);
      }).map(row => row.id) };
  }, chosen);
  assert.deepEqual(evidence.selection, [chosen]); assert.deepEqual(evidence.bad, [], 'Survey pucker filter disagrees with raw residue classes');
  assert(evidence.plottedRows > 0, 'Pucker survey test did not exercise any observations');
  record('Pucker filter propagates to residue survey', evidence);
  await page.selectOption('#puckerGroup', 'all'); await waitReady(page);
}
