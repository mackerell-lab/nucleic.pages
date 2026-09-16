import assert from 'node:assert/strict';
import { waitReady } from './helpers.mjs';

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
      return Number.isFinite(row.values.opening) && (policy === 'all' || policy === 'near' ? policy !== 'near' || near : !near && !alternative)
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
