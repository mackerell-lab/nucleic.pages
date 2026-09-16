import test from 'node:test';
import assert from 'node:assert/strict';
import { deepFreeze } from '../core/repository.js';
import { createPlotSnapshot, csv } from '../core/export.js';
import { distribution } from '../core/analysis.js';

const parameter = { id: 'chi', period: 360 };
test('Owned results freeze in place while captured control state remains independent', () => {
  const result = distribution([{ id: 'U1', comp_id: 'U', values: { chi: -173.1234567890123 } }], parameter, { groupBy: 'none' });
  const controls = { contexts: ['U'] };
  const baseline = createPlotSnapshot({ result, selectionSpec: controls, snapshot_id: 'same', build_id: 'full' });
  const snapshot = createPlotSnapshot({ result, selectionSpec: controls, snapshot_id: 'same', build_id: 'full', transferResult: true });
  assert.equal(snapshot.result, result);
  assert.equal(csv(snapshot), csv(baseline));
  controls.contexts.push('G');
  assert.deepEqual(snapshot.selection_spec.contexts, ['U']);
  assert.throws(() => { result.series[0].values[0] = 1; }, TypeError);
  assert.throws(() => { result.series[0].rows[0].values.chi = 2; }, TypeError);
});

test('Snapshots share certified immutable source data and copy mutable analysis arrays', () => {
  const table = deepFreeze({ rows: [{ id: 'U1', pdb_id: '1RNA', comp_id: 'U', values: { chi: -173.1234567890123 }, statuses: { chi: 'ok' } }] });
  const selected = [{ ...table.rows[0], functions: ['ribozyme'] }];
  const result = distribution(selected, parameter, { groupBy: 'none' });
  const snapshot = createPlotSnapshot({ result });
  const stored = snapshot.result.series[0].rows[0];
  assert.equal(stored.values, table.rows[0].values, 'Immutable scientific values were duplicated');
  assert.equal(stored.statuses, table.rows[0].statuses);
  assert.notEqual(stored, selected[0]); assert.notEqual(snapshot.result.series[0].values, result.series[0].values);
  const before = csv(snapshot);
  result.series[0].values[0] = 42; selected[0].functions.push('other');
  assert.equal(csv(snapshot), before); assert.deepEqual(stored.functions, ['ribozyme']);
  assert.throws(() => { stored.values.chi = 0; }, TypeError);
});

test('A shallow frozen caller object cannot leak mutable children into snapshots', () => {
  const values = { chi: 123.5 };
  const row = Object.freeze({ id: 'r', values });
  const snapshot = createPlotSnapshot({ result: { kind: 'distribution', parameter, series: [{ key: 'All', rows: [row], values: [values.chi] }] } });
  values.chi = 99;
  assert.equal(snapshot.result.series[0].rows[0].values.chi, 123.5);
  const parent = Object.freeze({ child: { id: 'nested' } });
  deepFreeze(parent);
  assert(Object.isFrozen(parent.child), 'Deep freeze skipped a shallow frozen root');
});

test('Mutable built-ins and getters are not certified as shareable data', () => {
  const row = deepFreeze({ when: new Date(0), get child() { return { value: 4 }; } });
  const snapshot = createPlotSnapshot({ result: { rows: [row] } });
  assert.notEqual(snapshot.result.rows[0], row);
  assert(Object.isFrozen(snapshot.result.rows[0].child));
});
