import test from 'node:test';
import assert from 'node:assert/strict';
import { deepFreeze, deepFreezeOwned, isImmutableData } from '../core/repository.js';
import { createOwnedPlotSnapshot, createPlotSnapshot, restyleTraceSnapshot, csv, provenance } from '../core/export.js';

const scheduling = { checkpoint: async () => true, maxOperations: 2, timeBudgetMs: 1000 };

test('Owned async snapshots preserve scientific aliases, signed zero, raw CSV and provenance', async () => {
  const row = { id: 'U1', pdb_id: 'TEST', values: { chi: -0 } };
  const result = { kind: 'distribution', parameter: { id: 'chi', period: 360 },
    displaySpec: { traceStyle: 'filled' }, series: [{ key: 'U', rows: [row, row], values: [-0, 2.123456789012345], weights: [1, 1], statistics: { n: 2 } }] };
  const options = { result, snapshot_id: 'review', created_at: 'fixed', buildId: 'test', selectionSpec: { contexts: ['U'] } };
  const baseline = createPlotSnapshot(options);
  const snapshot = await createOwnedPlotSnapshot(options, scheduling);
  assert.equal(snapshot.result, result);
  assert.equal(snapshot.result.series[0].rows[0], snapshot.result.series[0].rows[1]);
  assert(Object.is(snapshot.result.series[0].values[0], -0));
  assert(Object.is(row.values.chi, -0));
  assert.equal(csv(snapshot), csv(baseline));
  assert.equal(provenance(snapshot), provenance(baseline));
  const styled = restyleTraceSnapshot(snapshot, 'line');
  assert.equal(styled.result.series, snapshot.result.series);
  assert.equal(styled.result.displaySpec.traceStyle, 'line');
  assert.throws(() => { row.values.chi = 99; }, TypeError);
});

test('External snapshot metadata is captured before the first cooperative pause', async () => {
  let finish; const gate = new Promise(resolve => { finish = resolve; }); let pauses = 0;
  const selection = { contexts: ['U'] }, display = { traceStyle: 'filled' }, source = { labels: ['original'] };
  const promise = createOwnedPlotSnapshot({ result: { rows: [{ value: -0 }] }, selectionSpec: selection,
    displaySpec: display, provenance: source }, { ...scheduling, maxOperations: 1,
    checkpoint: async () => { if (++pauses === 1) await gate; } });
  assert.equal(pauses, 1);
  selection.contexts.push('A'); display.traceStyle = 'line'; source.labels[0] = 'changed';
  finish(); const snapshot = await promise;
  assert.deepEqual(snapshot.selection_spec.contexts, ['U']);
  assert.equal(snapshot.display_spec.traceStyle, 'filled');
  assert.deepEqual(snapshot.provenance.labels, ['original']);
});

test('Cooperative traversal matches enumerable getter order and preserves unsupported hidden edges', async () => {
  function graph(log) {
    const hidden = { value: 1 }, symbol = Symbol('hidden');
    const leaf = Object.assign(Object.create(null), { value: -0 });
    const array = []; array[3] = leaf; array.extra = leaf;
    const root = { get first() { log.push('first'); return { get nested() { log.push('nested'); return array; } }; },
      get second() { log.push('second'); return leaf; }, date: new Date(0) };
    root.self = root;
    Object.defineProperty(root, 'hidden', { value: hidden }); root[symbol] = hidden;
    return { root, leaf, array, hidden, symbol };
  }
  const syncLog = [], asyncLog = [], sync = graph(syncLog), owned = graph(asyncLog);
  deepFreeze(sync.root); await deepFreezeOwned(owned.root, scheduling);
  assert.deepEqual(asyncLog, syncLog); assert.deepEqual(asyncLog, ['first', 'second', 'nested']);
  for (const fixture of [sync, owned]) {
    assert.equal(fixture.root.self, fixture.root);
    assert.equal(fixture.array[3], fixture.array.extra);
    assert.equal(Object.hasOwn(fixture.array, 0), false);
    assert.equal(fixture.array.length, 4);
    assert.equal(Object.getPrototypeOf(fixture.leaf), null);
    assert(Object.is(fixture.leaf.value, -0));
    assert(Object.isFrozen(fixture.root)); assert(Object.isFrozen(fixture.array));
    assert.equal(Object.isFrozen(fixture.hidden), false);
    fixture.root.date.setTime(10); assert.equal(fixture.root.date.getTime(), 10);
    assert.equal(isImmutableData(fixture.root), false);
  }
});

test('Cancellation and checkpoint failure never return a restylable snapshot', async () => {
  for (const checkpoint of [async () => false, async () => { throw Error('checkpoint failed'); }]) {
    const result = { rows: [{ value: 1 }] };
    await assert.rejects(createOwnedPlotSnapshot({ result }, { ...scheduling, maxOperations: 1, checkpoint }));
    assert.equal(Object.isFrozen(result), false);
    assert.throws(() => restyleTraceSnapshot(Object.freeze({ result }), 'line'), /owned frozen/);
  }
  let checks = 0;
  await assert.rejects(createOwnedPlotSnapshot({ result: { value: 1 } }, {
    maxOperations: 10000, timeBudgetMs: 1000,
    current: () => ++checks < 3,
  }), { name: 'AbortError' });
  assert.equal(checks, 3, 'Final post-await capability check was bypassed');
});

test('A merely shallow-frozen root is traversed without false certification', async () => {
  const child = { nested: [-0] }, root = Object.freeze({ child });
  assert.equal(isImmutableData(root), false);
  await deepFreezeOwned(root, scheduling);
  assert(Object.isFrozen(child)); assert(Object.isFrozen(child.nested));
  assert(Object.is(child.nested[0], -0));
  assert.equal(isImmutableData(root), true);
  const typed = new Uint8Array([1]);
  assert.throws(() => deepFreeze(typed), TypeError);
  await assert.rejects(deepFreezeOwned(new Uint8Array([1]), scheduling), TypeError);
});
