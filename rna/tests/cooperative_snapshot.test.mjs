import test from 'node:test';
import assert from 'node:assert/strict';
import { createPlotSnapshot, createOwnedPlotSnapshot, restyleTraceSnapshot, csv } from '../core/export.js';
import { deepFreeze, deepFreezeOwned, isImmutableData } from '../core/repository.js';
import { distribution } from '../core/analysis.js';

const rows = () => [{ id: 'u', comp_id: 'U', values: { chi: -0 } }, { id: 'g', comp_id: 'G', values: { chi: -179.1234567890123 } }];
const options = result => ({ result, snapshot_id: 'snapshot', created_at: 'fixed', build_id: 'release',
  selectionSpec: { contexts: ['U', 'G'] }, parameterDefinitionIds: ['chi'], dataHashes: { family: 'fixed' },
  provenance: { molecule: 'RNA' }, coordinatePolicy: { model: 1 } });

test('Cooperative owned snapshot matches independent synchronous snapshot values and CSV', async () => {
  const input = rows(), result = distribution(input, { id: 'chi', period: 360 }, { groupBy: 'none' });
  const expected = createPlotSnapshot({ ...options(structuredClone(result)), transferResult: true });
  let yields = 0;
  const actual = await createOwnedPlotSnapshot(options(result), { maxOperations: 2, checkpoint: async () => { yields++; } });
  assert(yields > 1); assert.deepEqual(actual, expected); assert.equal(csv(actual), csv(expected));
  assert.equal(actual.result, result); assert.equal(actual.result.series[0].rows[0], input[0]);
  assert(Object.is(actual.result.series[0].values[0], -0));
  assert(Object.isFrozen(actual)); assert(Object.isFrozen(input[0].values));
  assert.equal(restyleTraceSnapshot(actual, 'line').result.series, result.series);
});

test('Metadata is copied before first yield while scientific ownership preserves aliases', async () => {
  const shared = { value: Number.MIN_VALUE }, result = { kind: 'distribution', series: [], left: shared, right: shared };
  const config = options(result); let changed = false;
  const snapshot = await createOwnedPlotSnapshot(config, { maxOperations: 1, checkpoint: async () => {
    if (!changed) { config.selectionSpec.contexts.push('A'); config.dataHashes.family = 'mutated'; changed = true; }
  } });
  assert.deepEqual(snapshot.selection_spec.contexts, ['U', 'G']); assert.equal(snapshot.data_hashes.family, 'fixed');
  assert.equal(snapshot.result.left, shared); assert.equal(snapshot.result.right, shared);
  assert.equal(shared.value, Number.MIN_VALUE); assert(Object.isFrozen(shared));
});

test('Revision cancellation before start and during work never returns a reusable snapshot', async () => {
  const result = { kind: 'distribution', series: [], rows: rows() };
  await assert.rejects(createOwnedPlotSnapshot(options(result), { current: () => false }), { name: 'AbortError' });
  assert.equal(Object.isFrozen(result), false);
  let current = true, yielded = 0;
  await assert.rejects(createOwnedPlotSnapshot(options(result), { maxOperations: 1, current: () => current,
    checkpoint: async () => { yielded++; current = false; } }), { name: 'AbortError' });
  assert.equal(yielded, 1);
  assert.throws(() => restyleTraceSnapshot({ result }, 'line'), /owned frozen snapshot/);
  await assert.rejects(createOwnedPlotSnapshot(options({ kind: 'distribution', series: [] }), {
    maxOperations: 1, checkpoint: async () => false,
  }), { name: 'AbortError' });
});

test('Final revision check closes cancellation after the last cooperative freeze microtask', async () => {
  let calls = 0;
  await assert.rejects(createOwnedPlotSnapshot(options({ kind: 'distribution', series: [] }), {
    current: () => ++calls < 3, maxOperations: 100000, timeBudgetMs: 100000,
  }), { name: 'AbortError' });
  assert.equal(calls, 3);
});

test('Traversal preserves synchronous getter order cycles sparse arrays and foreign frozen children', async () => {
  function graph(log) {
    const child = { get leaf() { log.push('child'); return { value: -0 }; } };
    const root = { get first() { log.push('first'); return child; }, get second() { log.push('second'); return child; } };
    root.self = root; root.array = []; root.array[3] = Object.freeze({ nested: { value: 5 } });
    root.array.extra = { value: 7 }; return root;
  }
  const expectedLog = [], actualLog = [], expected = graph(expectedLog), actual = graph(actualLog);
  deepFreeze(expected); await deepFreezeOwned(actual, { maxOperations: 1, checkpoint: async () => {} });
  assert.deepEqual(actualLog, expectedLog); assert.equal(actual.self, actual);
  assert.equal(actual.array.length, 4); assert.equal(Object.hasOwn(actual.array, 0), false);
  assert(Object.isFrozen(actual.array.extra)); assert(Object.isFrozen(actual.array[3].nested));
  assert.equal(isImmutableData(actual), false, 'Async freeze must not certify accessors or cycles');
});

test('Cooperative freeze skips only already certified data without trusting shallow freezes', async () => {
  let reads = 0;
  const certified = new Proxy({ value: { x: 1 } }, { get(target, key, receiver) { reads++; return Reflect.get(target, key, receiver); } });
  deepFreeze(certified); assert(isImmutableData(certified)); const before = reads;
  const child = { x: 2 }, shallow = Object.freeze({ child });
  await deepFreezeOwned({ certified, shallow }, { maxOperations: 1, checkpoint: async () => {} });
  assert.equal(reads, before); assert(Object.isFrozen(child));
});

test('Real default task checkpoints allow timers and large primitive arrays consume budget', async () => {
  let timerRan = false; const timer = new Promise(resolve => setTimeout(() => { timerRan = true; resolve(); }, 0));
  const values = Array.from({ length: 100 }, (_, i) => i), result = { kind: 'distribution', series: [], values };
  const pending = createOwnedPlotSnapshot(options(result), { maxOperations: 20 });
  await timer; assert(timerRan); const snapshot = await pending; assert(Object.isFrozen(snapshot.result.values));
  let checkpoints = 0;
  await deepFreezeOwned(Array(100).fill(0), { maxOperations: 10, checkpoint: async () => { checkpoints++; } });
  assert(checkpoints >= 10);
});

test('Scheduling validation and checkpoint failures reject without altering the sync API', async () => {
  for (const scheduling of [{ maxOperations: 0 }, { maxOperations: 0.5 }, { timeBudgetMs: NaN }, { timeBudgetMs: 0 }, { current: null }, { checkpoint: null }]) {
    await assert.rejects(createOwnedPlotSnapshot(options({}), scheduling), /budget|callbacks/);
  }
  await assert.rejects(createOwnedPlotSnapshot(null), /completed result/);
  await assert.rejects(createOwnedPlotSnapshot(options({}), { maxOperations: 1, checkpoint: async () => { throw new Error('scheduler failed'); } }), /scheduler failed/);
  const result = { kind: 'distribution', series: [] }, snapshot = createPlotSnapshot(options(result));
  assert.notEqual(snapshot.result, result); assert.equal(Object.isFrozen(result), false);
});
