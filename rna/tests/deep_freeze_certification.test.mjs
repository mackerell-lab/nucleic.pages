import test from 'node:test';
import assert from 'node:assert/strict';
import { deepFreeze, isImmutableData } from '../core/repository.js';
import { createPlotSnapshot, csv } from '../core/export.js';
import { distribution } from '../core/analysis.js';

test('External shallow freeze still descends to mutable children before certification', () => {
  const child = { nested: [1, 2] }, root = Object.freeze({ child });
  assert.equal(isImmutableData(root), false);
  assert.equal(deepFreeze(root), root);
  assert(Object.isFrozen(child)); assert(Object.isFrozen(child.nested));
  assert(isImmutableData(root));
  assert.throws(() => { child.nested.push(3); }, TypeError);
});

test('Cycles terminate safely without claiming immutable plain data certification', () => {
  const root = {}, child = { root }; root.child = child; root.self = root;
  assert.equal(deepFreeze(root), root);
  assert(Object.isFrozen(root)); assert(Object.isFrozen(child));
  assert.equal(isImmutableData(root), false); assert.equal(isImmutableData(child), false);
  assert.equal(deepFreeze(root), root);
  const shared = { value: 1 };
  const aliases = { left: shared, right: shared };
  assert.equal(deepFreeze(aliases), aliases); assert(isImmutableData(aliases));
});

test('Dates maps functions exotic prototypes and accessors never gain a plain data certificate', () => {
  const date = new Date(0), map = new Map([['a', 1]]), exotic = Object.create({ inherited: true }); exotic.value = 1;
  for (const value of [date, map, exotic, { callback() {} }]) {
    const root = { value };
    deepFreeze(root); assert.equal(isImmutableData(root), false);
    deepFreeze(root); assert.equal(isImmutableData(root), false);
  }
  date.setTime(10); map.set('b', 2);
  assert.equal(date.getTime(), 10); assert.equal(map.get('b'), 2);
  let reads = 0;
  const accessor = { get child() { reads++; return { generation: reads }; } };
  deepFreeze(accessor); assert.equal(reads, 1);
  assert.equal(isImmutableData(accessor), false); assert.equal(reads, 1, 'Certification must not call accessors');
  deepFreeze(accessor); assert.equal(reads, 2, 'Uncertified accessor keeps the established traversal semantics');
});

test('Hidden mutable edges cannot be certified, including symbol descriptors', () => {
  for (const key of ['hidden', Symbol('hidden')]) {
    const child = { value: 1 }, root = {};
    Object.defineProperty(root, key, { value: child, enumerable: false }); Object.freeze(root);
    assert.equal(isImmutableData(root), false, 'Regression: symbol descriptors were previously skipped');
    deepFreeze(root);
    assert.equal(Object.isFrozen(child), false, 'Preserve existing enumerable-only freeze traversal');
    assert.equal(isImmutableData(root), false);
    child.value = 2; assert.equal(root[key].value, 2);
  }
  const key = Symbol('hidden'), child = deepFreeze({ value: 1 }), root = Object.freeze({ [key]: child });
  deepFreeze(root); assert(isImmutableData(root));
  const getter = Object.freeze({ get [key]() { return { mutable: true }; } });
  assert.equal(isImmutableData(getter), false);
  const callback = Object.freeze({ [key]: () => 1 });
  assert.equal(isImmutableData(callback), false);
});

test('Null prototypes sparse arrays and aliases preserve original values and identities', () => {
  const record = Object.assign(Object.create(null), { label: 'RNA', nested: { value: null } });
  const rows = []; rows[2] = record;
  deepFreeze(rows);
  assert(isImmutableData(rows)); assert.equal(rows[2], record);
  assert.equal(Object.hasOwn(rows, 0), false); assert.equal(rows.length, 3);
  assert.equal(Object.getPrototypeOf(record), null);
});

test('Snapshot ownership still preserves raw measurements and rejects mutable hidden roots', () => {
  const source = deepFreeze([{ id: 'u1', comp_id: 'U', values: { chi: -179.1234567890123 } }]);
  const result = distribution(source, { id: 'chi', period: 360 }, { groupBy: 'none' });
  const options = { result, snapshot_id: 'fixed', build_id: 'test', created_at: 'fixed' };
  const copied = createPlotSnapshot(options), transferred = createPlotSnapshot({ ...options, transferResult: true });
  assert.equal(transferred.result, result); assert.equal(csv(copied), csv(transferred));
  assert.equal(copied.result.series[0].rows[0], source[0]);
  const symbol = Symbol('mutable'), hidden = { value: 1 }, row = Object.freeze({ id: 'r', [symbol]: hidden });
  const snapshot = createPlotSnapshot({ result: { rows: [row] } });
  assert.notEqual(snapshot.result.rows[0], row);
  hidden.value = 2; assert.equal(isImmutableData(row), false);
});
