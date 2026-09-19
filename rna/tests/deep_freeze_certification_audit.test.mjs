import test from 'node:test';
import assert from 'node:assert/strict';
import { deepFreeze, isImmutableData } from '../core/repository.js';
import { createPlotSnapshot } from '../core/export.js';

test('Certification includes symbol and hidden mutable data edges', () => {
  for (const key of ['hidden', Symbol('hidden')]) {
    const child = { value: 1 }, root = {};
    Object.defineProperty(root, key, { value: child, enumerable: false });
    Object.freeze(root);
    assert.equal(isImmutableData(root), false, 'Foreign shallow freeze cannot hide a mutable child');
    deepFreeze(root);
    if (!Object.isFrozen(child)) child.value = 2;
    assert.equal(isImmutableData(root), Object.isFrozen(child));
  }
  const key = Symbol('visible'), child = { value: 1 }, root = { [key]: child };
  Object.freeze(root);
  assert.equal(isImmutableData(root), false, 'Enumerable symbols must also be inspected');
  deepFreeze(root);
  assert.equal(isImmutableData(root), Object.isFrozen(child));
});

test('Accessors and setters remain untrusted even when frozen', () => {
  for (const key of ['field', Symbol('field')]) {
    let external = { value: 1 };
    const root = {};
    Object.defineProperty(root, key, {
      enumerable: true, configurable: true,
      get: () => external, set: value => { external = value; },
    });
    deepFreeze(root);
    assert.equal(isImmutableData(root), false);
    root[key] = { value: 2 };
    assert.equal(root[key].value, 2, 'Object.freeze does not prevent accessor-backed mutation');
    assert.equal(isImmutableData(root), false);
  }
});

test('Functions exotic objects and their containing graphs are never certified', () => {
  const exotic = [new Date(0), new Map([['x', 1]]), new Set([1]),
    new (class Custom { constructor() { this.value = 1; } })(),
    Object.assign(Object.create({ inherited: true }), { value: 1 }), () => 1];
  for (const value of exotic) {
    const root = { value };
    deepFreeze(root);
    assert.equal(isImmutableData(value), false);
    assert.equal(isImmutableData(root), false);
  }
});

test('Array custom properties follow the complete descriptor contract', () => {
  const child = { value: 1 }, valid = [child];
  valid.extra = child;
  deepFreeze(valid);
  assert.equal(isImmutableData(valid), true);
  assert.throws(() => { valid.extra.value = 2; }, TypeError);
  const hidden = [], hiddenChild = { value: 1 };
  Object.defineProperty(hidden, 'hidden', { value: hiddenChild });
  deepFreeze(hidden);
  assert.equal(isImmutableData(hidden), Object.isFrozen(hiddenChild));
  const accessor = [1];
  Object.defineProperty(accessor, 'extra', { get: () => 3 });
  deepFreeze(accessor);
  assert.equal(isImmutableData(accessor), false);
  const callable = [1]; callable.extra = () => 3;
  deepFreeze(callable);
  assert.equal(isImmutableData(callable), false);
  const subclass = new (class CustomArray extends Array {})();
  subclass.push(1);
  const inherited = [1];
  Object.setPrototypeOf(inherited, Object.create(Array.prototype));
  for (const exotic of [subclass, inherited]) {
    deepFreeze(exotic);
    assert.equal(isImmutableData(exotic), false, 'Exotic array prototypes are outside plain-data certification');
  }
});

test('Repeated aliases form a certifiable DAG while cycles remain unsupported', () => {
  const leaf = { value: 1 }, root = { first: leaf, second: [leaf] };
  deepFreeze(root);
  assert.equal(isImmutableData(root), true);
  assert.equal(root.first, root.second[0]);
  const self = {}; self.self = self;
  deepFreeze(self);
  assert(Object.isFrozen(self));
  assert.equal(isImmutableData(self), false);
  const first = {}, second = { first }; first.second = second;
  deepFreeze(first);
  assert(Object.isFrozen(first) && Object.isFrozen(second));
  assert.equal(isImmutableData(first), false);
  assert.equal(isImmutableData(second), false);
});

test('Foreign frozen roots with ordinary mutable children are traversed', () => {
  const child = { nested: [-0, null, Number.MIN_VALUE] }, root = Object.freeze({ child });
  assert.equal(isImmutableData(root), false);
  assert.equal(deepFreeze(root), root);
  assert.equal(isImmutableData(root), true);
  assert(Object.isFrozen(child) && Object.isFrozen(child.nested));
  assert(Object.is(child.nested[0], -0));
  assert.throws(() => { child.nested.push(2); }, TypeError);
  const nullPrototype = Object.assign(Object.create(null), { child: { x: 1 } });
  deepFreeze(nullPrototype);
  assert.equal(isImmutableData(nullPrototype), true);
});

test('Snapshot transfer freezes the owned graph and preserves aliases exactly', () => {
  const row = { id: 'r1', value: -0 }, result = { kind: 'audit', rows: [row], again: row };
  const controls = { contexts: ['U'] };
  const snapshot = createPlotSnapshot({ result, transferResult: true, selectionSpec: controls });
  assert.equal(snapshot.result, result);
  assert.equal(snapshot.result.rows[0], snapshot.result.again);
  assert.equal(isImmutableData(snapshot), true);
  assert(Object.is(snapshot.result.rows[0].value, -0));
  assert.throws(() => { result.rows[0].value = 1; }, TypeError);
  controls.contexts.push('A');
  assert.deepEqual(snapshot.selection_spec.contexts, ['U']);
});

test('Snapshot copy does not borrow an accessor-backed foreign frozen result', () => {
  let external = 1;
  const source = Object.freeze({ kind: 'audit', get value() { return external; } });
  assert.equal(isImmutableData(source), false);
  const snapshot = createPlotSnapshot({ result: source });
  external = 2;
  assert.equal(source.value, 2);
  assert.equal(snapshot.result.value, 1);
  assert.notEqual(snapshot.result, source);
  assert.equal(isImmutableData(snapshot), true);
});
