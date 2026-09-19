import test from 'node:test';
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { decodeFamilyRows, FAMILY_COLUMNAR_ENCODING } from '../core/survey-codec.js';
import { BUNDLED_FAMILY_ENCODING, FAMILY_BUNDLE_ENCODING, verifyFamilyBundle,
  expandBundledFamilyColumns } from '../core/bundled-family-codec.js';

const hash = value => createHash('sha256').update(JSON.stringify(value)).digest('hex');
function bundleFor(...columns) {
  const bundle = { encoding: FAMILY_BUNDLE_ENCODING,
    columns: Object.fromEntries(columns.map(column => [hash(column), column]).sort(([a], [b]) => a.localeCompare(b))) };
  return { bundle, reference: hash(bundle), ref: column => ({ bundle: hash(bundle), column: hash(column) }) };
}
const table = columns => ({ encoding: BUNDLED_FAMILY_ENCODING, row_count: 3, build_id: 'fixture', columns });

test('Whole family bundle verifies mixed dictionaries and nested arrays with UTF8 byte accounting', async () => {
  const fixture = bundleFor(['α', null], { dictionary: [{ alpha: 1.2345678901234567 }, null], indices: [0, 1, 0] });
  assert.deepEqual(await verifyFamilyBundle(fixture.reference, fixture.bundle), {
    bundle: fixture.bundle, byteLength: Buffer.byteLength(JSON.stringify(fixture.bundle))
  });
  const corrupted = structuredClone(fixture.bundle);
  const dict = Object.values(corrupted.columns).find(column => !Array.isArray(column));
  dict.indices[0] = 1;
  await assert.rejects(verifyFamilyBundle(fixture.reference, corrupted), /hash mismatch/);
  const renamed = { ...fixture.bundle, columns: { ['f'.repeat(64)]: Object.values(fixture.bundle.columns)[0] } };
  await assert.rejects(verifyFamilyBundle(fixture.reference, renamed), /hash mismatch/);
  await assert.rejects(verifyFamilyBundle('../escape', fixture.bundle), /content hash/);
  for (const malformed of [null, [], {}, { encoding: FAMILY_BUNDLE_ENCODING, columns: [] },
    { encoding: FAMILY_BUNDLE_ENCODING, columns: { wrong: [] } }]) {
    await assert.rejects(verifyFamilyBundle(hash(malformed), malformed), /Invalid family/);
  }
});

test('Family expansion preserves identity, absent versus null, nested values and dictionaries', async () => {
  const ids = ['a', 'b', 'c'];
  const annotation = { dictionary: [null, { type: 'stem', atoms: ['O2\u2032', null] }], indices: [0, 0, 1] };
  const fixture = bundleFor(ids, annotation);
  const values = [{ chi: 1.2345678901234567 }, { chi: null }, { chi: -0.0000000000000023 }];
  const input = { ...table({ id: fixture.ref(ids), residue_id: fixture.ref(ids), annotation: fixture.ref(annotation), values }),
    family: 'backbone', missing: { annotation: [0] }, provenance: { source: 'fixture' } };
  let reads = 0;
  const expanded = await expandBundledFamilyColumns(JSON.parse(JSON.stringify(input)), async reference => {
    reads++; return (await verifyFamilyBundle(reference, fixture.bundle)).bundle;
  });
  assert.equal(reads, 1);
  assert.equal(expanded.encoding, FAMILY_COLUMNAR_ENCODING);
  assert.equal(expanded.family, 'backbone');
  assert.equal(expanded.build_id, 'fixture');
  assert.deepEqual(expanded.missing, input.missing);
  assert.deepEqual(expanded.provenance, input.provenance);
  assert.deepEqual(expanded.columns.annotation, annotation);
  assert.deepEqual(decodeFamilyRows(expanded), [
    { id: 'a', residue_id: 'a', values: values[0] },
    { id: 'b', residue_id: 'b', annotation: null, values: values[1] },
    { id: 'c', residue_id: 'c', annotation: annotation.dictionary[1], values: values[2] },
  ]);
});

test('Malformed inline dictionaries and references fail before any network read', async () => {
  const fixture = bundleFor(['a', 'b', 'c']);
  let reads = 0;
  for (const invalid of [null, 'text', 0, {}, { dictionary: [0], indices: [0, 0] },
    { dictionary: [0], indices: [0, 0, 1] }, { dictionary: [0], indices: [0, 0, -1] },
    { dictionary: [0], indices: [0, 0, 0.5] }, { dictionary: [0], indices: [0, 0, null] },
    { dictionary: [0], indices: [0, 0, 0], extra: true },
    { bundle: fixture.reference }, { ...fixture.ref(['a', 'b', 'c']), extra: true },
    { bundle: '../escape', column: hash([]) }, { bundle: fixture.reference, column: '../escape' }, [1]]) {
    await assert.rejects(expandBundledFamilyColumns(table({ id: fixture.ref(['a', 'b', 'c']), hidden: invalid }),
      async () => { reads++; return fixture.bundle; }), /Invalid bundled family|length mismatch/);
  }
  assert.equal(reads, 0);
});

test('All referenced columns and missing-field maps are validated before rows are returned', async () => {
  const ids = ['a', 'b', 'c'];
  for (const column of [[1], { dictionary: [1], indices: [0, 0] }, { dictionary: [1], indices: [0, 0, 1] }]) {
    const fixture = bundleFor(ids, column);
    await assert.rejects(expandBundledFamilyColumns(table({ id: fixture.ref(ids), value: fixture.ref(column) }),
      async () => fixture.bundle), /length mismatch|dictionary index/);
    if (!Array.isArray(column) && column.indices.includes(1)) {
      await assert.rejects(verifyFamilyBundle(fixture.reference, fixture.bundle), /dictionary index/);
    }
  }
  const fixture = bundleFor(ids);
  let reads = 0;
  for (const missing of [null, [], { absent: [0] }, { id: null }, { id: [3] }, { id: [-1] }, { id: [0.5] }]) {
    await assert.rejects(expandBundledFamilyColumns({ ...table({ id: fixture.ref(ids) }), missing },
      async () => { reads++; return fixture.bundle; }), /missing-field/);
  }
  assert.equal(reads, 0);
  await assert.rejects(expandBundledFamilyColumns(table({ id: fixture.ref(['missing']) }),
    async () => fixture.bundle), /Missing family bundle column/);
  for (const row_count of [-1, 0.5, null, Number.MAX_SAFE_INTEGER + 1]) {
    await assert.rejects(expandBundledFamilyColumns({ ...table({}), row_count }, async () => {}), /row count/);
  }
  const empty = await expandBundledFamilyColumns({ ...table({ id: [], value: { dictionary: [], indices: [] } }), row_count: 0 }, async () => {});
  assert.deepEqual(decodeFamilyRows(empty), []);
});

test('Distinct family bundle reads overlap at most four and failed reads remain retryable', async () => {
  const fixtures = Array.from({ length: 9 }, (_, i) => bundleFor([i, null, i]));
  const assets = new Map(fixtures.map(item => [item.reference, item.bundle]));
  const input = table(Object.fromEntries(fixtures.map((item, i) => [`field${i}`, item.ref([i, null, i])])));
  let active = 0, maximum = 0, reads = 0;
  const read = async reference => {
    reads++; active++; maximum = Math.max(maximum, active);
    await new Promise(resolve => setTimeout(resolve, 1));
    active--; return assets.get(reference);
  };
  const expanded = await expandBundledFamilyColumns(input, read);
  assert.equal(maximum, 4); assert.equal(active, 0); assert.equal(reads, 9);
  assert.deepEqual(expanded.columns.field8, [8, null, 8]);
  await expandBundledFamilyColumns(input, read);
  assert.equal(reads, 18);
  await assert.rejects(expandBundledFamilyColumns(input, async () => { throw new Error('HTTP failure'); }), /HTTP failure/);
  assert.deepEqual(await expandBundledFamilyColumns(input, read), expanded);
});

test('Family transport retains prototype-like names as safe own columns', async () => {
  const ids = ['a', 'b', 'c'];
  const fixture = bundleFor(ids);
  const input = table(Object.fromEntries([['__proto__', fixture.ref(ids)], ['constructor', [null, null, null]]]));
  const expanded = await expandBundledFamilyColumns(input, async () => fixture.bundle);
  assert.equal(Object.getPrototypeOf(expanded.columns), Object.prototype);
  assert(Object.hasOwn(expanded.columns, '__proto__'));
  assert.deepEqual(expanded.columns.__proto__, ids);
});
