import test from 'node:test';
import assert from 'node:assert/strict';
import {createHash} from 'node:crypto';
import {gzipSync} from 'node:zlib';
import {RnaDataRepository} from '../core/repository.js';

const hash = value => createHash('sha256').update(JSON.stringify(value)).digest('hex');
// Independent IEEE writer and byte-lane loop; never invoke the production encoder.
function numeric(values) {
  const raw = Buffer.alloc(values.length * 8), lanes = Buffer.alloc(raw.length);
  values.forEach((value, index) => raw.writeDoubleLE(value, index * 8));
  for (let index = 0; index < values.length; index++) for (let byte = 0; byte < 8; byte++) lanes[byte * values.length + index] = raw[index * 8 + byte];
  return {encoding: 'float64-le-shuffled-base64-1', count: values.length, data: lanes.toString('base64')};
}
function setup() {
  const root = 'https://rna.test/release/', ids = ['a', 'b', 'c', 'd'];
  const column = hash(ids), bundle = {encoding: 'rna-family-column-bundle-1', columns: {[column]: ids}}, reference = hash(bundle);
  const descriptor = {path: `families/bundles/${reference}.json.gz`, content_sha256: reference};
  const packed = values => ({encoding: 'rna-family-float64-1', build_id: 'test', row_count: 4,
    columns: {id: {bundle: reference, column}, values: {encoding: 'rna-family-values-float64-1', parameters: {
      chi: {values: numeric(values), nulls: [2], missing: [3]},
    }}, status: ['ok', 'ok', 'missing', 'missing'], optional: [null, null, null, null]}, missing: {optional: [1]}});
  const manifest = {schema_version: 'rna-explorer-1', molecule_type: 'RNA', build_id: 'test',
    families: ['one', 'two'].map(id => ({id, path: `families/${id}.json.gz`, row_count: 4, encoding: 'rna-family-float64-1'})), family_bundles: {[reference]: descriptor}};
  const assets = new Map([['manifest.json', manifest], [descriptor.path, bundle],
    ['families/one.json.gz', packed([-0, Number.MIN_VALUE, 0, 0])], ['families/two.json.gz', packed([Number.MAX_VALUE, 1.2345678901234567, 0, 0])]]);
  const calls = [];
  const repository = new RnaDataRepository({manifestUrl: root + 'manifest.json', maxCachedFamilies: 1,
    fetchImpl: async url => { calls.push(url); const value = assets.get(url.slice(root.length));
      return value ? new Response(url.endsWith('.gz') ? gzipSync(JSON.stringify(value)) : JSON.stringify(value)) : new Response('unavailable', {status: 503}); }});
  return {repository, assets, manifest, calls, descriptor};
}

test('Packed families preserve exact values, field absence, lazy bundles and bounded caches', async () => {
  const {repository, calls, descriptor} = setup();
  await repository.loadManifest();
  assert.equal(calls.length, 1, 'Manifest must not preload family data');
  const [one, two] = await Promise.all([repository.loadFamily('one'), repository.loadFamily('two')]);
  assert(Object.is(one.rows[0].values.chi, -0));
  assert(Object.is(one.rows[1].values.chi, Number.MIN_VALUE));
  assert.equal(one.rows[2].values.chi, null);
  assert(!Object.hasOwn(one.rows[3].values, 'chi'));
  assert.equal(one.rows[0].optional, null);
  assert(!Object.hasOwn(one.rows[1], 'optional'));
  assert(Object.is(two.rows[0].values.chi, Number.MAX_VALUE));
  assert(Object.is(two.rows[1].values.chi, 1.2345678901234567));
  assert(Object.isFrozen(one.rows[0].values));
  assert(!Object.hasOwn(one, 'columns'));
  assert(!Object.hasOwn(one, 'missing'));
  assert.equal(calls.filter(url => url.endsWith(descriptor.path)).length, 1);
  assert(repository.resolvedFamilies.size <= 1);
});

test('Packed family corruption fails before returning any rows and can be retried', async () => {
  const {repository, assets} = setup(), data = assets.get('families/one.json.gz');
  const saved = structuredClone(data.columns.values.parameters.chi);
  data.columns.values.parameters.chi.values = numeric([-0, Infinity, 0, 0]);
  await assert.rejects(repository.loadFamily('one'), /finite/i);
  data.columns.values.parameters.chi = saved;
  assets.delete('families/one.json.gz');
  await assert.rejects(repository.loadFamily('one'), /503/);
  assets.set('families/one.json.gz', data);
  assert(Object.is((await repository.loadFamily('one')).rows[0].values.chi, -0));
});

test('Packed family identity, encoding, counts and masks are enforced', async () => {
  for (const mutate of [
    ({assets}) => { delete assets.get('families/one.json.gz').build_id; },
    ({assets}) => { assets.get('families/one.json.gz').build_id = 'other'; },
    ({manifest}) => { manifest.families[0].encoding = 'rna-family-columnar-1'; },
    ({assets}) => { assets.get('families/one.json.gz').encoding = 'rna-family-columnar-1'; },
    ({manifest}) => { manifest.families[0].row_count = 3; },
    ({assets}) => { assets.get('families/one.json.gz').columns.values.parameters.chi.nulls = [2, 2]; },
    ({assets}) => { assets.get('families/one.json.gz').columns.values.parameters.chi.missing = [2]; },
  ]) {
    const fixture = setup(); mutate(fixture);
    await assert.rejects(fixture.repository.loadFamily('one'), /build|encoding|count|mask|index|overlap|indices/i);
  }
});
