import test from 'node:test';
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { gzipSync } from 'node:zlib';
import { RnaDataRepository } from '../core/repository.js';

const hash = value => createHash('sha256').update(JSON.stringify(value)).digest('hex');
function bundle(encoding, column) {
  const columnHash = hash(column);
  const payload = { encoding, columns: { [columnHash]: column } };
  const reference = hash(payload), raw = JSON.stringify(payload), compressed = gzipSync(raw);
  return { payload, reference, columnHash, cost: Buffer.byteLength(raw),
    descriptor: { path: `families/bundles/${reference}.json.gz`, content_sha256: reference,
      sha256: createHash('sha256').update(compressed).digest('hex'), bytes: compressed.length, uncompressed_bytes: Buffer.byteLength(raw) } };
}
function setup(maxBundleCacheBytes) {
  const family = bundle('rna-family-column-bundle-1', { dictionary: ['a', 'b'], indices: [0, 1] });
  const survey = bundle('rna-survey-column-bundle-1', ['a', 'b']);
  const root = 'https://rna.test/releases/test/';
  const table = values => ({ encoding: 'rna-family-bundled-columns-1', build_id: 'test', row_count: 2,
    columns: { id: { bundle: family.reference, column: family.columnHash },
      values, comp_id: { dictionary: ['U'], indices: [0, 0] }, optional: [null, null] }, missing: { optional: [1] } });
  const descriptor = id => ({ id, path: `families/${id}.json.gz`, encoding: 'rna-family-bundled-columns-1', row_count: 2 });
  const manifest = { build_id: 'test', molecule_type: 'RNA', schema_version: 'rna-explorer-1',
    families: [descriptor('one'), descriptor('two')], family_bundles: { [family.reference]: family.descriptor } };
  const assets = new Map([
    [root + 'manifest.json', manifest],
    [root + 'families/one.json.gz', table([{ chi: 1.125 }, { chi: null }])],
    [root + 'families/two.json.gz', table([{ delta: 2.5 }, { delta: 3.5 }])],
    [root + family.descriptor.path, family.payload],
    [root + `survey/bundles/${survey.reference}.json.gz`, survey.payload],
  ]);
  const calls = [];
  const repository = new RnaDataRepository({ manifestUrl: root + 'manifest.json', maxCachedFamilies: 1, maxBundleCacheBytes,
    fetchImpl: async url => { calls.push(url); return assets.has(url)
      ? new Response(url.endsWith('.gz') ? gzipSync(JSON.stringify(assets.get(url))) : JSON.stringify(assets.get(url)))
      : new Response('missing', { status: 404 }); } });
  return { repository, family, survey, root, assets, calls, manifest };
}

test('Concurrent families reuse verified dictionaries without losing nulls or missing fields', async () => {
  const { repository, family, calls } = setup();
  const [one, two] = await Promise.all([repository.loadFamily('one'), repository.loadFamily('two')]);
  assert.deepEqual(one.rows, [{ id: 'a', values: { chi: 1.125 }, comp_id: 'U', optional: null },
    { id: 'b', values: { chi: null }, comp_id: 'U' }]);
  assert.equal(two.rows[1].values.delta, 3.5);
  assert.equal(calls.filter(url => url.includes('/families/bundles/')).length, 1);
  assert.equal(repository.bundleCacheBytes, family.cost);
  assert(Object.isFrozen(one.rows[0].values));
  assert(!Object.hasOwn(one, 'columns'));
  assert(!Object.hasOwn(one, 'missing'));
  await repository.loadFamily('one');
  assert.equal(calls.filter(url => url.includes('/families/bundles/')).length, 1);
  assert(repository.resolvedFamilies.size <= 1);
});

test('Family and Survey bundles share one bounded cache without invalidating returned rows', async () => {
  const cost = setup().family.cost;
  const { repository, family, survey, calls } = setup(cost);
  const one = await repository.loadFamily('one');
  await repository.loadSurveyBundle(survey.reference);
  assert.equal(repository.bundleCache.size, 1);
  assert.equal(repository.bundleCacheBytes, survey.cost);
  await repository.loadFamily('two');
  assert.equal(repository.bundleCache.size, 1);
  assert.equal(repository.bundleCacheBytes, family.cost);
  assert.equal(calls.filter(url => url.includes('/families/bundles/')).length, 2);
  assert.equal(one.rows[0].values.chi, 1.125);
  assert(Object.isFrozen(one.rows));
});

test('Failed family bundle requests retry after missing or corrupt transport', async () => {
  const { repository, family, assets, root } = setup(0);
  const url = root + family.descriptor.path;
  assets.delete(url);
  await assert.rejects(repository.loadFamily('one'), /404/);
  assets.set(url, { ...family.payload, columns: { [family.columnHash]: ['wrong', 'ids'] } });
  await assert.rejects(repository.loadFamily('one'), /hash mismatch/i);
  assets.set(url, family.payload);
  assert.equal((await repository.loadFamily('one')).rows[0].id, 'a');
  assert.equal(repository.bundleRequests.size, 0);
  assert.equal(repository.bundleCacheBytes, 0);
  assert.equal(repository.bundleCache.size, 0);
});

test('Family bundles reject unregistered, escaping, and cross-build transport', async () => {
  for (const mutate of [
    ({ manifest }) => { manifest.family_bundles = {}; },
    ({ family }) => { family.descriptor.path = '../escape.json.gz'; },
    ({ assets, root }) => { assets.get(root + 'families/one.json.gz').build_id = 'different'; },
    ({ manifest }) => { manifest.families[0].encoding = 'rna-family-columnar-1'; },
    ({ assets, root }) => { const table = assets.get(root + 'families/one.json.gz');
      table.encoding = 'rna-family-columnar-1'; table.columns.id = ['a', 'b']; },
  ]) {
    const fixture = setup(); mutate(fixture);
    await assert.rejects(fixture.repository.loadFamily('one'), /registry|registered|Cross-build|encoding/i);
  }
});
