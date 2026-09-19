import test from 'node:test';
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { gzipSync } from 'node:zlib';
import { RnaDataRepository } from '../core/repository.js';
import { BUNDLED_SURVEY_ENCODING, SURVEY_BUNDLE_ENCODING } from '../core/bundled-survey-codec.js';

const hash = data => createHash('sha256').update(JSON.stringify(data)).digest('hex');
function bundle(values) {
  const column = hash(values);
  const data = { encoding: SURVEY_BUNDLE_ENCODING, columns: { [column]: values } };
  return { data, column, reference: hash(data), bytes: Buffer.byteLength(JSON.stringify(data)) };
}
function setup(maxBundleCacheBytes) {
  const a = bundle(['a', 'b']), b = bundle(['c', 'd']);
  const term = (item, values) => ({ encoding: BUNDLED_SURVEY_ENCODING, build_id: 'test', row_count: 2,
    columns: { id: { bundle: item.reference, column: item.column }, value: values } });
  const root = 'https://rna.test/releases/test/';
  const assets = new Map([
    ['https://rna.test/current.json', { manifest: 'releases/test/manifest.json', build_id: 'test' }],
    [root + 'manifest.json', { build_id: 'test', molecule_type: 'RNA', schema_version: 'rna-explorer-1', survey: {
      scalars: { terms: { one: { path: 'one.json.gz', row_count: 2 }, two: { path: 'two.json.gz', row_count: 2 } } } } }],
    [root + 'one.json.gz', term(a, [1.125, null])], [root + 'two.json.gz', term(a, [2.5, 3.5])],
    [root + `survey/bundles/${a.reference}.json.gz`, a.data],
    [root + `survey/bundles/${b.reference}.json.gz`, b.data],
  ]);
  const calls = [];
  const repository = new RnaDataRepository({ manifestUrl: 'https://rna.test/current.json', maxBundleCacheBytes,
    fetchImpl: async url => {
      calls.push(url);
      if (!assets.has(url)) return new Response('missing', { status: 404 });
      const body = JSON.stringify(assets.get(url));
      return new Response(url.endsWith('.gz') ? gzipSync(body) : body);
    } });
  return { repository, assets, calls, root, a, b };
}

test('Gzip bundles follow release pointers and concurrent terms share one verified request', async () => {
  const { repository, calls, a } = setup();
  const [one, two] = await Promise.all([repository.loadSurveyScalars('one'), repository.loadSurveyScalars('two', { fields: ['value'] })]);
  assert.deepEqual(one.rows, [{ id: 'a', value: 1.125 }, { id: 'b', value: null }]);
  assert.deepEqual(two.rows, [{ value: 2.5 }, { value: 3.5 }]);
  assert.equal(calls.filter(url => url.includes('/survey/bundles/')).length, 1);
  assert.equal(repository.bundleCacheBytes, a.bytes);
  assert.equal(repository.bundleRequests.size, 0);
  assert(Object.isFrozen(repository.bundleCache.get(a.reference).bundle));
  assert(!Object.hasOwn(one, 'columns'));
  repository.releaseSurvey('scalars', 'one');
  await repository.loadSurveyScalars('one');
  assert.equal(calls.filter(url => url.includes('/survey/bundles/')).length, 1);
});

test('Bundle eviction respects serialized cost and leaves retained row results valid', async () => {
  const { a } = setup();
  const { repository, b, calls } = setup(a.bytes);
  const rows = (await repository.loadSurveyScalars('one')).rows;
  await repository.loadSurveyBundle(b.reference);
  assert.equal(repository.bundleCache.size, 1);
  assert(repository.bundleCacheBytes <= a.bytes);
  assert(!repository.bundleCache.has(a.reference));
  assert.deepEqual(rows, [{ id: 'a', value: 1.125 }, { id: 'b', value: null }]);
  await repository.loadSurveyBundle(a.reference);
  assert.equal(calls.filter(url => url.includes('/survey/bundles/')).length, 3);
  const disabled = setup(0);
  await disabled.repository.loadSurveyBundle(a.reference);
  assert.equal(disabled.repository.bundleCacheBytes, 0);
  assert.equal(disabled.repository.bundleCache.size, 0);
});

test('Missing and corrupt bundles fail projected loads and remain retryable', async () => {
  const { repository, assets, root, a } = setup(0);
  const url = root + `survey/bundles/${a.reference}.json.gz`;
  assets.delete(url);
  await assert.rejects(repository.loadSurveyScalars('one', { fields: ['value'] }), /404/);
  assets.set(url, { ...a.data, columns: { [a.column]: ['wrong', 'data'] } });
  await assert.rejects(repository.loadSurveyScalars('one', { fields: ['value'] }), /hash mismatch/i);
  assets.set(url, a.data);
  assert.deepEqual((await repository.loadSurveyScalars('one', { fields: ['value'] })).rows, [{ value: 1.125 }, { value: null }]);
  assert.equal(repository.bundleRequests.size, 0);
  await assert.rejects(repository.loadSurveyBundle('../escape'), /content hash/);
});

test('Verified bundle hits promote recency and oversized bundles are not retained', async () => {
  const cost = bundle(['a', 'b']).bytes;
  const { repository, assets, root, a, b } = setup(2 * cost);
  const c = bundle(['e', 'f']);
  assets.set(root + `survey/bundles/${c.reference}.json.gz`, c.data);
  await repository.loadSurveyBundle(a.reference);
  await repository.loadSurveyBundle(b.reference);
  await repository.loadSurveyBundle(a.reference);
  await repository.loadSurveyBundle(c.reference);
  assert.deepEqual([...repository.bundleCache.keys()], [a.reference, c.reference]);
  const large = bundle(['long'.repeat(500), 'large']);
  assets.set(root + `survey/bundles/${large.reference}.json.gz`, large.data);
  const result = await repository.loadSurveyBundle(large.reference);
  assert.equal(result.columns[large.column][0].length, 2000);
  assert(!repository.bundleCache.has(large.reference));
  assert.deepEqual([...repository.bundleCache.keys()], [a.reference, c.reference]);
  assert.equal(repository.bundleCacheBytes, 2 * cost);
});

test('Consumers share the same pending verified bundle operation', async () => {
  const { repository, a, calls } = setup();
  let release;
  const gate = new Promise(resolve => { release = resolve; });
  const fetch = repository.fetchImpl;
  repository.fetchImpl = async (...args) => {
    if (String(args[0]).includes('/survey/bundles/')) await gate;
    return fetch(...args);
  };
  const first = repository.loadSurveyBundle(a.reference);
  const second = repository.loadSurveyBundle(a.reference);
  assert.equal(first, second);
  assert.equal(repository.bundleRequests.size, 1);
  release();
  assert.equal(await first, await second);
  assert.equal(repository.bundleRequests.size, 0);
  assert.equal(calls.filter(url => url.includes('/survey/bundles/')).length, 1);
});
