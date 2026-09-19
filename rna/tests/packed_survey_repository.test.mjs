import test from 'node:test';
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { gzipSync } from 'node:zlib';
import { RnaDataRepository } from '../core/repository.js';

const hash = value => createHash('sha256').update(JSON.stringify(value)).digest('hex');
// Fixtures are independent of the production encoder.
function numeric(values) {
  const raw = Buffer.alloc(values.length * 8), lanes = Buffer.alloc(raw.length);
  values.forEach((value, index) => raw.writeDoubleLE(value, index * 8));
  for (let index = 0; index < values.length; index++) for (let byte = 0; byte < 8; byte++) lanes[byte * values.length + index] = raw[index * 8 + byte];
  return { encoding: 'float64-le-shuffled-base64-1', count: values.length, data: lanes.toString('base64') };
}
function setup(maxBundleCacheBytes = 32 * 1024 * 1024) {
  const root = 'https://rna.test/releases/test/', observations = ['residue:A', 'residue:B', 'pair:A:B'];
  const column = hash(observations), bundle = { encoding: 'rna-survey-column-bundle-1', columns: { [column]: observations } }, reference = hash(bundle);
  const packed = (term, values) => ({ encoding: 'rna-survey-float64-1', build_id: 'test', row_count: 3,
    columns: { id: { encoding: 'rna-survey-id-suffix-1', column: 'observation_id', suffix: `|survey|${term}` },
      observation_id: { bundle: reference, column }, term_id: { encoding: 'rna-survey-constant-1', value: term },
      value: { encoding: 'rna-survey-values-float64-1', values: numeric(values), nulls: [2] },
      status: ['ok', 'ok', 'missing_atom'], opening: [0, 5, null] } });
  const descriptor = { path: `survey/bundles/${reference}.json.gz`, content_sha256: reference };
  const manifest = { schema_version: 'rna-explorer-1', molecule_type: 'RNA', build_id: 'test', survey: {
    scalars: { terms: { one: { path: 'one.json.gz', row_count: 3, encoding: 'rna-survey-float64-1' },
      two: { path: 'two.json.gz', row_count: 3, encoding: 'rna-survey-float64-1' } } }, bundles: { [reference]: descriptor } } };
  const assets = new Map([['https://rna.test/current.json', { manifest: 'releases/test/manifest.json', build_id: 'test' }],
    [root + 'manifest.json', manifest], [root + descriptor.path, bundle],
    [root + 'one.json.gz', packed('one', [-0, Number.MIN_VALUE, 0])],
    [root + 'two.json.gz', packed('two', [Number.MAX_VALUE, 1.2345678901234567, 0])]]);
  const calls = [];
  const repository = new RnaDataRepository({ manifestUrl: 'https://rna.test/current.json', maxBundleCacheBytes,
    fetchImpl: async url => { calls.push(url); const value = assets.get(url);
      return value ? new Response(url.endsWith('.gz') ? gzipSync(JSON.stringify(value)) : JSON.stringify(value)) : new Response('unavailable', { status: 503 }); } });
  return { repository, assets, manifest, calls, descriptor, root, reference, bundle, column };
}

test('Packed Survey registry is lazy and preserves exact scoped projections and IDs', async () => {
  const { repository, calls, descriptor } = setup();
  const registry = await repository.loadSurveyScalars();
  assert.equal(registry.partitions.length, 2);
  assert.equal(calls.length, 2, 'Registry must only load pointer and manifest');
  const one = await repository.loadSurveyScalars('one');
  assert.equal(calls.length, 4);
  assert(!calls.some(url => url.endsWith('two.json.gz')), 'No unrelated scalar fetch');
  assert.deepEqual(one.rows.map(row => row.id), ['residue:A|survey|one', 'residue:B|survey|one', 'pair:A:B|survey|one']);
  assert(Object.is(one.rows[0].value, -0));
  assert(Object.is(one.rows[1].value, Number.MIN_VALUE));
  assert.equal(one.rows[2].value, null);
  assert(Object.isFrozen(one.rows[0]));
  assert(!Object.hasOwn(one, 'columns'));
  assert.equal(await repository.loadSurveyScalars('one'), one);
  const two = await repository.loadSurveyScalars('two', { fields: ['value'] });
  assert.deepEqual(two.rows, [{ value: Number.MAX_VALUE }, { value: 1.2345678901234567 }, { value: null }]);
  assert.equal(calls.filter(url => url.endsWith(descriptor.path)).length, 1);
  repository.releaseSurvey('scalars', 'one');
  assert.deepEqual((await repository.loadSurveyScalars('one')).rows, one.rows);
  assert.equal(calls.filter(url => url.endsWith(descriptor.path)).length, 1);
});

test('Packed Survey projected loads authenticate all bundles and retry missing or corrupt sources', async () => {
  const fixture = setup(0), { repository, assets, root, descriptor, bundle, column } = fixture;
  const bundleUrl = root + descriptor.path;
  assets.delete(bundleUrl);
  await assert.rejects(repository.loadSurveyScalars('one', { fields: ['value'] }), /503/);
  assets.set(bundleUrl, { ...bundle, columns: { [column]: ['bad', 'bad', 'bad'] } });
  await assert.rejects(repository.loadSurveyScalars('one', { fields: ['value'] }), /hash mismatch/);
  assets.set(bundleUrl, bundle);
  assert(Object.is((await repository.loadSurveyScalars('one', { fields: ['value'] })).rows[0].value, -0));
  assert.equal(repository.bundleCacheBytes, 0);
  assert.equal(repository.bundleRequests.size, 0);
  const original = assets.get(root + 'two.json.gz');
  assets.delete(root + 'two.json.gz');
  await assert.rejects(repository.loadSurveyScalars('two'), /503/);
  assets.set(root + 'two.json.gz', original);
  assert.equal((await repository.loadSurveyScalars('two')).rows.length, 3);
});

test('Packed Survey rejects build, encoding and count drift even for projected loads', async () => {
  for (const mutate of [
    f => { delete f.assets.get(f.root + 'one.json.gz').build_id; },
    f => { f.assets.get(f.root + 'one.json.gz').build_id = 'another'; },
    f => { f.manifest.survey.scalars.terms.one.encoding = 'rna-survey-bundled-columns-1'; },
    f => { f.assets.get(f.root + 'one.json.gz').encoding = 'rna-survey-bundled-columns-1'; },
    f => { f.manifest.survey.scalars.terms.one.row_count = 2; },
    f => { f.assets.get(f.root + 'one.json.gz').missing = { value: [0] }; },
    f => { f.assets.get(f.root + 'one.json.gz').columns.id.column = 'id'; },
    f => { f.assets.get(f.root + 'one.json.gz').columns.value.values = numeric([-0, Infinity, 0]); },
  ]) {
    const fixture = setup(); mutate(fixture);
    await assert.rejects(fixture.repository.loadSurveyScalars('one', { fields: ['status'] }), /build|encoding|count|missing|descriptor|finite/i);
  }
});

test('Concurrent packed terms share authentication while bundle cache cost remains bounded', async () => {
  const { repository, calls, descriptor, bundle } = setup();
  const [one, two] = await Promise.all([repository.loadSurveyScalars('one'), repository.loadSurveyScalars('two')]);
  assert.equal(one.rows[0].id, 'residue:A|survey|one');
  assert.equal(two.rows[0].id, 'residue:A|survey|two');
  assert.equal(calls.filter(url => url.endsWith(descriptor.path)).length, 1);
  assert.equal(repository.bundleCacheBytes, Buffer.byteLength(JSON.stringify(bundle)));
  assert(repository.bundleCacheBytes <= repository.maxBundleCacheBytes);
});
