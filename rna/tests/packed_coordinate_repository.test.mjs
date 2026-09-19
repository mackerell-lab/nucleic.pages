import test from 'node:test';
import assert from 'node:assert/strict';
import { gzipSync } from 'node:zlib';
import { RnaDataRepository } from '../core/repository.js';

const encoding = 'rna-coordinate-float64-1';
// Independent Node writer: the production codec uses browser DataView APIs.
function pack(values) {
  const raw = Buffer.alloc(values.length * 8), shuffled = Buffer.alloc(raw.length);
  values.forEach((value, i) => raw.writeDoubleLE(value, i * 8));
  for (let lane = 0; lane < 8; lane++) for (let i = 0; i < values.length; i++) shuffled[lane * values.length + i] = raw[i * 8 + lane];
  return { encoding: 'float64-le-shuffled-base64-1', count: values.length, data: shuffled.toString('base64') };
}
function fixture() {
  const table = (ids, pdbs, x, y, z) => ({ encoding, build_id: 'packed', row_count: ids.length,
    columns: { id: ids, pdb_id: pdbs, x: pack(x), y: pack(y), z: pack(z), status: ids.map(() => 'available') } });
  const first = table(['a1', 'b1'], ['1AAA', '1BBB'], [-0, Number.MIN_VALUE], [1.2345678901234567, -2.5], [Number.MAX_VALUE, 0]);
  first.columns.optional = [null, null]; first.missing = { optional: [1] };
  const second = table(['c1'], ['1CCC'], [0.25], [0.5], [0.75]);
  const manifest = { schema_version: 'rna-explorer-1', molecule_type: 'RNA', build_id: 'packed', survey: {
    coordinates: { groups: { G: { row_count: 3, partitions: [
      { path: 'ab.json.gz', encoding, row_count: 2, entry_ids: ['1AAA', '1BBB'] },
      { path: 'c.json.gz', encoding, row_count: 1, entry_ids: ['1CCC'] },
    ] } } } } };
  const assets = new Map([['manifest.json', manifest], ['ab.json.gz', first], ['c.json.gz', second]]), calls = [];
  const repository = new RnaDataRepository({ manifestUrl: 'https://rna.test/manifest.json', fetchImpl: async (url, options) => {
    const name = new URL(url).pathname.slice(1); calls.push({ name, signal: options?.signal });
    if (!assets.has(name)) return new Response('missing', { status: 404 });
    const json = JSON.stringify(assets.get(name));
    return new Response(name.endsWith('.gz') ? gzipSync(json) : json);
  } });
  return { repository, assets, calls, first, manifest };
}
async function collect(repository, options) {
  const chunks = [];
  for await (const chunk of repository.iterateSurveyCoordinates('G', options)) chunks.push(chunk);
  return chunks;
}

test('Packed coordinate streaming preserves exact numbers, absence and partition pruning', async () => {
  const { repository, calls } = fixture();
  const [first] = await collect(repository, { entryIds: ['1aaa', '1BBB'] });
  assert.deepEqual(first.rows.map(row => row.id), ['a1', 'b1']);
  assert(Object.is(first.rows[0].x, -0));
  assert(Object.is(first.rows[1].x, Number.MIN_VALUE));
  assert(Object.is(first.rows[0].y, 1.2345678901234567));
  assert(Object.is(first.rows[0].z, Number.MAX_VALUE));
  assert.equal(first.rows[0].optional, null);
  assert(!Object.hasOwn(first.rows[1], 'optional'));
  assert(Object.isFrozen(first.rows[0]));
  assert(!Object.hasOwn(first, 'columns'));
  assert(!Object.hasOwn(first, 'missing'));
  assert.deepEqual(calls.map(call => call.name), ['manifest.json', 'ab.json.gz']);
  assert.equal(repository.coordinateRequests.size, 0);
  assert.deepEqual([...repository.promises.keys()], ['manifest']);
  const before = calls.length;
  assert.deepEqual(await collect(repository, { entryIds: [] }), []);
  assert.equal(calls.length, before);
  const [again] = await collect(repository, { entryIds: ['1AAA'] });
  assert.equal(again.rows.length, 1);
  assert.equal(again.source_row_count, 2);
  assert.equal(calls.length, before + 1, 'Completed coordinate partitions must not remain cached');
});

test('Invalid packed data in an unselected row rejects before delivery and can retry', async () => {
  const { repository, assets, first } = fixture();
  const saved = first.columns.x;
  first.columns.x = pack([-0, Infinity]);
  await assert.rejects(collect(repository, { entryIds: ['1AAA'] }), /finite/i);
  assert.equal(repository.coordinateRequests.size, 0);
  first.columns.x = saved;
  assets.delete('ab.json.gz');
  await assert.rejects(collect(repository, { entryIds: ['1AAA'] }), /404/);
  assets.set('ab.json.gz', first);
  assert.equal((await collect(repository, { entryIds: ['1AAA'] }))[0].rows[0].id, 'a1');
});

test('Packed partitions require matching declared encoding, build, count and size limit', async () => {
  for (const mutate of [
    f => { f.first.build_id = 'another'; },
    f => { delete f.first.build_id; },
    f => { f.manifest.survey.coordinates.groups.G.partitions[0].encoding = 'rna-coordinate-columnar-1'; },
    f => { f.first.encoding = 'rna-coordinate-columnar-1'; f.first.columns.x = [0, 1]; f.first.columns.y = [1, 2]; f.first.columns.z = [2, 3]; },
    f => { f.manifest.survey.coordinates.groups.G.partitions[0].row_count = 1; },
    f => { f.first.row_count = 10001; },
  ]) {
    const f = fixture(); mutate(f);
    await assert.rejects(collect(f.repository), /Cross-build|encoding|row count|10000/i);
  }
});

test('Concurrent packed consumers share only in-flight requests; abort remains isolated', async () => {
  const { repository, calls } = fixture();
  let release, started;
  const gate = new Promise(resolve => { release = resolve; });
  const entered = new Promise(resolve => { started = resolve; });
  const fetch = repository.fetchImpl;
  repository.fetchImpl = async (...args) => {
    if (String(args[0]).endsWith('ab.json.gz')) { started(); await gate; }
    return fetch(...args);
  };
  const a = collect(repository, { entryIds: ['1AAA'] }), b = collect(repository, { entryIds: ['1BBB'] });
  await entered;
  assert.equal(repository.coordinateRequests.size, 1);
  release();
  const [left, right] = await Promise.all([a, b]);
  assert.equal(left[0].rows[0].id, 'a1'); assert.equal(right[0].rows[0].id, 'b1');
  assert.equal(calls.filter(call => call.name === 'ab.json.gz').length, 1);
  assert.equal(repository.coordinateRequests.size, 0);

  const controller = new AbortController();
  repository.fetchImpl = async (...args) => { const response = await fetch(...args);
    if (args[1]?.signal) controller.abort(new Error('superseded packed selection')); return response; };
  await assert.rejects(collect(repository, { signal: controller.signal }), /superseded packed selection/);
  assert.equal((await collect(repository, { entryIds: ['1CCC'] }))[0].rows[0].id, 'c1');
  assert.equal(repository.coordinateRequests.size, 0);
});
