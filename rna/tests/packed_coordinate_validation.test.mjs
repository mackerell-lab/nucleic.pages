import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import path from 'node:path';
import os from 'node:os';
import {gzipSync} from 'node:zlib';
import {sha256} from '../offline/output_scope.mjs';
import {validateRelease} from '../offline/assets.mjs';
import {encodeCoordinateRows} from '../core/survey-codec.js';
import {PACKED_COORDINATE_ENCODING, encodePackedCoordinates} from '../core/packed-coordinate-codec.js';

async function fixture(t) {
  const root = await fs.mkdtemp(path.join(os.tmpdir(), 'rna-packed-coordinates-'));
  t.after(() => fs.rm(root, {recursive: true, force: true}));
  const write = async (relative, value) => {
    const raw = Buffer.from(JSON.stringify(value)), bytes = gzipSync(raw), file = path.join(root, relative);
    await fs.mkdir(path.dirname(file), {recursive: true}); await fs.writeFile(file, bytes);
    return {path: relative, sha256: sha256(bytes), bytes: bytes.length, uncompressed_bytes: raw.length};
  };
  const rows = [0, 1].map(index => ({id: `c${index}`, pdb_id: 'TEST', entity_id: '1', model_id: '1',
    residue_id: 'r1', residue_ids: ['r1'], anchor_residue_id: 'r1', target_residue_id: 'r1',
    endpoint_entities: [{pdb_id: 'TEST', entity_id: '1'}], is_terminal_any: false,
    atom: index ? 'N2' : 'N1', x: 1.2345678901234567 + index, y: index ? 1e-12 : -5.6789,
    z: index ? -1e5 : 0, status: 'available', ...(index ? {} : {auth_seq_id: '42'})}));
  const packed = structuredClone(encodePackedCoordinates(encodeCoordinateRows(rows, 'packed-coordinate-release')));
  const manifest = {build_id: 'packed-coordinate-release', partial: false,
    counts: {entries: 1, entities: 1, residues: 1}, source: {candidate_count: 1},
    metadata: await write('metadata.json.gz', {entries: [{pdb_id: 'TEST', selected_model_id: '1'}], entities: [{pdb_id: 'TEST', entity_id: '1'}]}),
    families: [{...await write('families/backbone.json.gz', [{id: 'r1', pdb_id: 'TEST', entity_id: '1',
      is_terminal_any: false, values: {}, statuses: {}}]), id: 'backbone', level: 'residue', parameters: [], row_count: 1}],
    relations: {}, survey: {scalars: {terms: {}}, coordinates: {groups: {base: {row_count: 2, partitions: []}}}},
    provenance: {decisions: await write('decisions.json.gz', [{pdb_id: 'TEST', accepted: true}])}};
  const saveManifest = () => fs.writeFile(path.join(root, 'manifest.json'), JSON.stringify(manifest));
  const save = async (payload = packed) => {
    const descriptor = {...await write('survey/coordinates/base/00000.json.gz', payload),
      row_count: Array.isArray(payload) ? payload.length : payload.row_count, entry_ids: ['TEST'],
      ...(payload.encoding ? {encoding: payload.encoding} : {})};
    manifest.survey.coordinates.groups.base.partitions = [descriptor];
    manifest.survey.coordinates.groups.base.row_count = descriptor.row_count;
    await saveManifest();
  };
  await save();
  return {root, rows, packed, manifest, save, saveManifest,
    validate: () => validateRelease(path.join(root, 'manifest.json'))};
}

test('release validator expands exact packed coordinates and optional missing identity fields', async t => {
  const f = await fixture(t);
  assert.equal(f.packed.encoding, PACKED_COORDINATE_ENCODING);
  assert.equal(typeof f.packed.columns.x.data, 'string');
  assert.deepEqual(f.packed.missing.auth_seq_id, [1]);
  assert.equal((await f.validate()).ok, true);
});

test('packed coordinate release rejects malformed base64 descriptors and counts', async t => {
  const f = await fixture(t), original = structuredClone(f.packed.columns.x);
  const invalid = [
    {...original, data: '!invalid!'},
    {...original, data: Buffer.alloc(16).toString('base64').slice(0, -3) + 'B=='},
    {...original, data: Buffer.alloc(8).toString('base64')},
    {...original, data: undefined},
    {...original, count: 1},
    {...original, count: -1},
    {...original, encoding: 'other-numeric-format'},
  ];
  for (const descriptor of invalid) {
    f.packed.columns.x = descriptor; await f.save();
    await assert.rejects(f.validate(), /coordinate|base64|float|count|encoding/i);
  }
});

test('packed coordinate release rejects nonfinite binary values in every row', async t => {
  const f = await fixture(t);
  for (const invalid of [NaN, Infinity, -Infinity]) {
    // Encode independently to put the invalid value in the second row; checking
    // only the first row or the displayed subset would miss this corruption.
    const raw = Buffer.alloc(16), shuffled = Buffer.alloc(16);
    raw.writeDoubleLE(1.25, 0); raw.writeDoubleLE(invalid, 8);
    for (let byte = 0; byte < 8; byte++) for (let row = 0; row < 2; row++) shuffled[byte * 2 + row] = raw[row * 8 + byte];
    f.packed.columns.z.data = shuffled.toString('base64'); await f.save();
    await assert.rejects(f.validate(), /finite|coordinate/i);
  }
});

test('packed coordinates retain scientific checks on fallback missing and null numeric values', async t => {
  for (const missing of [false, true]) {
    const f = await fixture(t);
    if (missing) delete f.rows[1].x; else f.rows[1].x = null;
    const packed = encodePackedCoordinates(encodeCoordinateRows(f.rows, f.manifest.build_id));
    assert(Array.isArray(packed.columns.x), 'Missing or null numeric data must retain its original array');
    await f.save(packed);
    const result = await f.validate();
    assert.equal(result.ok, false);
    assert.ok(result.errors.includes('Coordinate value/status: c1'));
  }
});

test('rehashed packed coordinates retain endpoint and status scientific checks', async t => {
  for (const [field, value, expected] of [
    ['anchor_residue_id', 'missing-residue', 'Coordinate residue foreign key: c1'],
    ['status', 'missing_atoms', 'Coordinate value/status: c1'],
    ['endpoint_entities', [], 'Endpoint ownership: coordinate/c1'],
  ]) {
    const f = await fixture(t); f.rows[1][field] = value;
    await f.save(encodePackedCoordinates(encodeCoordinateRows(f.rows, f.manifest.build_id)));
    const result = await f.validate();
    assert.equal(result.ok, false); assert.ok(result.errors.includes(expected), JSON.stringify(result.errors));
  }
});

test('packed coordinate missing-field indices must remain valid', async t => {
  const f = await fixture(t); f.packed.missing.auth_seq_id = [2]; await f.save();
  await assert.rejects(f.validate(), /missing-field/i);
});

test('full packed release requires matching build and encoding; partial legacy data still validates', async t => {
  const f = await fixture(t);
  f.packed.build_id = 'another-release'; await f.save();
  await assert.rejects(f.validate(), /build ID mismatch/);
  delete f.packed.build_id; await f.save();
  await assert.rejects(f.validate(), /build ID mismatch/);
  f.manifest.partial = true; await f.saveManifest();
  assert.equal((await f.validate()).ok, true);
  f.manifest.survey.coordinates.groups.base.partitions[0].encoding = 'rna-coordinate-columnar-1'; await f.saveManifest();
  await assert.rejects(f.validate(), /encoding mismatch/);
  await f.save(f.rows);
  assert.equal((await f.validate()).ok, true);
  await f.save(encodeCoordinateRows(f.rows));
  assert.equal((await f.validate()).ok, true);
});

test('packed coordinate partitions keep the ten-thousand-row publication limit', async t => {
  const f = await fixture(t);
  const rows = Array.from({length: 10001}, (_, index) => ({...f.rows[0], id: `c${index}`}));
  await f.save(encodePackedCoordinates(encodeCoordinateRows(rows, f.manifest.build_id)));
  const result = await f.validate();
  assert.equal(result.ok, false); assert.ok(result.errors.includes('Coordinate partition row count'));
});
