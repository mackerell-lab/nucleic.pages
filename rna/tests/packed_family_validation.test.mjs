import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import path from 'node:path';
import os from 'node:os';
import {gzipSync} from 'node:zlib';
import {sha256} from '../offline/output_scope.mjs';
import {validateRelease} from '../offline/assets.mjs';
import {encodeFamilyRows} from '../core/survey-codec.js';
import {BUNDLED_FAMILY_ENCODING, FAMILY_BUNDLE_ENCODING} from '../core/bundled-family-codec.js';
import {PACKED_FAMILY_ENCODING, encodePackedFamily} from '../core/packed-family-codec.js';
import {encodeFloat64Column} from '../core/packed-coordinate-codec.js';

async function fixture(t) {
  const root = await fs.mkdtemp(path.join(os.tmpdir(), 'rna-packed-family-validation-'));
  t.after(() => fs.rm(root, {recursive: true, force: true}));
  const write = async (relative, value) => {
    const raw = Buffer.from(JSON.stringify(value)), bytes = gzipSync(raw), file = path.join(root, relative);
    await fs.mkdir(path.dirname(file), {recursive: true}); await fs.writeFile(file, bytes);
    return {path: relative, sha256: sha256(bytes), bytes: bytes.length, uncompressed_bytes: raw.length};
  };
  const rows = [{id: 'r1', pdb_id: 'TEST', entity_id: '1', is_terminal_any: false,
    values: {chi: 1.2345678901234567}, statuses: {chi: 'available'}},
  {id: 'r2', pdb_id: 'TEST', entity_id: '1', is_terminal_any: false, values: {chi: null}, statuses: {chi: 'missing_atoms'}}];
  const bundled = encodeFamilyRows(rows, 'packed-family-validation');
  bundled.encoding = BUNDLED_FAMILY_ENCODING;
  // Keep actual measurements inline so the new transport is exercised; identity
  // still passes through an authenticated content-addressed family bundle.
  bundled.columns.values = rows.map(row => row.values);
  const ids = bundled.columns.id, column = sha256(JSON.stringify(ids));
  const bundle = {encoding: FAMILY_BUNDLE_ENCODING, columns: {[column]: ids}}, reference = sha256(JSON.stringify(bundle));
  const bundleDescriptor = {...await write(`families/bundles/${reference}.json.gz`, bundle), content_sha256: reference};
  bundled.columns.id = {bundle: reference, column};
  const packed = structuredClone(encodePackedFamily(bundled));
  const manifest = {build_id: 'packed-family-validation', partial: false,
    counts: {entries: 1, entities: 1, residues: 2}, source: {candidate_count: 1},
    metadata: await write('metadata.json.gz', {entries: [{pdb_id: 'TEST'}], entities: [{pdb_id: 'TEST', entity_id: '1'}]}),
    families: [], family_bundles: {[reference]: bundleDescriptor}, relations: {},
    survey: {scalars: {terms: {}}, coordinates: {groups: {}}},
    provenance: {decisions: await write('decisions.json.gz', [{pdb_id: 'TEST', accepted: true}])}};
  const saveManifest = () => fs.writeFile(path.join(root, 'manifest.json'), JSON.stringify(manifest));
  const save = async (payload = packed) => {
    manifest.families = [{...await write('families/backbone.json.gz', payload), row_count: 2,
      id: 'backbone', level: 'residue', parameters: [{id: 'chi'}], encoding: payload.encoding}];
    await saveManifest();
  };
  await save();
  return {root, bundled, packed, rows, manifest, save, saveManifest,
    validate: () => validateRelease(path.join(root, 'manifest.json'))};
}

test('packed family release validates finite measurements nulls statuses and identity bundles', async t => {
  const f = await fixture(t);
  assert.equal(f.packed.encoding, PACKED_FAMILY_ENCODING);
  assert.equal(typeof f.packed.columns.values.parameters.chi.values.data, 'string');
  assert.deepEqual(f.packed.columns.values.parameters.chi.nulls, [1]);
  assert.equal((await f.validate()).ok, true);
});

test('packed family release rejects malformed numeric descriptors and nonfinite values', async t => {
  const f = await fixture(t), parameter = f.packed.columns.values.parameters.chi;
  const original = structuredClone(parameter.values);
  for (const invalid of [{...original, data: '!invalid'}, {...original, count: 1},
    {...original, encoding: 'wrong-float-format'}, {...original, data: Buffer.alloc(8).toString('base64')}]) {
    parameter.values = invalid; await f.save();
    await assert.rejects(f.validate(), /float|descriptor|base64|count|coordinate/i);
  }
  // Independent IEEE754 construction places infinity in the masked second row.
  // The mask must not hide invalid binary input.
  const raw = Buffer.alloc(16), shuffled = Buffer.alloc(16);
  raw.writeDoubleLE(1.25, 0); raw.writeDoubleLE(Infinity, 8);
  for (let byte = 0; byte < 8; byte++) for (let row = 0; row < 2; row++) shuffled[byte * 2 + row] = raw[row * 8 + byte];
  parameter.values = {...original, data: shuffled.toString('base64')}; await f.save();
  await assert.rejects(f.validate(), /finite/i);
});

test('packed family release validates all masks and hidden numeric slots', async t => {
  const f = await fixture(t), initial = structuredClone(f.packed.columns.values.parameters.chi);
  for (const invalid of [
    {...initial, nulls: [2]}, {...initial, nulls: [-1]}, {...initial, nulls: [1, 1]},
    {...initial, nulls: [1, 0]}, {...initial, missing: [1]}, {...initial, missing: '1'},
    {...initial, values: encodeFloat64Column([1.2345678901234567, 99])},
  ]) {
    f.packed.columns.values.parameters.chi = invalid; await f.save();
    await assert.rejects(f.validate(), /mask|slot/i);
  }
});

test('packed family release rejects rehashed mismatched scientific statuses and absent values', async t => {
  const f = await fixture(t);
  f.packed.columns.statuses = [{chi: 'available'}, {chi: 'available'}]; await f.save();
  let result = await f.validate();
  assert.equal(result.ok, false); assert.ok(result.errors.includes('Value/status: backbone/chi'));
  f.packed.columns.statuses = [{chi: 'available'}, {chi: 'missing_atoms'}];
  f.packed.columns.values.parameters.chi.nulls = [];
  f.packed.columns.values.parameters.chi.missing = [1]; await f.save();
  result = await f.validate();
  assert.equal(result.ok, false); assert.ok(result.errors.includes('Value/status: backbone/chi'));
});

test('packed family requires exact build and descriptor encoding in both directions', async t => {
  const f = await fixture(t);
  f.packed.build_id = 'wrong-build'; await f.save();
  await assert.rejects(f.validate(), /build ID mismatch/);
  delete f.packed.build_id; await f.save();
  await assert.rejects(f.validate(), /build ID mismatch/);
  f.packed.build_id = f.manifest.build_id; await f.save();
  f.manifest.families[0].encoding = BUNDLED_FAMILY_ENCODING; await f.saveManifest();
  await assert.rejects(f.validate(), /encoding mismatch/);
  await f.save(f.bundled);
  f.manifest.families[0].encoding = PACKED_FAMILY_ENCODING; await f.saveManifest();
  await assert.rejects(f.validate(), /encoding mismatch/);
});

test('packed family still enforces registered authenticated identity bundles', async t => {
  const f = await fixture(t), reference = f.packed.columns.id.bundle;
  const descriptor = f.manifest.family_bundles[reference];
  descriptor.sha256 = '0'.repeat(64); await f.saveManifest();
  await assert.rejects(f.validate(), /Asset hash mismatch/);
  delete f.manifest.family_bundles[reference]; await f.saveManifest();
  await assert.rejects(f.validate(), /Unregistered family bundle/);
});
