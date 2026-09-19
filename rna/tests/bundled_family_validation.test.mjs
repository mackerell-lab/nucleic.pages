import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import os from 'node:os';
import path from 'node:path';
import {gzipSync} from 'node:zlib';
import {createHash} from 'node:crypto';
import {validateRelease} from '../offline/assets.mjs';
import {verifyReleaseInventory} from '../offline/verify_release_inventory.mjs';
import {encodeFamilyRows} from '../core/survey-codec.js';
import {BUNDLED_FAMILY_ENCODING, FAMILY_BUNDLE_ENCODING} from '../core/bundled-family-codec.js';

const hash = bytes => createHash('sha256').update(bytes).digest('hex');

async function fixture(t) {
  const directory = await fs.mkdtemp(path.join(os.tmpdir(), 'rna-family-validation-'));
  const root = path.join(directory, 'release');
  await fs.mkdir(root);
  t.after(() => fs.rm(directory, {recursive: true, force: true}));
  const write = async (relative, value) => {
    const file = path.join(root, relative), raw = Buffer.from(JSON.stringify(value)), bytes = gzipSync(raw);
    await fs.mkdir(path.dirname(file), {recursive: true});
    await fs.writeFile(file, bytes);
    return {path: relative, sha256: hash(bytes), bytes: bytes.length, uncompressed_bytes: raw.length};
  };
  const packed = encodeFamilyRows(['r1', 'r2'].map(id => ({id, pdb_id: 'TEST', entity_id: '1',
    is_terminal_any: false, values: {chi: 1.2345678901234567}, statuses: {chi: 'available'}})), 'family-validation');
  packed.encoding = BUNDLED_FAMILY_ENCODING;
  const manifest = {
    schema_version: 'rna-explorer-1', molecule_type: 'RNA', build_id: 'family-validation', partial: false,
    counts: {entries: 1, entities: 1, residues: 2}, source: {candidate_count: 1},
    metadata: await write('metadata.json.gz', {entries: [{pdb_id: 'TEST'}], entities: [{pdb_id: 'TEST', entity_id: '1'}]}),
    families: [], family_bundles: {}, relations: {},
    survey: {scalars: {terms: {}}, coordinates: {groups: {}}},
    provenance: {decisions: await write('decisions.json.gz', [{pdb_id: 'TEST', accepted: true}])},
  };
  const bundleColumn = async values => {
    const column = hash(JSON.stringify(values));
    const payload = {encoding: FAMILY_BUNDLE_ENCODING, columns: {[column]: values}};
    const reference = hash(JSON.stringify(payload));
    const descriptor = {...await write(`families/bundles/${reference}.json.gz`, payload), content_sha256: reference};
    manifest.family_bundles[reference] = descriptor;
    return {bundle: reference, column};
  };
  packed.columns.values = await bundleColumn({dictionary: [{chi: 1.2345678901234567}], indices: [0, 0]});
  packed.columns.pdb_id = await bundleColumn(['TEST', 'TEST']);
  const saveManifest = () => fs.writeFile(path.join(root, 'manifest.json'), JSON.stringify(manifest));
  const save = async () => {
    manifest.families = [{...await write('families/backbone.json.gz', packed),
      id: 'backbone', level: 'residue', parameters: [{id: 'chi'}], row_count: 2, encoding: BUNDLED_FAMILY_ENCODING}];
    await saveManifest();
  };
  await save();
  return {directory, root, manifest, packed, write, bundleColumn, save, saveManifest,
    validate: () => validateRelease(path.join(root, 'manifest.json')),
    inventory: () => verifyReleaseInventory(path.join(root, 'manifest.json'))};
}

test('full family release validates dictionary and array bundles and exact inventory', async t => {
  const f = await fixture(t);
  assert.equal((await f.validate()).ok, true);
  const inventory = await f.inventory();
  assert.equal(inventory.ok, true);
  assert.equal(inventory.resource_count, 5);
  assert.equal(inventory.files.filter(item => item.path.startsWith('families/bundles/')).length, 2);
});

test('bundled family decoding retains scientific value/status validation', async t => {
  const f = await fixture(t);
  const previous = f.packed.columns.values.bundle;
  f.packed.columns.values = await f.bundleColumn({dictionary: [{chi: null}], indices: [0, 0]});
  delete f.manifest.family_bundles[previous];
  await f.save();
  const result = await f.validate();
  assert.equal(result.ok, false);
  assert.ok(result.errors.includes('Value/status: backbone/chi'));
});

test('family validator rejects malformed dictionary indices and array lengths', async t => {
  for (const values of [{dictionary: [{chi: 1}], indices: [0, 1]}, [{chi: 1}]]) {
    const f = await fixture(t), previous = f.packed.columns.values.bundle;
    f.packed.columns.values = await f.bundleColumn(values);
    delete f.manifest.family_bundles[previous];
    await f.save();
    await assert.rejects(f.validate(), /dictionary|length|column/i);
  }
});

test('full family release rejects missing registry entries and unreferenced bundles', async t => {
  const f = await fixture(t), reference = f.packed.columns.values.bundle;
  const descriptor = f.manifest.family_bundles[reference];
  delete f.manifest.family_bundles[reference];
  await f.saveManifest();
  await assert.rejects(f.validate(), /Unregistered family bundle/);
  await assert.rejects(f.inventory(), /missing or unreferenced/);
  f.manifest.family_bundles[reference] = descriptor;
  f.packed.columns.values = [{chi: 1}, {chi: 1}];
  await f.save();
  await assert.rejects(f.validate(), /Unreferenced family bundle/);
  descriptor.sha256 = '0'.repeat(64);
  await f.saveManifest();
  await assert.rejects(f.validate(), /Asset hash mismatch/);
});

test('family bundle registry authenticates paths hashes and byte sizes', async t => {
  const f = await fixture(t), descriptor = f.manifest.family_bundles[f.packed.columns.values.bundle];
  for (const [field, invalid, error] of [
    ['path', 'elsewhere.json.gz', /registry descriptor/],
    ['content_sha256', '0'.repeat(64), /registry descriptor/],
    ['sha256', '0'.repeat(64), /Asset hash mismatch/],
    ['bytes', descriptor.bytes + 1, /compressed byte size mismatch/],
    ['uncompressed_bytes', descriptor.uncompressed_bytes + 1, /uncompressed byte size mismatch/],
    ['uncompressed_bytes', undefined, /registry descriptor/],
  ]) {
    const original = descriptor[field];
    descriptor[field] = invalid;
    await f.saveManifest();
    await assert.rejects(f.validate(), error, field);
    descriptor[field] = original;
  }
});

test('rehashed family bundle bytes still require the declared content hash', async t => {
  const f = await fixture(t), {bundle: reference, column} = f.packed.columns.values;
  const relative = f.manifest.family_bundles[reference].path;
  f.manifest.family_bundles[reference] = {...await f.write(relative, {
    encoding: FAMILY_BUNDLE_ENCODING, columns: {[column]: {dictionary: [{chi: 99}], indices: [0, 0]}},
  }), content_sha256: reference};
  await f.saveManifest();
  await assert.rejects(f.validate(), /hash mismatch/i);
});

test('family transport rejects changed build identities and encoding declarations', async t => {
  const f = await fixture(t);
  f.packed.build_id = 'another-release';
  await f.save();
  await assert.rejects(f.validate(), /build ID mismatch/);
  f.packed.build_id = f.manifest.build_id;
  await f.save();
  f.manifest.families[0].encoding = 'rna-family-columnar-1';
  await f.saveManifest();
  await assert.rejects(f.validate(), /encoding mismatch/);
});

test('missing family bundles permit retry while external symlinks remain forbidden', async t => {
  const f = await fixture(t);
  const file = path.join(f.root, f.manifest.family_bundles[f.packed.columns.values.bundle].path);
  const moved = path.join(f.directory, 'saved-bundle.json.gz');
  await fs.rename(file, moved);
  await assert.rejects(f.validate(), {code: 'ENOENT'});
  await fs.symlink(moved, file);
  await assert.rejects(f.validate(), /escapes release root/);
  await fs.unlink(file);
  await fs.rename(moved, file);
  assert.equal((await f.validate()).ok, true);
});
