import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import os from 'node:os';
import path from 'node:path';
import { gzipSync } from 'node:zlib';
import { activateLocalRelease } from '../offline/activate_local_release.mjs';
import { sha256 } from '../offline/output_scope.mjs';

async function fixture(t) {
  const directory = await fs.mkdtemp(path.join(os.tmpdir(), 'rna-local-activation-'));
  t.after(() => fs.rm(directory, { recursive: true, force: true }));
  const staging = path.join(directory, 'staging'), assetsRoot = path.join(directory, 'pure_rna');
  await fs.mkdir(staging); await fs.mkdir(assetsRoot);
  const writeAsset = async (name, value) => {
    const raw = Buffer.from(JSON.stringify(value)), bytes = gzipSync(raw);
    await fs.writeFile(path.join(staging, name), bytes);
    return { path: name, bytes: bytes.length, uncompressed_bytes: raw.length, sha256: sha256(bytes) };
  };
  const manifest = {
    schema_version: 'rna-explorer-1', molecule_type: 'RNA', build_id: 'activation-test', partial: false,
    counts: { entries: 1, entities: 1, residues: 1 }, source: { candidate_count: 1, processed_candidate_count: 1 },
    metadata: await writeAsset('metadata.json.gz', {
      entries: [{ pdb_id: '1SDR', selected_model_id: '1' }],
      entities: [{ pdb_id: '1SDR', entity_id: '1' }],
    }),
    families: [{ id: 'backbone', level: 'residue', row_count: 1, parameters: [{ id: 'alpha' }],
      ...await writeAsset('backbone.json.gz', [{ id: 'residue-1', pdb_id: '1SDR', entity_id: '1', model_id: '1',
        values: { alpha: 30 }, statuses: { alpha: 'available' }, is_terminal_any: true }]) }],
    survey: { scalars: { terms: {} }, coordinates: { groups: {} } }, relations: {},
    provenance: { decisions: { row_count: 1,
      ...await writeAsset('decisions.json.gz', [{ pdb_id: '1SDR', accepted: true }]) } },
  };
  const manifestPath = path.join(staging, 'manifest.json');
  const save = () => fs.writeFile(manifestPath, JSON.stringify(manifest));
  await save();
  const pointer = { schema_version: 'rna-explorer-1', molecule_type: 'RNA', build_id: 'previous',
    manifest: 'releases/previous/manifest.json', partial: false };
  const pointerPath = path.join(assetsRoot, 'manifest.json');
  await fs.writeFile(pointerPath, JSON.stringify(pointer));
  await fs.mkdir(path.join(assetsRoot, 'releases', 'previous'), { recursive: true });
  await fs.writeFile(path.join(assetsRoot, 'releases', 'previous', 'sentinel'), 'retained');
  return { directory, staging, assetsRoot, manifest, manifestPath, save, pointer, pointerPath,
    activate: () => activateLocalRelease({ manifestPath, assetsRoot }) };
}

test('Activation validates, copies exclusively and atomically selects a complete RNA release', async t => {
  const f = await fixture(t);
  const result = await f.activate();
  assert.equal(result.ok, true); assert.equal(result.validation.ok, true);
  assert.equal(result.inventory.file_count, 4);
  assert.deepEqual(result.previous, f.pointer);
  assert.equal(JSON.parse(await fs.readFile(f.pointerPath)).build_id, f.manifest.build_id);
  assert.equal(await fs.readFile(path.join(f.assetsRoot, 'releases', 'previous', 'sentinel'), 'utf8'), 'retained');
  assert.deepEqual(await fs.readFile(result.manifest_path), await fs.readFile(f.manifestPath));
  await assert.rejects(f.activate(), /immutable/);
  await assert.rejects(fs.access(path.join(f.assetsRoot, '.activation.lock')), { code: 'ENOENT' });
});

test('Activation rejects invalid science, partial candidates and overwritten destinations', async t => {
  const f = await fixture(t), originalPointer = await fs.readFile(f.pointerPath);
  f.manifest.partial = true; await f.save(); await assert.rejects(f.activate(), /complete RNA/);
  f.manifest.partial = false; f.manifest.scalar_only = true;
  await f.save(); await assert.rejects(f.activate(), /complete RNA/);
  delete f.manifest.scalar_only; f.manifest.source.processed_candidate_count = 0;
  await f.save(); await assert.rejects(f.activate(), /complete candidate processing/);
  f.manifest.source.processed_candidate_count = 1; f.manifest.counts.residues = 2;
  await f.save(); await assert.rejects(f.activate(), /validation failed: Residue count/);
  f.manifest.counts.residues = 1; await f.save();
  const target = path.join(f.assetsRoot, 'releases', f.manifest.build_id);
  await fs.mkdir(target); await fs.writeFile(path.join(target, 'owned'), 'keep');
  await assert.rejects(f.activate(), { code: 'EEXIST' });
  assert.equal(await fs.readFile(path.join(target, 'owned'), 'utf8'), 'keep');
  assert.deepEqual(await fs.readFile(f.pointerPath), originalPointer);
});

test('Activation refuses unsafe roots, symlinks, staging overlap and concurrent ownership', async t => {
  const f = await fixture(t);
  await assert.rejects(activateLocalRelease({ manifestPath: f.manifestPath, assetsRoot: f.directory }), /pure_rna/);
  f.manifest.build_id = '../escape'; await f.save(); await assert.rejects(f.activate(), /safe build ID/);
  f.manifest.build_id = 'activation-test'; await f.save();
  const nested = path.join(f.assetsRoot, 'staging'); await fs.mkdir(nested);
  await fs.copyFile(f.manifestPath, path.join(nested, 'manifest.json'));
  await assert.rejects(activateLocalRelease({ manifestPath: path.join(nested, 'manifest.json'), assetsRoot: f.assetsRoot }), /disjoint/);
  await fs.symlink(path.join(f.staging, 'metadata.json.gz'), path.join(f.staging, 'alias'));
  await assert.rejects(f.activate(), /Symlink release resource/); await fs.unlink(path.join(f.staging, 'alias'));
  const lock = path.join(f.assetsRoot, '.activation.lock'); await fs.writeFile(lock, 'other owner');
  await assert.rejects(f.activate(), { code: 'EEXIST' });
  assert.equal(await fs.readFile(lock, 'utf8'), 'other owner');
});

test('A source mutation during copying cannot replace the active pointer', async t => {
  const f = await fixture(t), originalPointer = await fs.readFile(f.pointerPath);
  const copyFile = fs.copyFile;
  let changed = false;
  t.mock.method(fs, 'copyFile', async (...args) => {
    await copyFile(...args);
    if (!changed) {
      changed = true;
      await fs.appendFile(path.join(f.staging, 'metadata.json.gz'), 'corruption');
    }
  });
  await assert.rejects(f.activate(), /Compressed size mismatch/);
  assert.deepEqual(await fs.readFile(f.pointerPath), originalPointer);
  assert.equal(await fs.readFile(path.join(f.assetsRoot, 'releases', 'previous', 'sentinel'), 'utf8'), 'retained');
  await assert.rejects(fs.access(path.join(f.assetsRoot, '.activation.lock')), { code: 'ENOENT' });
});
