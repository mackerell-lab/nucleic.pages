import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import os from 'node:os';
import path from 'node:path';
import { gzipSync } from 'node:zlib';
import { sha256 } from '../offline/output_scope.mjs';
import { verifyReleaseInventory } from '../offline/verify_release_inventory.mjs';

test('Inventory verifies bytes and rejects unreferenced files, symlinks and incomplete sizes', async t => {
  const root = await fs.mkdtemp(path.join(os.tmpdir(), 'rna-inventory-'));
  t.after(() => fs.rm(root, { recursive: true, force: true }));
  const raw = Buffer.from(JSON.stringify({ entries: [] })), bytes = gzipSync(raw);
  const manifest = { schema_version: 'rna-explorer-1', molecule_type: 'RNA', build_id: 'test', metadata: {
    path: 'metadata.json.gz', bytes: bytes.length, uncompressed_bytes: raw.length, sha256: sha256(bytes) } };
  const file = path.join(root, 'manifest.json');
  const save = () => fs.writeFile(file, JSON.stringify(manifest));
  await save(); await fs.writeFile(path.join(root, manifest.metadata.path), bytes);
  const valid = await verifyReleaseInventory(file);
  assert.equal(valid.file_count, 2); assert.equal(valid.resource_bytes, bytes.length);
  await fs.writeFile(path.join(root, 'stray.txt'), 'unreferenced');
  await assert.rejects(verifyReleaseInventory(file), /unreferenced/);
  await fs.unlink(path.join(root, 'stray.txt'));
  await fs.symlink(file, path.join(root, 'alias'));
  await assert.rejects(verifyReleaseInventory(file), /Symlink/);
  await fs.unlink(path.join(root, 'alias'));
  manifest.metadata.uncompressed_bytes++;
  await save(); await assert.rejects(verifyReleaseInventory(file), /Uncompressed size/);
  manifest.metadata.uncompressed_bytes--;
  manifest.metadata.bytes++;
  await save(); await assert.rejects(verifyReleaseInventory(file), /Compressed size/);
  delete manifest.metadata.bytes;
  await save(); await assert.rejects(verifyReleaseInventory(file), /Missing inventory byte sizes/);
});
