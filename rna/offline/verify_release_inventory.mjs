/** Independent byte inventory for staging/publishing; does not activate releases. */
import fs from 'node:fs/promises';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { gunzipSync } from 'node:zlib';
import { sha256 } from './output_scope.mjs';

export function releaseDescriptors(manifest) {
  return [manifest.metadata, ...(manifest.families ?? []),
    ...Object.values(manifest.relations ?? {}), manifest.provenance?.decisions,
    ...Object.values(manifest.survey?.scalars?.terms ?? {}),
    ...Object.values(manifest.survey?.coordinates?.groups ?? {}).flatMap(group => group.partitions ?? [group]),
    ...Object.values(manifest.survey?.bundles ?? {}),
    ...Object.values(manifest.survey?.shared_columns ?? {})].filter(Boolean);
}

export async function verifyReleaseInventory(manifestPath) {
  const file = path.resolve(manifestPath), root = path.dirname(file);
  if (await fs.realpath(file) !== file) throw new Error('Symlink manifest path is forbidden');
  const manifestBytes = await fs.readFile(file), manifest = JSON.parse(manifestBytes);
  if (manifest.schema_version !== 'rna-explorer-1' || manifest.molecule_type !== 'RNA' || !manifest.build_id) throw new Error('Expected an RNA release manifest');
  const expected = new Map();
  for (const descriptor of releaseDescriptors(manifest)) {
    const relative = descriptor.path;
    if (typeof relative !== 'string' || !/^[a-zA-Z0-9_./-]+$/.test(relative)
        || relative.split('/').some(part => !part || part === '.' || part === '..')) throw new Error('Unsafe inventory resource path');
    if (relative === path.basename(file) || expected.has(relative)) throw new Error(`Duplicate inventory resource: ${relative}`);
    if (![descriptor.bytes, descriptor.uncompressed_bytes].every(value => Number.isSafeInteger(value) && value >= 0)) throw new Error(`Missing inventory byte sizes: ${relative}`);
    expected.set(relative, descriptor);
  }
  const actual = [];
  async function walk(relative = '') {
    for (const item of await fs.readdir(path.join(root, relative), { withFileTypes: true })) {
      const child = path.join(relative, item.name);
      if (item.isSymbolicLink()) throw new Error(`Symlink release resource: ${child}`);
      if (item.isDirectory()) await walk(child);
      else if (item.isFile()) actual.push(child.split(path.sep).join('/'));
      else throw new Error(`Unsupported release resource: ${child}`);
    }
  }
  await walk();
  const wanted = [...expected.keys(), path.basename(file)].sort(); actual.sort();
  if (JSON.stringify(actual) !== JSON.stringify(wanted)) throw new Error('Release inventory contains missing or unreferenced files');
  let resourceBytes = 0, uncompressedBytes = 0;
  const files = [];
  for (const [relative, descriptor] of expected) {
    const bytes = await fs.readFile(path.join(root, relative));
    if (bytes.length !== descriptor.bytes) throw new Error(`Compressed size mismatch: ${relative}`);
    if (sha256(bytes) !== descriptor.sha256) throw new Error(`Asset hash mismatch: ${relative}`);
    const raw = relative.endsWith('.gz') ? gunzipSync(bytes) : bytes;
    if (raw.length !== descriptor.uncompressed_bytes) throw new Error(`Uncompressed size mismatch: ${relative}`);
    resourceBytes += bytes.length; uncompressedBytes += raw.length;
    files.push({ path: relative, bytes: bytes.length, sha256: descriptor.sha256 });
  }
  return { ok: true, build_id: manifest.build_id, manifest_sha256: sha256(manifestBytes),
    resource_count: expected.size, file_count: actual.length, resource_bytes: resourceBytes,
    total_bytes: resourceBytes + manifestBytes.length, uncompressed_resource_bytes: uncompressedBytes,
    checked_at: new Date().toISOString(), files };
}

if (process.argv[1] && import.meta.url === pathToFileURL(path.resolve(process.argv[1])).href) {
  if (process.argv.length !== 3) throw new Error('Usage: node verify_release_inventory.mjs MANIFEST');
  console.log(JSON.stringify(await verifyReleaseInventory(process.argv[2]), null, 2));
}
