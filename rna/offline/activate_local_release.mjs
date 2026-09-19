/** Local RNA activation only. No Git operations or remote publication. */
import fs from 'node:fs/promises';
import { constants } from 'node:fs';
import path from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';
import { randomUUID } from 'node:crypto';
import { verifyReleaseInventory } from './verify_release_inventory.mjs';
import { validateRelease } from './assets.mjs';
import { sha256 } from './output_scope.mjs';

const inside = (root, target) => target === root || target.startsWith(`${root}${path.sep}`);
const validId = value => typeof value === 'string' && /^[A-Za-z0-9][A-Za-z0-9_-]*$/.test(value);
const inventoryIdentity = inventory => JSON.stringify({
  build_id: inventory.build_id, manifest_sha256: inventory.manifest_sha256,
  files: [...inventory.files].sort((a, b) => a.path.localeCompare(b.path)),
});

async function regularFileOrAbsent(file) {
  try {
    if (!(await fs.lstat(file)).isFile()) throw new Error(`Expected a regular file: ${file}`);
    return await fs.readFile(file);
  } catch (error) {
    if (error.code === 'ENOENT') return null;
    throw error;
  }
}

/** The caller supplies an existing RNA-owned root; the CLI fixes the real root.
 * A failed copy is retained for inspection and can never be overwritten here.
 * Scientific validation is always real and cannot be injected or bypassed.
 */
export async function activateLocalRelease({ manifestPath, assetsRoot }) {
  const sourceFile = path.resolve(manifestPath), root = path.resolve(assetsRoot);
  if (path.basename(root) !== 'pure_rna') throw new Error('Activation requires a pure_rna assets root');
  if (await fs.realpath(root) !== root || !(await fs.lstat(root)).isDirectory()) {
    throw new Error('Symlink or invalid RNA assets root');
  }
  if (path.basename(sourceFile) !== 'manifest.json' || await fs.realpath(sourceFile) !== sourceFile) {
    throw new Error('Activation requires a non-symlink manifest.json');
  }
  const sourceRoot = path.dirname(sourceFile);
  if (inside(root, sourceRoot) || inside(sourceRoot, root)) throw new Error('Staging and RNA assets must be disjoint');
  const originalManifestBytes = await fs.readFile(sourceFile);
  const manifest = JSON.parse(originalManifestBytes);
  if (manifest.schema_version !== 'rna-explorer-1' || manifest.molecule_type !== 'RNA'
      || manifest.partial !== false || manifest.scalar_only || !validId(manifest.build_id)) {
    throw new Error('Activation requires a complete RNA release with a safe build ID');
  }
  if (!['entries', 'entities', 'residues'].every(key => Number.isSafeInteger(manifest.counts?.[key]) && manifest.counts[key] > 0)
      || !Number.isSafeInteger(manifest.source?.candidate_count) || manifest.source.candidate_count < manifest.counts.entries
      || manifest.source.processed_candidate_count !== manifest.source.candidate_count) {
    throw new Error('Activation requires nonempty RNA counts and complete candidate processing');
  }

  const pointerFile = path.join(root, 'manifest.json');
  const lockFile = path.join(root, '.activation.lock');
  const lock = await fs.open(lockFile, 'wx');
  const temporaryPointer = path.join(root, `.manifest-${randomUUID()}.tmp`);
  let pointerCreated = false;
  try {
    await lock.writeFile(JSON.stringify({ pid: process.pid, build_id: manifest.build_id, started_at: new Date().toISOString() }));
    const previousBytes = await regularFileOrAbsent(pointerFile);
    const previous = previousBytes ? JSON.parse(previousBytes) : null;
    if (previous && (previous.schema_version !== 'rna-explorer-1' || previous.molecule_type !== 'RNA'
        || !validId(previous.build_id) || previous.manifest !== `releases/${previous.build_id}/manifest.json`)) {
      throw new Error('Existing pointer is not a valid RNA release pointer');
    }
    if (previous?.build_id === manifest.build_id) throw new Error('Active RNA releases are immutable');

    const before = await verifyReleaseInventory(sourceFile);
    if (before.manifest_sha256 !== sha256(originalManifestBytes)) throw new Error('Staging manifest changed before validation');
    const validation = await validateRelease(sourceFile);
    if (!validation.ok) throw new Error(`RNA release validation failed: ${validation.errors.join('; ')}`);
    const afterValidation = await verifyReleaseInventory(sourceFile);
    if (inventoryIdentity(before) !== inventoryIdentity(afterValidation)) throw new Error('Staging release changed during validation');

    const releasesRoot = path.join(root, 'releases');
    await fs.mkdir(releasesRoot, { recursive: true });
    if (await fs.realpath(releasesRoot) !== releasesRoot) throw new Error('Symlink releases directory is forbidden');
    const destination = path.join(releasesRoot, manifest.build_id);
    await fs.mkdir(destination); // Exclusive ownership; old releases are retained.
    for (const relative of [...before.files.map(file => file.path), 'manifest.json']) {
      const target = path.join(destination, relative);
      await fs.mkdir(path.dirname(target), { recursive: true });
      await fs.copyFile(path.join(sourceRoot, relative), target, constants.COPYFILE_EXCL);
    }
    const installedManifest = path.join(destination, 'manifest.json');
    const installed = await verifyReleaseInventory(installedManifest);
    const afterCopy = await verifyReleaseInventory(sourceFile);
    if (inventoryIdentity(before) !== inventoryIdentity(installed)
        || inventoryIdentity(before) !== inventoryIdentity(afterCopy)) {
      throw new Error('Release inventory changed during copy; active pointer preserved');
    }
    const currentPointerBytes = await regularFileOrAbsent(pointerFile);
    if (previousBytes?.toString() !== currentPointerBytes?.toString()) throw new Error('Active pointer changed during activation');
    const pointer = {
      schema_version: 'rna-explorer-1', molecule_type: 'RNA', build_id: manifest.build_id,
      manifest: `releases/${manifest.build_id}/manifest.json`, partial: false,
    };
    const handle = await fs.open(temporaryPointer, 'wx');
    pointerCreated = true;
    try {
      await handle.writeFile(JSON.stringify(pointer, null, 2) + '\n');
      await handle.sync();
    } finally { await handle.close(); }
    await fs.rename(temporaryPointer, pointerFile);
    pointerCreated = false;
    return { ok: true, activated_at: new Date().toISOString(), previous, pointer,
      manifest_path: installedManifest, validation, inventory: installed,
      source_inventory_sha256: sha256(inventoryIdentity(before)) };
  } finally {
    if (pointerCreated) await fs.unlink(temporaryPointer);
    await lock.close();
    await fs.unlink(lockFile);
  }
}

if (process.argv[1] && import.meta.url === pathToFileURL(path.resolve(process.argv[1])).href) {
  if (process.argv.length !== 3) throw new Error('Usage: node activate_local_release.mjs STAGING_MANIFEST');
  const assetsRoot = fileURLToPath(new URL('../../assets/pure_rna', import.meta.url));
  console.log(JSON.stringify(await activateLocalRelease({ manifestPath: process.argv[2], assetsRoot }), null, 2));
}
