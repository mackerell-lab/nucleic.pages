import fs from 'node:fs/promises';
import path from 'node:path';
import crypto from 'node:crypto';

export const sha256 = value => crypto.createHash('sha256').update(value).digest('hex');
const inside = (root, target) => target === root || target.startsWith(root + path.sep);

/** Every writer must use an explicit RNA-owned root; resolve symlink ancestors. */
export class OutputScope {
  constructor(roots) {
    if (!roots?.length) throw new Error('At least one RNA output root is required');
    this.roots = roots.map(root => path.resolve(root));
  }
  async resolve(target) {
    const absolute = path.resolve(target);
    const root = this.roots.find(candidate => inside(candidate, absolute));
    if (!root) throw new Error(`Output escapes RNA scope: ${absolute}`);
    for (const candidate of [root, absolute]) {
      let ancestor = candidate;
      while (true) {
        try {
          const real = await fs.realpath(ancestor);
          if (real !== ancestor) throw new Error(`Symlink output ancestor is forbidden: ${ancestor}`);
          break;
        } catch (error) {
          if (error.code !== 'ENOENT') throw error;
          const parent = path.dirname(ancestor);
          if (parent === ancestor) throw error;
          ancestor = parent;
        }
      }
    }
    return absolute;
  }
  async write(target, value) {
    const resolved = await this.resolve(target);
    await fs.mkdir(path.dirname(resolved), {recursive: true});
    // Recheck after mkdir so an existing symlink cannot hide in newly resolved parents.
    await this.resolve(resolved);
    const temporary = `${resolved}.tmp-${process.pid}-${crypto.randomBytes(6).toString('hex')}`;
    try { await fs.writeFile(temporary, value, {flag: 'wx'}); await fs.rename(temporary, resolved); }
    finally { await fs.rm(temporary, {force: true}); }
    return {path: resolved, sha256: sha256(value), bytes: Buffer.byteLength(value)};
  }
  json(target, value) { return this.write(target, JSON.stringify(value, null, 2) + '\n'); }
}

export async function readJson(file) { return JSON.parse(await fs.readFile(file, 'utf8')); }
export async function exists(file) { try { await fs.access(file); return true; } catch { return false; } }
