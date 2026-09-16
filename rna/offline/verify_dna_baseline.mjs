import { createHash } from 'node:crypto';
import { readFile } from 'node:fs/promises';
import { resolve, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';

/** Verify the protected working-tree bytes, including pre-existing user edits. */
export async function verifyDnaBaseline(baselinePath, repositoryRoot) {
  const baseline = JSON.parse(await readFile(baselinePath, 'utf8'));
  const failures = [];
  for (const [relativePath, expected] of Object.entries(baseline.files)) {
    try {
      const bytes = await readFile(resolve(repositoryRoot, relativePath));
      const actual = createHash('sha256').update(bytes).digest('hex');
      if (actual !== expected) failures.push({ path: relativePath, expected, actual });
    } catch (error) {
      failures.push({ path: relativePath, error: error.code || error.message });
    }
  }
  return { ok: failures.length === 0, checked: Object.keys(baseline.files).length, failures };
}

if (process.argv[1] && resolve(process.argv[1]) === fileURLToPath(import.meta.url)) {
  const here = dirname(fileURLToPath(import.meta.url));
  const baselinePath = process.argv[2] || resolve(here, '../../../data/pure_rna/dna_baseline.json');
  const result = await verifyDnaBaseline(baselinePath, resolve(here, '../..'));
  process.stdout.write(`${JSON.stringify(result, null, 2)}\n`);
  if (!result.ok) process.exitCode = 1;
}
