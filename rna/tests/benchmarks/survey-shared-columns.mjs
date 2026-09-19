/** Storage feasibility only: independent column streams are not release assets. */
import { readFile } from 'node:fs/promises';
import { createHash } from 'node:crypto';
import { gzipSync, gunzipSync } from 'node:zlib';
import path from 'node:path';
const manifestFile = process.argv[2];
if (!manifestFile) throw new Error('Provide a release manifest path');
const manifest = JSON.parse(await readFile(manifestFile));
const terms = Object.entries(manifest.survey.scalars.terms);
const seen = new Map(), byField = {};
let sourceBytes = 0, separateColumnBytes = 0, uniqueColumnBytes = 0;
for (const [term, descriptor] of terms) {
  const bytes = await readFile(path.resolve(path.dirname(manifestFile), descriptor.path));
  sourceBytes += bytes.length;
  const data = JSON.parse(gunzipSync(bytes));
  for (const [field, column] of Object.entries(data.columns)) {
    const serialized = JSON.stringify(column);
    const hash = createHash('sha256').update(serialized).digest('hex');
    const size = gzipSync(serialized, { level: 9 }).length;
    separateColumnBytes += size;
    const entry = byField[field] ??= { references: 0, repeatedBytes: 0 };
    entry.references++;
    if (seen.has(hash)) entry.repeatedBytes += size;
    else { seen.set(hash, { term, field, size }); uniqueColumnBytes += size; }
  }
}
console.log(JSON.stringify({ build: manifest.build_id, terms: terms.length, sourceBytes,
  separateColumnBytes, uniqueColumnBytes, repeatedColumnBytes: separateColumnBytes - uniqueColumnBytes,
  byField, limitation: 'Feasibility estimate excludes reference metadata and request overhead; not a deployable format or measured release reduction.' }, null, 2));
