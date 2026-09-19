/** Compare lossless transport candidates; never changes release assets. */
import { readFile } from 'node:fs/promises';
import { gzipSync, gunzipSync } from 'node:zlib';
import assert from 'node:assert/strict';
import { encodeFamilyRows, decodeFamilyRows, decodeSurveyRows } from '../../core/survey-codec.js';

if (!process.argv[2]) throw new Error('Provide one or more scalar .json.gz files');
for (const filename of process.argv.slice(2)) {
  const source = await readFile(filename);
  const data = JSON.parse(gunzipSync(source));
  const rows = decodeSurveyRows(data);
  const baseline = gzipSync(JSON.stringify(data), { level: 9 });
  // The existing family codec is a generic dictionary prototype. Its family
  // encoding tag is NOT suitable for publishing Survey data without migration.
  const candidate = encodeFamilyRows(rows, data.build_id);
  assert.deepEqual(decodeFamilyRows(candidate), rows);
  const compressed = gzipSync(JSON.stringify(candidate), { level: 9 });
  console.log(JSON.stringify({ filename, rows: rows.length, sourceBytes: source.length,
    baselineBytes: baseline.length, candidateBytes: compressed.length,
    savedBytes: baseline.length - compressed.length,
    savedPercent: 100 * (1 - compressed.length / baseline.length), semanticEquality: true }));
}
