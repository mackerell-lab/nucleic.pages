/** Isolated snapshot allocation measurement; run each mode in a fresh Node process. */
import { readFile, mkdir, writeFile } from 'node:fs/promises';
import path from 'node:path';
import { gunzipSync } from 'node:zlib';
import { createHash } from 'node:crypto';
import { performance } from 'node:perf_hooks';
import { deepFreeze } from '../../core/repository.js';
import { selectRows } from '../../core/selection.js';
import { familyParameters } from '../../core/registry.js';
import { distribution, histogram2D } from '../../core/analysis.js';
import { join } from '../../core/joints.js';
import { createPlotSnapshot } from '../../core/export.js';

const mode = process.argv[2];
if (!globalThis.gc || !['copy', 'transfer'].includes(mode)) throw new Error('Use node --expose-gc snapshot-memory.mjs copy|transfer');
const workspace = path.resolve(process.env.RNA_WORKSPACE || process.cwd());
const release = path.join(workspace, 'nucleic.pages/assets/pure_rna/releases/full_20260916');
const read = async relative => {
  const bytes = await readFile(path.join(release, relative));
  return JSON.parse((relative.endsWith('.gz') ? gunzipSync(bytes) : bytes).toString());
};
const manifest = await read('manifest.json');
const family = manifest.families.find(item => item.id === 'backbone');
const metadata = deepFreeze(await read(manifest.metadata.path));
const table = deepFreeze(await read(family.path));
const selectionSpec = { components: 'all', methods: [], resolutionMax: null, contexts: [], includeEnds: true };
const selected = selectRows(table, metadata, selectionSpec);
const parameters = familyParameters(manifest, 'backbone');
const chi = parameters.find(item => item.id === 'chi'), delta = parameters.find(item => item.id === 'delta');
const results = [distribution(selected.rows, chi, { groupBy: 'base' }), histogram2D(join(selected.rows, selected.rows).points, chi, delta)];
globalThis.gc();
const before = process.memoryUsage().heapUsed, started = performance.now();
const snapshots = results.map((result, i) => createPlotSnapshot({ result, snapshot_id: `allocation-${i}`, buildId: manifest.build_id,
  selectionSpec, transferResult: mode === 'transfer' }));
const elapsedMs = performance.now() - started;
globalThis.gc();
const after = process.memoryUsage().heapUsed;
const hash = createHash('sha256');
for (const series of snapshots[0].result.series) for (let i = 0; i < series.values.length; i++) hash.update(JSON.stringify([series.key, series.rowIds[i], series.values[i]]));
for (const point of snapshots[1].result.points) hash.update(JSON.stringify([point.left_id, point.right_id, point.x, point.y]));
const report = { mode, build: manifest.build_id, node: process.version, rows: selected.rows.length,
  finiteChi: snapshots[0].result.coverage.plottedRows, finiteJoint: snapshots[1].result.points.length,
  retainedInputRows: results[0].series.reduce((sum, series) => sum + series.rows.length, 0),
  beforeBytes: before, afterBytes: after, snapshotAdditionalBytes: after - before, elapsedMs,
  observationHash: hash.digest('hex'), measuredAt: new Date().toISOString(),
  scope: 'Snapshot allocation only, with source tables, selection and both original results retained; forced GC in a fresh Node process. Copy mode uses the current safe copying API, not historical browser code.' };
const output = path.join(workspace, 'data/pure_rna/validation/snapshot-allocation'); await mkdir(output, { recursive: true });
await writeFile(path.join(output, `${mode}.json`), JSON.stringify(report, null, 2));
console.log(JSON.stringify(report, null, 2));
