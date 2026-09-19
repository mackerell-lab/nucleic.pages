/** Actual streaming coordinate cancellation; no substituted cache or coordinate data. */
import assert from 'node:assert/strict';
import { mkdir, readFile, writeFile } from 'node:fs/promises';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { createHash } from 'node:crypto';
import { waitReady } from './helpers.mjs';

const workspace = path.resolve(process.env.RNA_WORKSPACE || process.cwd());
const output = process.env.RNA_BROWSER_OUTPUT || path.join(workspace, 'data/pure_rna/browser_validation/cooperative-coordinates');
const phase = process.env.RNA_COORDINATE_PHASE || 'after', source = process.env.RNA_HELD_SOURCE;
const base = process.env.RNA_URL || 'http://127.0.0.1:8767/nucleic.pages/rna/';
await mkdir(output, { recursive: true });
const { chromium } = await import(pathToFileURL(process.env.PLAYWRIGHT_MODULE || path.join(workspace, 'data/pure_rna/browser-tools/node_modules/playwright-core/index.mjs')).href);
const browser = await chromium.launch({ headless: true });
const page = await browser.newPage({ viewport: { width: 1440, height: 1000 } });
const report = { phase, started: new Date().toISOString(), errors: [], sourceSha256: {}, checks: [], limitation: 'One real streaming run per source. Browser HTTP may be warm, but parsed partitions are not cached or replaced. Existing network/decompression tasks may already allow cancellation. No blanket speed claim.' };
for (const file of ['app/PureRnaExplorer.js', 'app/NucleicAcidExplorer.js', 'core/repository.js', 'core/coordinates.js']) report.sourceSha256[file] = createHash('sha256').update(await readFile(source ? path.join(source, file) : new URL('../../' + file, import.meta.url))).digest('hex');
page.on('pageerror', error => report.errors.push(error.message));
page.on('console', message => { if (message.type() === 'error') report.errors.push(message.text()); });
if (source) await page.route(base + '**/*.js', async route => {
  const relative = new URL(route.request().url()).pathname.slice(new URL(base).pathname.length);
  await route.fulfill({ status: 200, contentType: 'application/javascript', body: await readFile(path.join(source, relative)) });
});
try {
  await page.goto(base); await waitReady(page);
  await page.click('#coordinatesLoad'); await waitReady(page);
  report.initial = await page.evaluate(() => { const a = window.rnaExplorer; return { build: a.manifest.build_id, group: a.state.survey.coordinateGroup,
    opening: a.state.survey.coordinateOpening, summaryGroups: a.coordinateSummary.length,
    partitions: a.manifest.survey.coordinates.groups[a.state.survey.coordinateGroup].partitions.length }; });
  assert(report.initial.summaryGroups > 0); assert.match(report.initial.group, /cytosine_standard_pair/);
  console.log(JSON.stringify({ phase, initial: report.initial }));
  const run = await page.evaluate(async () => {
    const a = window.rnaExplorer, original = a.repository.iterateSurveyCoordinates.bind(a.repository);
    const { CoordinateSummary } = await import('./core/coordinates.js');
    const originalAdd = CoordinateSummary.prototype.add;
    const started = performance.now(), events = [], heartbeats = [], longtasks = [];
    const obsoleteRevision = a.revision + 1; let queued = false, delivered = false, queuedAt, deliveredAt;
    let obsoleteYields = 0, obsoleteAdds = 0, pending = null;
    const observer = new PerformanceObserver(list => { for (const e of list.getEntries()) longtasks.push({ start: e.startTime - started, duration: e.duration }); });
    observer.observe({ type: 'longtask' }); const heartbeat = setInterval(() => heartbeats.push(performance.now() - started), 25);
    let resolveNew; const newFinished = new Promise(resolve => { resolveNew = resolve; });
    CoordinateSummary.prototype.add = function (row) { if (a.revision === obsoleteRevision) obsoleteAdds++; return originalAdd.call(this, row); };
    a.repository.iterateSurveyCoordinates = async function* (...args) {
      const owner = a.revision;
      for await (const chunk of original(...args)) {
        if (owner === obsoleteRevision) {
          obsoleteYields++; events.push({ stage: 'partition', path: chunk.partition_path, rows: chunk.rows.length, at: performance.now() - started });
          if (!queued) {
            queued = true; queuedAt = performance.now() - started;
            setTimeout(async () => {
              delivered = true; deliveredAt = performance.now() - started;
              const input = document.querySelector('#coordinateOpeningSelect'); input.value = 'small'; input.dispatchEvent(new Event('change', { bubbles: true }));
              pending = { loading: document.querySelector('#appStatus').dataset.state, completedCoordinateKey: a.completedCoordinateKey,
                mainExportDisabled: document.querySelector('#filteredCsvDownload').disabled, obsoleteAdds, obsoleteYields };
              while (document.querySelector('#appStatus').dataset.state === 'loading') await new Promise(resolve => setTimeout(resolve, 10));
              resolveNew();
            }, 0);
          }
        }
        yield chunk;
      }
    };
    try {
      const input = document.querySelector('#coordinateOpeningSelect'); input.value = 'large'; input.dispatchEvent(new Event('change', { bubbles: true }));
      await newFinished; await a.commitQueue; await new Promise(resolve => setTimeout(resolve, 40));
    } finally {
      a.repository.iterateSurveyCoordinates = original; CoordinateSummary.prototype.add = originalAdd;
      clearInterval(heartbeat); observer.disconnect();
    }
    return { obsoleteRevision, finalRevision: a.revision, queuedAt, deliveredAt, inputDelayMs: deliveredAt - queuedAt,
      obsoleteYields, obsoleteAdds, events, heartbeats, longtasks, pending, delivered,
      finalOpening: a.state.survey.coordinateOpening, finalKey: a.completedCoordinateKey, final: a.coordinateSummary,
      tableRows: document.querySelectorAll('#baseGeometryCoordBody tr').length, inflightCoordinates: a.repository.coordinateRequests.size };
  });
  report.run = run; await waitReady(page);
  assert(run.delivered); assert.equal(run.finalOpening, 'small'); assert.equal(JSON.parse(run.finalKey).opening, 'small');
  assert(run.finalRevision > run.obsoleteRevision); assert(run.obsoleteYields < report.initial.partitions);
  assert.equal(run.pending.loading, 'loading'); assert(run.pending.mainExportDisabled);
  if (phase === 'after') assert.equal(run.obsoleteAdds, 0, 'Queued input runs before obsolete first-chunk accumulation');
  assert.equal(run.inflightCoordinates, 0);
  report.checks.push({ name: 'Real partition delivery cancels obsolete opening selection', obsoleteYields: run.obsoleteYields, obsoleteAdds: run.obsoleteAdds, inputDelayMs: run.inputDelayMs });

  const independent = await page.evaluate(async () => {
    const a = window.rnaExplorer, { selectRows } = await import('./core/selection.js');
    const eligible = selectRows([], a.metadata, a.state.selection).entryIds, pairs = (await a.openingIndex(a.state)).pairs;
    const groups = new Map(); let partitions = 0, selectedAtoms = 0;
    for await (const chunk of a.repository.iterateSurveyCoordinates(a.state.survey.coordinateGroup, { entryIds: eligible })) {
      partitions++;
      for (const row of selectRows(chunk, a.metadata, { ...a.state.selection, contexts: [] }).rows) {
        if (!row.pair_id || !pairs.has(row.pair_id) || row.opening_bin !== 'small') continue;
        const context = row.context ?? row.sequence_context ?? row.base ?? row.base_code ?? '';
        const atom = row.atom_label ?? row.atom ?? row.atom_name ?? row.atom_id, xyz = row.xyz ?? [row.x, row.y, row.z];
        if (!atom || xyz.length !== 3 || !xyz.every(Number.isFinite) || (row.status && !['ok', 'available', 'computed', 'valid'].includes(row.status))) continue;
        const key = JSON.stringify([context, atom]);
        if (!groups.has(key)) groups.set(key, { context, atom_label: atom, n: 0, sums: [0, 0, 0], squared: 0, entries: new Set(), residues: new Set(), pairs: new Set() });
        const g = groups.get(key), pdb = String(row.pdb_id ?? row.accession ?? row.entry_id ?? '').toUpperCase();
        g.n++; selectedAtoms++; g.entries.add(pdb); const identity = id => JSON.stringify([pdb, row.model_id ?? null, id]);
        g.residues.add(identity(row.target_residue_id ?? row.residue_id)); g.pairs.add(identity(row.pair_id));
        xyz.forEach((value, axis) => { g.sums[axis] += value; g.squared += value * value; });
      }
    }
    return { partitions, selectedAtoms, groups: [...groups.values()].map(g => {
      const mean = g.sums.map(value => value / g.n);
      return { context: g.context, atom_label: g.atom_label, n: g.n, mean,
        rms: Math.sqrt(Math.max(0, g.squared / g.n - mean.reduce((sum, value) => sum + value * value, 0))),
        entries: g.entries.size, residues: g.residues.size, pairs: g.pairs.size };
    }) };
  });
  const actual = new Map(run.final.map(g => [JSON.stringify([g.context, g.atom_label]), g]));
  assert.equal(actual.size, independent.groups.length); assert(independent.selectedAtoms > 0);
  let maximumMeanError = 0, maximumRmsError = 0;
  for (const group of independent.groups) {
    const observed = actual.get(JSON.stringify([group.context, group.atom_label])); assert(observed);
    for (const field of ['n', 'entries', 'residues', 'pairs']) assert.equal(observed[field], group[field]);
    group.mean.forEach((value, axis) => { const error = Math.abs(value - observed.mean[axis]); maximumMeanError = Math.max(maximumMeanError, error); assert(error < 1e-10); });
    const error = Math.abs(group.rms - observed.rms); maximumRmsError = Math.max(maximumRmsError, error); assert(error < 1e-8);
  }
  report.independent = { ...independent, maximumMeanError, maximumRmsError };
  report.checks.push({ name: 'Independent source sums and identity sets match final small-opening aggregate', groups: actual.size, atoms: independent.selectedAtoms, maximumMeanError, maximumRmsError });
  assert.deepEqual(report.errors, []); report.passed = true;
  console.log(JSON.stringify({ phase, checks: report.checks.length, obsoleteYields: run.obsoleteYields, obsoleteAdds: run.obsoleteAdds, inputDelayMs: run.inputDelayMs, groups: actual.size, atoms: independent.selectedAtoms }));
} catch (error) { report.passed = false; report.failure = error.stack; throw error; }
finally { report.finished = new Date().toISOString(); await writeFile(path.join(output, `report-${phase}.json`), JSON.stringify(report, null, 2)); await browser.close(); }
