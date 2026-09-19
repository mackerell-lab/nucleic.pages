#!/usr/bin/env node
import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import {execFile} from 'node:child_process';
import {promisify} from 'node:util';
import {OutputScope, sha256, readJson, exists} from './output_scope.mjs';
import {discover, fetchCoordinates, parallelMap} from './discovery.mjs';
import {loadAnnotations, normalizeAnnotation} from './annotations.mjs';

const run = promisify(execFile);
const here = path.dirname(fileURLToPath(import.meta.url));
const website = path.resolve(here, '../..');
const workspace = path.dirname(website);
const dataRoot = path.join(workspace, 'data/pure_rna');
const assetsRoot = path.join(website, 'assets/pure_rna');
export const STAGES = ['discover', 'fetch', 'normalize', 'annotate', 'residue', 'geometry', 'survey', 'assets', 'validate-release'];
const HELP = `Pure RNA dataset builder (never writes DNA paths)
  --build-id ID                  Required stable build identifier
  --stage STAGE | --through STAGE
  --discovery-snapshot FILE      Replay an archived complete RCSB query
  --annotation-snapshot FILE     Replay archived NAKB batches
  --only-ids ID,ID --limit N      Permanently labels output partial
  --offline --resume            Require cache / reuse hash-matching checkpoints
  --concurrency N --timeout-ms N --retries N
  --python FILE --fr3d-path DIR --refresh-coordinates
Stages: ${STAGES.join(', ')}
`;

export function parseArgs(argv) {
  const options = {concurrency: 4, timeoutMs: 60000, retries: 3, python: path.join(dataRoot, 'venv/bin/python')};
  const values = new Map([['--build-id', 'buildId'], ['--stage', 'stage'], ['--through', 'through'],
    ['--discovery-snapshot', 'snapshot'], ['--annotation-snapshot', 'annotationSnapshot'], ['--only-ids', 'onlyIds'],
    ['--limit', 'limit'], ['--concurrency', 'concurrency'], ['--timeout-ms', 'timeoutMs'], ['--retries', 'retries'],
    ['--python', 'python'], ['--fr3d-path', 'fr3dPath']]);
  for (let i = 0; i < argv.length; i++) {
    const flag = argv[i];
    if (flag === '--help') { options.help = true; continue; }
    if (['--offline', '--resume', '--refresh-coordinates'].includes(flag)) { options[flag.slice(2).replace(/-([a-z])/g, (_, c) => c.toUpperCase())] = true; continue; }
    if (!values.has(flag)) throw new Error(`Unknown flag: ${flag}`);
    const value = argv[++i];
    if (!value || value.startsWith('--')) throw new Error(`Missing value for ${flag}`);
    options[values.get(flag)] = value;
  }
  if (options.help) return options;
  if (!options.buildId || !/^[A-Za-z0-9][A-Za-z0-9_.-]*$/.test(options.buildId)) throw new Error('A safe --build-id is required');
  if (options.stage && options.through) throw new Error('Choose --stage or --through');
  if (!options.stage && !options.through) options.through = 'validate-release';
  if (!STAGES.includes(options.stage ?? options.through)) throw new Error('Unknown build stage');
  for (const name of ['limit', 'concurrency', 'timeoutMs', 'retries']) if (options[name] != null) {
    options[name] = Number(options[name]);
    if (!Number.isInteger(options[name]) || options[name] < (name === 'retries' ? 0 : 1)) throw new Error(`Invalid ${name}`);
  }
  if (options.concurrency > 16) throw new Error('Concurrency must be at most 16');
  if (options.onlyIds) options.onlyIds = options.onlyIds.split(',').map(id => id.trim().toUpperCase()).filter(Boolean);
  if (options.refreshCoordinates && options.offline) throw new Error('Cannot refresh coordinates offline');
  return options;
}

async function fingerprint(files) {
  const records = [];
  for (const file of files) records.push([file, sha256(await fs.readFile(path.join(here, file)))]);
  return sha256(JSON.stringify(records));
}

async function stageFiles(stage) {
  const common = ['build_dataset.mjs','output_scope.mjs'];
  const numeric = ['numeric.mjs','base_frames.mjs','vendor/dna_geometry_core.mjs','vendor/base_templates.mjs','vendor/rna_base_templates.mjs'];
  const mapping = {
    discover:['discovery.mjs'], fetch:['discovery.mjs'], normalize:['normalize_mmcif.py','requirements.txt'],
    annotate:['annotations.mjs','../config/canonical_rna_v1.json','../config/components_v1.json'],
    residue:['residue_geometry.mjs','parameter_registry.mjs',...numeric],
    geometry:['interactions.mjs','geometry_adapter.mjs','stem_projection.mjs','fr3d_provider.py',...numeric],
    survey:['survey.mjs','survey_terms.mjs',...numeric],
    assets:['assets.mjs','../core/survey-codec.js','../core/shared-survey-codec.js','../core/bundled-survey-codec.js','../core/bundled-family-codec.js','../core/packed-coordinate-codec.js','parameter_registry.mjs','survey_terms.mjs','../config/geometry_parameters_v1.json'],
    'validate-release':['assets.mjs','../core/survey-codec.js','../core/shared-survey-codec.js','../core/bundled-survey-codec.js','../core/bundled-family-codec.js','../core/packed-coordinate-codec.js'],
  };
  return [...new Set([...common,...mapping[stage]])];
}

async function runtimeFingerprint(stage, options) {
  if (!['normalize','geometry'].includes(stage)) return null;
  const packages = stage === 'normalize' ? ['gemmi'] : ['numpy','scipy','biopython'];
  const code = 'import sys,json,importlib.metadata; print(json.dumps({"python":sys.version,"packages":{p:importlib.metadata.version(p) for p in sys.argv[1:]}}))';
  const result = JSON.parse((await run(options.python,['-c',code,...packages])).stdout);
  if (stage === 'geometry') {
    const directory = options.fr3dPath ?? path.join(dataRoot,'fr3d-python');
    result.fr3d_commit = (await run('git',['-C',directory,'rev-parse','HEAD'])).stdout.trim();
    result.fr3d_diff = sha256((await run('git',['-C',directory,'diff','HEAD','--','fr3d'])).stdout);
  }
  return result;
}

export class NucleicDatasetBuilder {
  constructor(options) {
    this.options = options;
    this.buildDir = path.join(dataRoot, 'builds', options.buildId);
    this.cacheDir = path.join(dataRoot, 'cache');
    this.scope = new OutputScope([dataRoot, assetsRoot]);
  }
  async initialize() {
    const file = path.join(this.buildDir, 'build.json');
    const requested = {limit: this.options.limit ?? null, only_ids: this.options.onlyIds ?? null};
    if (await exists(file)) {
      this.build = await readJson(file);
      // A resumed stage can omit population flags, but cannot change their meaning.
      if ((requested.limit || requested.only_ids) && JSON.stringify(requested) !== JSON.stringify(this.build.selection)) throw new Error('Build population is immutable; use a new build ID');
      this.options.limit = this.build.selection.limit;
      this.options.onlyIds = this.build.selection.only_ids;
    } else {
      this.build = {schema_version: 'rna-build-1', molecule_type: 'RNA', build_id: this.options.buildId,
        created_at: new Date().toISOString(), selection: requested, partial: Boolean(requested.limit || requested.only_ids),
        coordinate_cache_policy: 'immutable_content_hash_with_recorded_retrieval_time', stages: {}};
      await this.scope.json(file, this.build);
    }
  }
  async save() { await this.scope.json(path.join(this.buildDir, 'build.json'), this.build); }
  async read(relative) { return readJson(path.join(this.buildDir, relative)); }
  async output(relative, value) { return this.scope.json(path.join(this.buildDir, relative), value); }
  async accepted() { return (await this.read('eligibility/decisions.json')).filter(row => row.accepted); }

  async execute() {
    await this.initialize();
    const stages = this.options.stage ? [this.options.stage] : STAGES.slice(0, STAGES.indexOf(this.options.through) + 1);
    for (const stage of stages) await this.executeStage(stage);
    return this.build;
  }
  async executeStage(stage) {
    const index = STAGES.indexOf(stage);
    if (index && this.build.stages[STAGES[index - 1]]?.status !== 'complete') throw new Error(`Stage ${stage} requires completed ${STAGES[index - 1]}`);
    const files = await stageFiles(stage);
    const runtime = await runtimeFingerprint(stage,this.options);
    const previous = index ? this.build.stages[STAGES[index-1]] : null;
    const signature = sha256(JSON.stringify({source: await fingerprint(files), previous: previous ? {
      signature:previous.signature,outputs:previous.outputs?.map(output=>({path:output.path,sha256:output.sha256}))} : null,
      selection: this.build.selection, snapshot: this.options.snapshot ? sha256(await fs.readFile(this.options.snapshot)) : null,
      annotation_snapshot: stage === 'annotate' && this.options.annotationSnapshot ? sha256(await fs.readFile(this.options.annotationSnapshot)) : null,runtime}));
    const prior = this.build.stages[stage];
    if (this.options.resume && prior?.status === 'complete' && prior.signature === signature && !this.options.refreshCoordinates) {
      for (const output of prior.outputs ?? []) if (sha256(await fs.readFile(output.path)) !== output.sha256) throw new Error(`Checkpoint output changed: ${output.path}`);
      process.stderr.write(`resume ${stage}\n`); return;
    }
    for (const child of STAGES.slice(index)) delete this.build.stages[child];
    this.build.stages[stage] = {status: 'running', started_at: new Date().toISOString(), signature, runtime,
      source_files:files}; await this.save();
    process.stderr.write(`start ${stage}\n`);
    try {
      const outputs = await this[stage.replaceAll('-', '_')]();
      Object.assign(this.build.stages[stage], {status: 'complete', completed_at: new Date().toISOString(), outputs: outputs ?? []});
      await this.save(); process.stderr.write(`complete ${stage}\n`);
    } catch (error) {
      Object.assign(this.build.stages[stage], {status: 'failed', error: error.message, completed_at: new Date().toISOString()}); await this.save(); throw error;
    }
  }
  async discover() {
    await discover({...this.options, scope: this.scope, buildDir: this.buildDir});
    return this.hashOutputs(['discovery/source.json', 'discovery/candidates.json']);
  }
  async fetch() {
    const {ids} = await this.read('discovery/candidates.json');
    await fetchCoordinates(ids, {...this.options, scope: this.scope, cacheDir: this.cacheDir, buildDir: this.buildDir});
    return this.hashOutputs(['discovery/coordinates.json']);
  }
  async normalize() {
    const records = await this.read('discovery/coordinates.json');
    const normalizerHash = await fingerprint(['normalize_mmcif.py']);
    const normalized = await parallelMap(records, this.options.concurrency, async record => {
      const target = path.join(this.buildDir, 'identity', `${record.pdb_id}.json`);
      const checkpoint = target + '.checkpoint.json';
      const signature = sha256(record.sha256 + normalizerHash);
      try {
        if (await exists(checkpoint)) {
          const saved = await readJson(checkpoint);
          if (saved.signature === signature && sha256(await fs.readFile(target)) === saved.sha256) return saved;
        }
        await this.scope.resolve(target);
        await run(this.options.python, [path.join(here, 'normalize_mmcif.py'), record.path, '--out', target], {maxBuffer: 4 * 1024 * 1024});
        const entry = await readJson(target);
        if (entry.pdb_id !== record.pdb_id) throw new Error('Coordinate accession mismatch');
        const row = {pdb_id: record.pdb_id, path: target, status: 'available', signature, sha256: sha256(await fs.readFile(target)),
          canonical_rna: entry.eligibility.canonical_rna, reasons: entry.eligibility.reasons};
        await this.scope.json(checkpoint, row); return row;
      } catch (error) { return {pdb_id: record.pdb_id, status: 'normalization_failed', reasons: [error.message]}; }
    });
    await this.output('identity/index.json', normalized);
    if (normalized.some(row => row.status !== 'available')) throw new Error('Normalization failures retained in identity/index.json');
    return this.hashOutputs(['identity/index.json', ...normalized.map(row => path.relative(this.buildDir, row.path))]);
  }
  async annotate() {
    const {ids, all_ids} = await this.read('discovery/candidates.json');
    const annotations = await loadAnnotations(ids, {...this.options, snapshot: this.options.annotationSnapshot, scope: this.scope, buildDir: this.buildDir});
    const normalized = await this.read('identity/index.json');
    const decisions = [], entries = [], entities = [];
    for (const record of normalized) {
      const entry = await readJson(record.path), annotation = normalizeAnnotation(annotations.get(record.pdb_id), entry);
      const reasons = [...entry.eligibility.reasons];
      if (annotation.composition_conflict) reasons.push('nakb_composition_conflict');
      const accepted = reasons.length === 0;
      decisions.push({pdb_id: record.pdb_id, accepted, canonical_rna: entry.eligibility.canonical_rna, reasons,
        annotation_status: annotation.status, composition: annotation.composition, coordinate_status: record.status,
        source_sha256: entry.source_sha256, normalized_path: record.path, profiles: entry.eligibility.profiles});
      if (!accepted) continue;
      const profiles = {...entry.eligibility.profiles, all: true, relaxed: entry.eligibility.profiles.dna_compatible_relaxed_v1,
        conservative: entry.eligibility.profiles.dna_compatible_conservative_v1, mw100: entry.eligibility.profiles.dna_compatible_mw100_v1};
      entries.push({pdb_id: entry.pdb_id, ...entry.metadata, profiles, canonical_rna: true,
        model_count: entry.coordinate_policy.model_count, selected_model_id: entry.coordinate_policy.model_id,
        coordinate_policy: entry.coordinate_policy, residue_count: entry.residues.length, source_sha256: entry.source_sha256,
        annotation_status: annotation.status, rnaeq: annotation.rnaeq ?? null});
      for (const entity of entry.entities) entities.push({...entity, pdb_id: entry.pdb_id, sequence_length: entity.sequence?.length,
        ...annotation.entities.find(row => row.entity_id === entity.entity_id), functions: annotation.entities.find(row => row.entity_id === entity.entity_id)?.functions ?? [],
        subtypes: annotation.entities.find(row => row.entity_id === entity.entity_id)?.subtypes ?? [], structures: annotation.entities.find(row => row.entity_id === entity.entity_id)?.structures ?? []});
    }
    const selected = new Set(ids);
    for (const id of all_ids) if (!selected.has(id)) decisions.push({pdb_id: id, accepted: false, reasons: ['not_selected_partial_build'], coordinate_status: 'not_requested'});
    await this.output('eligibility/decisions.json', decisions);
    await this.output('tables/metadata.json', {entries, entities});
    return this.hashOutputs(['eligibility/decisions.json', 'tables/metadata.json', 'annotations/raw_batches.json', 'annotations/coverage.json']);
  }
  async computeEntries(stage, compute) {
    const records = await this.accepted(); let complete = 0;
    const outputs = await parallelMap(records, Math.min(this.options.concurrency,4), async record => {
      const target = `tables/${stage}/${record.pdb_id}.json`, checkpoint = `${target}.checkpoint.json`;
      try {
        const inputBytes = await fs.readFile(record.normalized_path);
        const signature = sha256(this.build.stages[stage].signature + sha256(inputBytes));
        if (await exists(path.join(this.buildDir,checkpoint))) {
          const saved = await this.read(checkpoint);
          if (saved.signature === signature && sha256(await fs.readFile(saved.output.path)) === saved.output.sha256) {
            complete++; if (complete % 50 === 0) process.stderr.write(`${stage} ${complete}/${records.length} (cached)\n`);
            return {...saved.output,pdb_id:record.pdb_id};
          }
        }
        const result = await compute(JSON.parse(inputBytes.toString()));
        const output = await this.output(target,result);
        await this.output(checkpoint,{signature,output});
        complete++; if (complete % 50 === 0) process.stderr.write(`${stage} ${complete}/${records.length}\n`);
        return {...output,pdb_id:record.pdb_id};
      } catch(error) { return {pdb_id:record.pdb_id,status:'failed',reason:error.message}; }
    });
    await this.output(`tables/${stage}/index.json`, outputs);
    if (outputs.some(row=>row.status==='failed')) throw new Error(`${stage} failures retained in tables/${stage}/index.json; retry reuses hash-matching entries`);
    return [...await this.hashOutputs([`tables/${stage}/index.json`]), ...outputs];
  }
  async residue() {
    const {computeResidueObservables} = await import('./residue_geometry.mjs');
    return this.computeEntries('residue', entry => computeResidueObservables(entry));
  }
  async geometry() {
    const {computeInteractionGeometry} = await import('./interactions.mjs');
    return this.computeEntries('geometry', entry => computeInteractionGeometry(entry, {python: this.options.python, fr3dPath: this.options.fr3dPath}));
  }
  async survey() {
    const {computeSurvey} = await import('./survey.mjs');
    return this.computeEntries('survey', async entry => {
      const geometry = await this.read(`tables/geometry/${entry.pdb_id}.json`);
      return computeSurvey(entry, {pairs: geometry.families?.base_pair ?? [], graph: geometry.graph});
    });
  }
  async assets() {
    const {buildAssets} = await import('./assets.mjs');
    return buildAssets({build: this.build, buildDir: this.buildDir, scope: this.scope, assetsRoot});
  }
  async validate_release() {
    const {validateRelease} = await import('./assets.mjs');
    const report = await validateRelease(path.join(assetsRoot, 'releases', this.build.build_id, 'manifest.json'));
    const output = await this.output('checks/release.json', report);
    if (!report.ok) throw new Error('Serialized release validation failed');
    // Activate only after checking the serialized files, and only within RNA assets.
    const descriptor = {schema_version: 'rna-explorer-1', molecule_type: 'RNA', build_id: this.build.build_id,
      manifest: `releases/${this.build.build_id}/manifest.json`, partial: this.build.partial};
    return [output, await this.scope.json(path.join(assetsRoot, 'manifest.json'), descriptor)];
  }
  async hashOutputs(relatives) {
    return Promise.all(relatives.map(async relative => { const file = path.join(this.buildDir, relative); const bytes = await fs.readFile(file); return {path: file, bytes: bytes.length, sha256: sha256(bytes)}; }));
  }
}

export async function main(argv = process.argv.slice(2)) {
  const options = parseArgs(argv);
  if (options.help) { process.stdout.write(HELP); return; }
  await new NucleicDatasetBuilder(options).execute();
}
if (process.argv[1] && path.resolve(process.argv[1]) === fileURLToPath(import.meta.url)) main().catch(error => { process.stderr.write(`${error.stack}\n`); process.exitCode = 1; });
