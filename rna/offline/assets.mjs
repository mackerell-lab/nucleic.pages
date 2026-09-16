import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import {gzipSync, gunzipSync} from 'node:zlib';
import {RESIDUE_PARAMETERS} from './parameter_registry.mjs';
import {TERM_REGISTRY,publicBaseGeometryConfig} from './survey_terms.mjs';
import {readJson, sha256} from './output_scope.mjs';

const here = path.dirname(fileURLToPath(import.meta.url));
const labels = {backbone:'Backbone Torsions',pseudo_torsion:'Pseudo Torsions',sugar_torsion:'Sugar Torsions',
  pucker:'Sugar Pucker',glycosidic_sugar_angles:'Glycosidic Sugar Angles',glycosidic_base_angles:'Glycosidic Base Angles',
  ribose_2oh:'Ribose 2′-OH Heavy Atoms',base_pair:'Base Pair',pair_quality:'Pair Quality',step:'Base-pair Step',
  helical:'Helical',step_position:'Step Position',same_strand:'Same Strand',helix_radius:'Helix Radius'};

function slimRow(row, parameters) {
  const allowed = ['id','pdb_id','entry_id','entity_id','label_asym_id','label_seq_id','auth_asym_id','auth_seq_id','comp_id',
    'model_id','altloc','pucker_class','sequence_context','pair_label','pair_type','interaction_family','family','residue1_id',
    'residue2_id','pair1_id','pair2_id','residue_ids','entity_ids','endpoint_entities','chain_ids','is_terminal','is_terminal_any','quality_flags'];
  const result = Object.fromEntries(allowed.filter(key => row[key] !== undefined).map(key => [key, row[key]]));
  result.pdb_id ??= row.entry_id;
  result.values = {}; result.statuses = {};
  for (const parameter of parameters) {
    const id = parameter.id ?? parameter.param_id, value = row.values?.[id];
    result.values[id] = Number.isFinite(value) ? value : null;
    result.statuses[id] = row.statuses?.[id] ?? (Number.isFinite(value) ? 'available' : 'missing');
  }
  return result;
}

/** Spool per-partition JSON so a worldwide survey is never retained as one array. */
class Partitions {
  constructor(scope, directory) { this.scope = scope; this.directory = directory; this.files = new Map(); this.coordinateParts=new Map(); }
  async append(key, rows) {
    if (!rows.length) return;
    if (!/^[a-zA-Z0-9_./-]+$/.test(key) || key.includes('..')) throw new Error('Unsafe partition key');
    let item = this.files.get(key);
    if (!item) {
      const file = await this.scope.resolve(path.join(this.directory, `${key}.json`));
      await fs.mkdir(path.dirname(file), {recursive:true});
      item = {file, count:0, entryIds:new Set(),handle:await fs.open(file, 'w')};
      await item.handle.write('['); this.files.set(key, item);
    }
    await item.handle.write((item.count ? ',' : '') + rows.map(row => JSON.stringify(row)).join(','));
    item.count += rows.length;
    for(const row of rows) if(row.pdb_id ?? row.entry_id) item.entryIds.add(row.pdb_id ?? row.entry_id);
  }
  async coordinates(group,rows) {
    let part=this.coordinateParts.get(group) ?? {index:0,count:0};
    while(rows.length) {
      const available=10000-part.count, chunk=rows.slice(0,available);
      const key=`survey/coordinates/${group}/${String(part.index).padStart(5,'0')}`;
      await this.append(key,chunk);part.count+=chunk.length;rows=rows.slice(chunk.length);
      if(part.count===10000) {
        const item=this.files.get(key);await item.handle.write(']\n');await item.handle.close();item.handle=null;
        part={index:part.index+1,count:0};
      }
    }
    this.coordinateParts.set(group,part);
  }
  async finish(releaseRoot) {
    const descriptors = new Map();
    for (const [key, item] of this.files) {
      if(item.handle) {await item.handle.write(']\n'); await item.handle.close();item.handle=null;}
      const raw = await fs.readFile(item.file), compressed = gzipSync(raw, {level:9});
      const relative = `${key}.json.gz`, output = await this.scope.write(path.join(releaseRoot, relative), compressed);
      descriptors.set(key, {path: relative,row_count:item.count,bytes:compressed.length,uncompressed_bytes:raw.length,sha256:output.sha256,
        ...(key.startsWith('survey/coordinates/') ? {entry_ids:[...item.entryIds].sort()} : {})});
    }
    return descriptors;
  }
}

export async function buildAssets({build, buildDir, scope, assetsRoot}) {
  const releaseRoot = path.join(assetsRoot, 'releases', build.build_id);
  const geometryParameters = await readJson(path.join(here, '../config/geometry_parameters_v1.json'));
  const allParameters = [...RESIDUE_PARAMETERS, ...geometryParameters];
  const definitions = new Map();
  for (const parameter of allParameters) {
    const id = parameter.family ?? parameter.family_id;
    if (!definitions.has(id)) definitions.set(id, []);
    definitions.get(id).push(parameter);
  }
  const metadata = await readJson(path.join(buildDir, 'tables/metadata.json'));
  // The browser entity facets describe RNA polymers; waters and associated
  // components remain in normalized evidence and the entry cleanliness policy.
  metadata.entities = metadata.entities.filter(entity=>entity.type==='polymer' && entity.polymer_type==='polyribonucleotide');
  const decisions = await readJson(path.join(buildDir, 'eligibility/decisions.json'));
  const discovery = await readJson(path.join(buildDir, 'discovery/candidates.json'));
  const partitions = new Partitions(scope, path.join(buildDir, 'pages_staging/partitions'));
  const capabilities = [], relationTypes = new Set();
  // Co-locate comparable method/profile entries to let browser metadata filters
  // skip coordinate partitions before network transfer.
  const assetEntries=[...metadata.entries].sort((a,b)=>`${a.method}|${a.profiles.relaxed}|${a.pdb_id}`.localeCompare(`${b.method}|${b.profiles.relaxed}|${b.pdb_id}`));
  for (const entry of assetEntries) {
    const residues = await readJson(path.join(buildDir, 'tables/residue', `${entry.pdb_id}.json`));
    const residueEntities = new Map(residues.map(row => [row.id,row.entity_id]));
    const withEndpoints = row => {
      const ids = row.residue_ids ?? [row.residue1_id,row.residue2_id].filter(Boolean);
      const entityIds = [...new Set(ids.map(id=>residueEntities.get(id)).filter(id=>id != null))];
      return {...row,pdb_id:entry.pdb_id,endpoint_entities:entityIds.map(entity_id=>({pdb_id:entry.pdb_id,entity_id}))};
    };
    for (const [family, parameters] of definitions) if (parameters[0].level === 'residue' || parameters[0].observation_level === 'residue')
      await partitions.append(`families/${family}`, residues.map(row => slimRow(row, parameters)));
    const geometry = await readJson(path.join(buildDir, 'tables/geometry', `${entry.pdb_id}.json`));
    capabilities.push(geometry.capabilities);
    for (const [family, rows] of Object.entries(geometry.families ?? {})) {
      const parameters = definitions.get(family);
      if (!parameters) throw new Error(`Unknown geometry family: ${family}`);
      await partitions.append(`families/${family}`, rows.map(row => slimRow(withEndpoints(row), parameters)));
    }
    const relations = geometry.relations ?? [];
    if (relations.length) { relationTypes.add('observations'); await partitions.append('relations/observations', relations); }
    await partitions.append('relations/interactions', geometry.interactions ?? []);
    if (geometry.interactions?.length) relationTypes.add('interactions');
    const survey = await readJson(path.join(buildDir, 'tables/survey', `${entry.pdb_id}.json`));
    const scalarGroups = Map.groupBy ? Map.groupBy(survey.scalars, row => row.term_id) : groupBy(survey.scalars, row => row.term_id);
    for (const [term, rows] of scalarGroups) await partitions.append(`survey/scalars/${term}`, rows.map(({term_label,survey_group,unit, ...row}) => row));
    const coordinateGroups = groupBy(survey.coordinates, row => `${row.anchor_frame}_${row.anchor_base}`);
    for (const [group, rows] of coordinateGroups) await partitions.coordinates(group, rows);
  }
  const descriptors = await partitions.finish(releaseRoot);
  const families = [];
  for (const [id, parameters] of definitions) {
    let descriptor = descriptors.get(`families/${id}`);
    if (!descriptor) {
      const bytes = gzipSync(Buffer.from('[]\n')), output = await scope.write(path.join(releaseRoot, `families/${id}.json.gz`), bytes);
      descriptor = {path:`families/${id}.json.gz`,row_count:0,sha256:output.sha256,bytes:bytes.length,uncompressed_bytes:3};
    }
    families.push({id,label:labels[id] ?? id,level:parameters[0].level ?? parameters[0].observation_level,parameters:parameters.map(parameter => ({
      ...parameter,id:parameter.id ?? parameter.param_id,label:parameter.label ?? parameter.display_name,
      period:parameter.period ?? null})),...descriptor});
  }
  const metadataOutput = await scope.json(path.join(releaseRoot, 'metadata.json'), metadata);
  const decisionBytes = gzipSync(Buffer.from(JSON.stringify(decisions))), decisionOutput = await scope.write(path.join(releaseRoot, 'provenance/decisions.json.gz'), decisionBytes);
  const scalarTerms = {}, coordinateGroups = {};
  for (const [key, descriptor] of descriptors) {
    if (key.startsWith('survey/scalars/')) scalarTerms[key.split('/').at(-1)] = descriptor;
    if (key.startsWith('survey/coordinates/')) {
      const group=key.split('/')[2];
      coordinateGroups[group] ??= {label:group.replaceAll('_',' '),row_count:0,partitions:[]};
      coordinateGroups[group].row_count+=descriptor.row_count;
      coordinateGroups[group].partitions.push(descriptor);
    }
  }
  const manifest = {schema_version:'rna-explorer-1',molecule_type:'RNA',build_id:build.build_id,generated_at:new Date().toISOString(),
    partial:build.partial,source:{name:'PDB experimental archive',retrieved_at:discovery.source_retrieved_at,
      candidate_count:discovery.source_count,processed_candidate_count:discovery.selected_count},
    counts:{entries:metadata.entries.length,entities:metadata.entities.length,residues:metadata.entries.reduce((n,row) => n + row.residue_count,0),
      excluded:decisions.filter(row => !row.accepted && !row.reasons.includes('not_selected_partial_build')).length},
    metadata:{path:'metadata.json',sha256:metadataOutput.sha256},families,
    relations:Object.fromEntries([...relationTypes].map(type => [type,descriptors.get(`relations/${type}`)])),
    survey:{terms:TERM_REGISTRY,opening_bins:publicBaseGeometryConfig().opening_bins,scalars:{terms:scalarTerms},coordinates:{groups:coordinateGroups}},
    coordinate_policy:metadata.entries[0]?.coordinate_policy ?? null,
    defaults:{family:'backbone',parameter:'chi',profile:'relaxed',method:'xray',max_resolution:3},
    profiles:{all:'all_associated_components',relaxed:'dna_compatible_relaxed_v1',conservative:'dna_compatible_conservative_v1',mw100:'dna_compatible_mw100_v1'},
    capabilities:{residue:true,geometry:capabilities.every(item => item?.base_pair === 'available'),survey:true,
      suites:false,representative_ife:false,modified_rna:false,assembly_expansion:false,hydrogen_orientation:false},
    limitations:['Canonical A/C/G/U only; modified polymers are retained in the exclusion ledger.',
      'One deposited model and a coherent residue conformer; no biological assembly expansion.',
      'NAKB missing annotations remain unknown; conflicting composition is excluded pending review.',
      'No DNA ABI or BI/BII/BIII labels; backbone suites and representative IFE selection are not implemented.',
      'Pair and step geometry describe supported FR3D interactions and selected stems, not every RNA contact.'],
    provenance:{decisions:{path:'provenance/decisions.json.gz',sha256:decisionOutput.sha256,row_count:decisions.length},
      build_stages:structuredClone(build.stages),coordinate_policy:metadata.entries[0]?.coordinate_policy ?? null}};
  await scope.json(path.join(releaseRoot, 'manifest.json'), manifest);
  const outputs = [...descriptors.values()].map(item => ({path:path.join(releaseRoot,item.path),sha256:item.sha256,bytes:item.bytes}));
  outputs.push(metadataOutput,decisionOutput,{path:path.join(releaseRoot,'manifest.json'),sha256:sha256(await fs.readFile(path.join(releaseRoot,'manifest.json')))});
  return outputs;
}

function groupBy(rows, key) { const groups = new Map(); for (const row of rows) { const id=key(row); if (!groups.has(id)) groups.set(id,[]); groups.get(id).push(row); } return groups; }

export async function validateRelease(manifestPath) {
  const manifest = await readJson(manifestPath), root = path.dirname(manifestPath), errors = [], checks = [];
  const load = async descriptor => {
    if (!descriptor?.path || path.isAbsolute(descriptor.path) || descriptor.path.split('/').includes('..')) throw new Error('Unsafe asset path');
    const bytes = await fs.readFile(path.join(root,descriptor.path));
    if (sha256(bytes) !== descriptor.sha256) throw new Error(`Asset hash mismatch: ${descriptor.path}`);
    return JSON.parse((descriptor.path.endsWith('.gz') ? gunzipSync(bytes) : bytes).toString());
  };
  const metadata = await load(manifest.metadata), entryIds = new Set(metadata.entries.map(row => row.pdb_id));
  if (entryIds.size !== metadata.entries.length || entryIds.size !== manifest.counts.entries) errors.push('Entry count or uniqueness');
  const residueIds = new Set();
  for (const family of manifest.families) {
    const rows = await load(family), ids = new Set();
    if (rows.length !== family.row_count) errors.push(`Row count: ${family.id}`);
    for (const row of rows) {
      if (!row.id || ids.has(row.id)) errors.push(`Duplicate identity: ${family.id}/${row.id}`);
      ids.add(row.id);
      if (!entryIds.has(row.pdb_id ?? row.entry_id)) errors.push(`Entry foreign key: ${family.id}/${row.id}`);
      for (const parameter of family.parameters) {
        const value = row.values?.[parameter.id], status = row.statuses?.[parameter.id];
        if (!(value === null || Number.isFinite(value)) || !status) errors.push(`Value/status: ${family.id}/${parameter.id}`);
        if (value === null && ['available','ok'].includes(status)) errors.push(`Null available value: ${family.id}/${parameter.id}`);
      }
      if (family.level === 'residue') residueIds.add(row.id);
    }
    checks.push({family:family.id,row_count:rows.length});
  }
  for (const [id, descriptor] of Object.entries(manifest.survey.scalars.terms)) {
    const rows = await load(descriptor), ids = new Set();
    if (rows.length !== descriptor.row_count) errors.push(`Survey count: ${id}`);
    for (const row of rows) {
      if (ids.has(row.id) || row.term_id !== id || !entryIds.has(row.pdb_id)) errors.push(`Survey identity: ${id}`);
      ids.add(row.id);
      if (row.residue_id && !residueIds.has(row.residue_id)) errors.push(`Survey residue key: ${id}`);
      if (!(row.value === null || Number.isFinite(row.value))) errors.push(`Survey value: ${id}`);
    }
  }
  for (const group of Object.values(manifest.survey.coordinates.groups)) {
    let count=0;
    for(const descriptor of group.partitions ?? [group]) {
      const rows=await load(descriptor);count+=rows.length;
      if(rows.length !== descriptor.row_count || (group.partitions && rows.length>10000)) errors.push('Coordinate partition row count');
      if(descriptor.entry_ids && rows.some(row=>!descriptor.entry_ids.includes(row.pdb_id))) errors.push('Coordinate partition entry index');
    }
    if(count!==group.row_count) errors.push('Coordinate group row count');
  }
  for (const descriptor of Object.values(manifest.relations)) await load(descriptor);
  const decisions = await load(manifest.provenance.decisions);
  if (decisions.length !== manifest.source.candidate_count || decisions.filter(row=>row.accepted).length !== entryIds.size) errors.push('Candidate ledger reconciliation');
  return {ok:!errors.length,checked_at:new Date().toISOString(),build_id:manifest.build_id,partial:manifest.partial,checks,errors};
}
