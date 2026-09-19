import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import {gzipSync, gunzipSync} from 'node:zlib';
import {RESIDUE_PARAMETERS} from './parameter_registry.mjs';
import {TERM_REGISTRY,publicBaseGeometryConfig} from './survey_terms.mjs';
import {readJson, sha256} from './output_scope.mjs';
import {encodeCoordinateRows, encodeFamilyRows, encodeInteractionRows, encodeSurveyRows, decodeCoordinateRows, decodeFamilyRows, decodeInteractionRows, decodeSurveyRows, COORDINATE_COLUMNAR_ENCODING, FAMILY_COLUMNAR_ENCODING, INTERACTION_COLUMNAR_ENCODING, SURVEY_COLUMNAR_ENCODING} from '../core/survey-codec.js';
import {SHARED_SURVEY_ENCODING, expandSharedSurveyColumns, verifySharedColumn} from '../core/shared-survey-codec.js';
import {BUNDLED_SURVEY_ENCODING, expandBundledSurveyColumns, verifySurveyBundle} from '../core/bundled-survey-codec.js';

const here = path.dirname(fileURLToPath(import.meta.url));
const labels = {backbone:'Backbone Torsions',pseudo_torsion:'Pseudo Torsions',sugar_torsion:'Sugar Torsions',
  pucker:'Sugar Pucker',glycosidic_sugar_angles:'Glycosidic Sugar Angles',glycosidic_base_angles:'Glycosidic Base Angles',
  ribose_2oh:'Ribose 2′-OH Heavy Atoms',base_pair:'Base Pair',pair_quality:'Pair Quality',step:'Base-pair Step',
  helical:'Helical',step_position:'Step Position',same_strand:'Same Strand',helix_radius:'Helix Radius'};

function slimRow(row, parameters) {
  const allowed = ['id','residue_id','pdb_id','entry_id','entity_id','label_asym_id','label_seq_id','auth_asym_id','auth_seq_id','insertion_code','comp_id',
    'model_id','altloc','pucker_class','pucker_classes','sequence_context','pair_label','pair_type','interaction_family','family','residue1_id',
    'residue2_id','pair1_id','pair2_id','residue_ids','entity_ids','endpoint_entities','chain_ids','is_terminal','is_terminal_any','quality_flags',
    'step_label','frame_convention','atom_roles','near','alternative','stem_eligible','topology','is_terminal_5p','is_terminal_3p'];
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

function publicCoordinatePolicy(policy) {
  if (!policy) return null;
  const {model_id,model_ids,model_count,...rules}=policy;
  return {...rules,selected_model_location:'metadata.entries[].selected_model_id',
    model_inventory_location:'metadata.entries[].coordinate_policy'};
}

function publicStages(stages) {
  return Object.fromEntries(Object.entries(stages ?? {}).map(([name,stage])=>[name,{
    status:stage.status,started_at:stage.started_at,completed_at:stage.completed_at,
    signature:stage.signature,runtime:stage.runtime ?? null,output_count:stage.outputs?.length ?? 0,
    // Hash content identities, never publish workstation paths or output lists.
    outputs_sha256:sha256(JSON.stringify((stage.outputs ?? []).map(({sha256,bytes})=>({sha256,bytes})))),
  }]));
}

/** Spool per-partition JSON so a worldwide survey is never retained as one array. */
class Partitions {
  constructor(scope, directory, {columnar = false, buildId = null} = {}) { this.scope = scope; this.directory = directory; this.files = new Map(); this.coordinateParts=new Map(); this.columnar = columnar; this.buildId = buildId; }
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
      const raw = await fs.readFile(item.file);
      let payload = raw, encoding = null;
      if (this.columnar && key.startsWith('survey/scalars/')) {
        payload = Buffer.from(JSON.stringify(encodeSurveyRows(JSON.parse(raw), this.buildId)));
        encoding = SURVEY_COLUMNAR_ENCODING;
      } else if (this.columnar && key.startsWith('survey/coordinates/')) {
        payload = Buffer.from(JSON.stringify(encodeCoordinateRows(JSON.parse(raw), this.buildId)));
        encoding = COORDINATE_COLUMNAR_ENCODING;
      } else if (this.columnar && key.startsWith('families/')) {
        payload = Buffer.from(JSON.stringify(encodeFamilyRows(JSON.parse(raw), this.buildId)));
        encoding = FAMILY_COLUMNAR_ENCODING;
      } else if (this.columnar && key === 'relations/interactions') {
        payload = Buffer.from(JSON.stringify(encodeInteractionRows(JSON.parse(raw), this.buildId)));
        encoding = INTERACTION_COLUMNAR_ENCODING;
      }
      const compressed = gzipSync(payload, {level:9});
      const relative = `${key}.json.gz`, output = await this.scope.write(path.join(releaseRoot, relative), compressed);
      descriptors.set(key, {path: relative,row_count:item.count,bytes:compressed.length,uncompressed_bytes:payload.length,sha256:output.sha256,
        ...(encoding ? {encoding} : {}),
        ...(key.startsWith('survey/coordinates/') ? {entry_ids:[...item.entryIds].sort()} : {})});
    }
    return descriptors;
  }
}

export async function buildAssets({build, buildDir, scope, assetsRoot}) {
  try {
    const active = await readJson(path.join(assetsRoot, 'manifest.json'));
    if (active.build_id === build.build_id) throw new Error('Active RNA releases are immutable; choose a new build ID');
  } catch (error) {
    if (error.code !== 'ENOENT') throw error;
  }
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
  const partitions = new Partitions(scope, path.join(buildDir, 'pages_staging/partitions'), {columnar: build.partial === false, buildId: build.build_id});
  const capabilities = [], relationTypes = new Set();
  // Co-locate comparable method/profile entries to let browser metadata filters
  // skip coordinate partitions before network transfer.
  const assetEntries=[...metadata.entries].sort((a,b)=>`${a.method}|${a.profiles.relaxed}|${a.pdb_id}`.localeCompare(`${b.method}|${b.profiles.relaxed}|${b.pdb_id}`));
  for (const entry of assetEntries) {
    const normalized=await readJson(path.join(buildDir,'identity',`${entry.pdb_id}.json`));
    const sourceResidues=new Map(normalized.residues.map(row=>[row.id,row]));
    const sourceIdentity=id=>{
      const source=sourceResidues.get(id);
      if(!source) throw new Error(`Missing normalized residue identity: ${id}`);
      return {residue_id:id,entity_id:source.entity_id,label_asym_id:source.label_asym_id,label_seq_id:source.label_seq_id,
        auth_asym_id:source.auth_asym_id,auth_seq_id:source.auth_seq_id,insertion_code:source.ins_code ?? null,
        altloc:source.altloc ?? '',model_id:source.model_id};
    };
    const residues = (await readJson(path.join(buildDir, 'tables/residue', `${entry.pdb_id}.json`)))
      .map(row=>({...row,...sourceIdentity(row.id)}));
    const residueIndex = new Map(residues.map(row => [row.id,row]));
    const geometry = await readJson(path.join(buildDir, 'tables/geometry', `${entry.pdb_id}.json`));
    const pairIndex = new Map((geometry.families?.base_pair ?? []).map(row=>[row.id,row]));
    const withEndpoints = row => {
      const pair=row.pair_id ? pairIndex.get(row.pair_id) : null;
      if(row.pair_id && !pair) throw new Error(`Unknown survey pair: ${row.pair_id}`);
      const ids = pair?.residue_ids ?? row.residue_ids ?? [row.residue1_id,row.residue2_id,
        row.residue_id ?? row.target_residue_id].filter(Boolean);
      if(!ids.length || ids.some(id=>!residueIndex.has(id))) throw new Error(`Unknown observation endpoints: ${row.id}`);
      const endpoints=ids.map(id=>residueIndex.get(id));
      const entityIds = [...new Set(endpoints.map(residue=>residue.entity_id))];
      return {...row,...(endpoints.length===1 ? sourceIdentity(ids[0]) : {}),
        pdb_id:entry.pdb_id,model_id:row.model_id ?? entry.selected_model_id,
        residue_ids:ids,is_terminal_any:endpoints.some(residue=>residue.is_terminal_any===true),
        pucker_classes:endpoints.map(residue=>residue.pucker_class ?? null),
        ...(endpoints.length===1 ? {pucker_class:endpoints[0].pucker_class ?? null} : {}),
        endpoint_entities:entityIds.map(entity_id=>({pdb_id:entry.pdb_id,entity_id}))};
    };
    for (const [family, parameters] of definitions) if (parameters[0].level === 'residue' || parameters[0].observation_level === 'residue')
      await partitions.append(`families/${family}`, residues.map(row => slimRow(row, parameters)));
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
    for (const [term, rows] of scalarGroups) await partitions.append(`survey/scalars/${term}`, rows.map(({term_label,survey_group,unit, ...row}) => withEndpoints(row)));
    const coordinateGroups = groupBy(survey.coordinates, row => `${row.anchor_frame}_${row.anchor_base}`);
    for (const [group, rows] of coordinateGroups) await partitions.coordinates(group, rows.map(withEndpoints));
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
  const metadataRaw=Buffer.from(JSON.stringify(metadata)),metadataBytes=gzipSync(metadataRaw,{level:9});
  const metadataOutput = await scope.write(path.join(releaseRoot, 'metadata.json.gz'), metadataBytes);
  const publicDecisions=decisions.map(({normalized_path,...record})=>record);
  const decisionBytes = gzipSync(Buffer.from(JSON.stringify(publicDecisions))), decisionOutput = await scope.write(path.join(releaseRoot, 'provenance/decisions.json.gz'), decisionBytes);
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
    metadata:{path:'metadata.json.gz',sha256:metadataOutput.sha256,bytes:metadataBytes.length,uncompressed_bytes:metadataRaw.length},families,
    relations:Object.fromEntries([...relationTypes].map(type => [type,descriptors.get(`relations/${type}`)])),
    survey:{terms:TERM_REGISTRY,opening_bins:publicBaseGeometryConfig().opening_bins,scalars:{terms:scalarTerms},coordinates:{groups:coordinateGroups}},
    coordinate_policy:publicCoordinatePolicy(metadata.entries[0]?.coordinate_policy),
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
      build_stages:publicStages(build.stages),coordinate_policy:publicCoordinatePolicy(metadata.entries[0]?.coordinate_policy)}};
  await scope.json(path.join(releaseRoot, 'manifest.json'), manifest);
  const outputs = [...descriptors.values()].map(item => ({path:path.join(releaseRoot,item.path),sha256:item.sha256,bytes:item.bytes}));
  outputs.push(metadataOutput,decisionOutput,{path:path.join(releaseRoot,'manifest.json'),sha256:sha256(await fs.readFile(path.join(releaseRoot,'manifest.json')))});
  return outputs;
}

function groupBy(rows, key) { const groups = new Map(); for (const row of rows) { const id=key(row); if (!groups.has(id)) groups.set(id,[]); groups.get(id).push(row); } return groups; }

export async function validateRelease(manifestPath) {
  const manifest = await readJson(manifestPath), root = path.dirname(manifestPath), errors = [], checks = [];
  const realRoot = await fs.realpath(root), fullRelease = manifest.partial === false;
  const usedBundles = new Set(), bundleRegistry = manifest.survey?.bundles;
  if (bundleRegistry !== undefined && (!bundleRegistry || typeof bundleRegistry !== 'object' || Array.isArray(bundleRegistry))) {
    throw new Error('Invalid Survey bundle registry');
  }
  const readAsset = async (descriptor, {requireHash = true} = {}) => {
    if (typeof descriptor?.path !== 'string' || !descriptor.path || path.isAbsolute(descriptor.path)
        || descriptor.path.includes('\\') || descriptor.path.split('/').some(part => !part || part === '.' || part === '..')) {
      throw new Error('Unsafe asset path');
    }
    const file = await fs.realpath(path.join(realRoot, descriptor.path));
    const relative = path.relative(realRoot, file);
    if (relative === '..' || relative.startsWith(`..${path.sep}`) || path.isAbsolute(relative)) throw new Error('Asset path escapes release root');
    const bytes = await fs.readFile(file);
    if ((requireHash || descriptor.sha256 !== undefined) && sha256(bytes) !== descriptor.sha256) throw new Error(`Asset hash mismatch: ${descriptor.path}`);
    if (descriptor.bytes !== undefined && (!Number.isSafeInteger(descriptor.bytes) || descriptor.bytes !== bytes.length)) {
      throw new Error(`Asset compressed byte size mismatch: ${descriptor.path}`);
    }
    const raw = descriptor.path.endsWith('.gz') ? gunzipSync(bytes) : bytes;
    if (descriptor.uncompressed_bytes !== undefined && (!Number.isSafeInteger(descriptor.uncompressed_bytes) || descriptor.uncompressed_bytes !== raw.length)) {
      throw new Error(`Asset uncompressed byte size mismatch: ${descriptor.path}`);
    }
    const data = JSON.parse(raw.toString());
    if (data && Object.hasOwn(data, 'build_id') && data.build_id !== manifest.build_id) throw new Error(`Asset build ID mismatch: ${descriptor.path}`);
    return data;
  };
  const readBundle = async reference => {
    if (!/^[a-f0-9]{64}$/.test(reference)) throw new Error('Invalid Survey bundle content hash');
    const fixedPath = `survey/bundles/${reference}.json.gz`;
    const descriptor = bundleRegistry && Object.hasOwn(bundleRegistry, reference) ? bundleRegistry[reference] : null;
    if (!descriptor && fullRelease) throw new Error(`Unregistered Survey bundle: ${reference}`);
    if (descriptor && (descriptor.path !== fixedPath || descriptor.content_sha256 !== reference
        || (fullRelease && (!Number.isSafeInteger(descriptor.bytes) || !Number.isSafeInteger(descriptor.uncompressed_bytes))))) {
      throw new Error(`Invalid Survey bundle registry descriptor: ${reference}`);
    }
    // Legacy partial fixtures may lack inventory descriptors; their content hash
    // still authenticates the bundle from the same fixed, contained path.
    const payload = await readAsset(descriptor ?? {path: fixedPath}, {requireHash: Boolean(descriptor)});
    return (await verifySurveyBundle(reference, payload)).bundle;
  };
  const encodings = new Set([BUNDLED_SURVEY_ENCODING, SHARED_SURVEY_ENCODING, SURVEY_COLUMNAR_ENCODING,
    COORDINATE_COLUMNAR_ENCODING, FAMILY_COLUMNAR_ENCODING, INTERACTION_COLUMNAR_ENCODING]);
  const load = async descriptor => {
    const data = await readAsset(descriptor);
    if (descriptor.encoding !== undefined || encodings.has(data?.encoding)) {
      if (descriptor.encoding !== data?.encoding) throw new Error(`Asset encoding mismatch: ${descriptor.path}`);
      if (!encodings.has(descriptor.encoding)) throw new Error(`Unsupported asset encoding: ${descriptor.encoding}`);
      if (fullRelease && data.build_id !== manifest.build_id) throw new Error(`Asset build ID mismatch: ${descriptor.path}`);
    }
    if (descriptor.encoding === BUNDLED_SURVEY_ENCODING) {
      const expanded = await expandBundledSurveyColumns(data, async reference => {
        usedBundles.add(reference);
        return readBundle(reference);
      });
      return decodeSurveyRows(expanded);
    }
    if (descriptor.encoding === SHARED_SURVEY_ENCODING) {
      const expanded = await expandSharedSurveyColumns(data, async reference => {
        if (!/^[a-f0-9]{64}$/.test(reference)) throw new Error('Invalid shared Survey content hash');
        const values = await readAsset({path:`survey/columns/${reference}.json.gz`}, {requireHash:false});
        return verifySharedColumn(reference, values);
      });
      return decodeSurveyRows(expanded);
    }
    if (descriptor.encoding === SURVEY_COLUMNAR_ENCODING) return decodeSurveyRows(data);
    if (descriptor.encoding === COORDINATE_COLUMNAR_ENCODING) return decodeCoordinateRows(data);
    if (descriptor.encoding === FAMILY_COLUMNAR_ENCODING) return decodeFamilyRows(data);
    if (descriptor.encoding === INTERACTION_COLUMNAR_ENCODING) return decodeInteractionRows(data);
    return data;
  };
  const metadata = await load(manifest.metadata), entryIds = new Set(metadata.entries.map(row => row.pdb_id));
  if (entryIds.size !== metadata.entries.length || entryIds.size !== manifest.counts.entries) errors.push('Entry count or uniqueness');
  const entryMap=new Map(metadata.entries.map(row=>[row.pdb_id,row]));
  const entityIds=new Set(metadata.entities.map(row=>`${row.pdb_id}|${row.entity_id}`));
  if(entityIds.size!==metadata.entities.length || entityIds.size!==manifest.counts.entities) errors.push('Entity count or uniqueness');
  if(metadata.entities.some(row=>!entryIds.has(row.pdb_id))) errors.push('Entity entry foreign key');
  const residueIds = new Set(),residueMap=new Map(),pairMap=new Map(),stepMap=new Map();
  const observationRefs=new Map();
  const validValue=(value,status)=> (value===null || Number.isFinite(value)) && typeof status==='string' && status.length>0 &&
    (['available','ok'].includes(status) ? Number.isFinite(value) : value===null);
  const ownership=(row,label)=>{
    const pdb=row.pdb_id ?? row.entry_id,entry=entryMap.get(pdb);
    if(!entry) errors.push(`Entry foreign key: ${label}`);
    if(row.entity_id!=null && !entityIds.has(`${pdb}|${row.entity_id}`)) errors.push(`Entity foreign key: ${label}`);
    for(const endpoint of row.endpoint_entities ?? []) if(endpoint.pdb_id!==pdb || !entityIds.has(`${endpoint.pdb_id}|${endpoint.entity_id}`)) errors.push(`Endpoint entity foreign key: ${label}`);
    if(row.model_id!=null && entry && String(row.model_id)!==String(entry.selected_model_id)) errors.push(`Selected model mismatch: ${label}`);
  };
  for (const family of manifest.families) {
    const rows = await load(family), ids = new Set();
    if (rows.length !== family.row_count) errors.push(`Row count: ${family.id}`);
    for (const row of rows) {
      if (!row.id || ids.has(row.id)) errors.push(`Duplicate identity: ${family.id}/${row.id}`);
      ids.add(row.id);
      ownership(row,`${family.id}/${row.id}`);
      for (const parameter of family.parameters) {
        const value = row.values?.[parameter.id], status = row.statuses?.[parameter.id];
        if (!validValue(value,status)) errors.push(`Value/status: ${family.id}/${parameter.id}`);
      }
      if (family.level === 'residue') {residueIds.add(row.id);residueMap.set(row.id,{pdb_id:row.pdb_id,entity_id:row.entity_id,is_terminal_any:row.is_terminal_any});}
      else {
        const reference={id:row.id,pdb_id:row.pdb_id,residue_ids:row.residue_ids,pair1_id:row.pair1_id,pair2_id:row.pair2_id,
          endpoint_entities:row.endpoint_entities,is_terminal_any:row.is_terminal_any};
        if(observationRefs.has(row.id) && JSON.stringify(observationRefs.get(row.id))!==JSON.stringify(reference)) errors.push(`Observation identity differs across families: ${row.id}`);
        else observationRefs.set(row.id,reference);
      }
      if(family.id==='base_pair') pairMap.set(row.id,row.residue_ids);
      if(family.id==='step') {stepMap.set(row.id,{pairs:[row.pair1_id,row.pair2_id],residue_ids:row.residue_ids});if(!row.step_label) errors.push(`Missing step context: ${row.id}`);}
    }
    checks.push({family:family.id,row_count:rows.length});
  }
  if(residueIds.size!==manifest.counts.residues) errors.push('Residue count');
  const endpoints=(row,label)=>{
    const ids=row.residue_ids ?? [row.residue_id ?? row.target_residue_id].filter(Boolean);
    const targets=ids.map(id=>residueMap.get(id));
    if(!ids.length || targets.some(residue=>!residue || residue.pdb_id!==row.pdb_id)) {errors.push(`Residue foreign key: ${label}`);return;}
    const expected=[...new Set(targets.map(residue=>`${row.pdb_id}|${residue.entity_id}`))].sort();
    const actual=(row.endpoint_entities ?? []).map(endpoint=>`${endpoint.pdb_id}|${endpoint.entity_id}`).sort();
    if(JSON.stringify(expected)!==JSON.stringify(actual)) errors.push(`Endpoint ownership: ${label}`);
    if(row.is_terminal_any!==targets.some(residue=>residue.is_terminal_any===true)) errors.push(`Terminal endpoint flag: ${label}`);
    if(row.pair_id && (!pairMap.has(row.pair_id) || JSON.stringify(ids)!==JSON.stringify(pairMap.get(row.pair_id)))) errors.push(`Pair foreign key or endpoint order: ${label}`);
  };
  for(const row of observationRefs.values()) {
    endpoints(row,row.id);
    for(const key of ['pair1_id','pair2_id']) if(row[key] && !pairMap.has(row[key])) errors.push(`Step pair foreign key: ${row.id}`);
  }
  for (const [id, descriptor] of Object.entries(manifest.survey.scalars.terms)) {
    const rows = await load(descriptor), ids = new Set();
    if (rows.length !== descriptor.row_count) errors.push(`Survey count: ${id}`);
    for (const row of rows) {
      if (ids.has(row.id) || row.term_id !== id || !entryIds.has(row.pdb_id)) errors.push(`Survey identity: ${id}`);
      ids.add(row.id);
      ownership(row,`survey/${id}`);endpoints(row,`survey/${id}`);
      if (row.residue_id && !residueIds.has(row.residue_id)) errors.push(`Survey residue key: ${id}`);
      if (!validValue(row.value,row.status)) errors.push(`Survey value/status: ${id}`);
    }
  }
  for (const reference of Object.keys(bundleRegistry ?? {})) {
    if (!usedBundles.has(reference)) {
      await readBundle(reference);
      if (fullRelease) throw new Error(`Unreferenced Survey bundle: ${reference}`);
    }
  }
  for (const group of Object.values(manifest.survey.coordinates.groups)) {
    let count=0;const ids=new Set();
    for(const descriptor of group.partitions ?? [group]) {
      const rows=await load(descriptor);count+=rows.length;
      if(rows.length !== descriptor.row_count || (group.partitions && rows.length>10000)) errors.push('Coordinate partition row count');
      if(descriptor.entry_ids && rows.some(row=>!descriptor.entry_ids.includes(row.pdb_id))) errors.push('Coordinate partition entry index');
      for(const row of rows) {
        if(!row.id || ids.has(row.id)) errors.push('Coordinate identity');ids.add(row.id);
        ownership(row,`coordinate/${row.id}`);endpoints(row,`coordinate/${row.id}`);
        if(!['x','y','z'].every(axis=>Number.isFinite(row[axis])) || !['available','ok'].includes(row.status)) errors.push(`Coordinate value/status: ${row.id}`);
        for(const key of ['anchor_residue_id','target_residue_id']) if(!residueMap.has(row[key]) || residueMap.get(row[key]).pdb_id!==row.pdb_id) errors.push(`Coordinate residue foreign key: ${row.id}`);
        if(row.residue_ids && (!row.residue_ids.includes(row.anchor_residue_id) || !row.residue_ids.includes(row.target_residue_id))) errors.push(`Coordinate endpoint membership: ${row.id}`);
      }
    }
    if(count!==group.row_count) errors.push('Coordinate group row count');
  }
  for (const [type,descriptor] of Object.entries(manifest.relations)) {
    const rows=await load(descriptor),ids=new Set();
    if(rows.length!==descriptor.row_count) errors.push(`Relation count: ${type}`);
    for(const row of rows) {
      if(!row.id || ids.has(row.id)) errors.push(`Relation identity: ${type}`);ids.add(row.id);
      if(type==='interactions') {
        const a=residueMap.get(row.residue1_id),b=residueMap.get(row.residue2_id);
        if(!a || !b || a.pdb_id!==b.pdb_id) errors.push(`Interaction residue foreign key: ${row.id}`);
      } else if(row.kind==='pair_residue') {
        if(pairMap.get(row.pair_id)?.[row.side-1]!==row.residue_id || !residueMap.has(row.residue_id)) errors.push(`Pair relation endpoint: ${row.id}`);
      } else if(row.kind==='step_pair') {
        if(stepMap.get(row.step_id)?.pairs[row.side-1]!==row.pair_id || !pairMap.has(row.pair_id)) errors.push(`Step pair relation endpoint: ${row.id}`);
      } else if(row.kind==='step_residue') {
        if(stepMap.get(row.step_id)?.residue_ids[row.side-1]!==row.residue_id || !residueMap.has(row.residue_id)) errors.push(`Step residue relation endpoint: ${row.id}`);
      } else errors.push(`Unsupported relation kind: ${row.kind}`);
    }
  }
  const decisions = await load(manifest.provenance.decisions);
  if (decisions.length !== manifest.source.candidate_count || decisions.filter(row=>row.accepted).length !== entryIds.size) errors.push('Candidate ledger reconciliation');
  if(new Set(decisions.map(row=>row.pdb_id)).size!==decisions.length || decisions.some(row=>row.accepted!==entryIds.has(row.pdb_id))) errors.push('Candidate ledger identities');
  return {ok:!errors.length,checked_at:new Date().toISOString(),build_id:manifest.build_id,partial:manifest.partial,checks,errors};
}
