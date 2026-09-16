import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import os from 'node:os';
import path from 'node:path';
import {gzipSync,gunzipSync} from 'node:zlib';
import {buildAssets,validateRelease} from '../offline/assets.mjs';
import {OutputScope,sha256,readJson} from '../offline/output_scope.mjs';
import {computeResidueObservables} from '../offline/residue_geometry.mjs';
import {computeGeometry,geometryFamilies} from '../offline/geometry_adapter.mjs';
import {computeSurvey} from '../offline/survey.mjs';

test('real RNA numerical output retains release identities and rejects rehashed scientific corruption',async t=>{
  const reference=await readJson(new URL('./reference/rna_x3dna_reference.json',import.meta.url));
  const entry=structuredClone(reference.entries.find(row=>row.pdb_id==='1sdr'));
  entry.pdb_id='1SDR';
  for(const residue of entry.residues) Object.assign(residue,{pdb_id:entry.pdb_id,
    entity_id:['A','C'].includes(residue.label_asym_id)?'1':'2',model_id:'30'});
  const residues=computeResidueObservables(entry);
  const geometry=computeGeometry(entry,reference.graphs['1sdr']);
  const survey=computeSurvey(entry,{pairs:geometry.pairs});
  assert.ok(geometry.steps.length>0 && survey.coordinates.length>0);
  const directory=await fs.mkdtemp(path.join(os.tmpdir(),'rna-release-integrity-'));
  t.after(()=>fs.rm(directory,{recursive:true,force:true}));
  const buildDir=path.join(directory,'build'),assetsRoot=path.join(directory,'assets'),scope=new OutputScope([directory]);
  const write=(relative,value)=>scope.json(path.join(buildDir,relative),value);
  const policy={id:'single_deposited_model_v1',scope:'deposited_asymmetric_unit',model_id:'30',model_ids:['30','31'],model_count:2};
  await write('tables/metadata.json',{entries:[{pdb_id:entry.pdb_id,method:'X-RAY DIFFRACTION',
    profiles:{all:true,relaxed:true},residue_count:residues.length,selected_model_id:'30',coordinate_policy:policy}],
    entities:['1','2'].map(entity_id=>({pdb_id:entry.pdb_id,entity_id,type:'polymer',polymer_type:'polyribonucleotide',functions:[`entity_${entity_id}`]}))});
  await write('eligibility/decisions.json',[{pdb_id:entry.pdb_id,accepted:true,reasons:[],normalized_path:'/home/private/identity.json'}]);
  await write('discovery/candidates.json',{source_count:1,selected_count:1});
  await write('tables/residue/1SDR.json',residues);
  await write('tables/geometry/1SDR.json',{families:geometryFamilies(geometry),relations:geometry.relations,
    interactions:reference.graphs['1sdr'].edges,capabilities:{base_pair:'available'}});
  await write('tables/survey/1SDR.json',survey);
  const build={build_id:'release-integrity',partial:true,stages:{geometry:{status:'complete',
    signature:'source-signature',runtime:{python:'3',packages:{numpy:'2'}},
    outputs:[{path:'/home/private/geometry.json',bytes:12,sha256:'f'.repeat(64)}]}}};
  await buildAssets({build,buildDir,scope,assetsRoot});
  const root=path.join(assetsRoot,'releases',build.build_id),manifestPath=path.join(root,'manifest.json');
  const manifest=await readJson(manifestPath);
  const load=async descriptor=>JSON.parse(gunzipSync(await fs.readFile(path.join(root,descriptor.path))));
  const valid=await validateRelease(manifestPath);
  assert.equal(valid.ok,true,JSON.stringify(valid.errors));
  assert.equal(manifest.coordinate_policy.model_id,undefined);
  assert.equal(manifest.coordinate_policy.model_ids,undefined);
  assert.equal(manifest.coordinate_policy.model_count,undefined);
  assert.equal(manifest.provenance.build_stages.geometry.output_count,1);
  assert.match(manifest.provenance.build_stages.geometry.outputs_sha256,/^[a-f0-9]{64}$/);
  assert.equal(JSON.stringify(manifest).includes('/home/'),false);
  assert.equal(JSON.stringify(await load(manifest.provenance.decisions)).includes('/home/'),false);
  assert.equal((await load(manifest.metadata)).entries[0].selected_model_id,'30');
  const stepRows=await load(manifest.families.find(row=>row.id==='step'));
  const pairRows=await load(manifest.families.find(row=>row.id==='base_pair'));
  const residueIndex=new Map(residues.map(row=>[row.id,row]));
  for(const row of [...stepRows,...pairRows]) {
    assert.equal(row.model_id,'30');
    assert.equal(row.is_terminal_any,row.residue_ids.some(id=>residueIndex.get(id).is_terminal_any));
    assert.ok(row.frame_convention);
    assert.deepEqual(row.endpoint_entities.map(entity=>entity.entity_id).sort(),['1','2']);
  }
  assert.ok(pairRows.some(row=>row.is_terminal_any));
  for(const row of stepRows) assert.equal(row.step_label,geometry.steps.find(step=>step.id===row.id).step_label);
  for(const row of pairRows) assert.deepEqual(row.atom_roles,geometry.pairs.find(pair=>pair.id===row.id).atom_roles);
  const pairDescriptor=manifest.survey.scalars.terms.same_pair_a_n6__u_o4;
  for(const row of await load(pairDescriptor)) {
    assert.deepEqual(row.endpoint_entities.map(entity=>entity.entity_id).sort(),['1','2']);
    assert.equal(row.is_terminal_any,row.residue_ids.some(id=>residueIndex.get(id).is_terminal_any));
  }
  const coordinateDescriptor=Object.values(manifest.survey.coordinates.groups)[0].partitions[0];
  const pairCoordinateGroup=Object.values(manifest.survey.coordinates.groups).find(group=>group.label.includes('cytosine'));
  assert.ok(pairCoordinateGroup);
  for(const row of await load(pairCoordinateGroup.partitions[0])) {
    assert.equal(row.residue_ids.length,2);
    assert.deepEqual(row.endpoint_entities.map(entity=>entity.entity_id).sort(),['1','2']);
    assert.ok(row.residue_ids.includes(row.target_residue_id) && row.residue_ids.includes(row.anchor_residue_id));
  }
  const familyDescriptor=manifest.families.find(row=>row.id==='step');
  const relationDescriptor=manifest.relations.observations;
  const cases=[
    ['nonfinite coordinate',coordinateDescriptor,rows=>{rows[0].x=null;},/Coordinate value\/status/],
    ['coordinate endpoint',coordinateDescriptor,rows=>{rows[0].anchor_residue_id='missing';},/Coordinate residue foreign key/],
    ['survey null available',pairDescriptor,rows=>{rows[0].value=null;rows[0].status='ok';},/Survey value\/status/],
    ['step pair reference',familyDescriptor,rows=>{rows[0].pair1_id='missing';},/Step pair foreign key/],
    ['relation reference',relationDescriptor,rows=>{rows[0].residue_id='missing';},/relation endpoint/],
    ['endpoint ownership',pairDescriptor,rows=>{rows[0].endpoint_entities=[];},/Endpoint ownership/],
  ];
  for(const [name,descriptor,mutate,pattern] of cases) await t.test(name,async()=>{
    const file=path.join(root,descriptor.path),original=await fs.readFile(file),hash=descriptor.sha256;
    const rows=JSON.parse(gunzipSync(original));mutate(rows);
    const changed=gzipSync(Buffer.from(JSON.stringify(rows)));
    await fs.writeFile(file,changed);descriptor.sha256=sha256(changed);
    await fs.writeFile(manifestPath,JSON.stringify(manifest));
    const result=await validateRelease(manifestPath);
    assert.equal(result.ok,false);assert.ok(result.errors.some(error=>pattern.test(error)),JSON.stringify(result.errors));
    await fs.writeFile(file,original);descriptor.sha256=hash;
    await fs.writeFile(manifestPath,JSON.stringify(manifest));
  });
});
