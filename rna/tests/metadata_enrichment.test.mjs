import test from 'node:test';
import assert from 'node:assert/strict';
import {enrichEntryMetadata, COVERAGE_POLICY} from '../offline/metadata_enrichment.mjs';
const sha='a'.repeat(64);
function fixture(lengths, copies, supported, models=1){
 const entities=lengths.map((length,index)=>({entity_id:String(index+1),type:'polymer',polymer_type:'polyribonucleotide',sequence:Array.from({length},(_,position)=>({seq_id:String(position+1),mon_id:'G'}))}));
 const struct_asym=copies.map((entity,index)=>({id:String.fromCharCode(65+index),entity_id:String(entity)}));
 const residues=struct_asym.flatMap(chain=>entities[Number(chain.entity_id)-1].sequence.map(position=>({id:`${chain.id}:${position.seq_id}`,entity_id:chain.entity_id,label_asym_id:chain.id,label_seq_id:position.seq_id,comp_id:'G',model_id:'1',atoms:{}})));
 for(let i=0;i<supported;i++)residues[i].atoms={P:[1,2,3]};
 return {entry:{source_sha256:sha,entities,residues,components:[],chemical_components:[],declared_components:{nonpolymer:[]},explicit_connections:[],coordinate_policy:{id:'single_deposited_model_v1',scope:'deposited_asymmetric_unit',model_id:'1',model_count:models}},registry:{source_sha256:sha,struct_asym}};
}
test('coverage distinguishes declared chain copies, NMR models and empty coordinate groups',()=>{
 for(const [label,lengths,copies,supported,models,expected]of [
  ['157D',[12],[1,1],24,1,24],['1AM0',[40],[1],16,8,40],['10ZT',[519],[1],96,1,519],['1K9W',[23],[1,1,1,1],80,1,92],['8ZZJ',[37,37,32,32],[1,2,3,4],0,10,138],
 ]){const f=fixture(lengths,copies,supported,models),before=structuredClone(f.entry),result=enrichEntryMetadata(f.entry,f.registry);assert.equal(result.coverage.status,'available',label);assert.equal(result.coverage.observed_position_count,supported,label);assert.equal(result.coverage.declared_position_count,expected,label);assert.equal(result.coverage.fraction,supported/expected,label);assert.equal(result.coverage.policy,COVERAGE_POLICY);assert.deepEqual(f.entry,before);}
});
test('entirely absent chains and duplicate alternatives keep declared denominator',()=>{
 const f=fixture([2],[1,1],2);f.entry.residues=f.entry.residues.filter(r=>r.label_asym_id==='A');f.entry.residues.push({...f.entry.residues[0],id:'alternate',altloc:'B'});f.entry.entities[0].sequence.push({seq_id:'1',mon_id:'A'});
 const result=enrichEntryMetadata(f.entry,f.registry);assert.equal(result.coverage.observed_position_count,2);assert.equal(result.coverage.declared_position_count,4);assert.equal(result.coverage.fraction,.5);
});
test('missing declarations or conflicting identities produce unknown coverage without clamping',()=>{
 const mutations=[f=>f.entry.entities[0].sequence=[],f=>f.registry.struct_asym=[],f=>f.entry.residues[0].label_asym_id='unknown',f=>f.entry.residues[0].model_id='2',f=>f.entry.residues[0].label_seq_id='999',f=>f.entry.coordinate_policy.scope='biological_assembly',f=>f.entry.entities[0].polymer_type='polydeoxyribonucleotide'];
 for(const mutate of mutations){const f=fixture([2],[1],2);mutate(f);const coverage=enrichEntryMetadata(f.entry,f.registry).coverage;assert.equal(coverage.status,'unknown');assert.equal(coverage.observed_position_count,null);assert.equal(coverage.declared_position_count,null);assert.equal(coverage.fraction,null);assert.ok(coverage.reasons.length);}
 const f=fixture([2],[1],2);assert.throws(()=>enrichEntryMetadata(f.entry,{...f.registry,source_sha256:'b'.repeat(64)}),/source identity/);
});
test('component counts separate instances atoms water declarations and covalent annotations',()=>{
 const f=fixture([2],[1],2);for(const [id,type]of [['2','non-polymer'],['3','water']])f.entry.entities.push({entity_id:id,type});f.registry.struct_asym.push({id:'B',entity_id:'2'},{id:'C',entity_id:'3'});
 f.entry.chemical_components=[{id:'BR',name:'BROMIDE ION',type:'non-polymer'},{id:'HOH',name:'WATER',type:'non-polymer'},{id:'MG',name:'MAGNESIUM ION'}];
 f.entry.declared_components.nonpolymer=[{entity_id:'2',comp_id:'BR'},{entity_id:'3',comp_id:'HOH'},{entity_id:'4',comp_id:'MG'}];
 f.entry.components=[{id:'br',comp_id:'BR',entity_id:'2',label_asym_id:'B',model_id:'1',is_water:false,atoms:{BR:[1,0,0]}},{id:'w1',comp_id:'HOH',entity_id:'3',label_asym_id:'C',model_id:'1',is_water:true,atoms:{O:[1,0,0],H1:[1,1,0],H2:[1,-1,0]}},{id:'w2',comp_id:'HOH',entity_id:'3',label_asym_id:'C',model_id:'1',is_water:true,atoms:{O:[2,0,0]}},{id:'empty',comp_id:'HOH',entity_id:'3',label_asym_id:'C',model_id:'1',is_water:true,atoms:{}}];
 f.entry.explicit_connections=[{conn_type_id:'covale',ptnr1_label_asym_id:'A',ptnr1_label_comp_id:'G',ptnr2_label_asym_id:'B',ptnr2_label_comp_id:'BR'},{conn_type_id:'metalc',ptnr1_label_asym_id:'A',ptnr1_label_comp_id:'G',ptnr2_label_asym_id:'B',ptnr2_label_comp_id:'BR'}];
 const result=enrichEntryMetadata(f.entry,f.registry);assert.equal(result.coverage.declared_position_count,2);assert.deepEqual(result.associated_components.declared_only,[{comp_id:'MG',name:'MAGNESIUM ION'}]);
 const br=result.associated_components.observed.find(c=>c.comp_id==='BR'),water=result.associated_components.observed.find(c=>c.comp_id==='HOH');
 assert.equal(br.name,'BROMIDE ION');assert.equal(br.observed_instance_count,1);assert.equal(br.declared_covalent_connection_count,1);assert.equal(br.declared_rna_covalent_connection_count,1);assert.equal(br.ion,undefined);
 assert.equal(water.observed_instance_count,2);assert.equal(water.observed_atom_count,4);assert.equal(water.water_instance_count,2);
});
test('component identity conflicts fail closed and unknown names stay null',()=>{
 const f=fixture([1],[1],1);f.entry.entities.push({entity_id:'2',type:'non-polymer'});f.registry.struct_asym.push({id:'B',entity_id:'2'});const component={id:'unknown',comp_id:'X',entity_id:'2',label_asym_id:'B',model_id:'1',atoms:{X:[0,0,0]}};f.entry.components=[component];assert.equal(enrichEntryMetadata(f.entry,f.registry).associated_components.observed[0].name,null);f.entry.components.push({...component});assert.throws(()=>enrichEntryMetadata(f.entry,f.registry),/component identity/);
});
