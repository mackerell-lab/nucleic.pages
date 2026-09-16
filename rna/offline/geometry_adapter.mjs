import { buildBaseFrame } from './base_frames.mjs';
import { projectStems, residueOrder } from './stem_projection.mjs';
import {
  computePairGeometry, computeStepGeometry, computeHelicalGeometry,
  computeStepPositionMetrics, computeSameStrandMetrics, computeHelixRadiusMetrics, reverseXZFrame,
} from './vendor/dna_geometry_core.mjs';

const PAIR = ['shear','stretch','stagger','buckle','propeller','opening'];
const QUALITY = ['lambda_1','lambda_2','c1c1','rn9_yn1','rc8_yc6'];
const STEP = ['shift','slide','rise','tilt','roll','twist'];
const HELICAL = ['x_disp','y_disp','h_rise','inclination','tip','h_twist'];
const POSITION = ['xp','yp','zp','xph','yph','zph'];
const SAME = ['strand_i_p_p','strand_i_c1_c1','strand_ii_p_p','strand_ii_c1_c1'];
const RADIUS = ['strand_i_p_radius','strand_i_o4_radius','strand_i_c1_radius','strand_ii_p_radius','strand_ii_o4_radius','strand_ii_c1_radius'];
export const GEOMETRY_FAMILIES = {base_pair:PAIR,pair_quality:QUALITY,step:STEP,helical:HELICAL,step_position:POSITION,same_strand:SAME,helix_radius:RADIUS};
const distance = (a,b) => a && b ? Math.hypot(...a.map((v,i) => v-b[i])) : null;
function angle(a,b,c) {
  if (!a || !b || !c) return null;
  const x=a.map((v,i)=>v-b[i]), y=c.map((v,i)=>v-b[i]);
  const norm=Math.hypot(...x)*Math.hypot(...y);
  return norm > 1e-12 ? Math.acos(Math.max(-1,Math.min(1,x.reduce((s,v,i)=>s+v*y[i],0)/norm)))*180/Math.PI : null;
}

// x3DNA ana_fncs.c uses pyrimidine N1/C6 and purine N9/C8. U is
// explicitly a pyrimidine here; no T alias and no DNA-only helper dispatch.
export function computeRnaLambda(frame1,frame2) {
  const a=frame1.residue.atoms, b=frame2.residue.atoms;
  const na=a['AG'.includes(frame1.base_code)?'N9':'N1'];
  const nb=b['AG'.includes(frame2.base_code)?'N9':'N1'];
  return {lambda_1:angle(na,a["C1'"],b["C1'"]), lambda_2:angle(nb,b["C1'"],a["C1'"]),
    c1c1:distance(a["C1'"],b["C1'"]), rn9_yn1:distance(na,nb),
    rc8_yc6:distance(a['AG'.includes(frame1.base_code)?'C8':'C6'],b['AG'.includes(frame2.base_code)?'C8':'C6'])};
}

export function assertStepTopology(pairA,pairB) {
  const a=pairA.nt1,b=pairA.nt2,c=pairB.nt1,d=pairB.nt2;
  if ([a,b,c,d].some(f => !Number.isInteger(f?.chain_pos) || !f?.chain_id)) throw new Error('Missing topology ordinal in RNA geometry adapter');
  const delta=d.chain_pos-b.chain_pos;
  if(a.chain_id!==c.chain_id || c.chain_pos-a.chain_pos!==1 || b.chain_id!==d.chain_id || ![-1,1].includes(delta)) throw new Error('Unsupported RNA step topology');
  return delta;
}

function measurements(keys,raw={},failure='missing_atoms') {
  const values={},statuses={};
  for(const key of keys) { const good=Number.isFinite(raw[key]); values[key]=good?raw[key]:null; statuses[key]=good?'available':failure; }
  return {values,statuses};
}

// x3DNA ana_fncs.c:1303-1330 orients a whole connected helix when every
// paired base normal points against the ordered pair progression. Apply the
// same geometric criterion per supported stem, never by an RNA type label.
function stemFrameConventions(projection,frames) {
  const neighbors=new Map(projection.pairs.map(p=>[p.id,new Set()]));
  for(const step of projection.steps) {
    neighbors.get(step.pair1_id).add(step.pair2_id);
    neighbors.get(step.pair2_id).add(step.pair1_id);
  }
  const conventions=new Map(),visited=new Set();
  for(const pair of projection.pairs) {
    if(visited.has(pair.id)) continue;
    const members=new Set(),queue=[pair.id];
    while(queue.length) {const id=queue.pop(); if(members.has(id)) continue; members.add(id); visited.add(id); queue.push(...neighbors.get(id));}
    const steps=projection.steps.filter(s=>members.has(s.pair1_id));
    const directions=[];
    for(const step of steps) {
      const fs=step.residue_ids.map(id=>frames.get(id));
      if(!fs.every(Boolean)) {directions.push(null);continue;}
      // nt2 normal is aligned with nt1 by computePairGeometry before the
      // whole-helix convention is evaluated in x3DNA.
      const dot=(a,b)=>a.reduce((s,v,i)=>s+v*b[i],0);
      for(const [first,second] of [[0,2],[1,3]]) {
        const displacement=fs[second].origin.map((x,i)=>x-fs[first].origin[i]);
        for(const index of [first,second]) {
          const reference=index<2?fs[0]:fs[2];
          const sign=dot(fs[index].z_axis,reference.z_axis)<0?-1:1;
          directions.push(sign*dot(displacement,fs[index].z_axis));
        }
      }
    }
    const negative=directions.length>0 && directions.every(x=>x!==null&&x<-1e-8);
    const positive=directions.length>0 && directions.every(x=>x!==null&&x>1e-8);
    const convention=negative?'x3dna_stem_reverse_xz':positive?'x3dna_stem_standard':steps.length?'mixed_or_unresolved_stem':'isolated_standard';
    for(const id of members) conventions.set(id,convention);
  }
  return conventions;
}

export function computeGeometry(entry,graph) {
  if(graph.status!=='available') return {pairs:[],steps:[],relations:[],diagnostics:[{status:'interaction_provider_unavailable'}]};
  const projection=projectStems(entry,graph), topology=projection.topology;
  const compare=residueOrder(entry);
  const residues=new Map(entry.residues.map(r=>[r.id,r]));
  const frames=new Map(entry.residues.map(r=> {const p=topology.positions.get(r.id); return [r.id,buildBaseFrame(r,{chain_id:p?.segment,chain_pos:p?.ordinal})];}));
  const conventions=stemFrameConventions(projection,frames);
  const eligible=new Set(projection.pairs.map(p=>p.edge_id));
  const pairs=[],steps=[],relations=[],seen=new Set();
  for(const edge of graph.edges) {
    if(edge.category!=='basepair') continue;
    const ids=[edge.residue1_id,edge.residue2_id].sort(compare);
    // Reciprocal FR3D directed families refer to one geometric pair.
    const family=ids[0]===edge.residue1_id?edge.family:reversePairFamily(edge.family);
    const key=ids.join('\u0000')+'\u0000'+family;
    if(seen.has(key)) continue;
    seen.add(key);
    const convention=conventions.get(edge.id)??'isolated_standard';
    let a=frames.get(ids[0]),b=frames.get(ids[1]);
    if(convention==='x3dna_stem_reverse_xz') {if(a)a=reverseXZFrame(a);if(b)b=reverseXZFrame(b);}
    const ra=residues.get(ids[0]),rb=residues.get(ids[1]);
    if(!ra||!rb) throw new Error('Interaction graph references an absent residue');
    const raw={...(a&&b?computePairGeometry(a,b).params_obj:{}),...computeRnaLambda({base_code:ra.comp_id,residue:ra},{base_code:rb.comp_id,residue:rb})};
    const measured=measurements([...PAIR,...QUALITY],raw,'missing_atoms');
    if(!a||!b) for(const key of PAIR) measured.statuses[key]='missing_base_frame';
    const row={id:edge.id,entry_id:entry.pdb_id,pdb_id:entry.pdb_id,residue1_id:ids[0],residue2_id:ids[1],residue_ids:ids,
      family,pair_label:ids.map(id=>residues.get(id)?.comp_id).join('-'),
      stem_eligible:eligible.has(edge.id),frame_convention:convention,near:edge.near,alternative:edge.alternative,
      atom_roles:{residue1:{glycosidic_n:'AG'.includes(ra.comp_id)?'N9':'N1',adjacent_c:'AG'.includes(ra.comp_id)?'C8':'C6'},
        residue2:{glycosidic_n:'AG'.includes(rb.comp_id)?'N9':'N1',adjacent_c:'AG'.includes(rb.comp_id)?'C8':'C6'}},
      ...measured};
    pairs.push(row);
    ids.forEach((residue_id,side)=>relations.push({id:`${row.id}/residue/${side+1}`,kind:'pair_residue',pair_id:row.id,residue_id,side:side+1}));
  }
  for(const candidate of projection.steps) {
    const convention=conventions.get(candidate.pair1_id);
    const fs=candidate.residue_ids.map(id=> {const f=frames.get(id);return f&&convention==='x3dna_stem_reverse_xz'?reverseXZFrame(f):f;});
    let raw={},status='missing_base_frame';
    const supported=convention!=='mixed_or_unresolved_stem' && candidate.topology==='antiparallel';
    if(!supported) status='unsupported_stem_orientation';
    if(fs.every(Boolean) && supported) {
      const first=computePairGeometry(fs[0],fs[1]),second=computePairGeometry(fs[2],fs[3]);
      const delta=assertStepTopology(first,second);
      if(delta!==candidate.strand2_delta) throw new Error('RNA step topology provider mismatch');
      // These are geometric frame conventions, never an A/B/Z assignment.
      const options={applyBzAdjustment:false,handedness:convention==='x3dna_stem_reverse_xz'?'left':'right'};
      const step=computeStepGeometry(first.frame,second.frame,options);
      const helix=computeHelicalGeometry(first.frame,second.frame,options);
      const position=computeStepPositionMetrics(first,second,step.frame,helix.frame);
      const radii=computeHelixRadiusMetrics(first,second,helix.frame);
      // Legacy means allow a single endpoint. RNA exposes the declared two-end
      // observable only when both endpoints exist.
      for(const [side,indices] of [['i',[0,2]],['ii',[1,3]]]) for(const [atom,label] of [["O4'",'o4'],["C1'",'c1']]) {
        if(indices.some(i=>!fs[i].residue.atoms[atom])) radii[`strand_${side}_${label}_radius`]=null;
      }
      // Average phosphate position needs both strand phosphates.
      if(!fs[2].residue.atoms.P || !(delta===1?fs[3]:fs[1]).residue.atoms.P) POSITION.forEach(k=>{position[k]=null;});
      raw={...step.params_obj,...helix.params_obj,...Object.fromEntries(POSITION.map(k=>[k,position[k]])),...computeSameStrandMetrics(first,second),...radii};
      status='missing_atoms';
    }
    // Same-strand distances require atoms and confirmed links, not fitted bases.
    const rs=candidate.residue_ids.map(id=>residues.get(id));
    for(const [side,indices] of [['i',[0,2]],['ii',[1,3]]]) for(const [atom,label] of [['P','p'],["C1'",'c1']]) {
      raw[`strand_${side}_${label}_${label}`]=distance(rs[indices[0]].atoms[atom],rs[indices[1]].atoms[atom]);
    }
    const measured=measurements([...STEP,...HELICAL,...POSITION,...SAME,...RADIUS],raw,status);
    for(const key of SAME) if(measured.values[key]===null) measured.statuses[key]='missing_atoms';
    const row={id:candidate.id,entry_id:entry.pdb_id,pdb_id:entry.pdb_id,pair1_id:candidate.pair1_id,pair2_id:candidate.pair2_id,
      residue_ids:candidate.residue_ids,topology:candidate.topology,frame_convention:convention,
      step_label:candidate.residue_ids.map(id=>residues.get(id)?.comp_id).join(''),
      ...measured};
    steps.push(row);
    [row.pair1_id,row.pair2_id].forEach((pair_id,index)=>relations.push({id:`${row.id}/pair/${index+1}`,kind:'step_pair',step_id:row.id,pair_id,side:index+1}));
    row.residue_ids.forEach((residue_id,index)=>relations.push({id:`${row.id}/residue/${index+1}`,kind:'step_residue',step_id:row.id,residue_id,side:index+1}));
  }
  return {pairs,steps,relations,diagnostics:[...topology.diagnostics,...projection.ambiguous_edge_ids.map(edge_id=>({edge_id,status:'ambiguous_stem_partner'}))],
    provenance:{frame_source:'x3dna_2.4_atomic_templates',pair_kernel:'vendored_dna_geometry_core',step_options:{applyBzAdjustment:false,handedness:'per_stem_geometry'},stem_policy:projection.policy,abi:'not_applicable'}};
}

function reversePairFamily(code) {
  return code.replace(/([ct])([WHSB])([WHSB])/i,(_,o,a,b)=>o+b+a);
}

export function geometryFamilies(result) {
  return Object.fromEntries(Object.entries(GEOMETRY_FAMILIES).map(([family,keys])=>[family,
    (['base_pair','pair_quality'].includes(family)?result.pairs:result.steps).map(row=>({...row,
      values:Object.fromEntries(keys.map(k=>[k,row.values[k]])),statuses:Object.fromEntries(keys.map(k=>[k,row.statuses[k]]))}))]));
}
