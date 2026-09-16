import test from 'node:test';
import assert from 'node:assert/strict';
import { computeResidueObservables, residueNeighbors } from '../offline/residue_geometry.mjs';
import { buildBaseFrame } from '../offline/base_frames.mjs';
import { RNA_BASE_TEMPLATES } from '../offline/vendor/rna_base_templates.mjs';
import { computePucker, dihedralSigned, classifyPucker } from '../offline/numeric.mjs';
import { RESIDUE_PARAMETERS } from '../offline/parameter_registry.mjs';

const near = (a, b, eps=1e-8) => assert.ok(Math.abs(a-b) < eps, `${a} != ${b}`);
const make = (id, comp_id='U') => ({ id, comp_id, label_asym_id: 'A', atoms: {} });

test('signed torsions retain chirality and reject coincident/collinear atoms', () => {
  near(dihedralSigned([1,0,0], [0,0,0], [0,1,0], [0,1,1]), -90);
  near(dihedralSigned([1,0,0], [0,0,0], [0,1,0], [0,1,-1]), 90);
  assert.equal(dihedralSigned([0,0,0], [0,0,0], [0,1,0], [0,1,1]), null);
  assert.equal(dihedralSigned([0,0,0], [0,1,0], [0,2,0], [0,3,0]), null);
});

test('Altona phase/amplitude recovers ideal five-angle series across the circle', () => {
  for (const phase of [0, 18, 90, 162, 270, 359.9]) {
    const ring = Array.from({length:5}, (_, i) => 35*Math.cos((phase+144*(i-2))*Math.PI/180));
    const p = computePucker(...ring); near(p.p, phase); near(p.tm, 35);
  }
  assert.deepEqual(computePucker(0,0,0,0,0), {p:null,tm:0});
  assert.equal(classifyPucker(18), "C3'-endo");
  assert.equal(classifyPucker(162), "C2'-endo");
});

test('native uracil frame recovers imposed proper rotation and translation', () => {
  const r = make('u');
  r.quality_flags = { missing_atoms: 2, zero_occupancy: 0 };
  // A known +90-degree z rotation followed by translation.
  r.atoms = Object.fromEntries(Object.entries(RNA_BASE_TEMPLATES.U.coords).map(([n,[x,y,z]]) => [n,[-y+8,x-3,z+7]]));
  const frame = buildBaseFrame(r, {chain_id:'segment:u',chain_pos:0});
  assert.equal(frame.base_code,'U'); assert.equal(frame.reference,'x3dna_2.4_Atomic_U');
  frame.origin.forEach((x,i)=>near(x,[8,-3,7][i]));
  frame.x_axis.forEach((x,i)=>near(x,[0,1,0][i]));
  frame.y_axis.forEach((x,i)=>near(x,[-1,0,0][i]));
  near(frame.rmsd,0); assert.equal(frame.chain_pos,0);
  const row = computeResidueObservables({residues:[r],links:[]})[0];
  assert.ok(row.quality_flags.includes('missing_atoms'));
  assert.ok(!row.quality_flags.includes('zero_occupancy'));
  assert.equal(buildBaseFrame({...r,comp_id:'T'}),null);
  assert.equal(buildBaseFrame({...r,comp_id:'PSU'}),null);
  const linear = Object.fromEntries(Object.keys(r.atoms).map((n,i)=>[n,[i,0,0]]));
  assert.equal(buildBaseFrame({...r,atoms:linear}),null);
});

test('U uses N1/C2, individual angle coverage, and distinguishes O2 from O2-prime', () => {
  const r = make('u');
  r.atoms = {"O4'":[1,0,0], "C1'":[0,0,0], N1:[0,1,0], C2:[0,1,1], "C2'":[0,0,1], O2:[10,10,10]};
  const row = computeResidueObservables({residues:[r],links:[]})[0];
  near(row.values.chi,-90); near(row.values.o4_c1_n,90);
  near(row.values.c1_n1_c2,90);
  assert.equal(row.statuses.c1_n9_c4,'not_applicable');
  assert.equal(row.statuses.c1_n1_c6,'missing_atoms');
  assert.equal(row.statuses.c2_o2_length,'missing_atoms');
  r.atoms["O2'"] = [0,1,1];
  const withO2 = computeResidueObservables({residues:[r],links:[]})[0];
  near(withO2.values.c2_o2_length,1);
  near(withO2.values.c1_c2_o2,90);
  delete r.atoms["C2'"];
  const partial = computeResidueObservables({residues:[r],links:[]})[0];
  near(partial.values.o4_c1_n,90);
  assert.equal(partial.values.c2_c1_n,null);
});

test('array order never bridges a broken RNA backbone; Dp does not require a frame', () => {
  const a = make('a'), b = make('b');
  a.atoms = {"C1'":[0,0,0], N1:[1,0,0], "O3'":[0,0,1], "C3'":[1,0,1], "C4'":[1,1,1]};
  b.atoms = {P:[1,2,0], "O5'":[2,2,0], "C5'":[2,2,1]};
  const entry = {residues:[b,a],links:[{from_id:'a',to_id:'b',status:'connected'}]};
  const rows = computeResidueObservables(entry), ar = rows.find(r=>r.id==='a');
  near(ar.values.dp,2); assert.equal(ar.base_frame,null);
  assert.equal(ar.statuses.sszp,'missing_atoms_or_frame');
  assert.equal(ar.chain_pos,0); assert.equal(rows.find(r=>r.id==='b').chain_pos,1);
  entry.links[0].status='broken_modeled_bond';
  assert.equal(computeResidueObservables(entry).find(r=>r.id==='a').statuses.dp,'missing_next_link');
});

test('unsupported cyclic/branched graphs cannot fabricate linear neighbor ordinals', () => {
  const residues=[make('a'),make('b'),make('c')];
  const edge=(a,b)=>({from_id:a,to_id:b,status:'connected'});
  assert.throws(()=>residueNeighbors({residues,links:[edge('a','b'),edge('a','c')]}),/Branching/);
  const rows=computeResidueObservables({residues,links:[edge('a','b'),edge('b','a')]});
  assert.equal(rows[0].statuses.alpha,'unsupported_cyclic_topology');
  assert.equal(rows[0].chain_pos,undefined);
});

test('registry contains 29 inherited RNA-compatible fields and four explicit O2-prime metrics', () => {
  assert.equal(RESIDUE_PARAMETERS.length,33);
  assert.equal(new Set(RESIDUE_PARAMETERS.map(p=>p.id)).size,33);
  assert.equal(RESIDUE_PARAMETERS.filter(p=>p.family==='ribose_2oh').length,4);
  assert.equal(RESIDUE_PARAMETERS.find(p=>p.id==='e_z').circular,false);
  assert.ok(!RESIDUE_PARAMETERS.some(p=>p.id==='abi'||p.id==='bi_bii'));
});
