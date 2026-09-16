import test from 'node:test';
import assert from 'node:assert/strict';
import { projectStems, buildTopology } from '../offline/stem_projection.mjs';
import { computeRnaLambda, computeGeometry, assertStepTopology } from '../offline/geometry_adapter.mjs';

const residue=(id,comp_id,atoms={})=>({id,comp_id,atoms,pdb_id:'TEST'});
const edge=(id,a,b,family='cWW')=>({id,residue1_id:a,residue2_id:b,family,category:'basepair',near:family.startsWith('n'),alternative:family.endsWith('a')});
const graph=edges=>({status:'available',edges});
const connected=(a,b)=>({from_id:a,to_id:b,status:'connected'});

test('RNA U lambda uses N1 and C6, independently of O2 prime and T methyl',()=>{
  const a=residue('a','A',{"C1'":[0,0,0],N9:[0,1,0],C8:[0,2,0]});
  const u=residue('u','U',{"C1'":[4,0,0],N1:[4,1,0],C6:[4,2,0],O2:[90,90,90],"O2'":[-90,-90,-90]});
  assert.deepEqual(computeRnaLambda({residue:a,base_code:'A'},{residue:u,base_code:'U'}),
    {lambda_1:90,lambda_2:90,c1c1:4,rn9_yn1:4,rc8_yc6:4});
  const result=computeGeometry({pdb_id:'TEST',residues:[a,u],links:[]},graph([edge('p','a','u')]));
  assert.equal(result.pairs[0].values.lambda_2,90);
  assert.equal(result.pairs[0].statuses.shear,'missing_base_frame');
});

test('noncanonical and near interactions survive graph while stem projection requires conventional cWW',()=>{
  const entry={residues:[residue('a','A'),residue('u','U'),residue('g','G'),residue('c','C')],links:[]};
  const input=graph([edge('p','a','u'),edge('n','a','g','tHS'),edge('q','g','c','ncWW')]);
  const before=structuredClone(input),projection=projectStems(entry,input);
  assert.deepEqual(input,before); assert.equal(projection.pairs.length,1); assert.equal(projection.steps.length,0);
});

test('competing cWW partners are ambiguous instead of greedy matching',()=>{
  const entry={residues:[residue('g','G'),residue('u','U'),residue('c','C')],links:[]};
  const result=projectStems(entry,graph([edge('gu','g','u'),edge('gc','g','c')]));
  assert.equal(result.pairs.length,0); assert.deepEqual(result.ambiguous_edge_ids,['gu','gc']);
});

test('step requires observed covalent adjacency on both strands',()=>{
  const entry={residues:[residue('a','A'),residue('b','G'),residue('c','C'),residue('d','U')],links:[connected('a','b'),connected('c','d')]};
  const input=graph([edge('p','a','d'),edge('q','b','c')]);
  const valid=projectStems(entry,input);
  assert.equal(valid.steps.length,1); assert.equal(valid.steps[0].strand2_delta,-1);
  entry.links[1].status='broken_modeled_bond';
  assert.equal(projectStems(entry,input).steps.length,0);
});

test('cyclic and branched topologies cannot fabricate linear ordinals',()=>{
  const entry={residues:[residue('a','A'),residue('b','U')],links:[connected('a','b'),connected('b','a')]};
  const topology=buildTopology(entry);
  assert.equal(topology.positions.get('a').ordinal,null);
  assert.equal(topology.diagnostics.length,2);
  entry.links.push(connected('a','a'));
  assert.throws(()=>buildTopology(entry),/Branched/);
});

test('finite coordinates with missing or corrupt ordinals fail adapter contract',()=>{
  const frame=(chain_id,chain_pos)=>({chain_id,chain_pos,origin:[1,2,3]});
  const first={nt1:frame('a',0),nt2:frame('b',1)},second={nt1:frame('a',1),nt2:frame('b',0)};
  assert.equal(assertStepTopology(first,second),-1);
  delete second.nt2.chain_pos;
  assert.throws(()=>assertStepTopology(first,second),/Missing topology/);
  second.nt2.chain_pos=8;
  assert.throws(()=>assertStepTopology(first,second),/Unsupported RNA/);
});

test('hairpin pair orientation follows numeric polymer order instead of opaque ID sorting',()=>{
  const first={...residue('A.G2','G'),label_asym_id:'A',label_seq_id:'2'};
  const second={...residue('A.C11','C'),label_asym_id:'A',label_seq_id:'11'};
  const entry={pdb_id:'TEST',residues:[second,first],links:[]};
  const input=graph([edge('p',second.id,first.id)]);
  const projection=projectStems(entry,input);
  assert.equal(projection.pairs[0].residue1_id,first.id);
  const calculated=computeGeometry(entry,input);
  assert.equal(calculated.pairs[0].residue1_id,first.id);
  assert.equal(calculated.pairs[0].pair_label,'G-C');
});
