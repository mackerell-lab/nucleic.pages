import { spawn } from 'node:child_process';
import { createHash } from 'node:crypto';
import { fileURLToPath } from 'node:url';
import path from 'node:path';
import { buildTopology } from './stem_projection.mjs';
import { computeGeometry, geometryFamilies } from './geometry_adapter.mjs';

const workspace = fileURLToPath(new URL('../../../', import.meta.url));
export const FR3D_COMMIT = '994e54ea8fcea1a0484ea8c082f4a41c5406191d';

export async function annotateInteractions(entry, options = {}) {
  const topology = buildTopology(entry);
  const view = {pdb_id:entry.pdb_id,source_sha256:entry.source_sha256,coordinate_policy:entry.coordinate_policy,
    residues:entry.residues,links:entry.links};
  const coordinate_view_hash=createHash('sha256').update(JSON.stringify(view)).digest('hex');
  const input=JSON.stringify({entry:view,coordinate_view_hash,topology:Object.fromEntries(topology.positions)});
  const python=options.python ?? path.join(workspace,'data/pure_rna/venv/bin/python');
  const fr3dPath=options.fr3dPath ?? path.join(workspace,'data/pure_rna/fr3d-python');
  const output=await new Promise((resolve,reject)=>{
    const child=spawn(python,[fileURLToPath(new URL('./fr3d_provider.py',import.meta.url)),'--fr3d-path',fr3dPath],
      {stdio:['pipe','pipe','pipe'],env:{...process.env,OPENBLAS_NUM_THREADS:'1',OMP_NUM_THREADS:'1',MKL_NUM_THREADS:'1'}});
    const stdout=[],stderr=[];
    child.stdout.on('data',b=>stdout.push(b)); child.stderr.on('data',b=>stderr.push(b));
    child.on('error',reject);
    child.on('close',code=>code===0?resolve(Buffer.concat(stdout).toString()):reject(new Error(`FR3D exited ${code}: ${Buffer.concat(stderr).toString().slice(-12000)}`)));
    child.stdin.on('error',()=>{}); child.stdin.end(input);
  });
  const graph=JSON.parse(output);
  if(graph.coordinate_view_hash!==coordinate_view_hash || graph.provider_input_sha256!==createHash('sha256').update(input).digest('hex')) throw new Error('FR3D selected-coordinate view mismatch');
  if(graph.provider?.commit!==FR3D_COMMIT) throw new Error('FR3D provider version mismatch');
  const residueIds=new Set(entry.residues.map(r=>r.id)),edgeIds=new Set();
  if(residueIds.size!==entry.residues.length) throw new Error('Duplicate selected RNA residue identity');
  for(const edge of graph.edges) {
    if(!residueIds.has(edge.residue1_id)||!residueIds.has(edge.residue2_id)) throw new Error('Unmapped FR3D interaction endpoint');
    if(edgeIds.has(edge.id)) throw new Error('Duplicate FR3D interaction identity');
    edgeIds.add(edge.id);
  }
  return graph;
}

export { computeGeometry };

export async function computeInteractionGeometry(entry,options={}) {
  const graph=await annotateInteractions(entry,options);
  const geometry=computeGeometry(entry,graph);
  return {...geometry,graph,families:geometryFamilies(geometry),interactions:graph.edges,
    capabilities:{base_pair:'available',pair_quality:'available',step:'available',helical:'available',
      step_position:'available_except_DNA_ABI',same_strand:'available',helix_radius:'available',
      interaction_graph:'FR3D_pinned',stem_projection:'unambiguous_conventional_cWW',
      symmetry_mates:'not_computed',modified_bases:'not_supported'}};
}
