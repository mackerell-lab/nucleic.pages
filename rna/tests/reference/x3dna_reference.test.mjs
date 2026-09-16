import test from 'node:test';
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import { computeResidueObservables } from '../../offline/residue_geometry.mjs';
const fixture=JSON.parse(await readFile(new URL('./rna_x3dna_reference.json',import.meta.url),'utf8'));
for(const entry of fixture.entries){
  test(`${entry.pdb_id}: independently generated 3DNA RNA residue reference`,()=>{
    const actual=computeResidueObservables(entry);
    assert.equal(actual.length,entry.residues.length);
    let compared=0;
    for(const row of actual){
      const expected=fixture.reference[entry.pdb_id][row.id];
      assert.ok(expected,`Missing independent reference ${row.id}`);
      for(const [key,raw] of Object.entries(expected)){
        const reference=key==='e_z' && raw!==null ? ((raw+180)%360+360)%360-180 : raw;
        const value=row.values[key];
        if(reference===null){assert.equal(value,null,`${row.id} ${key} should be missing`);continue;}
        assert.ok(Number.isFinite(value),`${row.id} ${key} must be finite`);
        const linear=['tm','sszp','dp','e_z'].includes(key);
        const error=linear ? Math.abs(value-reference) : Math.abs(((value-reference+180)%360+360)%360-180);
        const tolerance=['sszp','dp'].includes(key) ? 0.006 : 0.051;
        assert.ok(error<=tolerance,`${row.id} ${key}: actual ${value}, 3DNA ${reference}, error ${error}`);
        compared++;
      }
    }
    assert.ok(compared>100,'Fixture must exercise real finite RNA metrics');
  });
}

for(const entry of fixture.entries){
  test(`${entry.pdb_id}: independently evaluated glycosidic and ribose O2' atoms`,()=>{
    for(const row of computeResidueObservables(entry)){
      for(const [parameter,reference] of Object.entries(fixture.analytical_local_reference[entry.pdb_id][row.id])){
        const actual=row.values[parameter];
        if(reference===null)assert.equal(actual,null,`${row.id} ${parameter}`);
        else assert.ok(Number.isFinite(actual)&&Math.abs(actual-reference)<1e-8,`${row.id} ${parameter}: ${actual} versus ${reference}`);
      }
    }
  });
}

const { computeGeometry }=await import('../../offline/geometry_adapter.mjs');
for(const entry of fixture.entries){
  test(`${entry.pdb_id}: independently generated 3DNA pair and step reference`,()=>{
    const geometry=computeGeometry(entry,fixture.graphs[entry.pdb_id]);
    assert.ok(geometry.pairs.length>0);assert.ok(geometry.steps.length>0);
    let compared=0;
    for(const [level,groups] of [['pairs',['base_pairs','lambda']],['steps',['steps','helical','step_position','same_strand','helix_radius']]]){
      for(const row of geometry[level]){
        const key=row.residue_ids.join('|');
        for(const group of groups){
          const reference=fixture.geometry_reference[entry.pdb_id][group][key];
          assert.ok(reference,`Missing ${group} independent reference for ${key}`);
          for(const [parameter,expected] of Object.entries(reference)){
            const actual=row.values[parameter];
            if(expected===null){assert.equal(actual,null,`${key} ${parameter}`);continue;}
            assert.ok(Number.isFinite(actual),`${key} ${parameter} expected finite`);
            const circular=['buckle','propeller','opening','tilt','roll','twist','inclination','tip','h_twist'].includes(parameter);
            const error=circular ? Math.abs(((actual-expected+180)%360+360)%360-180) : Math.abs(actual-expected);
            const tolerance=group==='lambda' ? 0.051 : 0.006;
            assert.ok(error<=tolerance,`${key} ${parameter}: actual ${actual}, 3DNA ${expected}, error ${error}`);
            compared++;
          }
        }
      }
    }
    assert.ok(compared>100);
  });
}

test('Pair and step frame conventions do not depend on opaque residue IDs',()=>{
  for(const entry of fixture.entries){
    const graph=fixture.graphs[entry.pdb_id];
    const baseline=computeGeometry(entry,graph);
    const aliases=new Map(entry.residues.map((residue,index)=>[residue.id,`opaque:${entry.residues.length-index}`]));
    const changed={...entry,residues:entry.residues.map(residue=>({...residue,id:aliases.get(residue.id)})),
      links:entry.links.map(link=>({...link,from_id:aliases.get(link.from_id),to_id:aliases.get(link.to_id)}))};
    const changedGraph={...graph,edges:graph.edges.map(edge=>({...edge,residue1_id:aliases.get(edge.residue1_id),residue2_id:aliases.get(edge.residue2_id)}))};
    const actual=computeGeometry(changed,changedGraph);
    for(const level of ['pairs','steps']){
      assert.equal(actual[level].length,baseline[level].length,`${entry.pdb_id} ${level} count`);
      const original=new Map(baseline[level].map(row=>[row.id,row]));
      for(const row of actual[level])assert.deepEqual(row.values,original.get(row.id)?.values,`${entry.pdb_id} ${row.id} changed after ID relabeling`);
    }
  }
});
