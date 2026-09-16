import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import os from 'node:os';
import path from 'node:path';
import {OutputScope} from '../offline/output_scope.mjs';
import {parseArgs} from '../offline/build_dataset.mjs';
import {RNA_QUERY,reconcileDiscovery} from '../offline/discovery.mjs';
import {normalizeAnnotation} from '../offline/annotations.mjs';

test('output scope rejects DNA siblings, traversal and symlink escapes',async()=>{
  const root=await fs.mkdtemp(path.join(os.tmpdir(),'rna-scope-'));
  try {
    const owned=path.join(root,'rna'), outside=path.join(root,'dna');
    await fs.mkdir(owned);await fs.mkdir(outside);
    const scope=new OutputScope([owned]);
    await scope.json(path.join(owned,'nested/test.json'),{ok:true});
    await assert.rejects(scope.write(path.join(outside,'manifest.json'),'bad'),/escapes/);
    await assert.rejects(scope.write(path.join(owned,'../dna/manifest.json'),'bad'),/escapes/);
    await fs.symlink(outside,path.join(owned,'link'));
    await assert.rejects(scope.write(path.join(owned,'link/result.json'),'bad'),/Symlink/);
    assert.equal(await fs.readFile(path.join(owned,'nested/test.json'),'utf8'),'{\n  "ok": true\n}\n');
  } finally {await fs.rm(root,{recursive:true,force:true});}
});

test('CLI unknown flags and unsafe build scopes fail closed',()=>{
  assert.throws(()=>parseArgs(['--build-id','good','--out-dir','/tmp']),/Unknown flag/);
  assert.throws(()=>parseArgs(['--build-id','../dna']),/safe/);
  assert.throws(()=>parseArgs(['--build-id','good','--concurrency','0']),/Invalid/);
  assert.throws(()=>parseArgs(['--build-id','good','--stage','fetch','--through','assets']),/Choose/);
  assert.throws(()=>parseArgs(['--build-id','good','--offline','--refresh-coordinates']),/refresh/);
  assert.deepEqual(parseArgs(['--build-id','good','--only-ids','1rna,1sdr']).onlyIds,['1RNA','1SDR']);
});

test('discovery counts and uniqueness must reconcile; RNA includes CG-only sequences',()=>{
  assert.equal(RNA_QUERY.query.nodes[0].parameters.attribute,'rcsb_entry_info.polymer_entity_count_RNA');
  assert.deepEqual(reconcileDiscovery({total_count:2,result_set:[{identifier:'1T4X'},{identifier:'1RNA'}]}),['1RNA','1T4X']);
  assert.throws(()=>reconcileDiscovery({total_count:2,result_set:[{identifier:'1RNA'}]}),/incomplete/);
  assert.throws(()=>reconcileDiscovery({total_count:2,result_set:[{identifier:'1RNA'},{identifier:'1RNA'}]}),/duplicate/);
});

test('NAKB annotations preserve separate entities and retain composition conflicts',()=>{
  const entry={pdb_id:'TEST',entities:[{entity_id:'1'},{entity_id:'2'},{entity_id:'3'}]};
  const raw={pdbid:'TEST',polyclass:'Protein/RNA',NAKBna:['function > makesprotein > transferrna','function > switch > vitaminswitch','nastructure > double > aform'],
    'NAKBna.entityids':['1','2'],'NAKBna.entityannot':['transferrna, aform','vitaminswitch']};
  const normalized=normalizeAnnotation(raw,entry);
  assert.equal(normalized.composition_conflict,true);
  assert.deepEqual(normalized.entities[0].functions,['makesprotein']);
  assert.deepEqual(normalized.entities[1].functions,['switch']);
  assert.deepEqual(normalized.entities[1].structures,[]);
  assert.equal(normalized.entities[2].annotation_status,'unknown');
  assert.equal(normalizeAnnotation(null,entry).status,'unknown');
  assert.throws(()=>normalizeAnnotation({...raw,'NAKBna.entityids':['1']},entry),/array mismatch/);
});
