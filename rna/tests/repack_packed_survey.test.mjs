import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import path from 'node:path';
import os from 'node:os';
import {fileURLToPath} from 'node:url';
import {spawnSync} from 'node:child_process';
import {gzipSync, gunzipSync} from 'node:zlib';
import {sha256} from '../offline/output_scope.mjs';
import {verifyReleaseInventory} from '../offline/verify_release_inventory.mjs';
import {encodeFamilyRows, encodeSurveyRows, encodeCoordinateRows, decodeSurveyRows} from '../core/survey-codec.js';
import {BUNDLED_FAMILY_ENCODING} from '../core/bundled-family-codec.js';
import {BUNDLED_SURVEY_ENCODING, SURVEY_BUNDLE_ENCODING} from '../core/bundled-survey-codec.js';
import {encodePackedCoordinates} from '../core/packed-coordinate-codec.js';
import {encodePackedFamily} from '../core/packed-family-codec.js';
import {PACKED_SURVEY_ENCODING, encodePackedSurvey, expandPackedSurvey} from '../core/packed-survey-codec.js';
const script = fileURLToPath(new URL('../offline/repack_bundled_release.mjs', import.meta.url));
const read = async (root, descriptor) => JSON.parse(gunzipSync(await fs.readFile(path.join(root, descriptor.path))));
async function write(root, relative, data) {
  const raw = Buffer.from(JSON.stringify(data)), bytes = gzipSync(raw);
  await fs.mkdir(path.dirname(path.join(root, relative)), {recursive: true}); await fs.writeFile(path.join(root, relative), bytes);
  return {path: relative, sha256: sha256(bytes), bytes: bytes.length, uncompressed_bytes: raw.length};
}
async function fixture(t) {
  const root = await fs.mkdtemp(path.join(os.tmpdir(), 'rna-repack-packed-survey-'));
  t.after(() => fs.rm(root, {recursive:true, force:true}));
  const sourceRoot = path.join(root,'source'), candidateRoot = path.join(root,'candidate');
  await fs.mkdir(sourceRoot); await fs.mkdir(candidateRoot);
  const observations = ['r1','r2'];
  const column = sha256(JSON.stringify(observations));
  const bundle = {encoding:SURVEY_BUNDLE_ENCODING,columns:{[column]:observations}}, reference=sha256(JSON.stringify(bundle));
  const bundleDescriptor={...await write(sourceRoot,`survey/bundles/${reference}.json.gz`,bundle),content_sha256:reference,column_count:1};
  await write(candidateRoot,bundleDescriptor.path,bundle);
  const rows=observations.map((observation_id,index)=>({id:`${observation_id}|survey|angle`,observation_id,term_id:'angle',value:index?null:1.2345678901234567,status:index?'missing_atom':'ok',opening:index?null:-8,atom:"O2'"}));
  const scalar={...encodeSurveyRows(rows,'source-build'),encoding:BUNDLED_SURVEY_ENCODING};
  scalar.columns.observation_id={bundle:reference,column};
  const family={...encodeFamilyRows([{id:'r1',values:{chi:1.2345678901234567}}],'source-build'),encoding:BUNDLED_FAMILY_ENCODING};
  const packedFamily=encodePackedFamily(family), coordinates=encodePackedCoordinates(encodeCoordinateRows([{id:'atom',x:Number.MIN_VALUE,y:1.2345678901234567,z:-1e-30}],'source-build'));
  const source={schema_version:'rna-explorer-1',molecule_type:'RNA',build_id:'source-build',partial:false,counts:{entries:1},source:{provider:'unchanged'},
    metadata:await write(sourceRoot,'metadata.json.gz',{entries:[{pdb_id:'TEST'}]}),
    families:[{...await write(sourceRoot,'families/backbone.json.gz',packedFamily),encoding:packedFamily.encoding,id:'backbone',row_count:1,parameters:[{id:'chi',unit:'deg'}]}],family_bundles:{},
    relations:{edges:{...await write(sourceRoot,'relations/edges.json.gz',[{id:'edge',source:'r1',target:'r2'}]),row_count:1}},
    survey:{terms:[{term_id:'angle',unit:'deg'}],opening_bins:[{id:'middle',min:-8,max:2}],bundles:{[reference]:bundleDescriptor},
      scalars:{policy:'lazy',terms:{angle:{...await write(sourceRoot,'survey/scalars/angle.json.gz',scalar),encoding:scalar.encoding,row_count:2}}},
      coordinates:{unit:'angstrom',groups:{base:{row_count:1,partitions:[{...await write(sourceRoot,'survey/coordinates/base/0.json.gz',coordinates),encoding:coordinates.encoding,row_count:1,entry_ids:['TEST']}]}}}},
    provenance:{decisions:await write(sourceRoot,'decisions.json.gz',[{accepted:true}]),build_stages:{geometry:{status:'complete'}},repack:{operation:'previous'},source_release:{build_id:'initial'}}};
  const sourceFile=path.join(sourceRoot,'manifest.json'); await fs.writeFile(sourceFile,JSON.stringify(source));
  const packed=structuredClone(await encodePackedSurvey(scalar,async()=>bundle));
  const candidate={schema_version:'rna-survey-packed-candidate-1',survey_only:true,build_id:source.build_id,source_manifest:{sha256:sha256(await fs.readFile(sourceFile))},survey:structuredClone(source.survey)};
  Object.assign(candidate.survey.scalars.terms.angle,await write(candidateRoot,'survey/scalars/angle.json.gz',packed),{encoding:PACKED_SURVEY_ENCODING});
  const candidateFile=path.join(candidateRoot,'candidate.json');const save=()=>fs.writeFile(candidateFile,JSON.stringify(candidate));await save();
  const run=(name='output')=>spawnSync(process.execPath,[script,sourceFile,candidateFile,path.join(root,name),'new-build'],{encoding:'utf8',timeout:30000});
  return {root,sourceRoot,candidateRoot,sourceFile,candidateFile,source,candidate,packed,bundle,reference,rows,save,run};
}

test('packed Survey assembly preserves full mixed release and exact decoded rows',async t=>{
  const f=await fixture(t), sourceBytes=await fs.readFile(f.sourceFile),candidateBytes=await fs.readFile(f.candidateFile);
  const result=f.run();assert.equal(result.status,0,result.stderr);
  const output=path.join(f.root,'output'),manifest=JSON.parse(await fs.readFile(path.join(output,'manifest.json'))),report=JSON.parse(result.stdout);
  assert.equal(report.candidate_kind,'scalar');assert.ok(report.scalar_candidate);assert.equal(report.all_decoded_equal,true);
  assert.equal(manifest.provenance.repack.operation,'lossless_packed_survey_transport');assert.equal(manifest.provenance.repack.scalar_encoding,PACKED_SURVEY_ENCODING);
  assert.ok(manifest.provenance.repack.code_sha256['../core/packed-survey-codec.js']);
  assert.deepEqual(manifest.provenance.source_release.build_stages,f.source.provenance.build_stages);
  assert.deepEqual(manifest.provenance.source_release.repack,f.source.provenance.repack);assert.deepEqual(manifest.provenance.source_release.source_release,f.source.provenance.source_release);
  assert.deepEqual(manifest.survey.bundles,f.source.survey.bundles);assert.deepEqual(manifest.survey.terms,f.source.survey.terms);assert.deepEqual(manifest.survey.opening_bins,f.source.survey.opening_bins);
  assert.deepEqual(manifest.counts,f.source.counts);assert.equal(manifest.survey.scalars.policy,'lazy');
  const scalar=await read(output,manifest.survey.scalars.terms.angle);assert.equal(scalar.build_id,'new-build');
  const rows=decodeSurveyRows(await expandPackedSurvey(scalar,async()=>f.bundle));assert.deepEqual(rows,f.rows);assert.ok(Object.is(rows[0].value,f.rows[0].value));
  const pairs=[[manifest.families[0],f.source.families[0]],[manifest.survey.coordinates.groups.base.partitions[0],f.source.survey.coordinates.groups.base.partitions[0]]];
  for(const [actual,original]of pairs){assert.equal(actual.encoding,original.encoding);const a=await read(output,actual),b=await read(f.sourceRoot,original);assert.equal(a.build_id,'new-build');assert.deepEqual({...a,build_id:b.build_id},b);}
  for(const descriptor of [manifest.metadata,manifest.relations.edges,manifest.provenance.decisions,...Object.values(manifest.survey.bundles)])assert.deepEqual(await fs.readFile(path.join(output,descriptor.path)),await fs.readFile(path.join(f.sourceRoot,descriptor.path)));
  assert.ok((await verifyReleaseInventory(path.join(output,'manifest.json'))).ok);
  assert.ok((await fs.readFile(f.sourceFile)).equals(sourceBytes));assert.ok((await fs.readFile(f.candidateFile)).equals(candidateBytes));
  // A subsequent family migration must also read the now-packed Survey source.
  const nextRoot=path.join(f.root,'next-candidate');await fs.mkdir(nextRoot);
  await fs.cp(path.join(output,'families'),path.join(nextRoot,'families'),{recursive:true});
  const nextFile=path.join(nextRoot,'candidate.json'),nextSource=path.join(output,'manifest.json');
  await fs.writeFile(nextFile,JSON.stringify({schema_version:'rna-family-packed-candidate-1',family_only:true,build_id:manifest.build_id,
    source_manifest:{sha256:sha256(await fs.readFile(nextSource))},families:manifest.families,family_bundles:manifest.family_bundles}));
  const nextOutput=path.join(f.root,'next-output');const next=spawnSync(process.execPath,[script,nextSource,nextFile,nextOutput,'next-build'],{encoding:'utf8',timeout:30000});
  assert.equal(next.status,0,next.stderr);const nextManifest=JSON.parse(await fs.readFile(path.join(nextOutput,'manifest.json')));
  assert.equal(nextManifest.survey.scalars.terms.angle.encoding,PACKED_SURVEY_ENCODING);
  assert.deepEqual(decodeSurveyRows(await expandPackedSurvey(await read(nextOutput,nextManifest.survey.scalars.terms.angle),async()=>f.bundle)),f.rows);
  assert.ok((await verifyReleaseInventory(path.join(nextOutput,'manifest.json'))).ok);
});

test('packed Survey rejects identity, scope and metadata mutations before staging',async t=>{
 const cases=[
  [c=>c.build_id='wrong',/Candidate source identity/],[c=>c.source_manifest.sha256='0'.repeat(64),/Candidate source manifest identity/],
  [c=>c.scalar_only=true,/Exactly one/],[c=>{delete c.survey_only;c.scalar_only=true;},/Packed Survey candidate scope/],
  [c=>c.schema_version='rna-survey-bundled-candidate-1',/Packed Survey candidate scope/],
  [c=>delete c.survey.scalars.terms.angle,/Complete scalar registry/],
  [c=>c.survey.coordinates.unit='nm',/Unchanged packed Survey metadata/],
  [c=>c.survey.scalars.policy='eager',/Unchanged packed Survey scalar metadata/],
  [c=>c.survey.scalars.terms.angle.path='../outside.json.gz',/Unchanged packed Survey scalar path/],
  [c=>c.survey.scalars.terms.angle.row_count=1,/Unchanged packed Survey scientific descriptor/],
  [c=>c.survey.bundles={},/Unchanged packed Survey metadata/],
 ];
 for(const [mutate,expected]of cases){const f=await fixture(t);mutate(f.candidate);await f.save();const result=f.run();assert.notEqual(result.status,0);assert.match(result.stderr,expected);await assert.rejects(fs.access(path.join(f.root,'output')),{code:'ENOENT'});}
});

test('packed Survey rejects corrupt authenticated dependencies and numeric or status changes',async t=>{
 for(const mutation of ['reference','bundle','value','status','build','encoding','symlink']){
  const f=await fixture(t),descriptor=f.candidate.survey.scalars.terms.angle;
  if(mutation==='reference')f.packed.columns.observation_id.bundle='f'.repeat(64);
  if(mutation==='value'){f.packed.columns.value=[99,null];}
  if(mutation==='status')f.packed.columns.status[0]='missing_atom';
  if(mutation==='build')f.packed.build_id='wrong';
  if(mutation==='encoding')descriptor.encoding=BUNDLED_SURVEY_ENCODING;
  Object.assign(descriptor,await write(f.candidateRoot,descriptor.path,f.packed));await f.save();
  if(mutation==='bundle'){const file=path.join(f.candidateRoot,f.candidate.survey.bundles[f.reference].path),bytes=await fs.readFile(file);bytes[bytes.length-1]^=1;await fs.writeFile(file,bytes);}
  if(mutation==='symlink'){const file=path.join(f.candidateRoot,descriptor.path),target=path.join(f.root,'outside.gz');await fs.rename(file,target);await fs.symlink(target,file);}
  const result=f.run();assert.notEqual(result.status,0,mutation);assert.match(result.stderr,/Registered Survey bundle reference|Resource checksum|All original scalar transport metadata and columns|Candidate build identity|Resource encoding|Resource path or symlink escaped/);
  await assert.rejects(fs.access(path.join(f.root,'output-repack.json')),{code:'ENOENT'});
 }
});
