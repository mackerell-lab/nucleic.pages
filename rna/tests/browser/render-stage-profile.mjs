/** Instrumented production-data diagnosis; shared-host timings are not acceptance ceilings. */
import path from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';
import assert from 'node:assert/strict';
import { writeFile, readFile, mkdir } from 'node:fs/promises';
import { createHash } from 'node:crypto';
const workspace = path.resolve(process.env.RNA_WORKSPACE || fileURLToPath(new URL('../../../..', import.meta.url)));
const out = path.resolve(process.env.RNA_BROWSER_OUTPUT || path.join(workspace, 'data/pure_rna/profiling/render-stage-profile'));
const base = process.env.RNA_URL || 'http://127.0.0.1:8767/nucleic.pages/rna/';
const { chromium } = await import(pathToFileURL(process.env.PLAYWRIGHT_MODULE || path.join(workspace, 'data/pure_rna/browser-tools/node_modules/playwright-core/index.mjs')).href);
await mkdir(out, { recursive: true });
const pinnedRepository = process.env.RNA_PROFILE_REPOSITORY_SOURCE ? path.resolve(process.env.RNA_PROFILE_REPOSITORY_SOURCE) : null;
const pinnedBody = pinnedRepository ? await readFile(pinnedRepository, 'utf8') : null;
const overrideManifest = process.env.RNA_PROFILE_SOURCE_OVERRIDES ? path.resolve(process.env.RNA_PROFILE_SOURCE_OVERRIDES) : null;
const sourceOverrides = overrideManifest ? JSON.parse(await readFile(overrideManifest, 'utf8')) : {};
assert(sourceOverrides && !Array.isArray(sourceOverrides) && typeof sourceOverrides === 'object', 'Source overrides must be a relative-module-to-file object');
const overrideBodies = new Map();
for (const [file, source] of Object.entries(sourceOverrides)) {
  assert(/^[A-Za-z0-9_./-]+\.js$/.test(file) && !file.startsWith('/') && !file.split('/').includes('..') && typeof source === 'string', `Invalid source override: ${file}`);
  sourceOverrides[file] = path.resolve(path.dirname(overrideManifest), source);
  overrideBodies.set(file, await readFile(sourceOverrides[file], 'utf8'));
}
if (pinnedBody !== null) {
  assert(!overrideBodies.has('core/repository.js'), 'Repository source supplied twice');
  overrideBodies.set('core/repository.js', pinnedBody);
  sourceOverrides['core/repository.js'] = pinnedRepository;
}
const sha256 = value => createHash('sha256').update(value).digest('hex');
const report={repositoryOverride: pinnedRepository, sourceOverrides, loadedSourceSha256: {}, servedSourceSha256: {}, started:new Date().toISOString(),phases:[],errors:[],sourceSha256:{},limitations:'Single instrumented shared-host run; wall intervals inclusive and overlapping; CPU sampling approximate, no universal performance claim.'};
for(const file of ['app/PureRnaExplorer.js','core/selection.js','core/analysis.js','core/export.js','core/repository.js','math/numeric.js']) report.sourceSha256[file]=createHash('sha256').update(await readFile(new URL('../../'+file, import.meta.url))).digest('hex');
report.workspaceSourceSha256 = { ...report.sourceSha256 };
for (const [file, body] of overrideBodies) report.sourceSha256[file] = sha256(body);
const browser=await chromium.launch({headless:true});
const page=await browser.newPage({viewport:{width:1440,height:1000}});
page.on('pageerror', e=>report.errors.push(e.message));
await page.addInitScript(()=>{
 window.__intervals=[];window.__longtasks=[];
 new PerformanceObserver(list=>{for(const e of list.getEntries())window.__longtasks.push({start:e.startTime,duration:e.duration});}).observe({type:'longtask',buffered:true});
 window.__measurePrototype=(proto,prefix,methods)=>{for(const name of methods){const f=proto[name];if(typeof f!=='function')continue;proto[name]=function(...args){const start=performance.now(); const detail=name==='plot'?(args[0]?.id||args[0]?.className):['snapshot','snapshotAsync'].includes(name)?args[0]?.result?.kind:name==='loadFamily'?args[0]:'';let v;const done=()=>window.__intervals.push({name:prefix+'.'+name,detail,start,duration:performance.now()-start});try{v=f.apply(this,args);}catch(e){done();throw e;}if(v&&typeof v.then==='function')return v.then(r=>{done();return r;},e=>{done();throw e;});done();return v;};}};
});
const sourceBase = new URL(base);
await page.route(url => url.origin === sourceBase.origin && url.pathname.startsWith(sourceBase.pathname) && url.pathname.endsWith('.js'), async route => {
 const url = new URL(route.request().url());
 const file = url.pathname.slice(sourceBase.pathname.length);
 const response = await route.fetch();
 let body = overrideBodies.has(file) ? overrideBodies.get(file) : await response.text();
 report.loadedSourceSha256[file] = sha256(body);
 if (file === 'main.js') {
   assert(body.includes('const root ='), 'Bootstrap instrumentation anchor missing');
   body = body.replace('const root =', `window.__measurePrototype(PureRnaExplorer.prototype,'app',['start','render','snapshot','snapshotAsync','renderFamilyOverview','renderJoint','renderTable','updateContexts','updatePuckerControls','updateInteractionControls','requestJointOnly']);
 window.__measurePrototype(Object.getPrototypeOf(PureRnaExplorer.prototype),'base',['plot']);
 window.__measurePrototype(RnaDataRepository.prototype,'repository',['readJson','loadFamily','loadFamilyBundle','loadMetadata']);
 const root =`);
 }
 report.servedSourceSha256[file] = sha256(body);
 await route.fulfill({ response, body });
});

// Hash the complete result graph, including all rows and metadata, after timing.
// Sorted object keys and explicit reference tokens retain identity and aliasing;
// number tokens preserve signed zero, NaN and infinities. Snapshot metadata is
// outside result and is intentionally excluded. Chunk hashes bound allocations.
async function resultDigests() {
 return page.evaluate(async () => {
   const encoder = new TextEncoder();
   const hex = bytes => Array.from(new Uint8Array(bytes), byte => byte.toString(16).padStart(2, '0')).join('');
   async function digestGraph(root) {
     if (!root) return null;
     const seen = new Map(), chunks = []; let text = '', bytes = 0, objects = 0, numbers = 0;
     async function flush() {
       if (!text) return;
       const encoded = encoder.encode(text); bytes += encoded.byteLength;
       chunks.push(hex(await crypto.subtle.digest('SHA-256', encoded))); text = '';
     }
     const stack = [{ value: root }];
     while (stack.length) {
       const task = stack.pop();
       if (Object.hasOwn(task, 'token')) text += task.token;
       else {
         const value = task.value;
         if (value === null) text += 'null;';
         else if (typeof value === 'number') { numbers++; text += `n${Object.is(value, -0) ? '-0' : String(value)};`; }
         else if (typeof value !== 'object') text += `${typeof value}:${JSON.stringify(value)};`;
         else if (seen.has(value)) text += `ref:${seen.get(value)};`;
         else {
           const id = objects++; seen.set(value, id);
           const keys = Object.keys(value).sort();
           text += `${Array.isArray(value) ? 'array' : 'object'}:${id}:${Array.isArray(value) ? value.length : ''}{`;
           stack.push({ token: '};' });
           for (let index = keys.length - 1; index >= 0; index--) {
             const key = keys[index]; stack.push({ value: value[key] }); stack.push({ token: `${JSON.stringify(key)}:` });
           }
         }
       }
       if (text.length >= 1024 * 1024) await flush();
     }
     await flush();
     return { sha256: hex(await crypto.subtle.digest('SHA-256', encoder.encode(JSON.stringify(chunks)))),
       chunks: chunks.length, serializedBytes: bytes, uniqueObjects: objects, numericValues: numbers,
       encoding: 'canonical-result-graph-tokens-sha256-chunks-v1' };
   }
   const snapshots = window.rnaExplorer.snapshots;
   return { distribution: await digestGraph(snapshots.distribution?.result), joint: await digestGraph(snapshots.joint?.result) };
 });
}
const cdp=await page.context().newCDPSession(page);await cdp.send('Profiler.enable');await cdp.send('Profiler.setSamplingInterval',{interval:1000});
function summarize(profile){const nodes=new Map(profile.nodes.map(n=>[n.id,n]));const byFn=new Map(),byFile=new Map();for(let i=0;i<(profile.samples||[]).length;i++){const f=nodes.get(profile.samples[i])?.callFrame||{};const key=`${f.url}:${f.lineNumber+1} ${f.functionName||'(anonymous)'}`;const dt=(profile.timeDeltas[i]||0)/1000;byFn.set(key,(byFn.get(key)||0)+dt);byFile.set(f.url||f.functionName||'(native)',(byFile.get(f.url||f.functionName||'(native)')||0)+dt);}const sort=m=>[...m].map(([name,ms])=>({name,ms})).sort((a,b)=>b.ms-a.ms);return {topFunctions:sort(byFn).slice(0,35),files:sort(byFile),sampledMs:(profile.endTime-profile.startTime)/1000};}
async function phase(name, action){await cdp.send('Profiler.start');const wall=Date.now();const before=await page.evaluate(()=>({t:performance.now(),i:window.__intervals?.length||0,l:window.__longtasks?.length||0}));if(name==='initial-default') before.t=0;await action();await page.waitForFunction(()=>['ready','error'].includes(document.querySelector('#appStatus')?.dataset.state),null,{timeout:180000});await page.waitForTimeout(100);const snap=await page.evaluate(({i,l,t})=>({intervals:window.__intervals.slice(i),longtasks:window.__longtasks.slice(l).filter(task => task.start >= t),state:window.rnaExplorer.state,build:window.rnaExplorer.manifest.build_id,status:document.querySelector('#appStatus').dataset.state,distribution:window.rnaExplorer.snapshots.distribution?.result.coverage,joint:window.rnaExplorer.snapshots.joint?.result.coverage,heap:performance.memory?.usedJSHeapSize}),before);const wallMs = Date.now() - wall; const {profile}=await cdp.send('Profiler.stop');await writeFile(`${out}/${name}.cpuprofile`,JSON.stringify(profile));const digests = await resultDigests(); const p={name,wallMs,...snap,resultDigests:digests,cpu:summarize(profile)};report.phases.push(p);await writeFile(`${out}/report.json`,JSON.stringify(report,null,2));console.log(JSON.stringify({name,wallMs:p.wallMs,status:p.status,distribution:p.distribution,joint:p.joint,longtasks:p.longtasks,top:p.cpu.topFunctions.slice(0,8),intervals:p.intervals.filter(x=>['app.snapshot','app.snapshotAsync','app.renderFamilyOverview','app.renderJoint','repository.loadFamily'].includes(x.name))}));if(p.status!=='ready')throw Error('Render error');}
try{
 await phase('initial-default',()=>page.goto(base,{waitUntil:'domcontentloaded',timeout:120000}));
 await phase('broad-selection',()=>page.evaluate(async()=>{const a=window.rnaExplorer;await a.setSelection({components:'all',methods:[],resolutionMax:null,resolution:'any',contexts:[],functions:[],subtypes:[],structures:[],puckerStates:[],search:''});}));
 await phase('broad-same-family-joint',()=>page.evaluate(async()=>{const a=window.rnaExplorer;a.state.family2Id=a.state.familyId;a.state.parameter2Id=a.parameters().find(p=>p.id!==a.state.parameterId).id;a.updateSelectors();await a.requestJointOnly();}));
 await phase('joint-palette-only',()=>page.evaluate(async()=>{const a=window.rnaExplorer;a.state.joint.palette='viridis';await a.requestJointOnly();}));
 await phase('main-trace-style-only',()=>page.evaluate(async()=>{await window.rnaExplorer.setDisplay({traceStyle:'line'});}));
 assert.equal(report.errors.length, 0, 'Browser page errors');
 for (const phase of report.phases) if (process.env.RNA_EXPECT_BUILD_ID) assert.equal(phase.build, process.env.RNA_EXPECT_BUILD_ID);
 for (const [file, hash] of Object.entries(report.sourceSha256)) assert.equal(report.loadedSourceSha256[file], hash, `Loaded source differs from recorded source: ${file}`);
 report.finished=new Date().toISOString();
}finally{await writeFile(`${out}/report.json`,JSON.stringify(report,null,2));await browser.close();}
