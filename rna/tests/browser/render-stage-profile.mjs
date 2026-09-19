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
const report={started:new Date().toISOString(),phases:[],errors:[],sourceSha256:{},limitations:'Single instrumented shared-host run; wall intervals inclusive and overlapping; CPU sampling approximate, no universal performance claim.'};
for(const file of ['app/PureRnaExplorer.js','core/selection.js','core/analysis.js','core/export.js','core/repository.js','math/numeric.js']) report.sourceSha256[file]=createHash('sha256').update(await readFile(new URL('../../'+file, import.meta.url))).digest('hex');
const browser=await chromium.launch({headless:true});
const page=await browser.newPage({viewport:{width:1440,height:1000}});
page.on('pageerror', e=>report.errors.push(e.message));
await page.addInitScript(()=>{
 window.__intervals=[];window.__longtasks=[];
 new PerformanceObserver(list=>{for(const e of list.getEntries())window.__longtasks.push({start:e.startTime,duration:e.duration});}).observe({type:'longtask',buffered:true});
 window.__measurePrototype=(proto,prefix,methods)=>{for(const name of methods){const f=proto[name];if(typeof f!=='function')continue;proto[name]=function(...args){const start=performance.now(); const detail=name==='plot'?(args[0]?.id||args[0]?.className):name==='snapshot'?args[0]?.result?.kind:name==='loadFamily'?args[0]:'';let v;const done=()=>window.__intervals.push({name:prefix+'.'+name,detail,start,duration:performance.now()-start});try{v=f.apply(this,args);}catch(e){done();throw e;}if(v&&typeof v.then==='function')return v.then(r=>{done();return r;},e=>{done();throw e;});done();return v;};}};
});
await page.route(base+'main.js',async route=>{
 const response=await route.fetch();let body=await response.text();
 body=body.replace("const root =", `window.__measurePrototype(PureRnaExplorer.prototype,'app',['start','render','snapshot','renderFamilyOverview','renderJoint','renderTable','updateContexts','updatePuckerControls','updateInteractionControls','requestJointOnly']);
 window.__measurePrototype(Object.getPrototypeOf(PureRnaExplorer.prototype),'base',['plot']);
 window.__measurePrototype(RnaDataRepository.prototype,'repository',['readJson','loadFamily','loadFamilyBundle','loadMetadata']);
 const root =`);
 await route.fulfill({response,body});
});
const cdp=await page.context().newCDPSession(page);await cdp.send('Profiler.enable');await cdp.send('Profiler.setSamplingInterval',{interval:1000});
function summarize(profile){const nodes=new Map(profile.nodes.map(n=>[n.id,n]));const byFn=new Map(),byFile=new Map();for(let i=0;i<(profile.samples||[]).length;i++){const f=nodes.get(profile.samples[i])?.callFrame||{};const key=`${f.url}:${f.lineNumber+1} ${f.functionName||'(anonymous)'}`;const dt=(profile.timeDeltas[i]||0)/1000;byFn.set(key,(byFn.get(key)||0)+dt);byFile.set(f.url||f.functionName||'(native)',(byFile.get(f.url||f.functionName||'(native)')||0)+dt);}const sort=m=>[...m].map(([name,ms])=>({name,ms})).sort((a,b)=>b.ms-a.ms);return {topFunctions:sort(byFn).slice(0,35),files:sort(byFile),sampledMs:(profile.endTime-profile.startTime)/1000};}
async function phase(name, action){await cdp.send('Profiler.start');const wall=Date.now();const before=await page.evaluate(()=>({t:performance.now(),i:window.__intervals?.length||0,l:window.__longtasks?.length||0}));await action();await page.waitForFunction(()=>['ready','error'].includes(document.querySelector('#appStatus')?.dataset.state),null,{timeout:180000});await page.waitForTimeout(100);const snap=await page.evaluate(({i,l})=>({intervals:window.__intervals.slice(i),longtasks:window.__longtasks.slice(l),state:window.rnaExplorer.state,build:window.rnaExplorer.manifest.build_id,status:document.querySelector('#appStatus').dataset.state,distribution:window.rnaExplorer.snapshots.distribution?.result.coverage,joint:window.rnaExplorer.snapshots.joint?.result.coverage,heap:performance.memory?.usedJSHeapSize}),before);const {profile}=await cdp.send('Profiler.stop');await writeFile(`${out}/${name}.cpuprofile`,JSON.stringify(profile));const p={name,wallMs:Date.now()-wall,...snap,cpu:summarize(profile)};report.phases.push(p);await writeFile(`${out}/report.json`,JSON.stringify(report,null,2));console.log(JSON.stringify({name,wallMs:p.wallMs,status:p.status,distribution:p.distribution,joint:p.joint,longtasks:p.longtasks,top:p.cpu.topFunctions.slice(0,8),intervals:p.intervals.filter(x=>['app.snapshot','app.renderFamilyOverview','app.renderJoint','repository.loadFamily'].includes(x.name))}));if(p.status!=='ready')throw Error('Render error');}
try{
 await phase('initial-default',()=>page.goto(base,{waitUntil:'domcontentloaded',timeout:120000}));
 await phase('broad-selection',()=>page.evaluate(async()=>{const a=window.rnaExplorer;await a.setSelection({components:'all',methods:[],resolutionMax:null,resolution:'any',contexts:[],functions:[],subtypes:[],structures:[],puckerStates:[],search:''});}));
 await phase('broad-same-family-joint',()=>page.evaluate(async()=>{const a=window.rnaExplorer;a.state.family2Id=a.state.familyId;a.state.parameter2Id=a.parameters().find(p=>p.id!==a.state.parameterId).id;a.updateSelectors();await a.requestJointOnly();}));
 await phase('joint-palette-only',()=>page.evaluate(async()=>{const a=window.rnaExplorer;a.state.joint.palette='viridis';await a.requestJointOnly();}));
 await phase('main-trace-style-only',()=>page.evaluate(async()=>{await window.rnaExplorer.setDisplay({traceStyle:'line'});}));
 assert.equal(report.errors.length, 0, 'Browser page errors');
 for (const phase of report.phases) if (process.env.RNA_EXPECT_BUILD_ID) assert.equal(phase.build, process.env.RNA_EXPECT_BUILD_ID);
 report.finished=new Date().toISOString();
}finally{await writeFile(`${out}/report.json`,JSON.stringify(report,null,2));await browser.close();}
