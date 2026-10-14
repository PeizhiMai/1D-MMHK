import fs from 'node:fs/promises';
import path from 'node:path';
import crypto from 'node:crypto';
import {pathToFileURL} from 'node:url';
import {PresentationFile,FileBlob} from '@oai/artifact-tool';

const ROOT='/Users/cosdis/Desktop/projects/HK_orbital_cluster/1d';
const SKILL='/Users/cosdis/.codex/plugins/cache/openai-primary-runtime/presentations/26.1011.11113/skills/presentations';
const CAMPAIGN=ROOT+'/response/results/twist50_common_mu_U0to12_n0p7_0p8_0p9';
const WORK=CAMPAIGN+'/analysis_ba_step0p2_eta0p1';
const SOURCE=WORK+'/presentation/spin_charge_10MMHK_50twists_eta0p1_v2.pptx';
const BUILD=WORK+'/.ppt_build/seven_slides';
const FINAL=WORK+'/presentation/spin_charge_spectral_10MMHK_eta0p1_7slides.pptx';
const ASSETS=BUILD+'/spectral_assets';
process.env.RUNTIME_NODE_MODULES??='/Users/cosdis/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules';
process.env.RUNTIME_NODE??='/Users/cosdis/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/bin/node';
await fs.mkdir(BUILD,{recursive:true});
const {finalizePresentation,applyPresentationChartFont}=await import(pathToFileURL(SKILL+'/container_tools/artifact_tool_utils.mjs').href);
const sourceSha=crypto.createHash('sha256').update(await fs.readFile(SOURCE)).digest('hex');
const P=await PresentationFile.importPptx(await FileBlob.load(SOURCE));
const before=await P.inspect({select:'$.slides[*]',include:['id','index','layout','title','elements.id','elements.kind','elements.text','elements.title','elements.position','elements.style','notes.text']});
await fs.writeFile(BUILD+'/before.ndjson',before.ndjson);
const original=before.records;
if(original.length!==4)throw new Error('Expected the original four-slide deck');
const velocity=P.resolve(original[0].id);
const titleStyle=original[0].elements.find(e=>e.kind==='textbox').style;
const FONT=titleStyle.typeface;
const findText=value=>{
 const e=original[0].elements.find(e=>e.text===value);
 if(!e)throw new Error('Missing source text: '+value);
 return P.resolve(e.id);
};
findText('Charge: red     Spin: blue dashed     BA: lines     Peak/q: markers').text=
 'Charge: red     Spin: blue dashed     BA: lines     MMHK: markers';
findText('Velocity unit: t a₀/ℏ. Peak/q gives a finite-q, broadening-dependent estimate.').text=
 'Velocity unit: t a₀/ℏ';
velocity.speakerNotes.text+='\nThe MMHK markers are dominant-peak omega(q₁)/q₁ values at q₁=π/5 and eta/t=0.1. They are finite-q, broadening-dependent estimators, rather than a q→0 extrapolation.';
const velocityCharts=original[0].elements.filter(e=>e.kind==='chart');
for(const [i,n] of [.9,.8,.7].entries()){
 const e=velocityCharts.find(e=>e.title==='n = '+n);
 if(!e)throw new Error('Missing density chart');
 P.resolve(e.id).position={...e.position,left:36+i*413};
}

// Explicit repair of the imported chart/workbook lineage. The importer keeps
// all literal caches but loses their workbooks. Recreate only those chart
// objects from the exact cached arrays, retaining the original builder's
// axes, styles, titles and placement. The source workbooks contain no formulas.
const data=JSON.parse(await fs.readFile(WORK+'/presentation_data.json','utf8'));
const colors={charge:'#c43c64',spin:'#234f91'};
const repair=await P.inspect({select:'$.slides[*]',include:['id','index','elements.id','elements.kind','elements.title','elements.position']});
let repaired=0;
for(const slide of repair.records){
 for(const element of slide.elements.filter(e=>e.kind==='chart')){
  const old=P.resolve(element.id);
  const series=old.series.items.map(s=>{
   const ch=s.name.includes('spin')?'spin':'charge';
   const actual=s.name.startsWith('BA actual density');
   const peak=s.name.startsWith('Peak/q');
   const only=actual||peak;
   return {name:s.name,xValues:[...s.xValues],values:[...s.values],
    line:{style:ch==='spin'?'dashed':'solid',fill:only?'none':colors[ch],width:only?0:1.6},
    marker:{symbol:actual?'diamond':peak?(ch==='charge'?'circle':'square'):'none',size:5,fill:colors[ch]}};
  });
  P.delete(element.id);
  const big=slide.index===0;
  const chart=P.resolve(slide.id).charts.add('scatter',{position:element.position,title:element.title,
   titleTextStyle:{fontSize:big?24:18,fill:'#222222'},titlePlacement:'aboveChart',series,hasLegend:false,
   scatterOptions:{style:'lineWithMarkers'},chartFill:'#ffffff',plotAreaFill:'#ffffff',
   xAxis:{min:0,max:big?12:2,majorUnit:big?3:1,numberFormatCode:'0',line:{style:'solid',fill:'#333333',width:.8},textStyle:{fontSize:big?17:13},majorGridlines:null},
   yAxis:{min:0,max:big?2.55:Math.ceil(data.ymax/.05)*.05,majorUnit:big?.5:.2,numberFormatCode:'0.0',line:{style:'solid',fill:'#333333',width:.8},textStyle:{fontSize:big?17:13},majorGridlines:{style:'solid',fill:'#eeeeee',width:.5}}});
  applyPresentationChartFont(chart,{fontFamily:FONT});
  repaired++;
 }
}
if(repaired!==42)throw new Error('Expected exactly 42 native chart repairs');

function text(slide,value,x,y,w,h,size){
 const o=slide.shapes.add({geometry:'textbox',position:{left:x,top:y,width:w,height:h},
  fill:'none',line:{style:'solid',fill:'none',width:0}});
 o.text=value;o.text.style={typeface:FONT,fontSize:size,color:'#222222'};
 return o;
}
const order=[original[0].id];
for(const n of [.9,.8,.7]){
 const tag=String(n).replace('.','p');
 const response=original.find(s=>s.title.startsWith('Spin and charge response')&&s.title.endsWith(String(n)));
 if(!response)throw new Error('Missing response slide for n='+n);
 order.push(response.id);
 const s=P.slides.add({layout:'Blank',background:{fill:'#ffffff'}});
 const spectralTitle=`Single-particle spectra at n${n===.9?' ≈':' ='} ${n}`;
 text(s,spectralTitle,42,16,1196,57,43);
 text(s,'10-MMHK, η/t = 0.1, 101 momenta from 0 to π',44,75,1150,33,22);
 text(s,'U/t = 0–12     A(k,ω) per spin     Common color scale',44,112,1165,32,21);
 const bytes=await fs.readFile(ASSETS+`/spectral_slide_n${tag}.png`);
 s.images.add({blob:bytes,contentType:'image/png',fit:'contain',
  alt:`13 spectral heatmaps at target n=${n}, U/t=0..12, eta/t=0.1, k=0..pi, omega/t=-6..10`,
  position:{left:28,top:153,width:1224,height:498}});
 text(s,'Chemical potentials from the 50-twist density calculation',44,668,1190,32,19);
 s.speakerNotes.text=`Sources: ${CAMPAIGN}/single_particle_eta0p1_k101_depth600/U*/n${tag}/Akw.npz and metadata.json. The 13 panels show U/t=0..12, eta/t=0.1 and exactly 101 directly calculated physical momenta k_j=j*pi/100 without momentum interpolation. The 101-point grid selects20 folded phases, with11 inversion representatives. The common mu and reference density are those of the original50-twist campaign and are not retuned for this coarser spectral sampling. Per-spin A=-Im(Gup+Gdn)/(2*pi), with unit full-frequency weight. The display clips to omega/t=-6..10, while stored NPZ arrays cover -6..16. Addition Em(N+1)-E0-mu; removal E0-Em(N-1)-mu. Color range0..1.5, with an upper saturation arrow. ${n===.9?'Reference densities: U7 .896, U9 .902, U10 .898, otherwise .9. ':''}`;
 const added=await P.inspect({select:'$.slides[*]',include:['id','title']});
 const canonical=added.records.find(record=>record.title===spectralTitle);
 if(!canonical)throw new Error('Cannot resolve the new spectral slide');
 order.push(canonical.id);
}
P.slides.reorder(order);
const after=await P.inspect({select:'$.slides[*]',include:['id','index','layout','title','elements.id','elements.kind','elements.text','elements.title','elements.position','notes.text']});
await fs.writeFile(BUILD+'/after.ndjson',after.ndjson);
const candidate=BUILD+'/candidate.pptx';
await (await PresentationFile.exportPptx(P)).save(candidate);
await finalizePresentation({workspaceDir:ROOT,candidatePath:candidate,finalPath:FINAL,
 pythonExecutable:'/Users/cosdis/.venvs/myenv/bin/python',
 integrityValidatorPath:SKILL+'/container_tools/inspect_presentation_package_integrity.py',
 layoutValidatorPath:SKILL+'/container_tools/inspect_presentation_layout_geometry.py',
 layoutArgs:['--expected-slide-size-emu','12192000,6858000','--validate-heading-fit'],
 explicitTotalSlideCount:7,requiredNativeChartOwnerSlides:[1,2,4,6],requiredNativeTableOwnerSlides:[],
 // Explicit workbook repair: the importer drops externalData relationships.
 // All 42 complete literal caches exactly match the source; the original
 // 42 workbooks contain no formulas. See chart_workbook_repair.json.
 materializeLiteralChartWorkbooks:true,
 fontPolicy:{basis:'reference',families:[FONT],referencePath:SOURCE,referenceSha256:sourceSha},
 verifyArtifactToolImport:true,receiptPath:BUILD+'/validation.json'});
const result=await PresentationFile.importPptx(await FileBlob.load(FINAL));
for(const [i,slide] of result.slides.items.entries()){
 const png=await result.export({slide,format:'png',scale:1.5});
 await fs.writeFile(BUILD+`/slide${i+1}.png`,new Uint8Array(await png.arrayBuffer()));
 const layout=await slide.export({format:'layout'});
 await fs.writeFile(BUILD+`/slide${i+1}.layout.json`,await layout.text());
}
console.log('FINAL_PPTX',FINAL);
