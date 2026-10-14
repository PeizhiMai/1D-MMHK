import fs from 'node:fs/promises';
import path from 'node:path';
import {pathToFileURL} from 'node:url';
import {Presentation,PresentationFile,FileBlob} from '@oai/artifact-tool';
const ROOT='/Users/cosdis/Desktop/projects/HK_orbital_cluster/1d';
const SKILL='/Users/cosdis/.codex/plugins/cache/openai-primary-runtime/presentations/26.1011.11113/skills/presentations';
const WORK=ROOT+'/response/results/twist50_common_mu_U0to12_n0p7_0p8_0p9/analysis_ba_step0p2_eta0p1';
const BUILD=WORK+'/.ppt_build',OUT=WORK+'/presentation';
await fs.mkdir(BUILD,{recursive:true});await fs.mkdir(OUT,{recursive:true});
const {resolvePresentationFont,applyPresentationChartFont,finalizePresentation}=await import(pathToFileURL(SKILL+'/container_tools/artifact_tool_utils.mjs').href);
const FONT=resolvePresentationFont();
const d=JSON.parse(await fs.readFile(WORK+'/presentation_data.json','utf8'));
const P=Presentation.create({slideSize:{width:1280,height:720}});
const colors={charge:'#c43c64',spin:'#234f91'};
function text(s,value,x,y,w,h,size=22,color='#222222'){
 const o=s.shapes.add({geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',line:{style:'solid',fill:'none',width:0}});
 o.text=value;o.text.style={typeface:FONT,fontSize:size,color};return o;
}
function chart(s,title,series,position,xmax,ymax,big=false){
 const c=s.charts.add('scatter',{position,title,titleTextStyle:{fontSize:big?24:18,fill:'#222222'},
  titlePlacement:'aboveChart',series,hasLegend:false,scatterOptions:{style:'lineWithMarkers'},chartFill:'#ffffff',plotAreaFill:'#ffffff',
  xAxis:{min:0,max:xmax,majorUnit:xmax===12?3:1,numberFormatCode:'0',line:{style:'solid',fill:'#333333',width:.8},textStyle:{fontSize:big?17:13},majorGridlines:null},
  yAxis:{min:0,max:ymax,majorUnit:big?.5:.2,numberFormatCode:'0.0',line:{style:'solid',fill:'#333333',width:.8},textStyle:{fontSize:big?17:13},majorGridlines:{style:'solid',fill:'#eeeeee',width:.5}}});
 applyPresentationChartFont(c,{fontFamily:FONT});return c;
}
// Excel stores 15 significant digits. Retain 12 here; full-precision source
// arrays remain in presentation_data.json and the spectrum NPZ files.
function line(ch,name,x,y,marker='none',only=false){return {name,xValues:x.map(v=>Number(v.toPrecision(12))),values:y.map(v=>Number(v.toPrecision(12))),
 line:{style:ch==='spin'?'dashed':'solid',fill:only?'none':colors[ch],width:only?0:1.6},
 marker:{symbol:marker,size:5,fill:colors[ch]}};}
let s=P.slides.add();s.background.fill='#ffffff';
text(s,'Spin and charge velocities',42,20,1196,60,43);
text(s,'10-MMHK, 50 twists, η/t = 0.1, q₁ = π/5',44,83,1150,34,23);
text(s,'Charge: red     Spin: blue dashed     BA: lines     Peak/q: markers',44,125,1150,33,21);
for(let i=0;i<3;i++){
 const n=[.7,.8,.9][i],r=d.rows.filter(r=>r.target_density===n),b=d.ba.filter(r=>r.density===n),series=[];
 for(const ch of ['charge','spin']){
  series.push(line(ch,'BA '+ch,b.map(r=>r.U_over_t),b.map(r=>r['v_'+ch])));
  series.push(line(ch,'Peak/q '+ch,r.map(r=>r.U_over_t),r.map(r=>r[ch+'_dominant_velocity']),ch==='charge'?'circle':'square',true));
  const changed=r.filter(r=>Math.abs(r.actual_density-n)>1e-10);
  if(changed.length)series.push(line(ch,'BA actual density '+ch,changed.map(r=>r.U_over_t),changed.map(r=>r[ch+'_BA_velocity_actual_n']),'diamond',true));
 }
 chart(s,'n = '+n,series,{left:36+i*413,top:172,width:403,height:414},12,2.55,true);
}
text(s,'U/t',606,590,120,30,23);
text(s,'Velocity unit: t a₀/ℏ. Peak/q gives a finite-q, broadening-dependent estimate.',45,627,1190,30,20);
text(s,'BA grid: Δ(U/t) = 0.2. Diamonds use actual density at U/t = 7, 9, 10 near n = 0.9.',45,664,1190,30,18);
s.speakerNotes.text='BA source: Essler et al., The One-Dimensional Hubbard Model (Cambridge, 2005), Eqs.5.102–5.105 and6.44–6.47. BA uses fixed physical n=Ne/L, T=B=0. Actual near-n0.9 densities: U7 .896, U9 .902, U10 .898. Neutral source pole data: '+WORK+'/velocity_comparison.csv. No chemical-potential shift in neutral frequencies.';
for(const n of [.7,.8,.9]){
 s=P.slides.add();s.background.fill='#ffffff';
 text(s,`Spin and charge response at n${n===.9?' ≈':' ='} ${n}`,42,16,1196,57,43);
 text(s,'10-MMHK, 50 twists, η/t = 0.1, q₁ = π/5',44,75,1150,33,22);
 text(s,'Charge: red solid     Spin (2Sᶻ): blue dashed     Common absolute intensity scale',44,112,1165,32,21);
 for(let U=0;U<=12;U++){
  const series=['charge','spin'].map(ch=>{const z=d.responses.find(z=>z.n===n&&z.U===U&&z.channel===ch);return line(ch,ch,z.x,z.y);});
  const r=d.rows.find(r=>r.target_density===n&&r.U_over_t===U);
  const title=`U/t = ${U}`+(Math.abs(r.actual_density-n)>1e-10?`   n̄ = ${r.actual_density.toFixed(3)}`:'');
  chart(s,title,series,{left:28+(U%5)*249,top:153+Math.floor(U/5)*166,width:240,height:158},2,Math.ceil(d.ymax/.05)*.05);
 }
 text(s,'Frequency axis: ω/t     Intensity: S(q₁,ω). This slide shows ω/t = 0–2.',44,668,1190,32,19);
 s.speakerNotes.text='13 panels cover U/t=0..12. Physical n=Ne/L. Coherent densities rho_c/s=sum_j exp(iqj)(nup +/- ndn)/sqrt(L); spin is 2Sz, divide spin intensity by4 for Sz. Same common chemical potentials and complete spin/translation ground ensembles as original50-twist campaign. Eta0.1 rebroadened original poles and re-extracted peak positions. Full range figures (omega/t0..65) accompany the deck. Literal source values: '+WORK+'/presentation_data.json';
}
const candidate=BUILD+'/candidate.pptx';await (await PresentationFile.exportPptx(P)).save(candidate);
const final=OUT+'/spin_charge_10MMHK_50twists_eta0p1_v2.pptx';
await finalizePresentation({workspaceDir:ROOT,candidatePath:candidate,finalPath:final,
 pythonExecutable:'/Users/cosdis/.venvs/myenv/bin/python',
 integrityValidatorPath:SKILL+'/container_tools/inspect_presentation_package_integrity.py',
 layoutValidatorPath:SKILL+'/container_tools/inspect_presentation_layout_geometry.py',
 layoutArgs:['--expected-slide-size-emu','12192000,6858000','--validate-heading-fit'],
 explicitTotalSlideCount:4,requiredNativeChartOwnerSlides:[1,2,3,4],requiredNativeTableOwnerSlides:[],
 materializeLiteralChartWorkbooks:true,fontPolicy:{basis:'design',families:[FONT]},
 verifyArtifactToolImport:true,receiptPath:BUILD+'/validation_v2.json'});
const FINAL=await PresentationFile.importPptx(await FileBlob.load(final));
for(const [i,sl] of FINAL.slides.items.entries()){
 const png=await FINAL.export({slide:sl,format:'png',scale:1.5});
 await fs.writeFile(BUILD+`/slide${i+1}.png`,new Uint8Array(await png.arrayBuffer()));
}
console.log('PPTX:',final,'Font:',FONT);
