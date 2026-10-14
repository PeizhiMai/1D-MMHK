// Public, laptop-portable exporter; no Codex runtime or workstation paths.
import fs from 'node:fs';
import path from 'node:path';
import {createRequire} from 'node:module';
const require = createRequire(process.env.MMHK_SLIDES_ROOT
  ? path.join(process.env.MMHK_SLIDES_ROOT, 'package.json') : import.meta.url);
const PptxGenJS = require('pptxgenjs');
const [input, output] = process.argv.slice(2);
if (!input || !output) throw new Error('Usage: node export.mjs presentation_data.json OUTPUT');
const d = JSON.parse(fs.readFileSync(input, 'utf8'));
const pptx = new PptxGenJS();
pptx.layout = 'LAYOUT_WIDE';
pptx.author = '1D-MMHK reproducible pipeline';
pptx.subject = 'Physical Hubbard TBC spin/charge responses and single-particle spectra';
pptx.title = `${d.config.L}-MMHK spin, charge and spectra`;
pptx.lang = 'en-US';
pptx.theme = {headFontFace:'Arial', bodyFontFace:'Arial', lang:'en-US'};
const tag = n => String(n).replace('.', 'p');
function slide(title, image, note) {
  const s = pptx.addSlide(); s.background = {color:'FFFFFF'};
  s.addText(title, {x:.45,y:.18,w:12.4,h:.46,fontSize:24,bold:true,color:'17283A',margin:0});
  const p = path.resolve(output,image);
  if (!fs.existsSync(p)) throw new Error(`Missing plot ${p}`);
  s.addImage({path:p,...pptxgenImageSizingContain(p,.38,.78,12.55,6.10)});
  s.addText(`${d.config.L} orbitals · ${d.config.twists} density twists · η/t = ${d.config.eta}`, {
    x:.45,y:7.12,w:11,h:.18,fontSize:10,color:'596675',margin:0});
  s.addText(String(pptx._slides.length), {x:12.4,y:7.1,w:.4,h:.2,fontSize:10,align:'right',margin:0});
  s.addNotes(note);
}
// Public helper exports differ between CommonJS versions. Use instance API.
function pptxgenImageSizingContain(p,x,y,w,h) {return pptxgenSizing(p,x,y,w,h);}
const pptxgenSizing = pptx.constructor === undefined ? null : (p,x,y,w,h) => {
  // PNG IHDR gives dimensions without a private runtime or image package.
  const b=fs.readFileSync(p); if(b.toString('hex',0,8)!=='89504e470d0a1a0a')throw new Error('Expected PNG');
  const iw=b.readUInt32BE(16),ih=b.readUInt32BE(20),scale=Math.min(w/iw,h/ih);
  return {x:x+(w-iw*scale)/2,y:y+(h-ih*scale)/2,w:iw*scale,h:ih*scale};
};
slide('Spin and charge velocities','spin_charge_velocity_slide.png',
 'BA curves solve thermodynamic Hubbard integral equations at nominal density. MMHK markers are dominant positive-frequency peak divided by q_min=2*pi/L: a finite-q, broadening-dependent estimator, not an extrapolated sound velocity. Diamonds show BA at the actual density when different. No intensity normalization.');
for(const n of [...d.config.physical_densities].sort((a,b)=>b-a)) {
  slide(`Spin and charge response · target n = ${n}`, `response_slide_n${tag(n)}.png`,
   'Coherent charge and longitudinal spin (2Sz) positive inelastic Lehmann spectra. One common chemical potential, complete ground-state trace per twist, not a fixed number or Sz. Elastic weight is excluded. The CSV files report actual densities and peak diagnostics; broadening eta/t='+d.config.eta+'.');
  slide(`Single-particle spectrum · target n = ${n}`, `spectral_slide_n${tag(n)}.png`,
   `Per-spin A(k,omega)=-Im(Gup+Gdn)/(2*pi). ${d.config.k_points} directly calculated k points from 0 to pi, no interpolation. Common mu fixed by the density twist grid. The display is clipped to [${d.config.plot_y_min},${d.config.plot_y_max}], not the stored frequency range. Colormap saturation is 1.5.`);
}
fs.mkdirSync(output,{recursive:true});
const name='spin_charge_spectral.pptx';
await pptx.writeFile({fileName:path.resolve(output,name)});
fs.writeFileSync(path.resolve(output,'presentation_validation.json'),JSON.stringify({
  slide_count:1+2*d.config.physical_densities.length, exporter:'PptxGenJS 4.0.1 (public)',
  plot_assets:'unmodified scientific PNG plots; editable slide text and speaker notes',
  source_data:'presentation_data.json and plot_inputs.json', file:name},null,2)+'\n');
console.log(`Wrote ${path.resolve(output,name)}`);
