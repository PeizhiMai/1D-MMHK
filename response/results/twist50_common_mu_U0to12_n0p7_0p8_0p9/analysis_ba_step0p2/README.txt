Replot of the existing L=10, 50-twist U/t=0..12 spin/charge response campaign.

BA: 61 directly calculated U points per density, Delta(U/t)=0.2.
The U grid is five times finer than the original step-1 grid when step=0.2.
Lines: fixed nominal n=Ne/L=0.7,0.8,0.9. Not splines of the old 13 points.
Open diamonds: BA at actual density when different (U7:0.896, U9:0.902, U10:0.898).
No density interpolation and no arbitrary ensemble mixing introduced.
BA mesh convergence checked at every point, using 400/600 base quadrature nodes.
Source: Essler et al., The One-Dimensional Hubbard Model, Eqs.5.102--5.105 and6.44--6.47.
Thermodynamic T=B=0; u=U/(4t); t=a0=hbar=1.

Response: one 13-panel figure for each density; spin and charge overlaid in every panel.
U/t=0 through12, eta/t=0.2, q_min=pi/5. Spin operator is longitudinal2Sz.
Identical absolute-intensity scales across densities and U; no per-panel normalization.
response_panel_n*: low-energy omega/t=0..2 view with dominant-peak dots.
response_panel_n*_full_range: full saved omega/t=0..65, log-intensity supplemental view.
The log view displays intensities down to1e-6; smaller stored tails remain in source NPZs.
Peak/q estimates retain their original finite-q/eta interpretation and values.

Original analysis, native solver, raw poles, spectra and open LaTeX note were not changed.
No Hubbard or student solver was run. All cached spectra were checked against original poles.
Files: PNG/SVG figures; BA_velocities_fine_grid.csv; unchanged velocity_comparison.csv;
validation.json and metadata.json. Source spectrum NPZs remain in ../analysis.
