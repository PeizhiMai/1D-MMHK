#!/usr/bin/env python3
"""Slide-sized spectral panels from the completed 101-k arrays; no solver."""
import argparse
import hashlib
import json
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
CAMPAIGN = HERE / 'results/twist50_common_mu_U0to12_n0p7_0p8_0p9'
SPECTRA = CAMPAIGN / 'single_particle_eta0p1_k101_depth600'
BUILD = CAMPAIGN / 'analysis_ba_step0p2_eta0p1/.ppt_build/seven_slides'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=BUILD / 'spectral_assets')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault('MPLCONFIGDIR', str(HERE / '.mplconfig'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.font_manager import findfont

    findfont('Helvetica Neue', fallback_to_default=False)
    plt.rcParams.update({'font.family': 'Helvetica Neue', 'font.size': 10})
    inputs = {}
    panels = []
    for n in [.9, .8, .7]:
        tag = f'{n:g}'.replace('.', 'p')
        fig, axes = plt.subplots(3, 5, figsize=(12.24, 4.98), sharex=True, sharey=True)
        fig.subplots_adjust(left=.052, right=.987, bottom=.10, top=.965,
                            wspace=.16, hspace=.28)
        for U, ax in zip(range(13), axes.flat):
            path = SPECTRA / f'U{U}' / f'n{tag}' / 'Akw.npz'
            inputs[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
            with np.load(path, allow_pickle=False) as d:
                k, omega, A = d['k'], d['omega'], d['A_per_spin']
                np.testing.assert_array_equal(k, np.linspace(0, np.pi, 101))
                assert A.shape == (101, 2201)
                np.testing.assert_allclose(d['eta'], .1, atol=1e-12)
                assert np.isfinite(A).all() and np.all(A >= 0)
                nr = float(d['actual_density'])
                im = ax.pcolormesh(k, omega, A.T, shading='nearest',
                                  cmap='rainbow', vmin=0, vmax=1.5, rasterized=True)
            ax.set(xlim=(0, np.pi), ylim=(-6, 10),
                   xticks=[0, np.pi/2, np.pi], xticklabels=['0', 'π/2', 'π'],
                   yticks=[-5, 0, 5, 10])
            ax.tick_params(labelsize=10, length=3, width=.7, pad=3)
            ax.axhline(0, color='white', lw=.5, ls='--')
            ax.text(.035, .945, f'U/t = {U}', transform=ax.transAxes,
                    color='white', fontsize=11.5, va='top')
            if abs(nr-n) > 1e-10:
                ax.text(.96, .945, f'n̄ = {nr:.3f}', transform=ax.transAxes,
                        color='white', fontsize=9.5, va='top', ha='right')
        for ax in list(axes.flat)[13:]:
            ax.set_visible(False)
        # U8/U9 have no live axes directly below them.
        for ax in axes[1, 3:]:
            ax.tick_params(labelbottom=True)
        fig.supylabel('ω/t', x=.003, fontsize=12)
        fig.supxlabel('k', y=.005, fontsize=12)
        # Use the two unused panel positions, keeping the heatmaps tall.
        cax = fig.add_axes([.654, .22, .309, .032])
        cb = fig.colorbar(im, cax=cax, orientation='horizontal', extend='max',
                         ticks=[0, .5, 1, 1.5])
        cb.set_label('A(k,ω) per spin', fontsize=11, labelpad=4)
        cb.ax.tick_params(labelsize=10, length=3)
        stem = args.output / f'spectral_slide_n{tag}'
        fig.savefig(stem.with_suffix('.png'), dpi=300, facecolor='white')
        fig.savefig(stem.with_suffix('.svg'), facecolor='white')
        plt.close(fig)
        panels.append(dict(target_density=n, U=list(range(13)), k_points=101,
                           eta=.1, displayed_omega_range=[-6, 10],
                           stored_omega_range=[-6, 16], color_range=[0, 1.5],
                           png=str(stem.with_suffix('.png'))))
    (args.output / 'inputs.json').write_text(json.dumps(dict(
        input_sha256=inputs, panels=panels, many_body_solver_launched=False,
        interpolated=False), indent=2)+'\n')
    print(f'Created three slide-sized spectral panels: {args.output}')


if __name__ == '__main__':
    main()
