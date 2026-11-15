"""Command-line entry point for synthetic demo or explicitly selected ROOT input."""
from __future__ import annotations
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import numpy as np
from .simulation import toy_pseudodata
from .statistics import fit_local_signal, background_bin_counts, sideband_mask


def plot_mass(counts, edges, out, *, title, source, fit=None, blinded=False):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    centers = (edges[:-1] + edges[1:]) / 2.0
    errs = np.sqrt(counts)
    fig, ax = plt.subplots(figsize=(11, 6))
    valid = sideband_mask(edges) if blinded else np.ones_like(counts, dtype=bool)
    ax.errorbar(centers[valid], counts[valid], yerr=errs[valid], fmt='o', color='#16213b',
                capsize=2, ms=4, label='Observed counts (sidebands only)' if blinded else 'Observed counts')
    if blinded:
        ax.axvspan(120,130, color='#ef4444', alpha=0.13, label='Blinded 120–130 GeV')
    if fit is not None:
        bkg = background_bin_counts(edges, fit.log_background_rate_at_125, fit.background_slope_per_GeV)
        ax.stairs(bkg, edges, color='#0ea5a4', lw=2, label='Fitted background (illustrative)')
    ax.set(xlabel=r'Diphoton invariant mass $m_{\gamma\gamma}$ [GeV]', ylabel='Events / 1 GeV', title=title)
    ax.legend(loc='upper right'); ax.grid(alpha=0.17)
    fig.text(0.01, 0.015, source, size=9, color='#64748b')
    fig.tight_layout(rect=(0, .025, 1, 1))
    fig.savefig(out, dpi=180, bbox_inches='tight')
    plt.close(fig)


def main(argv=None):
    p = argparse.ArgumentParser(description='Reproducible H→γγ educational toolkit; no ATLAS ROOT data are bundled.')
    p.add_argument('--mode', choices=['demo','root'], default='demo')
    p.add_argument('--input', nargs='*', help='Local ATLAS ROOT input files (mode=root)')
    p.add_argument('--tree', default='mini')
    p.add_argument('--units', choices=['MeV','GeV'], default='MeV')
    p.add_argument('--output', default='examples/generated')
    p.add_argument('--unblind', action='store_true', help='Explicitly release 120–130 GeV in real-data plots; exploratory only')
    p.add_argument('--seed', type=int, default=42)
    a = p.parse_args(argv)
    out = Path(a.output); out.mkdir(parents=True,exist_ok=True)
    if a.mode == 'demo':
        masses = toy_pseudodata(seed=a.seed)
        source = 'SYNTHETIC pseudodata — no ATLAS events or detector simulation'
        unblind = True
    else:
        if not a.input:
            p.error('--input ROOT files required in --mode root')
        from .root_io import iter_diphoton_masses
        batches = list(iter_diphoton_masses(a.input, tree=a.tree, units=a.units))
        masses = np.concatenate(batches) if batches else np.empty(0, dtype=float)
        source = 'ATLAS-style ROOT input (user supplied); exploratory analysis — see docs'
        unblind = a.unblind
    edges = np.linspace(100,160,61)
    counts, _ = np.histogram(masses,bins=edges)
    fit = fit_local_signal(counts, edges) if unblind else None
    plot_mass(counts,edges,out/'mass_spectrum.png',title=('Synthetic H→γγ mass demonstration' if a.mode=='demo' else 'ATLAS-style H→γγ analysis'),source=source,fit=fit,blinded=not unblind)
    # Don't publish unblinded real counts by default; outputs are explicitly labeled.
    with (out/'mass_histogram.csv').open('w',encoding='utf-8') as f:
        f.write('mass_low_GeV,mass_high_GeV,counts,blinded\n')
        for i,n in enumerate(counts):
            blind = a.mode=='root' and not unblind and not sideband_mask(edges)[i]
            f.write(f'{edges[i]:.1f},{edges[i+1]:.1f},{"" if blind else int(n)},{int(blind)}\n')
    report = {'source': a.mode, 'input_files': [str(x) for x in a.input] if a.input else [],
              'data_status': 'synthetic-demo' if a.mode=='demo' else 'user-provided-root',
              'unblinded': bool(unblind), 'hypothesis_mass_GeV': 125.09,
              'fit_status': ('exploratory-fixed-mass' if fit else 'not-run-blinded'),
              'model_limitations': 'Fixed Gaussian shape; exponential background; no systematics, trials factor, or official ATLAS likelihood.'}
    if fit:
        report['fit'] = asdict(fit)
    (out/'analysis_metadata.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print('Wrote', out)
    if fit:
        print(f"Exploratory local asymptotic Z={fit.z_local_asymptotic:.3f} (not a discovery significance)")

if __name__ == '__main__':
    main()
