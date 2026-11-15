# Reproduce the new educational analysis

## Environment

Python >= 3.10. On Windows, open PowerShell in repository root. On Linux/WSL run equivalent commands.

```bash
python -m venv .venv
# Windows: .\.venv\Scripts\Activate.ps1
# Linux/WSL: source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e '.[test]'
python -m pytest -q
higgs-review --mode demo --seed 42 --output examples/generated
```

This creates `mass_spectrum.png`, `mass_histogram.csv`, `analysis_metadata.json`. They are **SYNTHETIC**; expected Z values vary with seed, signal injection and model choices, not with ATLAS data.

## ROOT inputs (optional)

```bash
python -m pip install -e '.[root]'
higgs-review --mode root --input /absolute/path/to/data_A.GamGam.root --tree mini --units MeV --output results/real_blinded
```

The default omits 120–130 GeV from output CSVs and graphics. The CLI expects `photon_pt`, `photon_eta`, `photon_phi`, `photon_isTightID`, `photon_etcone20` and `trigP`. It raises an error when critical branches are missing rather than silently dropping the corresponding cut. The real-data path has not been integration-tested on the original user's ROOT inputs in this build.

## BDT holdout check (optional)

```bash
python -m pip install -e '.[ml]'
higgs-bdt-check --csv /path/to/your_labeled_features.csv --out results/holdout.json
```

Required CSV columns: `aeta1,aeta2,adeta,adphi,costh,iso1,iso2,label` (`label=1` signal, `0` background). No direct mass / pT-over-mass features by default. Carefully document class normalization, original files and any event-level correlations. This is a diagnostic tool, not an official TMVA model or significance estimate.

## Legacy ROOT macros

Original `higgs_diphoton_analysis.C` and TMVA XML models have not been modified. To reproduce the old results, obtain the original ROOT inputs and check the exact tree names, paths and features in those files before executing the ROOT macros. They are not tested in this upgrade environment because ROOT and collision data are not installed here.
