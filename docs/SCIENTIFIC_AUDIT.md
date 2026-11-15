# Scientific audit of the historical Higgs diphoton project

Reviewed source: the **current GitHub `main` snapshot** of `BurhanKARADUMAN/Higgs-Diphoton-Analysis`, especially `higgs_diphoton_analysis.C`, `analysis_results.txt`, `draw_pvalue_master.C`, `higgs_analysis.C`, and the TMVA weights directory. The original macros and compiled report remain untouched. This document describes source-visible issues, not a re-run of the original ATLAS ROOT data.

## Evidence available

- `analysis_results.txt`: log with data/MC entry counts, TMVA training details, AUC, threshold scan, fit results and performance table.
- `higgs_diphoton_analysis.C`: ATLAS-like event selections, C++ ROOT/TMVA training, BDT scan and sideband / signal fits.
- `dataset/weights/TMVAClassification_BDT.weights.xml`: pretrained weights for one historical experiment.
- Other `*.C` macros: multiple historical approaches and input schemas; they are **not interchangeable** without checking variable names, units and training procedure.
- `higgs_diphoton_report_compiled (2).pdf`: report present in repo; **not independently replicated** with raw ROOT inputs in this upgrade.

## Important observations

1. **The reported `3.49` is S/sqrt(B), not a calibrated discovery significance.** The analysis uses a fitted signal plus background to find the BDT cut maximizing S/sqrt(B) on an overlapping data sample. Cut scanning changes the null distribution, and model uncertainty / look-elsewhere effects are not quantified.
2. **Pre-BDT signal baseline `S0=0.0`.** This makes BDT efficiency and purity improvement ratios based on `S0` ill-defined or meaningless. The log's very large printed BDT efficiency is an evident symptom; it must not be cited as performance.
3. **Fitted `125.54 GeV` is a fitted peak location.** The logged `0.93 GeV` is the Gaussian *shape sigma*, **not** an uncertainty on the mass measurement.
4. **TMVA AUC=0.703** is a historical classifier diagnostic, not equivalent to a particle signal p-value. The log says negative-weight training events were present and TMVA raised a warning. Confirm negative-weight treatment and normalization before claiming physics performance.
5. **The original macro uses `pt1_norm` / `pt2_norm`** (photon pT divided by mγγ) as classifier features. These may sculpt the mass spectrum. The new BDT validation demo rejects direct mass / pT-over-mass features by default; further decorrelation studies remain necessary.
6. **Potential training/cut-selection overlap.** The root macro's scan loops over scored events and refits the same mass distribution to choose `bestT`; it does not establish an independent unblinded holdout for the quoted 3.49.
7. **`FullData` was bypassed in a recorded run.** The log explicitly states that a format-corrupted `FullData` sample was not used and another data sample was substituted. Capture exact data provenance, tree names and hashes before research publication.
8. **Legacy macro/schema differences.** For example, `draw_pvalue_master.C` reads variables `pt1,eta1,phi1,pt2,eta2,phi2` and its own XML path, whereas the integrated macro registers a different feature set. Do not reuse XML weights across these scripts without auditing variable order, transforms and units.
9. **Approximate sliding-window p-values.** The legacy local p-value macro uses a fitted polynomial, a normal Z approximation and a scan; these p-values are neither a complete likelihood ratio nor look-elsewhere corrected. Fit uncertainty and binning choices matter.
10. **Data not packaged.** The repository contains no verified complete ROOT event inputs for reproducing all recorded figures. Some compiled `.so`, `.d`, `.pcm` byproducts are present; they should not be treated as portable analysis source.

## What can safely be reported now?

The repository demonstrates practical experience with ROOT, TMVA, diphoton kinematics, BDT methods, invariant-mass spectroscopy and exploratory fits. The old run log has **historical reported** peak parameters and AUC. **No claim of a new Higgs discovery, official mass measurement, or calibrated discovery significance** should be made from this archive without the follow-up validation below.

## Minimum validation before a numerical Higgs inference

1. Verify the exact ATLAS Open Data release, ROOT files and event weights.
2. Predefine trigger, ID, isolation and η/pT cuts; implement transparent cut flow and stable MeV-to-GeV conversions.
3. Isolate signal training MC from background training data; avoid signal-region data for optimization; use independent test samples.
4. Validate mass-feature correlations and post-BDT background shape.
5. Specify a fixed signal mass/resolution hypothesis, a physically defensible background model and nuisance constraints.
6. Check fit convergence, goodness-of-fit, alternative background families, spurious signal and toys.
7. Compute local test statistics under a specified null; treat scan-dependent trials / global p-value separately.
8. Report confidence intervals with statistical and systematic components; preserve full commands, software versions and output hashes.

## Where new code lives

The new `higgs_review/` package is a **transparent educational review toolkit**. It does not silently replace the C++ ROOT publication path or use the historical XML model. It includes tests and obvious safeguards. ROOT file analysis and model comparison must still be validated with real data.
