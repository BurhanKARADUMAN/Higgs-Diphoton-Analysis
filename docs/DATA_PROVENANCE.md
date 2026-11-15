# Data provenance and reproducibility

## Current status

- Historical ROOT/TMVA code and weights are present in the original GitHub repository.
- The full collision-event ROOT files used for the archived analysis are **not checked into GitHub**.
- `analysis_results.txt` is a historical textual log; the numerical entries are reported, not independently reproduced by this upgrade.
- `examples/demo/` is **generated synthetic pseudodata** and must never be relabeled as real ATLAS events.

## Where to get real input

The [ROOT 2020 ATLAS H→γγ tutorial](https://root.cern.ch/doc/v632/df104__HiggsToTwoPhotons_8py.html) documents an available ATLAS Open Data mini-tree workflow, including 2016 13 TeV `GamGam/Data` ROOT samples and simulated samples. Make sure your local dataset release and TTree schema match the source code; do not assume all files share the same branches.

## Input metadata recommended for each real run

Record each input URL / source DOI, exact filename and SHA-256, release year, run period, integrated luminosity, cross sections, generator labels, tree schema, selection cuts, negative MC weight treatment, random seeds, ROOT/Python versions, and fit configuration. Keep raw `.root` files out of Git unless licenses and storage limits permit it.

## Units

The original integrated macro explicitly divides `photon_pt` and `photon_E` by 1000 to convert MeV to GeV. The new reader uses an explicit `--units MeV` switch instead of guessing from the first event.

## Reproducibility levels

1. **Synthetic demo:** runs today and produces repeatable labeled figures.
2. **User-supplied ROOT selection:** conditional on local ROOT files, matching branches and optional `uproot` dependencies; not independently tested here on real input.
3. **Historical TMVA reproduction:** pending exact ROOT datasets, XML provenance, preprocessing, weighting and a comparable ROOT build.
4. **Validated Higgs inference:** pending full statistical and systematic treatment.
