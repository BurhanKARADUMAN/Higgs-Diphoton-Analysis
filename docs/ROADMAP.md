# Research roadmap

## Phase A — documentation & code hygiene (this upgrade)
- [x] English portfolio-facing README
- [x] Independent Python educational analysis with explicit source metadata
- [x] Analytic kinematic checks and reproducible toy data
- [x] Basic Poisson likelihood demonstration and default blinding of real-data signal window
- [x] Optional BDT holdout diagnostics and mass-feature safeguards
- [x] Scientific audit of the historical ROOT/TMVA log and macros
- [x] CI tests and GitHub issue templates

## Phase B — reproduce real open-data selections
- [ ] Fetch and hash exact ATLAS Open Data release ROOT inputs
- [ ] Produce data and MC cut-flow tables (event weights and uncertainties)
- [ ] Capture signal production and background component normalizations
- [ ] Compare baseline histograms to documented official tutorial output
- [ ] Build predeclared train/validation/test partitions by source file/event
- [ ] Add calibrated BDT threshold and sideband shape tests

## Phase C — inference that can withstand review
- [ ] Mass-shape validation after BDT selection
- [ ] Nuisance parameter and alternative background shape model
- [ ] Spurious-signal and injection tests; blinded likelihood validation
- [ ] Local and global p-value study with toys / look-elsewhere effect
- [ ] Full environment, sample metadata and public reproduction instructions
