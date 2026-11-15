# Statistical inference and blinding

## Distinguish these quantities

- **ROC AUC:** binary ranking quality of a classifier on a specified validation/test sample.
- **S/√B:** a rough expected sensitivity diagnostic with assumptions; not a calibrated p-value, especially when S and B are fitted on the same data used to select cuts.
- **Peak center:** fit parameter for a signal component; not an independent mass measurement until calibration and uncertainties are established.
- **Peak width (Gaussian σ):** shape / detector-resolution proxy, not standard error on the peak center.
- **Profile-likelihood q₀:** a likelihood-ratio statistic for a fixed preselected hypothesis, with nuisances treated under the model. `sqrt(q0)` can have a local asymptotic interpretation under qualifying assumptions; **global significance is different**.

## Implementation in this project

The educational `fit_local_signal` uses a binned Poisson likelihood on 100–160 GeV in 1 GeV bins, a falling exponential background, and a Gaussian signal template at a fixed 125.09 GeV with σ=1.7 GeV. The signal yield is constrained to be nonnegative. This is a simple pedagogical analysis with **no background-shape alternatives, nuisance profiling, calibration or trials correction**. It is not a replacement for an ATLAS fit.

## Real-data blinding

- Default ROOT mode does **not** show, export or fit counts in 120–130 GeV.
- The analyst must consciously pass `--unblind` to use that region.
- Blinding is **not** a substitute for defining selections, model and optimizer settings before examining the signal region.

## Why tuning is dangerous

Selecting the BDT threshold, mass window, binning and background fit order by maximizing an observed excess can inflate the apparent significance. Test cut thresholds using separate simulation / control data, freeze the analysis plan, then perform a single independent inference on unseen signal-region data. Alternatively use a full statistically calibrated selection procedure with toy experiments.

## Further reading

Cowan, Cranmer, Gross and Vitells (2011), [Asymptotic formulae for likelihood-based tests of new physics](https://doi.org/10.1140/epjc/s10052-011-1554-0).
