# Simulation 1 — scalar series

Scalar Beta-GAM HMM Monte Carlo study used in the revised manuscript.

Paper mode uses 200 repeated state-number selection replications for:

- baseline;
- short series;
- overlapping states.

BIC is the primary criterion and ICL is secondary. The conditional-on-correct-`K` estimation study uses 100 replications. Additional outputs include RMSE, transition error, Viterbi confusion matrices, bootstrap coverage, pseudo-residual diagnostics, basis sensitivity, a two-covariate experiment, smoothing sensitivity, GEM monotonicity, computational scaling, and an optional `hmmTMB` benchmark.

Run from the repository root with:

```r
source("RUN_ALL_CHECKS_QUICK.R")  # quick checks of all components
source("RUN_FINAL_SCALAR.R")       # full scalar study
```
