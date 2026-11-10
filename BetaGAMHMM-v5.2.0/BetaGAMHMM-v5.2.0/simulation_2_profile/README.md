# Simulation 2 — repeated age profiles

Year-level profile Beta-GAM HMM Monte Carlo study used in the revised manuscript (the "vector-valued" simulation of the paper).

Paper mode uses 200 repeated state-number selection replications for:

- baseline, `T = 55`;
- baseline, `T = 250`;
- overlap + rare state, `T = 55`.

Each year contains a 41-component age profile and shares one latent annual state. BIC is primary and ICL secondary. The simulation records whether all generating states were realized and reports both unconditional recovery and recovery conditional on all states being visited.

The bootstrap coverage study evaluates percentile, basic, and bias-corrected intervals from the same refits.

Run from the repository root with:

```r
source("RUN_ALL_CHECKS_QUICK.R")  # quick checks of all components
source("RUN_FINAL_PROFILE.R")      # full vector-valued study
```
