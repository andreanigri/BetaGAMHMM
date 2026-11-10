# Reproducibility workflow

## 1. Software

Use R 4.4.2 for the closest reproduction of the manuscript environment. Package versions are recorded in `DEPENDENCIES.md`.

Run:

```r
source("CHECK_CODE.R")
source("RUN_TESTS.R")
```

before expensive analyses.

## 2. Data

Place the 27 Aug 2026 HMD Death Counts archive at `real_data/data/HMD_Deaths.zip`. See `docs/DATA.md`.

## 3. Quick smoke test

```r
source("RUN_ALL_CHECKS_QUICK.R")
```

This is only a runtime/reproducibility check; quick Monte Carlo results are not manuscript results.

## 4. Full manuscript analyses

Run one stage at a time:

```r
source("RUN_FINAL_SCALAR.R")
source("RUN_FINAL_PROFILE.R")
source("RUN_FINAL_REAL_DATA.R")
source("RUN_FINAL_INFERENCE_REAL_DATA.R")
```

or:

```r
source("RUN_ALL_FINAL.R")
```

The expensive final real-data inference is intentionally separated from real-data model selection, but `RUN_ALL_FINAL.R` executes both in the correct order and then redraws the published state-profile figures (`redraw_profile_figures.R`, output in `redrawn_figures/`).

## 5. Resume behavior

Each project uses checkpoints. Rerunning the same entry point loads completed stages unless `FORCE_RERUN=true` is set.

A fresh Git clone has empty output folders, so the first official run starts from scratch.

## 6. Clean reset

To remove generated outputs from all three projects:

```r
Sys.setenv(BGHMM_CONFIRM_RESET = "YES")
source("RESET_OUTPUTS.R")
```

The confirmation variable is required to prevent accidental deletion.
