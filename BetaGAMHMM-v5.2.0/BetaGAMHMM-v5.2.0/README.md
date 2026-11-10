# Beta-GAM HMM — reproducibility code

Reproducibility code for Nigri, Shang and Bonetti, *A Beta-GAM Hidden Markov Model for Scalar and Vector-Valued Proportion Time Series* (submitted to *Statistics and Computing*). Repository: <https://github.com/andreuccio88/BetaGAMHMM>.

## Quick start

1. Place the HMD archive in `real_data/data/HMD_Deaths.zip` (see "HMD data" below).
2. Open `BetaGAMHMM.Rproj` and run:

   ```r
   source("RUN_ALL_FINAL.R")   # all simulations, real-data analyses and published profile figures
   ```

All other scripts are optional checks or run individual parts of the workflow.

## Contents

The repository contains three independent components:

- `simulation_1_scalar/`: scalar Monte Carlo study (Section 6.1–6.2 of the paper);
- `simulation_2_profile/`: repeated age-profile Monte Carlo study — the **vector-valued** simulation of the paper (Section 6.3);
- `real_data/`: Human Mortality Database applications for Russia and Italy (Section 7).

`redraw_profile_figures.R` redraws the published state-profile figures (manuscript Figures 3 and 5, Supplement Figure S11) from the saved final-inference objects, without refitting.

The public release contains **code only**. Generated outputs, checkpoints, caches, and the HMD data archive are intentionally excluded from version control.

## Reference environment

The manuscript runs were produced with R 4.4.2 on Windows 11. Reference package versions are recorded in `DEPENDENCIES.md` and `environment/REFERENCE_SESSIONINFO.txt`.

For the exact package versions used in the manuscript, run:

```r
source("INSTALL_PACKAGES.R")
```

This installer is optional if the required versions are already available.

## HMD data — required for the real-data analyses

HMD data are not redistributed in this repository. For exact manuscript reproduction:

1. Go to <https://www.mortality.org/Data/ZippedDataFiles>.
2. Download the **Death Counts** archive from the HMD snapshot dated **27 Aug 2026** (the file metadata must report `Last modified: 27 Aug 2026`).
3. Rename the downloaded archive to `HMD_Deaths.zip`.
4. Place it in:

```text
real_data/data/HMD_Deaths.zip
```

Alternatively, keep the archive anywhere and set:

```r
Sys.setenv(BGHMM_DATA_FILE = "C:/path/to/HMD_Deaths.zip")
```

The code validates the snapshot date before fitting. A deliberately updated-data analysis can override the check with `BGHMM_ALLOW_HMD_VERSION_MISMATCH=true`, but that is **not** an exact reproduction of the manuscript.

See `docs/DATA.md` for details and HMD citation guidance.

## First check

Open `BetaGAMHMM.Rproj` and run:

```r
source("CHECK_CODE.R")
source("RUN_TESTS.R")
```

For quick end-to-end smoke runs of all three components:

```r
source("RUN_ALL_CHECKS_QUICK.R")
```

## Full reproduction

The full analyses can be run separately:

```r
source("RUN_FINAL_SCALAR.R")
source("RUN_FINAL_PROFILE.R")
source("RUN_FINAL_REAL_DATA.R")
source("RUN_FINAL_INFERENCE_REAL_DATA.R")
```

or sequentially with one command:

```r
source("RUN_ALL_FINAL.R")
```

`RUN_ALL_FINAL.R` performs the two simulation studies, the three real-data selections/robustness analyses, the final real-data bootstrap inference, and finally runs `redraw_profile_figures.R`, which writes the published state-profile figures to `redrawn_figures/`. The figure step can also be run on its own once the final inference outputs exist:

```r
source("redraw_profile_figures.R")
```

### Final real-data cases

- Russia, 1960–2014, ages 0–40;
- Italy, 1960–2014, ages 0–40;
- Italy, 1950–2019, ages 0–40 (extended-period sensitivity analysis).

## Bootstrap implementation

The conditional real-data bootstrap uses `B = 200` as pre-specified. Each replicate initially uses 6 random starts. Only failed refits are retried **once**, with 8 starts; successful first-pass refits are never recomputed. The output reports first-pass success, number retried, number recovered, and final success rate.

The manuscript interval convention is pre-specified:

- precision parameters: bias-corrected bootstrap intervals;
- transition probabilities: percentile intervals;
- initial probabilities: estimated and counted in the criteria, but not tabulated (a single sequence carries very little information about them);
- smooth age-profile curves: pointwise percentile intervals.

All interval formulas are retained from the same bootstrap refits.

## Parallel execution

Parallelism is across Monte Carlo replicates/grid jobs. BLAS/OpenMP threads are forced to one per worker to avoid nested oversubscription. Automatic worker selection is capped at 12; this is appropriate for the i9-class machine used for the manuscript runs. Override only when necessary:

```r
Sys.setenv(BGHMM_N_CORES = "10")
```

On Windows, extracting/cloning to a short path such as `C:/BGHMM/BetaGAMHMM` is recommended.

## Outputs and checkpoints

All generated material is written under each subproject's `outputs/` directory and is ignored by Git. Runs are checkpointed and can resume after interruption. To deliberately delete all generated outputs before a clean rerun:

```r
Sys.setenv(BGHMM_CONFIRM_RESET = "YES")
source("RESET_OUTPUTS.R")
```

## Reproducibility safeguards

- frozen deterministic seeds;
- deterministic candidate-level multi-start seeds;
- BIC as the primary state-number criterion and ICL as secondary;
- no occupancy-based state-number filtering;
- explicit support for `K = 1`;
- Beta-specific EDF in information criteria;
- GEM monotonicity checks/backtracking;
- profile HMM time unit is **year**, with ages as vector-emission components;
- basis, precision-bound, and denominator-aware sensitivity analyses;
- selection-stability bootstrap separated from conditional parameter uncertainty;
- HMD snapshot validation;
- session information and run configuration written by official runs.

## License

Code is released under the MIT License. HMD data are **not** covered by this software license and are not included in this repository.
