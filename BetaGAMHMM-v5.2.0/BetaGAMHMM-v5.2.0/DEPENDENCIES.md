# Reference software environment

Manuscript environment:

- R 4.4.2
- Windows 11
- ggplot2 4.0.2
- pbapply 1.7-4
- hmmTMB 1.1.2
- TMB 1.9.17
- tmbstan 1.0.91
- mgcv 1.9-3

Core model code otherwise uses base/recommended R packages including `stats`, `splines`, `parallel`, `utils`, `graphics`, and `grDevices`.

`hmmTMB` is used for the external benchmark in Simulation 1; if it is unavailable, that benchmark is explicitly recorded as skipped rather than affecting the core analyses.

`INSTALL_PACKAGES.R` can install/check the reference versions. Every official analysis writes `sessionInfo.txt` in its output directory.
