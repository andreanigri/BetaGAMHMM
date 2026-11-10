# Final scalar simulation runner. Source from this directory.
options(stringsAsFactors = FALSE)
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
if (!dir.exists("R") || !file.exists("config/study_config.R")) stop("Run from simulation_1_scalar root.")
for (f in sort(list.files("R", pattern="\\.R$", full.names=TRUE))) source(f, local=.GlobalEnv)
source("config/study_config.R", local=.GlobalEnv)
source("scripts/run_all_studies.R", local=.GlobalEnv)
