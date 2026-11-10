# Final real-data selection/robustness pipeline for one country.
# Run from the real_data project root after setting BGHMM_COUNTRY and the period variables.

options(stringsAsFactors = FALSE)
# Avoid nested BLAS/OpenMP oversubscription. Parallelism is controlled at the
# Monte Carlo/grid level, which is faster and reproducible on multi-core machines.
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
           VECLIB_MAXIMUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
if (!file.exists(file.path("config", "real_data_config.R")) || !dir.exists("R")) {
  stop("RUN_COUNTRY.R must be sourced from the project root.")
}
PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path("config", "real_data_config.R"), local = .GlobalEnv)
for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f, local = .GlobalEnv)

safe_dir_create(OUTPUT_ROOT)
CHECKPOINT_DIR <- safe_dir_create(file.path(OUTPUT_ROOT, "checkpoints"))
FORCE_RERUN <- identical(tolower(Sys.getenv("FORCE_RERUN", "false")), "true")
load_or_run <- function(path, label, code) {
  if (file.exists(path) && !FORCE_RERUN) { message(label, ": loading checkpoint"); return(readRDS(path)) }
  message(label); value <- force(code); saveRDS(value, path); value
}

message("[0/10] Preflight checks and hardware diagnostics")
run_internal_checks()
hw <- choose_auto_workers(); N_CORES <- hw$selected_workers
write_hardware_diagnostics(hw, file.path(OUTPUT_ROOT, "hardware"))
if (N_CORES > 1L && !test_parallel_backend(N_CORES, project_root = PROJECT_ROOT, seed = SEED + 11L)) N_CORES <- 1L
message("      CPU: ", hw$cpu_name)
message("      Physical/logical cores: ", hw$physical_cores, "/", hw$logical_cores)
message("      Effective PSOCK workers: ", N_CORES)

message("============================================================")
message("Final unified HMD application: ", APPLICATION_LABEL, " [", COUNTRY_CODE, "]")
message("Years: ", YEAR_RANGE[1], "-", YEAR_RANGE[2], "; ages: ", AGE_RANGE[1], "-", AGE_RANGE[2])
message("State selection: lambda by min AIC_EDF within K; K by min BIC_EDF = -2 logL + nu log(T)")
message("K grid: ", paste(K_GRID, collapse = ","), "; lambda grid: ", paste(LAMBDA_GRID, collapse = ","))
message("============================================================")

message("[1/10] Reading and auditing HMD death-count data")
panel <- prepare_hmd_panel(DATA_FILE, COUNTRY_CODE, APPLICATION_LABEL, YEAR_RANGE, AGE_RANGE, BETA_BOUNDARY_ACTION)
write_data_diagnostics(panel, file.path(OUTPUT_ROOT, "data_diagnostics"))
message("      HMM time points: ", panel$n_years, "; ages/year: ", panel$n_ages, "; cells: ", panel$n_cells)
message("      Same-data max |Y - Female/Total|: ", format(panel$same_data_max_abs_error, scientific = TRUE))
basis_spec <- build_age_basis(panel$ages, degree = SPLINE_DEGREE, n_inner_knots = N_INNER_KNOTS)

# Model-level Windows/PSOCK smoke test. If it fails, continue sequentially.
if (N_CORES > 1L) {
  smoke_ok <- tryCatch({
    invisible(fit_beta_year_hmm(panel, basis_spec, K = 1L, lambda = 1, n_start = 2L,
      n_cores = min(2L, N_CORES), max_iter = 12L, tol = FIT_TOL,
      phi_lower = PHI_LOWER, phi_upper = PHI_UPPER, optimizer_maxit = min(50L, OPTIMIZER_MAXIT),
      precision_saturation_fraction = PRECISION_SATURATION_FRACTION, seed = SEED + 777L)); TRUE
  }, error = function(e) { warning("Model-level PSOCK smoke test failed: ", conditionMessage(e),
                                  ". Falling back to sequential execution."); FALSE })
  if (!smoke_ok) N_CORES <- 1L
}
message("      Effective workers after model smoke test: ", N_CORES)

message("[2/10] Scalar annual aggregate benchmark")
scalar_series <- prepare_scalar_aggregate(panel)
scalar_grid <- load_or_run(file.path(CHECKPOINT_DIR, "scalar_grid.rds"), "      Fitting scalar annual Beta-HMM grid",
  run_scalar_beta_grid(scalar_series, K_GRID, SCALAR_N_START, N_CORES, MAX_ITER, FIT_TOL,
                       PHI_LOWER, PHI_UPPER, OPTIMIZER_MAXIT, SEED + 1000L, TRUE))
scalar_selection <- select_scalar_beta_model(scalar_grid)
message("      Scalar selected K=", scalar_selection$fit$K, " by standard year-level BIC.")

message("[3/10] Repeated-profile Beta-GAM HMM grid")
profile_grid <- load_or_run(file.path(CHECKPOINT_DIR, "profile_beta_grid.rds"), "      Fitting profile grid",
  run_beta_grid(panel, basis_spec, K_GRID, LAMBDA_GRID, PROFILE_N_START, N_CORES,
                MAX_ITER, FIT_TOL, PHI_LOWER, PHI_UPPER, OPTIMIZER_MAXIT,
                PRECISION_SATURATION_FRACTION, SEED + 2000L, TRUE))
profile_selection <- select_beta_model(profile_grid)
if (profile_selection$fit$K == max(K_GRID) && max(K_GRID) < FINAL_K_AUTO_MAX && REAL_MODE == "paper") {
  message("      BIC optimum is at K-grid boundary; extending to K=", FINAL_K_AUTO_MAX, ".")
  add <- load_or_run(file.path(CHECKPOINT_DIR, "profile_beta_grid_extension.rds"), "      Fitting boundary extension",
    run_beta_grid(panel, basis_spec, seq.int(max(K_GRID)+1L, FINAL_K_AUTO_MAX), LAMBDA_GRID,
                  PROFILE_N_START, N_CORES, MAX_ITER, FIT_TOL, PHI_LOWER, PHI_UPPER,
                  OPTIMIZER_MAXIT, PRECISION_SATURATION_FRACTION, SEED + 2100L, TRUE))
  profile_grid <- append_grid_results(profile_grid, add)
  profile_selection <- select_beta_model(profile_grid)
}
message("      Profile selected K=", profile_selection$fit$K, ", lambda=", profile_selection$fit$lambda[1L], ".")
message("      Selected-fit precision-bound flag: ", profile_selection$fit$precision_diagnostic$flagged)

K_FINAL_GRID <- sort(unique(summarize_beta_grid(profile_grid)$K))
message("[4/10] Basis-dimension sensitivity (M_inner = 4,6,8)")
basis_sensitivity <- if (RUN_BASIS_SENSITIVITY) load_or_run(
  file.path(CHECKPOINT_DIR, "basis_sensitivity.rds"), "      Re-selecting K/lambda for each basis size",
  run_basis_selection_sensitivity(panel, profile_grid, K_FINAL_GRID, LAMBDA_GRID,
    BASIS_INNER_KNOTS_GRID, N_INNER_KNOTS, SPLINE_DEGREE, BASIS_SENSITIVITY_N_START,
    N_CORES, MAX_ITER, FIT_TOL, PHI_LOWER, PHI_UPPER, OPTIMIZER_MAXIT,
    PRECISION_SATURATION_FRACTION, SEED + 3000L)) else NULL

message("[5/10] Same-data denominator-aware sensitivity")
binomial_grid <- binomial_selection <- NULL
if (RUN_DENOMINATOR_SENSITIVITY) {
  binomial_grid <- load_or_run(file.path(CHECKPOINT_DIR, "denominator_grid.rds"), "      Fitting denominator-aware grid",
    run_binomial_grid(panel, basis_spec, K_GRID, LAMBDA_GRID, DENOMINATOR_N_START, N_CORES,
                      MAX_ITER, FIT_TOL, OPTIMIZER_MAXIT, SEED + 4000L, TRUE))
  binomial_selection <- select_binomial_model(binomial_grid)
  if (binomial_selection$fit$K == max(K_GRID) && max(K_GRID) < FINAL_K_AUTO_MAX && REAL_MODE == "paper") {
    addb <- load_or_run(file.path(CHECKPOINT_DIR, "denominator_grid_extension.rds"), "      Extending denominator-aware grid",
      run_binomial_grid(panel, basis_spec, seq.int(max(K_GRID)+1L, FINAL_K_AUTO_MAX), LAMBDA_GRID,
                        DENOMINATOR_N_START, N_CORES, MAX_ITER, FIT_TOL, OPTIMIZER_MAXIT,
                        SEED + 4100L, TRUE))
    binomial_grid <- append_grid_results(binomial_grid, addb)
    binomial_selection <- select_binomial_model(binomial_grid)
  }
  message("      Denominator-aware selected K=", binomial_selection$fit$K,
          ", lambda=", binomial_selection$fit$lambda[1L], ".")
}

message("[6/10] Targeted precision-bound sensitivity")
phi_bound_sensitivity <- if (RUN_PHI_BOUND_SENSITIVITY) load_or_run(
  file.path(CHECKPOINT_DIR, "phi_bound_sensitivity.rds"),
  "      Re-fitting one primary lambda per K under relaxed phi bounds",
  run_phi_bound_sensitivity(panel, basis_spec, profile_selection, PHI_BOUND_SENSITIVITY,
    PHI_BOUND_N_START, N_CORES, MAX_ITER, FIT_TOL, OPTIMIZER_MAXIT,
    PRECISION_SATURATION_FRACTION, SEED + 5000L, primary_grid = profile_grid)) else NULL
if (!is.null(phi_bound_sensitivity)) {
  message("      Phi-bound selected K: ",
          paste(phi_bound_sensitivity$summary$setting,
                phi_bound_sensitivity$summary$selected_K, sep="=", collapse="; "))
}

message("[7/10] Diagnostics and persistence")
beta_diagnostics <- compute_beta_pseudo_residuals(profile_selection$fit)
beta_persistence <- summarize_state_persistence(profile_selection$fit, panel$years)
binomial_diagnostics <- if (!is.null(binomial_selection)) compute_binomial_kernel_residuals(binomial_selection$fit) else NULL

message("[8/10] Writing final-selection outputs")
analysis_dir <- file.path(OUTPUT_ROOT, "analysis_outputs")
write_final_selection_outputs(panel, scalar_selection, profile_grid, profile_selection,
  basis_sensitivity, beta_diagnostics, beta_persistence,
  binomial_grid, binomial_selection, binomial_diagnostics, phi_bound_sensitivity, analysis_dir)

message("[9/10] Reproducibility snapshot")
config_snapshot <- data.frame(parameter = c("REAL_MODE","COUNTRY_CODE","APPLICATION_LABEL","START_YEAR","END_YEAR","MIN_AGE","MAX_AGE",
  "K_GRID_INITIAL","K_GRID_FINAL","LAMBDA_GRID","N_INNER_KNOTS","PROFILE_N_START","SCALAR_N_START",
  "BASIS_N_START","DENOMINATOR_N_START","PHI_BOUND_N_START","HMD_EXPECTED_LAST_MODIFIED","HMD_EXPECTED_DEATHS_MD5",
  "BOOTSTRAP_N_START","BOOTSTRAP_RETRY_FAILED","BOOTSTRAP_RETRY_N_START","N_CORES_EFFECTIVE","SEED"),
  value = c(REAL_MODE,COUNTRY_CODE,APPLICATION_LABEL,YEAR_RANGE[1],YEAR_RANGE[2],AGE_RANGE[1],AGE_RANGE[2],
    paste(K_GRID,collapse=";"),paste(K_FINAL_GRID,collapse=";"),paste(LAMBDA_GRID,collapse=";"),N_INNER_KNOTS,
    PROFILE_N_START,SCALAR_N_START,BASIS_SENSITIVITY_N_START,DENOMINATOR_N_START,PHI_BOUND_N_START,
    HMD_EXPECTED_LAST_MODIFIED,HMD_EXPECTED_DEATHS_MD5,BOOTSTRAP_N_START,BOOTSTRAP_RETRY_FAILED,BOOTSTRAP_RETRY_N_START,
    N_CORES,SEED), stringsAsFactors=FALSE)
utils::write.csv(config_snapshot, file.path(OUTPUT_ROOT, "run_configuration.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUTPUT_ROOT, "sessionInfo.txt"))
saveRDS(list(panel=panel,basis_spec=basis_spec,hardware=hw,scalar_selection=scalar_selection,
  profile_grid=profile_grid,profile_selection=profile_selection,basis_sensitivity=basis_sensitivity,
  beta_diagnostics=beta_diagnostics,beta_persistence=beta_persistence,
  phi_bound_sensitivity=phi_bound_sensitivity,
  binomial_grid=binomial_grid,binomial_selection=binomial_selection,
  binomial_diagnostics=binomial_diagnostics,configuration=config_snapshot),
  file.path(OUTPUT_ROOT,"analysis_object.rds"))

message("[10/10] Country analysis complete")
message("Outputs: ", normalizePath(OUTPUT_ROOT, winslash="/", mustWork=FALSE))
