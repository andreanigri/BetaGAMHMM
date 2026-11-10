# Avoid nested BLAS/OpenMP oversubscription. Parallelism is controlled at the
# Monte Carlo/grid level, which is faster and reproducible on multi-core machines.
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
           VECLIB_MAXIMUM_THREADS="1", NUMEXPR_NUM_THREADS="1")

cfg <- make_profile_config()
safe_dir_create(cfg$output_root)
safe_dir_create(file.path(cfg$output_root, "metadata"))
force_rerun <- tolower(Sys.getenv("FORCE_RERUN", unset = "false")) %in% c("1", "true", "yes")

stage_done <- function(path) !force_rerun && file.exists(path)
mark_done <- function(path) writeLines(format(Sys.time(), tz = "UTC"), path)

saveRDS(cfg, file.path(cfg$output_root, "metadata", "profile_config.rds"))
writeLines(capture.output(str(cfg, max.level = 5)), file.path(cfg$output_root, "metadata", "profile_config.txt"))
utils::write.csv(profile_scenario_design_summary(cfg$scenarios),
                 file.path(cfg$output_root, "metadata", "scenario_transition_design.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 1. Repeated state-number selection under the profile formulation
# -----------------------------------------------------------------------------
sel_dir <- file.path(cfg$output_root, "01_state_number_selection")
safe_dir_create(sel_dir); safe_dir_create(file.path(sel_dir, "chunks"))
sel_done <- file.path(sel_dir, ".complete")
selection_results <- NULL

if (!stage_done(sel_done)) {
  message("[1/5] Repeated-profile state-number selection")
  all_rows <- list()
  rr <- 0L
  for (sc_name in names(cfg$scenarios)) {
    sc <- cfg$scenarios[[sc_name]]
    message("  Scenario ", sc_name, ": ", sc$n_rep, " replicates")
    ids <- seq_len(sc$n_rep)
    chunks <- split(ids, ceiling(seq_along(ids) / cfg$selection$chunk_size))
    for (cc in seq_along(chunks)) {
      chunk_file <- file.path(sel_dir, "chunks", sprintf("%s_chunk_%03d.rds", sc_name, cc))
      if (!force_rerun && file.exists(chunk_file)) {
        chunk_rows <- readRDS(chunk_file)
      } else {
        jobs <- lapply(chunks[[cc]], function(r) list(
          scenario = sc, replicate = r,
          seed = cfg$seed + match(sc_name, names(cfg$scenarios)) * 1000000L + r,
          degree = cfg$basis$degree, n_inner_knots = cfg$basis$n_inner_knots,
          K_grid = cfg$fit$K_grid, lambda_grid = cfg$fit$lambda_grid,
          n_start = cfg$fit$n_start, max_iter = cfg$fit$max_iter,
          tol = cfg$fit$tol, phi_lower = cfg$fit$phi_lower, phi_upper = cfg$fit$phi_upper,
          optimizer_maxit = cfg$fit$optimizer_maxit,
          precision_saturation_fraction = cfg$fit$precision_saturation_fraction
        ))
        ans <- parallel_lapply(jobs, profile_selection_rep_worker,
          n_cores = cfg$n_cores, project_root = getwd(),
          seed = cfg$seed + 7000000L + cc)
        chunk_rows <- do.call(rbind, lapply(ans, `[[`, "results"))
        saveRDS(chunk_rows, chunk_file)
      }
      rr <- rr + 1L
      all_rows[[rr]] <- chunk_rows
      message(sprintf("    chunk %d/%d complete", cc, length(chunks)))
    }
  }
  selection_results <- do.call(rbind, all_rows)
  rownames(selection_results) <- NULL
  utils::write.csv(selection_results, file.path(sel_dir, "selection_replicates.csv"), row.names = FALSE)
  utils::write.csv(summarize_profile_selection(selection_results), file.path(sel_dir, "selection_recovery_summary.csv"), row.names = FALSE)
  utils::write.csv(summarize_profile_estimation(selection_results), file.path(sel_dir, "estimation_summary_correct_K.csv"), row.names = FALSE)
  utils::write.csv(profile_selected_k_distribution(selection_results), file.path(sel_dir, "selected_K_distribution.csv"), row.names = FALSE)
  mark_done(sel_done)
} else {
  message("[1/5] Repeated-profile state-number selection: already complete")
  selection_results <- utils::read.csv(file.path(sel_dir, "selection_replicates.csv"), stringsAsFactors = FALSE)
}

# -----------------------------------------------------------------------------
# 2. Conditional bootstrap coverage in the application-scale setting
# -----------------------------------------------------------------------------
cov_dir <- file.path(cfg$output_root, "02_bootstrap_coverage")
safe_dir_create(cov_dir)
cov_done <- file.path(cov_dir, ".complete")
coverage_results <- NULL
if (!stage_done(cov_done)) {
  message("[2/5] Conditional bootstrap coverage (T=55, J=41, K=4)")
  sc <- cfg$scenarios$baseline_T55
  jobs <- lapply(seq_len(cfg$coverage$n_outer), function(r) list(
    scenario = sc, replicate = r, seed = cfg$seed + 20000000L + r,
    degree = cfg$basis$degree, n_inner_knots = cfg$basis$n_inner_knots,
    lambda_grid = cfg$fit$lambda_grid,
    n_start_fit = cfg$coverage$n_start_fit, n_start_boot = cfg$coverage$n_start_boot,
    B_boot = cfg$coverage$B_boot, level = cfg$coverage$level,
    max_iter = cfg$fit$max_iter, tol = cfg$fit$tol,
    phi_lower = cfg$fit$phi_lower, phi_upper = cfg$fit$phi_upper,
    optimizer_maxit = cfg$fit$optimizer_maxit
  ))
  res <- parallel_lapply(jobs, profile_coverage_worker,
    n_cores = cfg$n_cores, project_root = getwd(), seed = cfg$seed + 21000000L)
  coverage_results <- do.call(rbind, Filter(Negate(is.null), res))
  if (is.null(coverage_results)) coverage_results <- data.frame()
  utils::write.csv(coverage_results, file.path(cov_dir, "coverage_replicates.csv"), row.names = FALSE)
  mark_done(cov_done)
} else {
  message("[2/5] Conditional bootstrap coverage: already complete")
  coverage_results <- utils::read.csv(file.path(cov_dir, "coverage_replicates.csv"), stringsAsFactors = FALSE)
}

# -----------------------------------------------------------------------------
# 3. Basis-dimension sensitivity with full K/lambda reselection
# -----------------------------------------------------------------------------
basis_dir <- file.path(cfg$output_root, "03_basis_sensitivity")
safe_dir_create(basis_dir)
basis_done <- file.path(basis_dir, ".complete")
basis_results <- NULL
if (!stage_done(basis_done)) {
  message("[3/5] Profile basis-dimension sensitivity")
  sc <- cfg$scenarios$baseline_T55
  jobs <- list(); jj <- 0L
  for (r in seq_len(cfg$basis_sensitivity$n_rep)) for (m in cfg$basis_sensitivity$n_inner_knots) {
    jj <- jj + 1L
    jobs[[jj]] <- list(
      scenario = sc, replicate = r, n_inner_knots = m,
      seed = cfg$seed + 30000000L + r,
      degree = cfg$basis$degree, K_grid = cfg$fit$K_grid, lambda_grid = cfg$fit$lambda_grid,
      n_start = cfg$basis_sensitivity$n_start,
      max_iter = cfg$fit$max_iter, tol = cfg$fit$tol,
      phi_lower = cfg$fit$phi_lower, phi_upper = cfg$fit$phi_upper,
      optimizer_maxit = cfg$fit$optimizer_maxit
    )
  }
  ans <- parallel_lapply(jobs, profile_basis_worker,
    n_cores = cfg$n_cores, project_root = getwd(), seed = cfg$seed + 31000000L)
  basis_results <- do.call(rbind, Filter(Negate(is.null), ans))
  utils::write.csv(basis_results, file.path(basis_dir, "basis_sensitivity_replicates.csv"), row.names = FALSE)
  mark_done(basis_done)
} else {
  message("[3/5] Profile basis-dimension sensitivity: already complete")
  basis_results <- utils::read.csv(file.path(basis_dir, "basis_sensitivity_replicates.csv"), stringsAsFactors = FALSE)
}

# -----------------------------------------------------------------------------
# 4. Representative fit, diagnostics, and computational scaling
# -----------------------------------------------------------------------------
rep_dir <- file.path(cfg$output_root, "04_diagnostics_and_computation")
safe_dir_create(rep_dir)
rep_done <- file.path(rep_dir, ".complete")
example <- NULL
benchmark <- NULL
if (!stage_done(rep_done)) {
  message("[4/5] Representative profile fit, diagnostics, and timing")
  sim <- simulate_profile_scenario(cfg$scenarios$baseline_T55, cfg$seed + 40000000L)
  basis <- build_age_basis(sim$panel$ages, cfg$basis$degree, cfg$basis$n_inner_knots)
  grid <- run_profile_grid(
    sim$panel, basis, K_grid = sim$scenario$K, lambda_grid = cfg$fit$lambda_grid,
    n_start = max(4L, cfg$fit$n_start), n_cores = cfg$n_cores,
    max_iter = cfg$fit$max_iter, tol = cfg$fit$tol,
    phi_lower = cfg$fit$phi_lower, phi_upper = cfg$fit$phi_upper,
    optimizer_maxit = cfg$fit$optimizer_maxit,
    seed = cfg$seed + 40100000L, progress = TRUE
  )
  tab <- summarize_profile_grid(grid)
  ok <- tab[is.na(tab$error) & tab$converged & tab$monotonicity_ok, , drop = FALSE]
  best <- ok[order(ok$aic_edf), , drop = FALSE][1L, ]
  idx <- which(vapply(grid, function(f) !inherits(f, "profile_grid_error") && abs(f$lambda[1L] - best$lambda) < 1e-12, logical(1)))
  fit <- align_profile_fit_to_truth(grid[[idx[1L]]], sim)
  diagnostics <- compute_beta_pseudo_residuals(fit)
  utils::write.csv(diagnostics$summary, file.path(rep_dir, "pseudo_residual_diagnostics.csv"), row.names = FALSE)
  saveRDS(list(sim = sim, fit = fit), file.path(rep_dir, "representative_fit.rds"))
  benchmark <- run_profile_computational_benchmark(cfg)
  utils::write.csv(benchmark, file.path(rep_dir, "computational_benchmark.csv"), row.names = FALSE)
  example <- list(sim = sim, fit = fit, diagnostics = diagnostics)
  mark_done(rep_done)
} else {
  message("[4/5] Representative profile fit and timing: already complete")
  tmp <- readRDS(file.path(rep_dir, "representative_fit.rds"))
  example <- list(sim = tmp$sim, fit = tmp$fit, diagnostics = compute_beta_pseudo_residuals(tmp$fit))
  benchmark <- utils::read.csv(file.path(rep_dir, "computational_benchmark.csv"), stringsAsFactors = FALSE)
}

# -----------------------------------------------------------------------------
# 5. Paper-ready outputs
# -----------------------------------------------------------------------------
out_done <- file.path(cfg$output_root, "05_paper", ".complete")
if (!stage_done(out_done)) {
  message("[5/5] Paper-ready profile tables and figures")
  write_profile_paper_outputs(cfg, selection_results, coverage_results, basis_results, benchmark, example)
  safe_dir_create(dirname(out_done)); mark_done(out_done)
} else message("[5/5] Paper-ready outputs: already complete")

writeLines(capture.output(sessionInfo()), file.path(cfg$output_root, "metadata", "sessionInfo.txt"))
message("Profile simulation study complete: ", normalizePath(cfg$output_root, winslash = "/", mustWork = FALSE))
