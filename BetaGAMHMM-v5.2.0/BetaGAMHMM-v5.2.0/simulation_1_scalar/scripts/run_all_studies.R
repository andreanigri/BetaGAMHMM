# Final reproducibility workflow for scalar simulations.
# Avoid nested BLAS/OpenMP oversubscription. Parallelism is controlled at the
# Monte Carlo/grid level, which is faster and reproducible on multi-core machines.
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
           VECLIB_MAXIMUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
# Each expensive stage is checkpointed. Set FORCE_RERUN=true to ignore checkpoints.

cfg <- make_study_config()
output_root <- cfg$output_root
ensure_dir(output_root)
ensure_dir(file.path(output_root, "metadata"))
force_rerun <- tolower(Sys.getenv("FORCE_RERUN", unset = "false")) %in% c("1", "true", "yes")
stage_complete <- function(path) !force_rerun && file.exists(path)
mark_stage_complete <- function(path) writeLines(format(Sys.time(), tz = "UTC"), path)

save_rds_safe(cfg, file.path(output_root, "metadata", "study_config.rds"))
writeLines(capture.output(str(cfg, max.level = 5)), file.path(output_root, "metadata", "study_config.txt"))
write_csv_safe(system_information_table(), file.path(output_root, "metadata", "system_information.csv"))

# 1. Repeated K selection: BIC primary, ICL secondary.
main_dir <- file.path(output_root, "01_state_selection"); ensure_dir(main_dir); ensure_dir(file.path(main_dir, "tables"))
main_done <- file.path(main_dir, ".complete")
if (!stage_complete(main_done)) {
  message("[1/8] Repeated scalar K-selection study")
  cfg_main <- cfg; cfg_main$output_dir <- main_dir; cfg_main$monte_carlo <- cfg$selection_study
  selection_study <- run_revision_selection_study(cfg_main, resume = TRUE)
  write_csv_safe(selection_study$rule_results, file.path(main_dir, "tables", "selection_results_BIC_ICL.csv"))
  write_csv_safe(selection_study$performance, file.path(main_dir, "tables", "selection_performance_BIC_ICL.csv"))
  write_csv_safe(selection_study$fit_metrics, file.path(main_dir, "tables", "BIC_estimation_metrics.csv"))
  write_csv_safe(selection_study$grid_summary, file.path(main_dir, "tables", "grid_fit_summary.csv"))
  kdist <- as.data.frame(table(selection_study$rule_results$scenario, selection_study$rule_results$rule,
                               selection_study$rule_results$selected_K, useNA = "ifany"), stringsAsFactors = FALSE)
  names(kdist) <- c("scenario", "criterion", "selected_K", "frequency")
  write_csv_safe(kdist, file.path(main_dir, "tables", "selected_K_distribution.csv"))
  mono <- stats::aggregate(monotonicity_ok ~ scenario + K + lambda, selection_study$grid_summary,
                           function(x) mean(x %in% TRUE, na.rm = TRUE))
  names(mono)[names(mono)=="monotonicity_ok"] <- "monotonicity_pass_rate"
  write_csv_safe(mono, file.path(main_dir, "tables", "GEM_monotonicity_summary.csv"))
  mark_stage_complete(main_done)
} else message("[1/8] Scalar K-selection: already complete")

baseline_sim <- simulate_revision_scenario(cfg$scenarios$baseline, cfg$basis, cfg$seed + 500000L)

# 2. Representative residual diagnostics.
diag_dir <- file.path(output_root, "02_residual_diagnostics"); ensure_dir(diag_dir); ensure_dir(file.path(diag_dir, "figures"))
diag_done <- file.path(diag_dir, ".complete")
if (!stage_complete(diag_done)) {
  message("[2/8] Scalar pseudo-residual diagnostics")
  grid <- run_grid_search_hmm(
    baseline_sim$y, baseline_sim$basis_spec, cfg$fit$K_grid, cfg$fit$lambda_grid,
    n_start = cfg$fit$n_start, max_iter = cfg$fit$max_iter,
    optimizer_maxit = cfg$fit$optimizer_maxit, diagnostic = cfg$diagnostic,
    n_cores = cfg$parallel$n_cores, n_start_cores = 1L, progress = TRUE,
    seed = cfg$seed + 505000L
  )
  sel <- select_with_rule(grid, cfg$selection_rules$BIC)
  pr <- compute_pseudo_residuals(baseline_sim$y, sel$selected$fit, baseline_sim$basis_spec)
  write_csv_safe(pseudo_residual_diagnostics(pr), file.path(diag_dir, "pseudo_residual_diagnostics.csv"))
  pp <- plot_pseudo_residuals(pr)
  save_plot(pp$acf, file.path(diag_dir, "figures", "pseudo_residual_acf.png"), width = 7, height = 4.5)
  save_plot(pp$qq, file.path(diag_dir, "figures", "pseudo_residual_qq.png"), width = 6, height = 5)
  mark_stage_complete(diag_done)
} else message("[2/8] Residual diagnostics: already complete")

# 3. Basis dimension sensitivity.
basis_dir <- file.path(output_root, "03_basis_sensitivity"); ensure_dir(basis_dir); basis_done <- file.path(basis_dir, ".complete")
if (!stage_complete(basis_done)) {
  message("[3/8] Scalar basis-dimension sensitivity")
  ans <- run_basis_dimension_sensitivity(
    baseline_sim, n_inner_knots = cfg$basis_sensitivity$n_inner_knots,
    lambda_grid = cfg$fit$lambda_grid, n_start = cfg$basis_sensitivity$n_start,
    max_iter = cfg$fit$max_iter, n_cores = 1L, diagnostic = cfg$diagnostic,
    seed = cfg$seed + 510000L
  )
  write_csv_safe(ans, file.path(basis_dir, "basis_dimension_sensitivity.csv"))
  mark_stage_complete(basis_done)
} else message("[3/8] Basis sensitivity: already complete")

# 4. Two-covariate feasibility.
multicov_dir <- file.path(output_root, "04_multicovariate"); ensure_dir(multicov_dir); multicov_done <- file.path(multicov_dir, ".complete")
if (!stage_complete(multicov_done)) {
  message("[4/8] Two-covariate feasibility study")
  mc <- run_multicovariate_study(
    n_rep = cfg$multicovariate_followup$n_rep, T = cfg$multicovariate_followup$T,
    lambda_grid = cfg$fit$lambda_grid, n_start = cfg$multicovariate_followup$n_start,
    max_iter = cfg$fit$max_iter, n_cores = cfg$parallel$n_cores,
    diagnostic = cfg$diagnostic, seed = cfg$seed + 1500000L
  )
  write_csv_safe(mc, file.path(multicov_dir, "multicovariate_results.csv"))
  write_csv_safe(
    summarize_replicate_metrics(mc, metrics = c("phi_rmse", "A_rmse", "mean_curve_rmse", "viterbi_accuracy", "elapsed_seconds"), group_cols = "p"),
    file.path(multicov_dir, "multicovariate_summary.csv")
  )
  mark_stage_complete(multicov_done)
} else message("[4/8] Two-covariate study: already complete")

# 5. Common vs state-specific smoothing.
smooth_dir <- file.path(output_root, "05_smoothing_sensitivity"); ensure_dir(smooth_dir); smooth_done <- file.path(smooth_dir, ".complete")
if (!stage_complete(smooth_done)) {
  message("[5/8] Common versus state-specific smoothing")
  ss <- run_smoothing_penalty_followup(
    T = cfg$smoothing_followup$T, K = cfg$smoothing_followup$K,
    lambda_grid = cfg$fit$lambda_grid, n_start = cfg$smoothing_followup$n_start,
    n_cycles = cfg$smoothing_followup$n_cycles, seed = cfg$seed + 1600000L
  )
  write_csv_safe(ss$common_grid, file.path(smooth_dir, "smoothing_common_grid.csv"))
  write_csv_safe(ss$comparison, file.path(smooth_dir, "smoothing_penalty_comparison.csv"))
  write_csv_safe(ss$trace, file.path(smooth_dir, "smoothing_state_specific_trace.csv"))
  mark_stage_complete(smooth_done)
} else message("[5/8] Smoothing sensitivity: already complete")

# 6. Bootstrap coverage.
cov_dir <- file.path(output_root, "06_bootstrap_coverage"); ensure_dir(cov_dir); cov_done <- file.path(cov_dir, ".complete")
if (!stage_complete(cov_done)) {
  message("[6/8] Scalar bootstrap coverage")
  bc <- run_bootstrap_coverage_followup(
    n_rep = cfg$coverage_followup$n_rep, B_boot = cfg$coverage_followup$B_boot,
    T = cfg$coverage_followup$T, K = cfg$coverage_followup$K,
    lambda = cfg$coverage_followup$lambda,
    n_start_fit = cfg$coverage_followup$n_start_fit,
    n_start_boot = cfg$coverage_followup$n_start_boot,
    n_cores = cfg$parallel$n_cores, seed = cfg$seed + 1700000L
  )
  write_csv_safe(bc, file.path(cov_dir, "bootstrap_coverage.csv"))
  write_csv_safe(summarize_bootstrap_coverage_followup(bc), file.path(cov_dir, "bootstrap_coverage_summary.csv"))
  mark_stage_complete(cov_done)
} else message("[6/8] Bootstrap coverage: already complete")

# 7. External comparison with hmmTMB.
hmmtmb_dir <- file.path(output_root, "07_hmmtmb_benchmark"); ensure_dir(hmmtmb_dir); hmmtmb_done <- file.path(hmmtmb_dir, ".complete")
if (!stage_complete(hmmtmb_done)) {
  message("[7/8] External comparison with hmmTMB")
  cmp <- run_hmmtmb_comparison_study(
    n_rep = cfg$hmmtmb_benchmark$n_rep, T = cfg$hmmtmb_benchmark$T,
    K = cfg$hmmtmb_benchmark$K, basis_k = cfg$hmmtmb_benchmark$basis_k,
    lambda_grid = cfg$hmmtmb_benchmark$lambda_grid,
    n_start_ours = cfg$hmmtmb_benchmark$n_start_ours,
    n_cores = cfg$parallel$n_cores, seed = cfg$seed + 1800000L
  )
  cs <- summarize_hmmtmb_comparison(cmp)
  write_csv_safe(cmp, file.path(hmmtmb_dir, "hmmtmb_comparison.csv"))
  write_csv_safe(cs$method_summary, file.path(hmmtmb_dir, "hmmtmb_comparison_summary.csv"))
  write_csv_safe(cs$paired_summary, file.path(hmmtmb_dir, "hmmtmb_paired_summary.csv"))
  mark_stage_complete(hmmtmb_done)
} else message("[7/8] hmmTMB benchmark: already complete")

# 8. Computational benchmark + paper-compatible outputs.
bench_dir <- file.path(output_root, "08_paper"); ensure_dir(bench_dir); paper_done <- file.path(bench_dir, ".complete")
if (!stage_complete(paper_done)) {
  message("[8/8] Computational benchmark and paper-ready outputs")
  scaling <- run_computational_scaling_benchmark(
    n_start = cfg$computational_benchmark$n_start,
    max_iter = cfg$computational_benchmark$max_iter,
    seed = cfg$seed + 800000L
  )
  write_csv_safe(scaling, file.path(bench_dir, "computational_scaling.csv"))

  paper_root <- file.path(bench_dir, "paper")
  main_table_dir <- file.path(paper_root, "main_text", "tables")
  main_figure_dir <- file.path(paper_root, "main_text", "figures")
  supp_dir <- file.path(paper_root, "supplement")
  ensure_dir(main_table_dir); ensure_dir(main_figure_dir); ensure_dir(supp_dir)

  selection_perf_file <- file.path(output_root, "01_state_selection", "tables", "selection_performance_BIC_ICL.csv")
  if (file.exists(selection_perf_file)) {
    sel_perf <- utils::read.csv(selection_perf_file, stringsAsFactors = FALSE)
    write_csv_safe(sel_perf, file.path(main_table_dir, "TableS1_scalar_K_selection_recovery.csv"))
    write_simple_latex_table(
      sel_perf, file.path(main_table_dir, "TableS1_scalar_K_selection_recovery.tex"),
      caption = "Repeated recovery of the data-generating number of states in the scalar simulation study. BIC is primary; ICL is a secondary classification-oriented comparison.",
      label = "tab:scalar-k-recovery", digits = 3
    )
  }

  representative_grid <- run_grid_search_hmm(
    baseline_sim$y, baseline_sim$basis_spec,
    K_grid = cfg$fit$K_grid, lambda_grid = cfg$fit$lambda_grid,
    n_start = cfg$fit$n_start, max_iter = cfg$fit$max_iter,
    optimizer_maxit = cfg$fit$optimizer_maxit, diagnostic = cfg$diagnostic,
    n_cores = cfg$parallel$n_cores, n_start_cores = 1L, progress = TRUE,
    seed = cfg$seed + 2650000L
  )
  tab1 <- paper_candidate_k_table(representative_grid, cfg$selection_rules$BIC, cfg$fit$K_grid)
  write_csv_safe(tab1, file.path(main_table_dir, "Table1_model_comparison.csv"))
  write_simple_latex_table(tab1, file.path(main_table_dir, "Table1_model_comparison.tex"),
    caption = "Candidate-state comparison for a representative scalar data set. Lambda is selected by AIC with corrected effective degrees of freedom; K is selected by BIC. ICL is reported as a secondary classification-oriented criterion.",
    label = "tab:simulation-model-comparison")
  rep_sel <- select_with_rule(representative_grid, cfg$selection_rules$BIC)
  rep_fit <- rep_sel$selected$fit
  if (rep_fit$K == baseline_sim$K_true) {
    rep_fit <- align_fit_to_reference(rep_fit,
      list(beta=baseline_sim$beta_true, phi=baseline_sim$phi_true, pi=baseline_sim$pi_true, A=baseline_sim$A_true),
      basis_spec=baseline_sim$basis_spec)
  }
  write_csv_safe(data.frame(state=seq_len(rep_fit$K), phi_hat=rep_fit$phi),
                 file.path(main_table_dir, "representative_phi_estimates.csv"))
  A_df <- as.data.frame(as.table(rep_fit$A), stringsAsFactors=FALSE)
  names(A_df) <- c("from_state","to_state","probability")
  write_csv_safe(A_df, file.path(main_table_dir, "representative_transition_matrix.csv"))
  rep_v <- viterbi_decode(baseline_sim$y, baseline_sim$basis_spec, rep_fit)
  rep_summary <- data.frame(selected_K=rep_fit$K, selected_lambda=rep_fit$lambda[1L],
    loglik=rep_fit$loglik_obs, BIC_EDF=rep_fit$criteria$bic_edf, ICL_EDF=rep_fit$criteria$icl_edf,
    viterbi_accuracy=if(rep_fit$K==baseline_sim$K_true) mean(rep_v==baseline_sim$z) else NA_real_,
    stringsAsFactors=FALSE)
  write_csv_safe(rep_summary, file.path(main_table_dir, "representative_fit_summary.csv"))

  conditional_dir <- file.path(paper_root, "conditional")
  conditional <- run_conditional_estimation_study(cfg, conditional_dir, resume = TRUE)
  write_csv_safe(conditional$metrics, file.path(conditional_dir, "conditional_estimation_metrics.csv"))
  write_csv_safe(conditional$curves, file.path(conditional_dir, "conditional_estimation_curves.csv"))
  lambda_freq <- paper_lambda_frequency_table(conditional$metrics)
  write_csv_safe(lambda_freq, file.path(conditional_dir, "lambda_selection_frequencies.csv"))
  cm <- aggregate_confusion_matrix(conditional$confusion, cfg$scenarios$baseline$K_true)
  write_csv_safe(cm$counts, file.path(conditional_dir, "viterbi_confusion_matrix_counts.csv"))
  write_csv_safe(cm$row_normalized, file.path(conditional_dir, "viterbi_confusion_matrix_row_normalized.csv"))
  tab2 <- paper_monte_carlo_table(conditional$metrics)
  write_csv_safe(tab2, file.path(main_table_dir, "Table2_monte_carlo_performance.csv"))
  write_simple_latex_table(tab2, file.path(main_table_dir, "Table2_monte_carlo_performance.tex"),
    caption = sprintf("Monte Carlo performance over %d replications conditional on the data-generating K=%d.", cfg$paper_reporting$n_rep, cfg$scenarios$baseline$K_true),
    label = "tab:simulation-performance")
  make_paper_figure1(conditional$curves, baseline_sim,
    file.path(main_figure_dir, "Figure1_monte_carlo_curves.png"),
    file.path(main_figure_dir, "Figure1_monte_carlo_curves.pdf"),
    include_observed_points = cfg$paper_reporting$include_observed_points)
  export_supplementary_revision_tables(output_root, supp_dir)
  mark_stage_complete(paper_done)
} else message("[8/8] Paper outputs: already complete")

writeLines(capture.output(sessionInfo()), file.path(output_root, "metadata", "sessionInfo.txt"))
message("Final scalar simulation study complete: ", normalizePath(output_root, winslash = "/", mustWork = FALSE))
