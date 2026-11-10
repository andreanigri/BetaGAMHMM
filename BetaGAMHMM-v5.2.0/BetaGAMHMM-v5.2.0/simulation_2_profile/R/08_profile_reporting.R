# Tables and figures for the final profile simulation study.

profile_selected_k_distribution <- function(results) {
  d <- results[is.finite(results$selected_K), , drop = FALSE]
  tab <- as.data.frame(table(d$scenario, d$criterion, d$selected_K), stringsAsFactors = FALSE)
  names(tab) <- c("scenario", "criterion", "selected_K", "frequency")
  totals <- ave(tab$frequency, tab$scenario, tab$criterion, FUN = sum)
  tab$proportion <- ifelse(totals > 0, tab$frequency / totals, NA_real_)
  tab
}

save_profile_selection_figures <- function(recovery, kdist, out_dir) {
  safe_dir_create(out_dir)
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))

  rec <- recovery
  rec$exact_percent <- 100 * rec$exact_rate
  p1 <- ggplot2::ggplot(rec, ggplot2::aes(x = criterion, y = exact_percent, fill = criterion)) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::facet_wrap(~ scenario, ncol = 2) +
    ggplot2::coord_cartesian(ylim = c(0, 100)) +
    ggplot2::labs(x = NULL, y = "Correct K selection (%)") +
    ggplot2::theme_minimal(base_size = 11)
  ggplot2::ggsave(file.path(out_dir, "FigureP1_K_recovery_BIC_ICL.png"), p1, width = 9, height = 6.5, dpi = 180)
  ggplot2::ggsave(file.path(out_dir, "FigureP1_K_recovery_BIC_ICL.pdf"), p1, width = 9, height = 6.5)

  kd <- kdist
  kd$selected_K <- factor(kd$selected_K)
  p2 <- ggplot2::ggplot(kd, ggplot2::aes(x = selected_K, y = proportion, fill = criterion)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::facet_wrap(~ scenario, ncol = 2) +
    ggplot2::labs(x = "Selected K", y = "Selection proportion", fill = "Criterion") +
    ggplot2::theme_minimal(base_size = 11)
  ggplot2::ggsave(file.path(out_dir, "FigureP2_selected_K_distribution.png"), p2, width = 9, height = 6.5, dpi = 180)
  ggplot2::ggsave(file.path(out_dir, "FigureP2_selected_K_distribution.pdf"), p2, width = 9, height = 6.5)

  scaling <- rec[rec$scenario %in% c("baseline_T55", "baseline_T250"), , drop = FALSE]
  scaling$T <- as.integer(sub("baseline_T", "", scaling$scenario))
  p3 <- ggplot2::ggplot(scaling, ggplot2::aes(x = T, y = exact_percent, linetype = criterion, shape = criterion)) +
    ggplot2::geom_line() + ggplot2::geom_point(size = 2.5) +
    ggplot2::coord_cartesian(ylim = c(0, 100)) +
    ggplot2::labs(x = "Number of HMM time points (T)", y = "Correct K selection (%)", linetype = "Criterion", shape = "Criterion") +
    ggplot2::theme_minimal(base_size = 11)
  ggplot2::ggsave(file.path(out_dir, "FigureP3_recovery_vs_T.png"), p3, width = 7.5, height = 5, dpi = 180)
  ggplot2::ggsave(file.path(out_dir, "FigureP3_recovery_vs_T.pdf"), p3, width = 7.5, height = 5)
  invisible(list(recovery = p1, distribution = p2, scaling = p3))
}

save_profile_example_curves <- function(sim, fit, path_png, path_pdf) {
  fit <- align_profile_fit_to_truth(fit, sim)
  K <- sim$scenario$K
  ages <- sim$panel$ages
  draw <- function() {
    op <- graphics::par(mfrow = c(ceiling(K / 2), 2), mar = c(4, 4, 2, 1))
    on.exit(graphics::par(op), add = TRUE)
    for (k in seq_len(K)) {
      idx <- which(fit$viterbi == k)
      yy <- c(sim$mu_true[, k], fit$mu_age[, k])
      if (length(idx)) yy <- c(yy, sim$panel$Y[idx, , drop=FALSE])
      graphics::plot(ages, sim$mu_true[, k], type = "l", lty = 2, lwd = 2,
        ylim = range(yy, finite=TRUE),
        xlab = "Profile coordinate (age)", ylab = expression(mu(x)), main = paste("State", k))
      if (length(idx)) {
        for (t in idx) graphics::points(ages, sim$panel$Y[t, ], pch=16, cex=0.22,
                                        col=grDevices::adjustcolor("grey45", alpha.f=0.15))
      }
      graphics::lines(ages, fit$mu_age[, k], lwd = 2)
      graphics::lines(ages, sim$mu_true[, k], lty=2, lwd=2)
      graphics::legend("bottomright", legend = c("Truth", "Estimate", "Observed (decoded state)"),
                       lty = c(2, 1, NA), pch=c(NA,NA,16), bty = "n", cex = 0.72)
    }
  }
  grDevices::png(path_png, width = 1500, height = 1100, res = 180); draw(); grDevices::dev.off()
  grDevices::pdf(path_pdf, width = 8.5, height = 6.2); draw(); grDevices::dev.off()
}

profile_representative_outputs <- function(sim, fit) {
  fit <- align_profile_fit_to_truth(fit, sim)
  cm <- as.data.frame(table(
    true_state=factor(sim$z, levels=seq_len(sim$scenario$K)),
    decoded_state=factor(fit$viterbi, levels=seq_len(sim$scenario$K))), stringsAsFactors=FALSE)
  names(cm) <- c("true_state","decoded_state","count")
  cm$true_state <- as.integer(as.character(cm$true_state))
  cm$decoded_state <- as.integer(as.character(cm$decoded_state))
  denom <- ave(cm$count, cm$true_state, FUN=sum)
  cmn <- cm; cmn$proportion <- ifelse(denom>0, cmn$count/denom, NA_real_)
  A <- as.data.frame(as.table(fit$A), stringsAsFactors=FALSE)
  names(A) <- c("from_state","to_state","probability")
  list(
    summary=data.frame(K=fit$K, lambda=fit$lambda[1L], loglik=fit$loglik_obs,
      BIC_EDF=fit$criteria$bic_edf, ICL_EDF=fit$criteria$icl_edf,
      viterbi_accuracy=mean(fit$viterbi==sim$z), ari=adjusted_rand_index(fit$viterbi,sim$z),
      curve_rmse=sqrt(mean((fit$mu_age-sim$mu_true)^2)), stringsAsFactors=FALSE),
    phi=data.frame(state=seq_len(fit$K), phi_hat=fit$phi, phi_true=sim$phi_true),
    transition=A, confusion_counts=cm, confusion_row_normalized=cmn
  )
}


write_profile_paper_outputs <- function(cfg, selection_results, coverage, basis_results, benchmark, example) {
  root <- file.path(cfg$output_root, "05_paper")
  main_tables <- file.path(root, "main_text", "tables")
  main_figs <- file.path(root, "main_text", "figures")
  supp_tables <- file.path(root, "supplement", "tables")
  supp_figs <- file.path(root, "supplement", "figures")
  lapply(c(main_tables, main_figs, supp_tables, supp_figs), safe_dir_create)

  recovery <- summarize_profile_selection(selection_results)
  estimation <- summarize_profile_estimation(selection_results)
  kdist <- profile_selected_k_distribution(selection_results)

  utils::write.csv(recovery, file.path(main_tables, "TableP1_state_number_recovery_BIC_ICL.csv"), row.names = FALSE)
  write_latex_table(recovery, file.path(main_tables, "TableP1_state_number_recovery_BIC_ICL.tex"),
    caption = "Recovery of the data-generating number of states. BIC is the proposed state-number criterion; ICL is reported as a classification-oriented secondary comparison.",
    label = "tab:profile-k-recovery", digits = 3)

  utils::write.csv(estimation, file.path(main_tables, "TableP2_profile_estimation_performance.csv"), row.names = FALSE)
  write_latex_table(estimation, file.path(main_tables, "TableP2_profile_estimation_performance.tex"),
    caption = "Profile-model estimation performance conditional on selecting the correct number of states.",
    label = "tab:profile-estimation", digits = 3)

  if (!is.null(coverage) && nrow(coverage)) {
    methods <- c("percentile", "basic", "bc")
    cov_sum <- do.call(rbind, lapply(methods, function(m) data.frame(
      method = m,
      n_outer = nrow(coverage),
      bootstrap_success_rate = mean(coverage$bootstrap_success_rate, na.rm = TRUE),
      phi_coverage = mean(coverage[[paste0("phi_coverage_", m)]], na.rm = TRUE),
      A_coverage_all = mean(coverage[[paste0("A_coverage_", m, "_all")]], na.rm = TRUE),
      A_coverage_diagonal = mean(coverage[[paste0("A_coverage_", m, "_diagonal")]], na.rm = TRUE),
      curve_coverage = mean(coverage[[paste0("curve_coverage_", m)]], na.rm = TRUE),
      mean_phi_width = mean(coverage[[paste0("mean_phi_width_", m)]], na.rm = TRUE),
      mean_A_width = mean(coverage[[paste0("mean_A_width_", m)]], na.rm = TRUE),
      mean_curve_width = mean(coverage[[paste0("mean_curve_width_", m)]], na.rm = TRUE),
      stringsAsFactors = FALSE
    )))
    utils::write.csv(cov_sum, file.path(main_tables, "TableP3_profile_bootstrap_coverage.csv"), row.names = FALSE)
    write_latex_table(cov_sum, file.path(main_tables, "TableP3_profile_bootstrap_coverage.tex"),
      caption = "Empirical coverage of the conditional parametric bootstrap in the application-scale profile setting.",
      label = "tab:profile-coverage", digits = 3)
  }

  utils::write.csv(kdist, file.path(supp_tables, "selected_K_distribution.csv"), row.names = FALSE)
  utils::write.csv(selection_results, file.path(supp_tables, "profile_selection_replicates.csv"), row.names = FALSE)
  if (!is.null(coverage)) utils::write.csv(coverage, file.path(supp_tables, "profile_coverage_replicates.csv"), row.names = FALSE)
  if (!is.null(basis_results) && nrow(basis_results)) {
    utils::write.csv(basis_results, file.path(supp_tables, "profile_basis_sensitivity.csv"), row.names = FALSE)
    bg <- split(basis_results, interaction(basis_results$n_inner_knots, basis_results$criterion, drop = TRUE))
    basis_summary <- do.call(rbind, lapply(bg, function(d) data.frame(
      n_inner_knots = d$n_inner_knots[1L], criterion = d$criterion[1L],
      n = nrow(d), exact_rate = mean(d$exact, na.rm = TRUE),
      mean_selected_K = mean(d$selected_K, na.rm = TRUE),
      mean_viterbi_accuracy = mean(d$viterbi_accuracy, na.rm = TRUE),
      mean_curve_rmse = mean(d$curve_rmse, na.rm = TRUE),
      precision_flag_rate = mean(d$precision_flag %in% TRUE, na.rm = TRUE),
      stringsAsFactors = FALSE
    )))
    rownames(basis_summary) <- NULL
    utils::write.csv(basis_summary, file.path(supp_tables, "profile_basis_sensitivity_summary.csv"), row.names = FALSE)
    write_latex_table(basis_summary, file.path(supp_tables, "TableSP1_basis_sensitivity.tex"),
      caption = "Sensitivity of profile state-number selection to the spline basis dimension.",
      label = "tab:profile-basis", digits = 3)
  }
  if (!is.null(benchmark) && nrow(benchmark)) {
    utils::write.csv(benchmark, file.path(supp_tables, "profile_computational_benchmark.csv"), row.names = FALSE)
    bm <- split(benchmark, interaction(benchmark$T, benchmark$J, drop = TRUE))
    benchmark_summary <- do.call(rbind, lapply(bm, function(d) data.frame(
      T = d$T[1L], J = d$J[1L], n_cells = d$n_cells[1L],
      mean_seconds = mean(d$elapsed_seconds), sd_seconds = stats::sd(d$elapsed_seconds),
      convergence_rate = mean(d$converged), stringsAsFactors = FALSE
    )))
    rownames(benchmark_summary) <- NULL
    utils::write.csv(benchmark_summary, file.path(supp_tables, "profile_computational_benchmark_summary.csv"), row.names = FALSE)
    write_latex_table(benchmark_summary, file.path(supp_tables, "TableSP2_computational_scaling.tex"),
      caption = "Computational scaling of the profile Beta-GAM HMM.",
      label = "tab:profile-computation", digits = 3)
  }

  save_profile_selection_figures(recovery, kdist, main_figs)
  if (!is.null(example)) {
    save_profile_example_curves(example$sim, example$fit,
      file.path(main_figs, "FigureP4_true_vs_estimated_profiles.png"),
      file.path(main_figs, "FigureP4_true_vs_estimated_profiles.pdf"))
    rex <- profile_representative_outputs(example$sim, example$fit)
    utils::write.csv(rex$summary, file.path(main_tables, "representative_fit_summary.csv"), row.names=FALSE)
    utils::write.csv(rex$phi, file.path(main_tables, "representative_phi_estimates.csv"), row.names=FALSE)
    utils::write.csv(rex$transition, file.path(main_tables, "representative_transition_matrix.csv"), row.names=FALSE)
    utils::write.csv(rex$confusion_counts, file.path(supp_tables, "viterbi_confusion_matrix_counts.csv"), row.names=FALSE)
    utils::write.csv(rex$confusion_row_normalized, file.path(supp_tables, "viterbi_confusion_matrix_row_normalized.csv"), row.names=FALSE)
    if (!is.null(example$diagnostics)) save_beta_diagnostic_plots(example$diagnostics, example$sim$panel$years, supp_figs)
  }
  invisible(root)
}
