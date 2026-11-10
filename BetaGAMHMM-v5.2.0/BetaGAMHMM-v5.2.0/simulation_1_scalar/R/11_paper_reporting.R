# Paper-compatible reporting utilities.
#
# These functions intentionally preserve the presentation used in the submitted
# manuscript: a candidate-K comparison table, a compact Monte Carlo performance
# table, and a four-panel figure of estimated state-dependent mean curves.
# Additional revision analyses are exported separately as supplementary tables.

format_mean_sd <- function(x, digits = 3L, multiplier = 1) {
  x <- x[is.finite(x)] * multiplier
  if (!length(x)) return("--")
  sprintf(paste0("%.", digits, "f (%.", digits, "f)"), mean(x), stats::sd(x))
}

latex_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([%&#_$])", "\\\\\\1", x)
  x <- gsub("~", "\\\\textasciitilde{}", x, fixed = TRUE)
  x <- gsub("\\^", "\\\\textasciicircum{}", x)
  x
}

write_simple_latex_table <- function(df, file, caption = NULL, label = NULL,
                                     align = NULL, digits = NULL) {
  ensure_dir(dirname(file))
  x <- df
  if (!is.null(digits)) {
    num <- vapply(x, is.numeric, logical(1))
    x[num] <- lapply(x[num], function(z) ifelse(is.na(z), NA, round(z, digits)))
  }
  x[] <- lapply(x, function(z) {
    z <- ifelse(is.na(z), "--", as.character(z))
    latex_escape(z)
  })
  if (is.null(align)) align <- paste0("l", paste(rep("c", ncol(x) - 1L), collapse = ""))
  lines <- c("\\begin{table}[htbp]", "\\centering")
  if (!is.null(caption)) lines <- c(lines, paste0("\\caption{", caption, "}"))
  if (!is.null(label)) lines <- c(lines, paste0("\\label{", label, "}"))
  lines <- c(lines, paste0("\\begin{tabular}{", align, "}"), "\\toprule")
  lines <- c(lines, paste(latex_escape(names(df)), collapse = " & "), "\\\\")
  lines <- c(lines, "\\midrule")
  for (i in seq_len(nrow(x))) {
    lines <- c(lines, paste(unlist(x[i, , drop = FALSE]), collapse = " & "), "\\\\")
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(lines, file)
  invisible(file)
}

paper_candidate_k_table <- function(grid_results, selection_rule, K_grid) {
  sel <- select_with_rule(grid_results, selection_rule)
  by_k <- sel$candidate_by_K
  rows <- lapply(K_grid, function(K) {
    dd <- by_k[by_k$K == K, , drop = FALSE]
    if (!nrow(dd)) {
      return(data.frame(
        K = as.character(K),
        lambda = "--",
        AIC_EDF = "--",
        BIC_EDF = "--",
        ICL_EDF = "--",
        stringsAsFactors = FALSE
      ))
    }
    data.frame(
      K = as.character(K),
      lambda = format(dd$lambda[1L], trim = TRUE, scientific = FALSE),
      AIC_EDF = sprintf("%.1f", dd$aic_edf[1L]),
      BIC_EDF = sprintf("%.1f", dd$bic_edf[1L]),
      ICL_EDF = sprintf("%.1f", dd$icl_edf[1L]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

run_one_conditional_estimation_rep <- function(rep_id, cfg) {
  spec <- cfg$scenarios$baseline
  seed_i <- cfg$seed + 2500000L + rep_id * 100L
  sim <- simulate_revision_scenario(spec, cfg$basis, seed_i)

  grid <- run_grid_search_hmm(
    y = sim$y,
    basis_spec = sim$basis_spec,
    K_grid = spec$K_true,
    lambda_grid = cfg$fit$lambda_grid,
    n_start = cfg$paper_reporting$n_start,
    max_iter = cfg$paper_reporting$max_iter,
    tol = cfg$fit$tol,
    phi_lower = cfg$fit$phi_lower,
    phi_upper = cfg$fit$phi_upper,
    optimizer_maxit = cfg$fit$optimizer_maxit,
    diagnostic = cfg$diagnostic,
    n_cores = 1L,
    n_start_cores = 1L,
    progress = FALSE,
    seed = seed_i + 1L
  )

  sel <- tryCatch(
    select_best_model_hierarchical(
      grid,
      criterion_primary = "bic_edf",
      criterion_secondary = "icl_edf",
      lambda_criterion_primary = "aic_edf",
      lambda_criterion_secondary = "bic_edf",
      exclude_degenerate = FALSE,
      require_converged = TRUE,
      require_monotone = TRUE,
      enforce_occupancy = FALSE,
      equivalence_margin = 0,
      prefer_richer_K_within_margin = FALSE
    ),
    error = function(e) NULL
  )

  if (is.null(sel)) {
    return(list(
      metrics = data.frame(
        rep = rep_id, success = FALSE, selected_lambda = NA_real_,
        curve_rmse = NA_real_, phi_rmse = NA_real_, A_rmse = NA_real_,
        viterbi_accuracy = NA_real_, stringsAsFactors = FALSE
      ),
      curves = NULL,
      confusion = NULL
    ))
  }

  fit <- align_fit_to_reference(
    sel$selected$fit,
    list(beta = sim$beta_true, phi = sim$phi_true, pi = sim$pi_true, A = sim$A_true),
    basis_spec = sim$basis_spec
  )
  comp <- compare_fit_to_truth(fit, sim)

  xg <- seq(min(sim$x), max(sim$x), length.out = cfg$paper_reporting$x_grid_n)
  Bg <- predict_basis(xg, sim$basis_spec)
  eta_funs <- default_eta_functions(spec$K_true, scale = spec$eta_scale %||% 1)

  curves <- do.call(rbind, lapply(seq_len(spec$K_true), function(k) {
    data.frame(
      rep = rep_id,
      state = factor(k, levels = seq_len(spec$K_true)),
      x = xg,
      mu_hat = inv_logit(drop(Bg %*% fit$beta[[k]])),
      mu_true = inv_logit(eta_funs[[k]](xg)),
      stringsAsFactors = FALSE
    )
  }))

  zhat <- viterbi_decode(sim$y, sim$basis_spec, fit)
  cm <- table(
    true_state = factor(sim$z, levels = seq_len(spec$K_true)),
    decoded_state = factor(zhat, levels = seq_len(spec$K_true))
  )
  confusion <- as.data.frame(cm, stringsAsFactors = FALSE)
  names(confusion) <- c("true_state", "decoded_state", "count")
  confusion$rep <- rep_id

  list(
    metrics = data.frame(
      rep = rep_id,
      success = TRUE,
      selected_lambda = as.numeric(sel$winner_row$lambda[1L]),
      curve_rmse = comp$curve_rmse,
      phi_rmse = comp$phi_rmse,
      A_rmse = comp$A_rmse,
      viterbi_accuracy = comp$viterbi_accuracy,
      stringsAsFactors = FALSE
    ),
    curves = curves,
    confusion = confusion
  )
}

run_conditional_estimation_study <- function(cfg, output_dir, resume = TRUE) {
  ensure_dir(output_dir)
  cp_dir <- file.path(output_dir, "checkpoints")
  ensure_dir(cp_dir)
  n_rep <- cfg$paper_reporting$n_rep
  chunks <- chunk_indices(n_rep, cfg$paper_reporting$chunk_size)
  metrics_all <- list()
  curves_all <- list()
  confusion_all <- list()

  for (chunk_id in seq_along(chunks)) {
    reps <- chunks[[chunk_id]]
    cp <- file.path(cp_dir, sprintf("conditional_chunk_%03d.rds", chunk_id))
    if (resume && file.exists(cp)) {
      obj <- readRDS(cp)
    } else {
      worker <- function(r) run_one_conditional_estimation_rep(r, cfg)
      ans <- parallel_lapply_crossplatform(
        reps, worker,
        n_cores = min(cfg$parallel$n_cores, length(reps)),
        exports = c("worker", "cfg"),
        project_root = getwd(),
        seed = cfg$seed + 2600000L + chunk_id,
        progress = TRUE,
        progress_label = sprintf("Conditional estimation chunk %d/%d", chunk_id, length(chunks))
      )
      obj <- list(
        metrics = merge_data_frames(lapply(ans, `[[`, "metrics")),
        curves = merge_data_frames(Filter(Negate(is.null), lapply(ans, `[[`, "curves"))),
        confusion = merge_data_frames(Filter(Negate(is.null), lapply(ans, `[[`, "confusion")))
      )
      save_rds_safe(obj, cp)
    }
    metrics_all[[length(metrics_all) + 1L]] <- obj$metrics
    curves_all[[length(curves_all) + 1L]] <- obj$curves
    confusion_all[[length(confusion_all) + 1L]] <- obj$confusion %||% data.frame()
  }

  list(
    metrics = merge_data_frames(metrics_all),
    curves = merge_data_frames(curves_all),
    confusion = merge_data_frames(confusion_all)
  )
}

aggregate_confusion_matrix <- function(confusion, K) {
  if (is.null(confusion) || !nrow(confusion)) return(list(counts=data.frame(), row_normalized=data.frame()))
  agg <- stats::aggregate(count ~ true_state + decoded_state, confusion, sum)
  full <- expand.grid(true_state=seq_len(K), decoded_state=seq_len(K))
  agg$true_state <- as.integer(as.character(agg$true_state))
  agg$decoded_state <- as.integer(as.character(agg$decoded_state))
  counts <- merge(full, agg, by=c("true_state","decoded_state"), all.x=TRUE)
  counts$count[is.na(counts$count)] <- 0
  denom <- ave(counts$count, counts$true_state, FUN=sum)
  norm <- counts
  norm$proportion <- ifelse(denom > 0, norm$count / denom, NA_real_)
  list(counts=counts, row_normalized=norm)
}

paper_monte_carlo_table <- function(metrics) {
  ok <- metrics[metrics$success %in% TRUE, , drop = FALSE]
  data.frame(
    Metric = c(
      "Mean curve RMSE",
      "Mean phi RMSE",
      "Mean A RMSE",
      "Mean Viterbi accuracy (%)"
    ),
    `Mean (SD)` = c(
      format_mean_sd(ok$curve_rmse, digits = 3L),
      format_mean_sd(ok$phi_rmse, digits = 2L),
      format_mean_sd(ok$A_rmse, digits = 3L),
      format_mean_sd(ok$viterbi_accuracy, digits = 1L, multiplier = 100)
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

paper_lambda_frequency_table <- function(metrics) {
  ok <- metrics[metrics$success %in% TRUE & is.finite(metrics$selected_lambda), , drop = FALSE]
  if (!nrow(ok)) return(data.frame())
  tab <- as.data.frame(table(ok$selected_lambda), stringsAsFactors = FALSE)
  names(tab) <- c("lambda", "count")
  tab$percentage <- 100 * tab$count / sum(tab$count)
  tab
}

make_paper_figure1 <- function(curves, representative_sim, file_png, file_pdf,
                               include_observed_points = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required for paper figures.")
  ensure_dir(dirname(file_png))

  truth <- unique(curves[, c("state", "x", "mu_true")])
  p <- ggplot2::ggplot()

  if (isTRUE(include_observed_points)) {
    pts <- data.frame(
      state = factor(representative_sim$z, levels = seq_len(representative_sim$K_true)),
      x = representative_sim$x,
      y = representative_sim$y
    )
    # A light subsample avoids obscuring the Monte Carlo curves while satisfying
    # the request to display observed responses in the corresponding state panel.
    if (nrow(pts) > 1200L) {
      set.seed(20260915L)
      pts <- pts[sample.int(nrow(pts), 1200L), , drop = FALSE]
    }
    p <- p + ggplot2::geom_point(
      data = pts,
      ggplot2::aes(x = x, y = y),
      size = 0.45, alpha = 0.16
    )
  }

  p <- p +
    ggplot2::geom_line(
      data = curves,
      ggplot2::aes(x = x, y = mu_hat, group = interaction(rep, state)),
      linewidth = 0.25, alpha = 0.18
    ) +
    ggplot2::geom_line(
      data = truth,
      ggplot2::aes(x = x, y = mu_true),
      linewidth = 0.8, linetype = "dashed"
    ) +
    ggplot2::facet_wrap(~ state, ncol = 2, labeller = ggplot2::labeller(state = function(x) paste("State", x))) +
    ggplot2::labs(x = "Covariate", y = "Conditional mean") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "none"
    )

  ggplot2::ggsave(file_png, p, width = 7.2, height = 5.8, dpi = 300)
  ggplot2::ggsave(file_pdf, p, width = 7.2, height = 5.8, device = "pdf")
  invisible(p)
}

copy_table_if_exists <- function(src, dst) {
  if (file.exists(src)) {
    ensure_dir(dirname(dst))
    file.copy(src, dst, overwrite = TRUE)
  }
}

export_supplementary_revision_tables <- function(output_root, supplement_dir) {
  ensure_dir(supplement_dir)
  files <- c(
    "01_state_selection/tables/selection_performance_BIC_ICL.csv",
    "01_state_selection/tables/BIC_estimation_metrics.csv",
    "02_residual_diagnostics/pseudo_residual_diagnostics.csv",
    "03_basis_sensitivity/basis_dimension_sensitivity.csv",
    "04_multicovariate/multicovariate_summary.csv",
    "05_smoothing_sensitivity/smoothing_penalty_comparison.csv",
    "06_bootstrap_coverage/bootstrap_coverage_summary.csv",
    "07_hmmtmb_benchmark/hmmtmb_comparison_summary.csv",
    "07_hmmtmb_benchmark/hmmtmb_paired_summary.csv",
    "08_paper/computational_scaling.csv"
  )
  for (rel in files) {
    src <- file.path(output_root, rel)
    dst <- file.path(supplement_dir, basename(rel))
    copy_table_if_exists(src, dst)
  }
  for (nm in c("pseudo_residual_acf.png", "pseudo_residual_qq.png")) {
    src <- file.path(output_root, "02_residual_diagnostics", "figures", nm)
    if (file.exists(src)) file.copy(src, file.path(supplement_dir, nm), overwrite = TRUE)
  }
  invisible(TRUE)
}
