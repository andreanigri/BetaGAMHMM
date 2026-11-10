# Final helper analyses retained because they directly address reviewer requests:
# (i) bootstrap interval coverage and (ii) common versus state-specific smoothing.

bootstrap_basic_interval <- function(samples, theta_hat, level = 0.95) {
  samples <- samples[is.finite(samples)]
  if (!length(samples)) return(c(NA_real_, NA_real_))
  alpha <- 1 - level
  q <- stats::quantile(samples, probs = c(1 - alpha / 2, alpha / 2), names = FALSE, type = 7)
  sort(2 * theta_hat - q)
}

bootstrap_bc_interval <- function(samples, theta_hat, level = 0.95) {
  samples <- samples[is.finite(samples)]
  B <- length(samples)
  if (B < 5L) return(c(NA_real_, NA_real_))

  # Bias-corrected percentile interval (BC; acceleration fixed at zero).
  # This is distinct from BCa, which would additionally require a jackknife
  # estimate of the acceleration term.
  p_less <- mean(samples < theta_hat)
  p_less <- min(max(p_less, 0.5 / B), 1 - 0.5 / B)
  z0 <- stats::qnorm(p_less)
  alpha <- 1 - level
  z_alpha <- stats::qnorm(c(alpha / 2, 1 - alpha / 2))
  probs <- stats::pnorm(2 * z0 + z_alpha)
  probs <- pmin(pmax(probs, 0.5 / B), 1 - 0.5 / B)
  stats::quantile(samples, probs = probs, names = FALSE, type = 7)
}

run_bootstrap_coverage_followup <- function(n_rep = 30L,
                                            B_boot = 100L,
                                            T = 700L,
                                            K = 4L,
                                            lambda = 0.1,
                                            n_start_fit = 4L,
                                            n_start_boot = 2L,
                                            n_cores = 1L,
                                            seed = 16300L) {
  one <- function(r) {
    sim_seed <- seed + r * 1000L
    sim <- simulate_default_scenario(T = T, K = K, diagonal = 0.90, seed = sim_seed)
    fit <- fit_beta_gam_hmm(
      sim$y, sim$basis_spec, K = K, lambda = lambda,
      n_start = n_start_fit, max_iter = 160L, optimizer_maxit = 80L,
      seed = sim_seed + 1L
    )
    fit <- align_fit_to_reference(
      fit,
      list(beta = sim$beta_true, phi = sim$phi_true, pi = sim$pi_true, A = sim$A_true),
      sim$basis_spec
    )

    boot <- bootstrap_beta_gam_hmm(
      fit, B_boot = B_boot, n_start = n_start_boot, max_iter = 120L,
      optimizer_maxit = 60L, n_cores = 1L, seed = sim_seed + 100L
    )

    phi_cover_pct <- phi_cover_basic <- phi_cover_bc <- logical(K)
    phi_width_pct <- phi_width_basic <- phi_width_bc <- numeric(K)
    phi_bias <- numeric(K)

    for (k in seq_len(K)) {
      samp <- boot$phi_samples[, k]
      ci_pct <- boot$phi_ci[k, ]
      ci_basic <- bootstrap_basic_interval(samp, fit$phi[k], level = boot$level)
      ci_bc <- bootstrap_bc_interval(samp, fit$phi[k], level = boot$level)

      phi_cover_pct[k] <- sim$phi_true[k] >= ci_pct[1L] && sim$phi_true[k] <= ci_pct[2L]
      phi_cover_basic[k] <- sim$phi_true[k] >= ci_basic[1L] && sim$phi_true[k] <= ci_basic[2L]
      phi_cover_bc[k] <- sim$phi_true[k] >= ci_bc[1L] && sim$phi_true[k] <= ci_bc[2L]

      phi_width_pct[k] <- diff(ci_pct)
      phi_width_basic[k] <- diff(ci_basic)
      phi_width_bc[k] <- diff(ci_bc)
      phi_bias[k] <- mean(samp, na.rm = TRUE) - fit$phi[k]
    }

    A_cover <- matrix(FALSE, K, K)
    for (i in seq_len(K)) for (j in seq_len(K)) {
      A_cover[i, j] <- sim$A_true[i, j] >= boot$A_ci[i, j, 1L] &&
        sim$A_true[i, j] <= boot$A_ci[i, j, 2L]
    }

    # Pointwise coverage of the state-specific mean curves.
    mu_ci <- boot$mu_curve_ci
    eta_funs <- default_eta_functions(K, scale = sim$scenario$eta_scale %||% 1)
    curve_cover <- logical(nrow(mu_ci))
    for (ii in seq_len(nrow(mu_ci))) {
      kk <- as.integer(sub("State ", "", mu_ci$state[ii], fixed = TRUE))
      mu_true <- inv_logit(eta_funs[[kk]](mu_ci$x[ii]))
      curve_cover[ii] <- mu_true >= mu_ci$lower[ii] && mu_true <= mu_ci$upper[ii]
    }

    data.frame(
      rep = r,
      phi_coverage_percentile = mean(phi_cover_pct),
      phi_coverage_basic = mean(phi_cover_basic),
      phi_coverage_bc = mean(phi_cover_bc),
      phi_width_percentile = mean(phi_width_pct),
      phi_width_basic = mean(phi_width_basic),
      phi_width_bc = mean(phi_width_bc),
      phi_bootstrap_bias_rmse = sqrt(mean(phi_bias^2)),
      transition_coverage_percentile = mean(A_cover),
      mean_curve_pointwise_coverage = mean(curve_cover),
      bootstrap_success_rate = boot$B_eff / B_boot,
      stringsAsFactors = FALSE
    )
  }

  ans <- parallel_lapply_crossplatform(
    seq_len(n_rep), one,
    n_cores = min(as.integer(n_cores %||% 1L), n_rep),
    exports = c(
      "one", "B_boot", "T", "K", "lambda", "n_start_fit", "n_start_boot", "seed",
      "bootstrap_basic_interval", "bootstrap_bc_interval"
    ),
    project_root = getwd(), seed = seed, progress = TRUE,
    progress_label = sprintf("Bootstrap follow-up: %d outer replicates", n_rep)
  )
  do.call(rbind, ans)
}

summarize_bootstrap_coverage_followup <- function(x) {
  metrics <- setdiff(names(x), "rep")
  do.call(rbind, lapply(metrics, function(m) {
    v <- x[[m]]
    data.frame(
      metric = m,
      n = sum(is.finite(v)),
      mean = mean(v, na.rm = TRUE),
      sd = stats::sd(v, na.rm = TRUE),
      median = stats::median(v, na.rm = TRUE),
      q25 = stats::quantile(v, 0.25, na.rm = TRUE, names = FALSE),
      q75 = stats::quantile(v, 0.75, na.rm = TRUE, names = FALSE),
      stringsAsFactors = FALSE
    )
  }))
}

run_smoothing_penalty_followup <- function(T = 1200L,
                                           K = 4L,
                                           lambda_grid = c(0.01, 0.05, 0.1, 0.5, 1),
                                           n_start = 4L,
                                           n_cycles = 2L,
                                           seed = 19100L) {
  sim <- simulate_default_scenario(T = T, K = K, diagonal = 0.92, seed = seed)

  common_grid <- run_grid_search_hmm(
    sim$y, sim$basis_spec,
    K_grid = K, lambda_grid = lambda_grid,
    n_start = n_start, max_iter = 150L, optimizer_maxit = 80L,
    n_cores = 1L, n_start_cores = 1L, progress = FALSE,
    seed = seed + 1L
  )
  common_sel <- select_best_model_hierarchical(
    common_grid,
    criterion_primary = "bic_edf", criterion_secondary = "icl_edf",
    lambda_criterion_primary = "aic_edf", lambda_criterion_secondary = "bic_edf",
    exclude_degenerate = FALSE, enforce_occupancy = FALSE,
    equivalence_margin = 0, prefer_richer_K_within_margin = FALSE
  )
  lambda_common <- as.numeric(common_sel$winner_row$lambda[1L])

  ss <- run_state_specific_lambda_sensitivity(
    sim,
    common_lambda = lambda_common,
    candidate_grid = lambda_grid,
    n_cycles = n_cycles,
    n_start = n_start,
    max_iter = 150L,
    seed = seed + 1000L
  )

  ss$comparison$lambda_specification <- c(
    sprintf("common=%g", lambda_common),
    paste0("state-specific=", paste(format(ss$selected_lambda_by_state, trim = TRUE), collapse = " | "))
  )

  list(
    common_grid = summarize_grid_results(common_grid),
    common_selected_lambda = lambda_common,
    selected_lambda_by_state = ss$selected_lambda_by_state,
    comparison = ss$comparison,
    trace = ss$trace
  )
}
