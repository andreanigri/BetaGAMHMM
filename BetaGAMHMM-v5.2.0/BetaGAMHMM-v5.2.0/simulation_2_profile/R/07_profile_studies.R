# Targeted studies for the final repeated-profile formulation.

profile_mc_binomial_summary <- function(success, conf_level = 0.95) {
  success <- success[!is.na(success)]
  n <- length(success)
  if (!n) return(c(rate = NA_real_, mc_se = NA_real_, ci_low = NA_real_, ci_high = NA_real_))
  x <- sum(success %in% TRUE)
  rate <- x / n
  mc_se <- sqrt(rate * (1 - rate) / n)
  ci <- stats::binom.test(x, n, conf.level = conf_level)$conf.int
  c(rate = rate, mc_se = mc_se, ci_low = ci[1L], ci_high = ci[2L])
}

summarize_profile_selection <- function(results) {
  groups <- split(results, interaction(results$scenario, results$criterion, drop = TRUE))
  out <- lapply(groups, function(d) {
    valid <- is.finite(d$selected_K)
    ex <- if (any(valid)) profile_mc_binomial_summary(d$exact[valid]) else c(rate=NA,mc_se=NA,ci_low=NA,ci_high=NA)
    data.frame(
      scenario = d$scenario[1L], criterion = d$criterion[1L],
      n = nrow(d), failure_rate = mean(!valid), K_true = d$K_true[1L],
      exact_rate = unname(ex["rate"]),
      exact_mc_se = unname(ex["mc_se"]),
      exact_ci_low = unname(ex["ci_low"]),
      exact_ci_high = unname(ex["ci_high"]),
      underfit_rate = if (any(valid)) mean(d$underfit[valid], na.rm = TRUE) else NA_real_,
      overfit_rate = if (any(valid)) mean(d$overfit[valid], na.rm = TRUE) else NA_real_,
      mean_selected_K = if (any(valid)) mean(d$selected_K[valid], na.rm = TRUE) else NA_real_,
      mean_abs_K_error = if (any(valid)) mean(abs(d$selected_K[valid] - d$K_true[1L]), na.rm = TRUE) else NA_real_,
      precision_flag_rate = if (any(valid)) mean(d$precision_flag[valid] %in% TRUE, na.rm = TRUE) else NA_real_,
      all_states_visited_rate = mean(d$all_states_visited %in% TRUE, na.rm = TRUE),
      n_all_states_visited = sum(d$all_states_visited %in% TRUE),
      exact_rate_given_all_states_visited = if (any(valid & d$all_states_visited %in% TRUE))
        mean(d$exact[valid & d$all_states_visited %in% TRUE], na.rm = TRUE) else NA_real_,
      mean_min_true_state_share = mean(d$min_true_state_share, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  ans <- do.call(rbind, out)
  rownames(ans) <- NULL
  ans[order(ans$scenario, ans$criterion), , drop = FALSE]
}

summarize_profile_estimation <- function(results) {
  d <- results[results$exact %in% TRUE, , drop = FALSE]
  if (!nrow(d)) return(data.frame())
  groups <- split(d, interaction(d$scenario, d$criterion, drop = TRUE))
  out <- lapply(groups, function(x) data.frame(
    scenario = x$scenario[1L], criterion = x$criterion[1L], n_correct_K = nrow(x),
    viterbi_accuracy_mean = mean(x$viterbi_accuracy, na.rm = TRUE),
    viterbi_accuracy_sd = stats::sd(x$viterbi_accuracy, na.rm = TRUE),
    ari_mean = mean(x$ari, na.rm = TRUE),
    curve_rmse_mean = mean(x$curve_rmse, na.rm = TRUE),
    phi_rmse_mean = mean(x$phi_rmse, na.rm = TRUE),
    A_rmse_mean = mean(x$A_rmse, na.rm = TRUE),
    min_state_share_mean = mean(x$min_state_share, na.rm = TRUE),
    stringsAsFactors = FALSE
  ))
  ans <- do.call(rbind, out)
  rownames(ans) <- NULL
  ans[order(ans$scenario, ans$criterion), , drop = FALSE]
}

profile_coverage_worker <- function(job) {
  sim <- simulate_profile_scenario(job$scenario, job$seed)
  basis <- build_age_basis(sim$panel$ages, degree = job$degree, n_inner_knots = job$n_inner_knots)
  grid <- run_profile_grid(
    sim$panel, basis,
    K_grid = job$scenario$K, lambda_grid = job$lambda_grid,
    n_start = job$n_start_fit, n_cores = 1L,
    max_iter = job$max_iter, tol = job$tol,
    phi_lower = job$phi_lower, phi_upper = job$phi_upper,
    optimizer_maxit = job$optimizer_maxit,
    seed = job$seed + 10000L, progress = FALSE
  )
  tab <- summarize_profile_grid(grid)
  ok <- tab[is.na(tab$error) & tab$converged & tab$monotonicity_ok, , drop = FALSE]
  if (!nrow(ok)) return(NULL)
  best <- ok[order(ok$aic_edf, ok$lambda), , drop = FALSE][1L, ]
  idx <- which(vapply(grid, function(f) !inherits(f, "profile_grid_error") &&
    f$K == job$scenario$K && abs(f$lambda[1L] - best$lambda) < 1e-12, logical(1)))
  if (!length(idx)) return(NULL)
  fit <- align_profile_fit_to_truth(grid[[idx[1L]]], sim)

  boot <- tryCatch(run_beta_parametric_bootstrap(
    fit = fit, B = job$B_boot,
    n_start = job$n_start_boot, n_cores = 1L,
    max_iter = job$max_iter, tol = job$tol,
    phi_lower = job$phi_lower, phi_upper = job$phi_upper,
    optimizer_maxit = job$optimizer_maxit,
    level = job$level, seed = job$seed + 20000L
  ), error = function(e) NULL)
  if (is.null(boot)) return(NULL)

  K <- job$scenario$K
  methods <- c("percentile", "basic", "bc")
  out <- list(
    replicate = job$replicate,
    selected_lambda = best$lambda,
    bootstrap_success_rate = boot$success_rate
  )
  for (m in methods) {
    phi_ci <- boot$phi_ci_all[[m]]
    phi_cov <- vapply(seq_len(K), function(k) {
      sim$phi_true[k] >= phi_ci[k, 1L] && sim$phi_true[k] <= phi_ci[k, 2L]
    }, logical(1))
    A_ci <- boot$A_ci_all[[m]]
    A_cov <- matrix(FALSE, K, K)
    for (i in seq_len(K)) for (j in seq_len(K)) {
      A_cov[i, j] <- sim$A_true[i, j] >= A_ci$low[i, j] && sim$A_true[i, j] <= A_ci$high[i, j]
    }
    mu_ci <- boot$mu_ci_all[[m]]
    curve_cov <- (sim$mu_true >= mu_ci$low) & (sim$mu_true <= mu_ci$high)
    out[[paste0("phi_coverage_", m)]] <- mean(phi_cov)
    out[[paste0("A_coverage_", m, "_all")]] <- mean(A_cov)
    out[[paste0("A_coverage_", m, "_diagonal")]] <- mean(diag(A_cov))
    out[[paste0("curve_coverage_", m)]] <- mean(curve_cov, na.rm = TRUE)
    out[[paste0("mean_phi_width_", m)]] <- mean(phi_ci[, 2L] - phi_ci[, 1L], na.rm = TRUE)
    out[[paste0("mean_A_width_", m)]] <- mean(A_ci$high - A_ci$low, na.rm = TRUE)
    out[[paste0("mean_curve_width_", m)]] <- mean(mu_ci$high - mu_ci$low, na.rm = TRUE)
  }
  as.data.frame(out, stringsAsFactors = FALSE)
}

profile_basis_worker <- function(job) {
  sim <- simulate_profile_scenario(job$scenario, job$seed)
  basis <- build_age_basis(sim$panel$ages, degree = job$degree, n_inner_knots = job$n_inner_knots)
  grid <- run_profile_grid(
    sim$panel, basis, job$K_grid, job$lambda_grid,
    n_start = job$n_start, n_cores = 1L,
    max_iter = job$max_iter, tol = job$tol,
    phi_lower = job$phi_lower, phi_upper = job$phi_upper,
    optimizer_maxit = job$optimizer_maxit,
    seed = job$seed + 10000L, progress = FALSE
  )
  cmap <- profile_criterion_map()
  out <- lapply(names(cmap), function(nm) {
    sel <- tryCatch(select_profile_model(grid, cmap[[nm]]), error = function(e) NULL)
    if (is.null(sel)) return(NULL)
    met <- profile_fit_metrics(sel$fit, sim)
    data.frame(
      replicate = job$replicate, n_inner_knots = job$n_inner_knots,
      criterion = nm, K_true = job$scenario$K, selected_K = sel$fit$K,
      exact = sel$fit$K == job$scenario$K,
      selected_lambda = sel$fit$lambda[1L],
      viterbi_accuracy = met$viterbi_accuracy,
      curve_rmse = met$curve_rmse,
      precision_flag = isTRUE(sel$fit$precision_diagnostic$flagged),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, Filter(Negate(is.null), out))
}

run_profile_computational_benchmark <- function(cfg) {
  designs <- data.frame(
    T = c(55L, 55L, 55L, 100L, 250L),
    J = c(21L, 41L, 81L, 41L, 41L),
    stringsAsFactors = FALSE
  )
  rows <- vector("list", nrow(designs) * cfg$benchmark$n_repeat)
  rr <- 0L
  for (i in seq_len(nrow(designs))) {
    for (r in seq_len(cfg$benchmark$n_repeat)) {
      rr <- rr + 1L
      sc <- cfg$scenarios$baseline_T55
      sc$T <- designs$T[i]
      sc$J <- designs$J[i]
      sc$name <- paste0("benchmark_T", sc$T, "_J", sc$J)
      sim <- simulate_profile_scenario(sc, cfg$seed + 800000L + 1000L * i + r)
      basis <- build_age_basis(sim$panel$ages, cfg$basis$degree, cfg$basis$n_inner_knots)
      tm <- system.time({
        fit <- fit_beta_year_hmm(
          sim$panel, basis, K = 4L, lambda = 1,
          n_start = cfg$benchmark$n_start, n_cores = 1L,
          max_iter = cfg$fit$max_iter, tol = cfg$fit$tol,
          phi_lower = cfg$fit$phi_lower, phi_upper = cfg$fit$phi_upper,
          optimizer_maxit = cfg$fit$optimizer_maxit,
          seed = cfg$seed + 810000L + 1000L * i + r
        )
      })[["elapsed"]]
      rows[[rr]] <- data.frame(
        T = designs$T[i], J = designs$J[i], n_cells = designs$T[i] * designs$J[i],
        replicate_id = r, elapsed_seconds = tm, converged = fit$converged, n_iter = fit$n_iter
      )
    }
  }
  do.call(rbind, rows)
}
