# Revision-focused simulation studies.
#
# These routines separate (i) parameter estimation conditional on a model from
# (ii) recovery of the number of states.  The latter is evaluated repeatedly and
# under several model-selection/diagnostic rules using the same fitted grid.

revision_selection_rules <- function(cfg = NULL) {
  if (!is.null(cfg$selection_rules)) return(cfg$selection_rules)
  common <- list(
    lambda_criterion_primary = "aic_edf",
    lambda_criterion_secondary = "bic_edf",
    require_converged = TRUE,
    require_monotone = TRUE,
    exclude_degenerate = FALSE,
    enforce_occupancy = FALSE,
    equivalence_margin = 0,
    prefer_richer_K_within_margin = FALSE
  )
  list(
    BIC = modifyList(common, list(
      criterion_primary = "bic_edf", criterion_secondary = "icl_edf"
    )),
    ICL = modifyList(common, list(
      criterion_primary = "icl_edf", criterion_secondary = "bic_edf"
    ))
  )
}

select_with_rule <- function(results, rule) {
  select_best_model_hierarchical(
    results,
    criterion_primary = rule$criterion_primary %||% "bic_edf",
    criterion_secondary = rule$criterion_secondary %||% "icl_edf",
    lambda_criterion_primary = rule$lambda_criterion_primary %||% "aic_edf",
    lambda_criterion_secondary = rule$lambda_criterion_secondary %||% "bic_edf",
    exclude_degenerate = rule$exclude_degenerate %||% FALSE,
    require_converged = rule$require_converged %||% TRUE,
    require_monotone = rule$require_monotone %||% TRUE,
    enforce_occupancy = rule$enforce_occupancy %||% FALSE,
    min_state_share = rule$min_state_share %||% 0,
    min_state_mass = rule$min_state_mass %||% 0,
    max_start_loglik_gap = rule$max_start_loglik_gap %||% Inf,
    equivalence_margin = rule$equivalence_margin %||% 0,
    prefer_richer_K_within_margin = rule$prefer_richer_K_within_margin %||% FALSE,
    diagnostic_override = rule$diagnostic %||% NULL
  )
}

simulate_revision_scenario <- function(spec, basis_cfg, seed) {
  simulate_default_scenario(
    T = spec$T,
    K = spec$K_true,
    degree = basis_cfg$degree,
    n_inner_knots = basis_cfg$n_inner_knots,
    diagonal = if (is.finite(spec$diagonal %||% NA_real_)) spec$diagonal else 0.90,
    phi_true = spec$phi_true,
    seed = seed,
    boundary_padding_frac = basis_cfg$boundary_padding_frac %||% 0,
    knot_quantile_range = basis_cfg$knot_quantile_range %||% c(0, 1),
    eta_scale = spec$eta_scale %||% 1,
    A_true = spec$A_true %||% NULL
  )
}

run_one_revision_selection_rep <- function(rep_id, scenario_name, spec, cfg) {
  seed_i <- cfg$seed + match(scenario_name, names(cfg$scenarios)) * 100000L + rep_id * 100L
  sim <- simulate_revision_scenario(spec, cfg$basis, seed_i)

  started <- proc.time()[[3L]]
  results <- run_grid_search_hmm(
    y = sim$y,
    basis_spec = sim$basis_spec,
    K_grid = cfg$fit$K_grid,
    lambda_grid = cfg$fit$lambda_grid,
    n_start = cfg$fit$n_start,
    max_iter = cfg$fit$max_iter,
    tol = cfg$fit$tol,
    phi_lower = cfg$fit$phi_lower,
    phi_upper = cfg$fit$phi_upper,
    optimizer_maxit = cfg$fit$optimizer_maxit %||% 80L,
    diagnostic = cfg$diagnostic,
    n_cores = 1L,
    n_start_cores = 1L,
    progress = FALSE,
    seed = seed_i + 1L
  )
  elapsed <- proc.time()[[3L]] - started

  rules <- revision_selection_rules(cfg)
  rule_results <- evaluate_selection_rules(results, rules, K_true = spec$K_true)
  rule_results$rep <- rep_id
  rule_results$scenario <- scenario_name

  primary_name <- cfg$primary_rule %||% "BIC"
  primary <- tryCatch(select_with_rule(results, rules[[primary_name]]), error = function(e) NULL)
  fit_metrics <- data.frame(
    rep = rep_id,
    scenario = scenario_name,
    T = spec$T,
    K_true = spec$K_true,
    selected_K = NA_integer_, selected_lambda = NA_real_,
    selected_true_K = NA, viterbi_accuracy = NA_real_,
    phi_rmse = NA_real_, A_rmse = NA_real_, curve_rmse = NA_real_, pi_l1 = NA_real_,
    elapsed_seconds = elapsed,
    stringsAsFactors = FALSE
  )

  if (!is.null(primary)) {
    fit <- primary$selected$fit
    fit_metrics$selected_K <- fit$K
    fit_metrics$selected_lambda <- as.numeric(primary$winner_row$lambda[1L])
    fit_metrics$selected_true_K <- fit$K == spec$K_true
    if (fit$K == spec$K_true) {
      fit <- align_fit_to_reference(
        fit,
        list(beta = sim$beta_true, phi = sim$phi_true, pi = sim$pi_true, A = sim$A_true),
        basis_spec = sim$basis_spec
      )
      comp <- compare_fit_to_truth(fit, sim)
      fit_metrics$viterbi_accuracy <- comp$viterbi_accuracy
      fit_metrics$phi_rmse <- comp$phi_rmse
      fit_metrics$A_rmse <- comp$A_rmse
      fit_metrics$curve_rmse <- comp$curve_rmse
      fit_metrics$pi_l1 <- comp$pi_l1
    }
  }

  grid_summary <- summarize_grid_results(results)
  grid_summary$rep <- rep_id
  grid_summary$scenario <- scenario_name

  list(rule_results = rule_results, fit_metrics = fit_metrics, grid_summary = grid_summary)
}

mc_binomial_summary <- function(success, conf_level = 0.95) {
  success <- success[!is.na(success)]
  n <- length(success)
  if (!n) return(c(rate = NA_real_, mc_se = NA_real_, ci_low = NA_real_, ci_high = NA_real_))
  x <- sum(success %in% TRUE)
  rate <- x / n
  mc_se <- sqrt(rate * (1 - rate) / n)
  ci <- stats::binom.test(x, n, conf.level = conf_level)$conf.int
  c(rate = rate, mc_se = mc_se, ci_low = ci[1L], ci_high = ci[2L])
}

summarize_revision_selection <- function(rule_results) {
  if (nrow(rule_results) == 0L) return(data.frame())
  keys <- unique(rule_results[, c("scenario", "rule")])
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    dd <- rule_results[
      rule_results$scenario == keys$scenario[i] & rule_results$rule == keys$rule[i],
      , drop = FALSE
    ]
    valid <- !is.na(dd$selected_K)
    true_K <- unique(dd$K_true %||% NA_integer_)
    if (!"K_true" %in% names(dd)) true_K <- NA_integer_
    k0 <- if (length(true_K)) true_K[1L] else NA_integer_
    ex <- if (any(valid)) mc_binomial_summary(dd$selected_true_K[valid]) else c(rate=NA,mc_se=NA,ci_low=NA,ci_high=NA)
    data.frame(
      scenario = keys$scenario[i],
      rule = keys$rule[i],
      n_rep = nrow(dd),
      failure_rate = mean(!valid),
      correct_K_rate = unname(ex["rate"]),
      correct_K_mc_se = unname(ex["mc_se"]),
      correct_K_ci_low = unname(ex["ci_low"]),
      correct_K_ci_high = unname(ex["ci_high"]),
      underfit_rate = if (any(valid) && is.finite(k0)) mean(dd$selected_K[valid] < k0) else NA_real_,
      overfit_rate = if (any(valid) && is.finite(k0)) mean(dd$selected_K[valid] > k0) else NA_real_,
      mean_selected_K = if (any(valid)) mean(dd$selected_K[valid]) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

run_revision_selection_study <- function(cfg, resume = TRUE) {
  ensure_dir(cfg$output_dir)
  checkpoint_dir <- file.path(cfg$output_dir, "checkpoints", "selection")
  ensure_dir(checkpoint_dir)
  scenarios <- cfg$scenarios
  n_rep <- cfg$monte_carlo$n_rep
  chunk_size <- cfg$monte_carlo$chunk_size %||% 5L

  all_rules <- list()
  all_metrics <- list()
  all_grids <- list()

  for (scenario_name in cfg$monte_carlo$scenarios) {
    spec <- scenarios[[scenario_name]]
    if (is.null(spec)) stop("Unknown revision scenario: ", scenario_name)
    chunks <- chunk_indices(n_rep, chunk_size)

    for (chunk_id in seq_along(chunks)) {
      reps <- chunks[[chunk_id]]
      cp <- file.path(checkpoint_dir, sprintf("%s_chunk_%03d.rds", scenario_name, chunk_id))
      if (resume && file.exists(cp)) {
        obj <- readRDS(cp)
      } else {
        worker <- function(r) run_one_revision_selection_rep(r, scenario_name, spec, cfg)
        ans <- parallel_lapply_crossplatform(
          reps, worker,
          n_cores = min(cfg$parallel$n_cores %||% 1L, length(reps)),
          exports = c("worker", "scenario_name", "spec", "cfg"),
          project_root = getwd(),
          seed = cfg$seed + chunk_id,
          progress = TRUE,
          progress_label = sprintf("%s: Monte Carlo chunk %d/%d", scenario_name, chunk_id, length(chunks))
        )
        obj <- list(
          rule_results = merge_data_frames(lapply(ans, `[[`, "rule_results")),
          fit_metrics = merge_data_frames(lapply(ans, `[[`, "fit_metrics")),
          grid_summary = merge_data_frames(lapply(ans, `[[`, "grid_summary"))
        )
        save_rds_safe(obj, cp)
      }
      all_rules[[length(all_rules) + 1L]] <- obj$rule_results
      all_metrics[[length(all_metrics) + 1L]] <- obj$fit_metrics
      all_grids[[length(all_grids) + 1L]] <- obj$grid_summary
    }
  }

  rules <- merge_data_frames(all_rules)
  metrics <- merge_data_frames(all_metrics)
  grids <- merge_data_frames(all_grids)
  # Store K_true explicitly in the rule-level table for under/over-selection rates.
  k_map <- vapply(cfg$scenarios, `[[`, integer(1), "K_true")
  rules$K_true <- unname(k_map[rules$scenario])
  performance <- summarize_revision_selection(rules)

  list(rule_results = rules, fit_metrics = metrics, grid_summary = grids, performance = performance)
}

# Sensitivity to the basis dimension M.  The data set is held fixed and the
# smoothing parameter is re-selected for every candidate basis dimension.
run_basis_dimension_sensitivity <- function(sim,
                                            n_inner_knots = c(4L, 6L, 8L, 10L),
                                            lambda_grid = c(0.01, 0.05, 0.1, 0.5, 1),
                                            n_start = 5L,
                                            max_iter = 150L,
                                            n_cores = 1L,
                                            diagnostic = list(
                                              phi_max = 500, eps = 1e-6, s_thresh = 2L,
                                              delta_abs = 50, delta_sum = 100
                                            ),
                                            seed = 4100L) {
  if (is.data.frame(sim$x)) stop("This basis-dimension sensitivity helper expects one covariate.")
  rows <- lapply(seq_along(n_inner_knots), function(i) {
    m <- n_inner_knots[i]
    basis <- build_bspline_basis(
      sim$x,
      degree = sim$basis_spec$degree,
      n_inner_knots = m,
      boundary_padding_frac = sim$basis_spec$boundary_padding_frac %||% 0
    )
    grid <- run_grid_search_hmm(
      sim$y, basis,
      K_grid = sim$K_true, lambda_grid = lambda_grid,
      n_start = n_start, max_iter = max_iter,
      optimizer_maxit = 80L,
      diagnostic = diagnostic,
      n_cores = n_cores, n_start_cores = 1L,
      progress = FALSE, seed = seed + i * 100L
    )
    sel <- select_best_model_hierarchical(
      grid,
      criterion_primary = "bic_edf", criterion_secondary = "icl_edf",
      lambda_criterion_primary = "aic_edf", lambda_criterion_secondary = "bic_edf",
      exclude_degenerate = FALSE, enforce_occupancy = FALSE,
      equivalence_margin = 0, prefer_richer_K_within_margin = FALSE
    )
    fit <- sel$selected$fit
    eta_funs <- default_eta_functions(sim$K_true, scale = sim$scenario$eta_scale %||% 1)
    beta_ref <- project_eta_functions_to_basis(sim$x, basis, eta_funs)
    fit <- align_fit_to_reference(
      fit,
      list(beta = beta_ref, phi = sim$phi_true, pi = sim$pi_true, A = sim$A_true),
      basis_spec = basis
    )
    xg <- seq(min(sim$x), max(sim$x), length.out = 300L)
    Bg <- predict_basis(xg, basis)
    mu_hat <- vapply(seq_len(sim$K_true), function(k) inv_logit(drop(Bg %*% fit$beta[[k]])), numeric(length(xg)))
    mu_true <- vapply(eta_funs, function(f) inv_logit(f(xg)), numeric(length(xg)))
    data.frame(
      n_inner_knots = m,
      basis_dimension = ncol(basis$B),
      selected_lambda = sel$winner_row$lambda[1L],
      edf_total = sum(fit$edf_spline),
      bic_edf = fit$criteria$bic_edf,
      curve_rmse = sqrt(mean((mu_hat - mu_true)^2)),
      timing_seconds = sum(summarize_grid_results(grid)$timing_seconds, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# Repeated two-covariate experiment.  K is fixed at its generating value so that
# the study isolates estimation/scalability as p increases from one to two.
run_multicovariate_study <- function(n_rep = 30L,
                                     T = 800L,
                                     lambda_grid = c(0.01, 0.05, 0.1, 0.5, 1),
                                     n_start = 4L,
                                     max_iter = 150L,
                                     n_cores = 1L,
                                     seed = 5200L,
                                     diagnostic = list(
                                       phi_max = 500, eps = 1e-6, s_thresh = 2L,
                                       delta_abs = 50, delta_sum = 100
                                     )) {
  one <- function(r) {
    sim <- simulate_multicovariate_scenario(T = T, seed = seed + r * 100L)
    started <- proc.time()[[3L]]
    grid <- run_grid_search_hmm(
      sim$y, sim$basis_spec,
      K_grid = sim$K_true, lambda_grid = lambda_grid,
      n_start = n_start, max_iter = max_iter,
      optimizer_maxit = 80L,
      diagnostic = diagnostic,
      n_cores = 1L, n_start_cores = 1L, progress = FALSE,
      seed = seed + r * 100L + 1L
    )
    sel <- select_best_model_hierarchical(
      grid,
      criterion_primary = "bic_edf", criterion_secondary = "icl_edf",
      lambda_criterion_primary = "aic_edf", lambda_criterion_secondary = "bic_edf",
      exclude_degenerate = FALSE, enforce_occupancy = FALSE,
      equivalence_margin = 0, prefer_richer_K_within_margin = FALSE
    )
    fit <- sel$selected$fit
    fit <- align_fit_to_reference(
      fit,
      list(beta = sim$beta_true, phi = sim$phi_true, pi = sim$pi_true, A = sim$A_true),
      basis_spec = sim$basis_spec
    )
    mu_hat <- vapply(
      seq_len(sim$K_true),
      function(k) inv_logit(drop(sim$basis_spec$B %*% fit$beta[[k]])),
      numeric(T)
    )
    mu_true <- sim$mu_by_state
    data.frame(
      rep = r,
      T = T,
      p = 2L,
      selected_lambda = sel$winner_row$lambda[1L],
      phi_rmse = sqrt(mean((fit$phi - sim$phi_true)^2)),
      A_rmse = sqrt(mean((fit$A - sim$A_true)^2)),
      mean_curve_rmse = sqrt(mean((mu_hat - mu_true)^2)),
      viterbi_accuracy = mean(viterbi_decode(sim$y, sim$basis_spec, fit) == sim$z),
      converged = fit$converged,
      elapsed_seconds = proc.time()[[3L]] - started,
      stringsAsFactors = FALSE
    )
  }
  ans <- parallel_lapply_crossplatform(
    seq_len(n_rep), one,
    n_cores = min(n_cores, n_rep),
    exports = c("one", "T", "lambda_grid", "n_start", "max_iter", "diagnostic", "seed"),
    project_root = getwd(), seed = seed, progress = TRUE,
    progress_label = sprintf("Two-covariate study: %d replicates", n_rep)
  )
  do.call(rbind, ans)
}

# Optional nested bootstrap study for empirical coverage.  This is intentionally
# separated from the main selection study because it is computationally costly.
run_bootstrap_coverage_study <- function(n_rep = 30L,
                                         B_boot = 50L,
                                         T = 700L,
                                         K = 4L,
                                         lambda = 0.1,
                                         n_start_fit = 4L,
                                         n_start_boot = 2L,
                                         n_cores = 1L,
                                         seed = 6300L) {
  one <- function(r) {
    sim <- simulate_default_scenario(T = T, K = K, diagonal = 0.90, seed = seed + r * 1000L)
    fit <- fit_beta_gam_hmm(
      sim$y, sim$basis_spec, K = K, lambda = lambda,
      n_start = n_start_fit, max_iter = 160L, optimizer_maxit = 80L,
      seed = seed + r * 1000L + 1L
    )
    fit <- align_fit_to_reference(
      fit,
      list(beta = sim$beta_true, phi = sim$phi_true, pi = sim$pi_true, A = sim$A_true),
      sim$basis_spec
    )
    boot <- bootstrap_beta_gam_hmm(
      fit, B_boot = B_boot, n_start = n_start_boot, max_iter = 120L,
      optimizer_maxit = 60L, n_cores = 1L, seed = seed + r * 1000L + 100L
    )
    phi_cover <- vapply(seq_len(K), function(k) {
      sim$phi_true[k] >= boot$phi_ci[k, 1L] && sim$phi_true[k] <= boot$phi_ci[k, 2L]
    }, logical(1))
    A_cover <- matrix(FALSE, K, K)
    for (i in seq_len(K)) for (j in seq_len(K)) {
      A_cover[i, j] <- sim$A_true[i, j] >= boot$A_ci[i, j, 1L] &&
        sim$A_true[i, j] <= boot$A_ci[i, j, 2L]
    }
    data.frame(
      rep = r,
      phi_coverage = mean(phi_cover),
      transition_coverage = mean(A_cover),
      bootstrap_success_rate = boot$B_eff / B_boot,
      stringsAsFactors = FALSE
    )
  }
  ans <- parallel_lapply_crossplatform(
    seq_len(n_rep), one,
    n_cores = min(n_cores, n_rep),
    exports = c(
      "one", "B_boot", "T", "K", "lambda", "n_start_fit", "n_start_boot", "seed"
    ),
    project_root = getwd(), seed = seed, progress = TRUE,
    progress_label = sprintf("Bootstrap coverage study: %d outer replicates", n_rep)
  )
  do.call(rbind, ans)
}

# Coordinate search used as a sensitivity analysis for the common-lambda
# restriction.  It is not presented as a replacement for REML/QREML; its purpose
# is to show whether allowing different state-specific penalties materially
# changes the fitted model in the study at hand.
run_state_specific_lambda_sensitivity <- function(sim,
                                                   common_lambda = 0.1,
                                                   candidate_grid = c(0.01, 0.05, 0.1, 0.5, 1),
                                                   n_cycles = 2L,
                                                   n_start = 4L,
                                                   max_iter = 150L,
                                                   seed = 9100L) {
  K <- sim$K_true
  lambda_vec <- rep(common_lambda, K)
  trace <- list()
  row_id <- 0L

  for (cycle in seq_len(n_cycles)) {
    changed <- FALSE
    for (k in seq_len(K)) {
      fits <- lapply(seq_along(candidate_grid), function(j) {
        lam <- lambda_vec
        lam[k] <- candidate_grid[j]
        fit <- tryCatch(
          fit_beta_gam_hmm(
            sim$y, sim$basis_spec,
            K = K, lambda = lam,
            n_start = n_start, max_iter = max_iter, optimizer_maxit = 80L,
            seed = seed + cycle * 10000L + k * 1000L + j * 10L
          ),
          error = function(e) NULL
        )
        if (is.null(fit)) return(NULL)
        row_id <<- row_id + 1L
        trace[[row_id]] <<- data.frame(
          cycle = cycle, state_updated = k, candidate_lambda = candidate_grid[j],
          lambda_vector = format_lambda(lam),
          aic_edf = fit$criteria$aic_edf, bic_edf = fit$criteria$bic_edf,
          converged = fit$converged, stringsAsFactors = FALSE
        )
        list(lambda = lam, fit = fit)
      })
      fits <- Filter(Negate(is.null), fits)
      if (!length(fits)) next
      scores <- vapply(fits, function(x) x$fit$criteria$aic_edf, numeric(1))
      best <- fits[[which.min(scores)]]
      if (!isTRUE(all.equal(lambda_vec[k], best$lambda[k]))) changed <- TRUE
      lambda_vec <- best$lambda
    }
    if (!changed) break
  }

  final_fit <- fit_beta_gam_hmm(
    sim$y, sim$basis_spec,
    K = K, lambda = lambda_vec,
    n_start = max(n_start, 6L), max_iter = max_iter, optimizer_maxit = 100L,
    seed = seed + 999999L
  )
  common_fit <- fit_beta_gam_hmm(
    sim$y, sim$basis_spec,
    K = K, lambda = common_lambda,
    n_start = max(n_start, 6L), max_iter = max_iter, optimizer_maxit = 100L,
    seed = seed + 888888L
  )

  list(
    selected_lambda_by_state = lambda_vec,
    trace = merge_data_frames(trace),
    comparison = data.frame(
      model = c("common_lambda", "state_specific_lambda"),
      aic_edf = c(common_fit$criteria$aic_edf, final_fit$criteria$aic_edf),
      bic_edf = c(common_fit$criteria$bic_edf, final_fit$criteria$bic_edf),
      edf_total = c(sum(common_fit$edf_spline), sum(final_fit$edf_spline)),
      stringsAsFactors = FALSE
    ),
    common_fit = common_fit,
    state_specific_fit = final_fit
  )
}
