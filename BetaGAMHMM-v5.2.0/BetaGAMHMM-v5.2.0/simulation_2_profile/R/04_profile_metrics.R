# Alignment and performance metrics against known simulation truth.

align_profile_fit_to_truth <- function(fit, sim) {
  if (fit$K != sim$scenario$K) return(fit)
  K <- fit$K
  if (K == 1L) return(fit)
  cost <- matrix(0, K, K)
  for (truth_k in seq_len(K)) {
    for (fit_k in seq_len(K)) {
      curve_cost <- mean((fit$mu_age[, fit_k] - sim$mu_true[, truth_k])^2)
      scale_phi <- max(stats::sd(sim$phi_true), 1)
      phi_cost <- ((fit$phi[fit_k] - sim$phi_true[truth_k]) / scale_phi)^2
      cost[truth_k, fit_k] <- curve_cost + 0.01 * phi_cost
    }
  }
  relabel_fit(fit, min_cost_assignment_dp(cost))
}

profile_fit_metrics <- function(fit, sim) {
  if (fit$K != sim$scenario$K) {
    return(data.frame(
      correct_K = FALSE, viterbi_accuracy = NA_real_, ari = NA_real_,
      curve_rmse = NA_real_, phi_rmse = NA_real_, A_rmse = NA_real_,
      min_state_share = fit$occupancy$min_share
    ))
  }
  fit <- align_profile_fit_to_truth(fit, sim)
  data.frame(
    correct_K = TRUE,
    viterbi_accuracy = mean(fit$viterbi == sim$z),
    ari = adjusted_rand_index(fit$viterbi, sim$z),
    curve_rmse = sqrt(mean((fit$mu_age - sim$mu_true)^2)),
    phi_rmse = sqrt(mean((fit$phi - sim$phi_true)^2)),
    A_rmse = sqrt(mean((fit$A - sim$A_true)^2)),
    min_state_share = fit$occupancy$min_share
  )
}

profile_selection_rep_worker <- function(job) {
  started <- proc.time()[[3L]]
  sim <- simulate_profile_scenario(job$scenario, job$seed)
  basis <- build_age_basis(sim$panel$ages,
    degree = job$degree, n_inner_knots = job$n_inner_knots
  )
  grid <- run_profile_grid(
    panel = sim$panel, basis_spec = basis,
    K_grid = job$K_grid, lambda_grid = job$lambda_grid,
    n_start = job$n_start, n_cores = 1L,
    max_iter = job$max_iter, tol = job$tol,
    phi_lower = job$phi_lower, phi_upper = job$phi_upper,
    optimizer_maxit = job$optimizer_maxit,
    precision_saturation_fraction = job$precision_saturation_fraction,
    seed = job$seed + 100000L, progress = FALSE
  )

  gtab <- summarize_profile_grid(grid)
  grid_monotonicity_pass_rate <- mean(gtab$monotonicity_ok %in% TRUE, na.rm=TRUE)
  grid_convergence_rate <- mean(gtab$converged %in% TRUE, na.rm=TRUE)
  true_counts_txt <- paste(sim$true_state_counts, collapse = " | ")
  true_min_count <- sim$min_true_state_count
  true_min_share <- sim$min_true_state_share
  realized_states <- sim$n_states_realized
  all_visited <- sim$all_states_visited
  cmap <- profile_criterion_map()
  selections <- lapply(names(cmap), function(nm) {
    sel <- tryCatch(select_profile_model(grid, cmap[[nm]]), error = function(e) NULL)
    if (is.null(sel)) {
      return(data.frame(
        scenario = job$scenario$name, replicate = job$replicate,
        criterion = nm, K_true = job$scenario$K, selected_K = NA_integer_,
        n_states_realized = realized_states, all_states_visited = all_visited,
        min_true_state_count = true_min_count, min_true_state_share = true_min_share,
        true_state_counts = true_counts_txt,
        selected_lambda = NA_real_, exact = FALSE, underfit = NA, overfit = NA,
        viterbi_accuracy = NA_real_, ari = NA_real_, curve_rmse = NA_real_,
        phi_rmse = NA_real_, A_rmse = NA_real_, min_state_share = NA_real_,
        precision_flag = NA,
        grid_convergence_rate = grid_convergence_rate,
        grid_monotonicity_pass_rate = grid_monotonicity_pass_rate
      ))
    }
    met <- profile_fit_metrics(sel$fit, sim)
    data.frame(
      scenario = job$scenario$name, replicate = job$replicate,
      criterion = nm, K_true = job$scenario$K, selected_K = sel$fit$K,
      n_states_realized = realized_states, all_states_visited = all_visited,
      min_true_state_count = true_min_count, min_true_state_share = true_min_share,
      true_state_counts = true_counts_txt,
      selected_lambda = sel$fit$lambda[1L], exact = sel$fit$K == job$scenario$K,
      underfit = sel$fit$K < job$scenario$K, overfit = sel$fit$K > job$scenario$K,
      viterbi_accuracy = met$viterbi_accuracy, ari = met$ari,
      curve_rmse = met$curve_rmse, phi_rmse = met$phi_rmse,
      A_rmse = met$A_rmse, min_state_share = met$min_state_share,
      precision_flag = isTRUE(sel$fit$precision_diagnostic$flagged),
      grid_convergence_rate = grid_convergence_rate,
      grid_monotonicity_pass_rate = grid_monotonicity_pass_rate
    )
  })

  list(
    results = do.call(rbind, selections),
    elapsed_seconds = proc.time()[[3L]] - started
  )
}

