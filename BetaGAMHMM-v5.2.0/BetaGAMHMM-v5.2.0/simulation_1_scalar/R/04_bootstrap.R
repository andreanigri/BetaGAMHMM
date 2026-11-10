# Parametric bootstrap utilities.
#
# `bootstrap_beta_gam_hmm()` is conditional on the selected model structure
# (K and lambda).  This is intentional and is now explicit in both the code and
# returned metadata.  `bootstrap_selection_stability()` is provided separately
# to quantify uncertainty induced by re-selecting K and lambda.

simulate_from_fit <- function(fit, basis_spec = fit$basis_spec, seed = NULL) {
  simulate_from_parameters(
    x = basis_spec$x,
    basis_spec = basis_spec,
    beta_list = fit$beta,
    phi = fit$phi,
    A = fit$A,
    pi = fit$pi,
    seed = seed
  )
}

bootstrap_beta_gam_hmm <- function(fit,
                                   B_boot = 200L,
                                   level = 0.95,
                                   n_start = 4L,
                                   max_iter = 150L,
                                   tol = 1e-6,
                                   optimizer_maxit = 80L,
                                   n_cores = 1L,
                                   seed = 1L,
                                   align_to_reference = TRUE,
                                   x_grid = NULL,
                                   diagnostic = list(
                                     phi_max = 500,
                                     eps = 1e-6,
                                     s_thresh = 2L,
                                     delta_abs = 50,
                                     delta_sum = 100
                                   )) {
  K <- fit$K
  basis_spec <- fit$basis_spec
  if (identical(basis_spec$type, "additive_bs") && is.null(x_grid)) {
    stop("For an additive multi-covariate basis, provide `x_grid` as a data frame.")
  }
  if (is.null(x_grid)) x_grid <- seq(min(basis_spec$x), max(basis_spec$x), length.out = 250L)
  B_grid <- predict_basis(x_grid, basis_spec)
  alpha <- 1 - level
  probs <- c(alpha / 2, 1 - alpha / 2)

  run_one <- function(b) {
    sim_b <- simulate_from_fit(fit, basis_spec = basis_spec, seed = seed + b)
    fit_b <- tryCatch(
      fit_beta_gam_hmm(
        y = sim_b$y,
        basis_spec = basis_spec,
        K = K,
        lambda = fit$lambda,
        n_start = n_start,
        n_start_cores = 1L,
        max_iter = max_iter,
        tol = tol,
        phi_lower = max(1e-3, min(fit$phi) / 10),
        phi_upper = max(500, 2 * max(fit$phi)),
        optimizer_maxit = optimizer_maxit,
        boundary_handling = "error",
        store_all_starts = FALSE,
        verbose = FALSE,
        diagnostic = diagnostic,
        seed = seed + 10000L + b
      ),
      error = function(e) NULL
    )
    if (is.null(fit_b)) return(NULL)
    if (align_to_reference) {
      fit_b <- align_fit_to_reference(fit_b, fit, basis_spec = basis_spec, x_grid = x_grid)
    }
    mu_grid <- vapply(seq_len(K), function(k) drop(inv_logit(B_grid %*% fit_b$beta[[k]])), numeric(nrow(B_grid)))
    if (K == 1L) mu_grid <- matrix(mu_grid, ncol = 1L)
    list(
      pi = fit_b$pi,
      phi = fit_b$phi,
      A = fit_b$A,
      mu_grid = mu_grid,
      converged = fit_b$converged,
      monotonicity_ok = fit_b$monotonicity_ok %||% TRUE,
      degenerate = fit_b$diagnostic$flagged,
      loglik_obs = fit_b$loglik_obs
    )
  }

  idx <- seq_len(B_boot)
  raw <- parallel_lapply_crossplatform(
    X = idx,
    FUN = run_one,
    n_cores = min(as.integer(n_cores %||% 1L), length(idx)),
    exports = c(
      "run_one", "fit", "basis_spec", "x_grid", "B_grid", "alpha", "probs", "K", "level",
      "n_start", "max_iter", "tol", "optimizer_maxit", "seed", "align_to_reference", "diagnostic"
    ),
    project_root = getwd(),
    seed = seed,
    progress = TRUE,
    progress_label = sprintf("Conditional parametric bootstrap: %d replicates", length(idx))
  )
  raw <- Filter(Negate(is.null), raw)
  if (length(raw) == 0L) stop("All bootstrap replicates failed.")

  B_eff <- length(raw)
  pi_mat <- do.call(rbind, lapply(raw, `[[`, "pi"))
  phi_mat <- do.call(rbind, lapply(raw, `[[`, "phi"))
  A_arr <- array(NA_real_, dim = c(B_eff, K, K))
  mu_arr <- array(NA_real_, dim = c(B_eff, nrow(B_grid), K))
  meta_df <- data.frame(
    rep = seq_len(B_eff),
    converged = vapply(raw, `[[`, logical(1), "converged"),
    monotonicity_ok = vapply(raw, `[[`, logical(1), "monotonicity_ok"),
    degenerate = vapply(raw, `[[`, logical(1), "degenerate"),
    loglik_obs = vapply(raw, `[[`, numeric(1), "loglik_obs")
  )

  for (b in seq_len(B_eff)) {
    A_arr[b, , ] <- raw[[b]]$A
    mu_arr[b, , ] <- raw[[b]]$mu_grid
  }

  pi_ci <- t(apply(pi_mat, 2L, safe_quantile, probs = probs))
  phi_ci <- t(apply(phi_mat, 2L, safe_quantile, probs = probs))
  A_ci <- array(NA_real_, dim = c(K, K, 2L))
  for (i in seq_len(K)) {
    for (j in seq_len(K)) A_ci[i, j, ] <- safe_quantile(A_arr[, i, j], probs = probs)
  }

  mu_est <- vapply(seq_len(K), function(k) drop(inv_logit(B_grid %*% fit$beta[[k]])), numeric(nrow(B_grid)))
  if (K == 1L) mu_est <- matrix(mu_est, ncol = 1L)
  mu_lower <- apply(mu_arr, c(2, 3), function(x) safe_quantile(x, probs = probs)[1L])
  mu_upper <- apply(mu_arr, c(2, 3), function(x) safe_quantile(x, probs = probs)[2L])

  mu_curve_ci <- do.call(rbind, lapply(seq_len(K), function(k) data.frame(
    grid_index = seq_len(nrow(B_grid)),
    x = if (is.data.frame(x_grid)) NA_real_ else as.numeric(x_grid),
    state = paste0("State ", k),
    estimate = mu_est[, k],
    lower = mu_lower[, k],
    upper = mu_upper[, k],
    stringsAsFactors = FALSE
  )))

  list(
    bootstrap_type = "conditional_parametric",
    conditional_on = list(K = K, lambda = fit$lambda),
    B_requested = B_boot,
    B_eff = B_eff,
    level = level,
    probs = probs,
    x_grid = x_grid,
    pi_ci = pi_ci,
    phi_ci = phi_ci,
    A_ci = A_ci,
    mu_curve_ci = mu_curve_ci,
    pi_samples = pi_mat,
    phi_samples = phi_mat,
    A_samples = A_arr,
    mu_samples = mu_arr,
    curve_draws = mu_arr,
    bootstrap_meta = meta_df,
    reference_fit = fit
  )
}

# Bootstrap selection stability: K and lambda are re-selected on every bootstrap
# data set.  This is more expensive than conditional inference, so it is intended
# as a sensitivity analysis rather than as the default CI engine.
bootstrap_selection_stability <- function(fit,
                                          K_grid = 1:5,
                                          lambda_grid = c(0.01, 0.05, 0.1, 0.5, 1),
                                          B_boot = 50L,
                                          n_start = 3L,
                                          max_iter = 120L,
                                          optimizer_maxit = 70L,
                                          n_cores = 1L,
                                          selection_cfg = list(),
                                          diagnostic = list(
                                            phi_max = 500,
                                            eps = 1e-6,
                                            s_thresh = 2L,
                                            delta_abs = 50,
                                            delta_sum = 100
                                          ),
                                          seed = 9001L) {
  basis_spec <- fit$basis_spec
  run_one <- function(b) {
    sim_b <- simulate_from_fit(fit, basis_spec = basis_spec, seed = seed + b)
    grid <- tryCatch(
      run_grid_search_hmm(
        y = sim_b$y,
        basis_spec = basis_spec,
        K_grid = K_grid,
        lambda_grid = lambda_grid,
        n_start = n_start,
        max_iter = max_iter,
        optimizer_maxit = optimizer_maxit,
        phi_lower = 1e-3,
        phi_upper = max(500, 2 * max(fit$phi)),
        diagnostic = diagnostic,
        n_cores = 1L,
        n_start_cores = 1L,
        progress = FALSE,
        seed = seed + 10000L + b
      ),
      error = function(e) NULL
    )
    if (is.null(grid)) return(NULL)
    sel <- tryCatch(
      select_best_model_hierarchical(
        grid,
        criterion_primary = selection_cfg$criterion_primary %||% "bic_edf",
        criterion_secondary = selection_cfg$criterion_secondary %||% "icl_edf",
        lambda_criterion_primary = selection_cfg$lambda_criterion_primary %||% "aic_edf",
        lambda_criterion_secondary = selection_cfg$lambda_criterion_secondary %||% "bic_edf",
        exclude_degenerate = selection_cfg$exclude_degenerate %||% TRUE,
        enforce_occupancy = selection_cfg$enforce_occupancy %||% FALSE,
        min_state_share = selection_cfg$min_state_share %||% 0,
        min_state_mass = selection_cfg$min_state_mass %||% 0,
        equivalence_margin = selection_cfg$equivalence_margin %||% 0,
        prefer_richer_K_within_margin = selection_cfg$prefer_richer_K_within_margin %||% FALSE
      ),
      error = function(e) NULL
    )
    if (is.null(sel)) return(NULL)
    data.frame(
      rep = b,
      selected_K = sel$winner_row$K[1L],
      selected_lambda = sel$winner_row$lambda[1L],
      same_K = sel$winner_row$K[1L] == fit$K,
      stringsAsFactors = FALSE
    )
  }

  raw <- parallel_lapply_crossplatform(
    seq_len(B_boot), run_one,
    n_cores = min(n_cores, B_boot),
    exports = c(
      "run_one", "fit", "basis_spec", "K_grid", "lambda_grid", "n_start", "max_iter",
      "optimizer_maxit", "selection_cfg", "diagnostic", "seed"
    ),
    project_root = getwd(), seed = seed, progress = TRUE,
    progress_label = sprintf("Selection bootstrap: %d replicates", B_boot)
  )
  raw <- Filter(Negate(is.null), raw)
  tab <- if (length(raw)) do.call(rbind, raw) else data.frame()
  list(
    bootstrap_type = "model_selection_stability",
    B_requested = B_boot,
    B_eff = nrow(tab),
    results = tab,
    same_K_rate = if (nrow(tab)) mean(tab$same_K) else NA_real_
  )
}
