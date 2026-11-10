# Conditional parametric bootstrap for a fitted profile Beta-GAM HMM.
#
# K and lambda are fixed at their selected values during each conditional
# bootstrap analysis. Bias-corrected intervals are used for phi; percentile
# intervals are used for transition probabilities. For smooth mean curves the
# simulation study reports percentile, basic and bias-corrected pointwise
# intervals so their finite-sample coverage can be assessed transparently.

simulate_state_path <- function(Tn, pi, A) {
  K <- length(pi)
  z <- integer(Tn)
  z[1L] <- sample.int(K, 1L, prob = pi)
  if (Tn >= 2L) {
    for (t in 2:Tn) z[t] <- sample.int(K, 1L, prob = A[z[t - 1L], ])
  }
  z
}

simulate_beta_panel_from_fit <- function(fit, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  panel <- fit$panel
  z <- simulate_state_path(panel$n_years, fit$pi, fit$A)
  Y <- matrix(NA_real_, panel$n_years, panel$n_ages)
  for (t in seq_len(panel$n_years)) {
    k <- z[t]
    mu <- fit$mu_age[, k]
    Y[t, ] <- stats::rbeta(panel$n_ages, mu * fit$phi[k], (1 - mu) * fit$phi[k])
  }
  sim_panel <- panel
  sim_panel$Y <- Y
  sim_panel$Y_raw <- Y
  sim_panel$boundary_cells <- matrix(integer(0), ncol = 2L)
  list(panel = sim_panel, z = z)
}

beta_bootstrap_worker <- function(job) {
  fit <- job$fit
  b <- job$b
  seed <- job$seed

  sim <- simulate_beta_panel_from_fit(fit, seed = seed + b)
  fb <- tryCatch(
    fit_beta_year_hmm(
      panel = sim$panel,
      basis_spec = fit$basis_spec,
      K = fit$K,
      lambda = fit$lambda,
      n_start = job$n_start,
      n_cores = 1L,
      max_iter = job$max_iter,
      tol = job$tol,
      phi_lower = job$phi_lower,
      phi_upper = job$phi_upper,
      optimizer_maxit = job$optimizer_maxit,
      precision_saturation_fraction = job$precision_saturation_fraction,
      seed = seed + 10000L + 100L * b
    ),
    error = function(e) NULL
  )
  if (is.null(fb)) return(NULL)
  if (!isTRUE(fb$converged) || !isTRUE(fb$monotonicity_ok) || isTRUE(fb$precision_diagnostic$flagged)) return(NULL)
  fb <- align_fit_to_reference(fb, fit)
  list(
    phi = fb$phi,
    A = fb$A,
    pi = fb$pi,
    mu_age = fb$mu_age,
    converged = fb$converged,
    monotonicity_ok = fb$monotonicity_ok,
    precision_flag = fb$precision_diagnostic$flagged
  )
}

run_beta_parametric_bootstrap <- function(fit,
                                          B = 200L,
                                          n_start = 6L,
                                          n_cores = 1L,
                                          max_iter = 300L,
                                          tol = 1e-6,
                                          phi_lower = 1,
                                          phi_upper = 2000,
                                          optimizer_maxit = 150L,
                                          precision_saturation_fraction = 0.98,
                                          level = 0.95,
                                          seed = 20262001L) {
  K <- fit$K
  age_grid <- fit$panel$ages

  jobs <- lapply(seq_len(B), function(b) list(
    fit = fit,
    b = b,
    seed = seed,
    n_start = n_start,
    max_iter = max_iter,
    tol = tol,
    phi_lower = phi_lower,
    phi_upper = phi_upper,
    optimizer_maxit = optimizer_maxit,
    precision_saturation_fraction = precision_saturation_fraction
  ))

  res <- if (n_cores > 1L) {
    parallel_lapply(
      jobs,
      beta_bootstrap_worker,
      n_cores = n_cores,
      project_root = getwd(),
      seed = seed + 900000L
    )
  } else {
    lapply(jobs, beta_bootstrap_worker)
  }

  success <- !vapply(res, is.null, logical(1))
  res_ok <- res[success]
  if (!length(res_ok)) stop("All parametric bootstrap refits failed.")

  phi_boot <- do.call(rbind, lapply(res_ok, `[[`, "phi"))
  pi_boot <- do.call(rbind, lapply(res_ok, `[[`, "pi"))
  A_boot <- array(NA_real_, dim = c(length(res_ok), K, K))
  mu_boot <- array(NA_real_, dim = c(length(res_ok), length(age_grid), K))
  for (b in seq_along(res_ok)) {
    A_boot[b, , ] <- res_ok[[b]]$A
    mu_boot[b, , ] <- res_ok[[b]]$mu_age
  }

  methods <- c("percentile", "basic", "bc")
  interval_fun <- function(theta_hat, boot, method) {
    if (method == "percentile") return(percentile_interval(boot, level = level))
    if (method == "basic") return(basic_interval(theta_hat, boot, level = level))
    bc_interval(theta_hat, boot, level = level)
  }
  phi_ci_all <- setNames(lapply(methods, function(m) {
    t(vapply(seq_len(K), function(k) interval_fun(fit$phi[k], phi_boot[, k], m), numeric(2)))
  }), methods)
  pi_ci_all <- setNames(lapply(methods, function(m) {
    t(vapply(seq_len(K), function(k) interval_fun(fit$pi[k], pi_boot[, k], m), numeric(2)))
  }), methods)
  A_ci_all <- setNames(lapply(methods, function(m) {
    lo <- hi <- matrix(NA_real_, K, K)
    for (i in seq_len(K)) for (j in seq_len(K)) {
      ci <- interval_fun(fit$A[i, j], A_boot[, i, j], m)
      lo[i, j] <- ci[1L]; hi[i, j] <- ci[2L]
    }
    list(low = lo, high = hi)
  }), methods)
  mu_ci_all <- setNames(lapply(methods, function(m) {
    lo <- hi <- matrix(NA_real_, nrow = length(age_grid), ncol = K)
    for (j in seq_along(age_grid)) for (k in seq_len(K)) {
      ci <- interval_fun(fit$mu_age[j, k], mu_boot[, j, k], m)
      lo[j, k] <- ci[1L]; hi[j, k] <- ci[2L]
    }
    list(low = lo, high = hi)
  }), methods)

  # Frozen manuscript convention. Coverage validates this choice; it is not
  # selected after inspecting coverage results.
  phi_ci <- phi_ci_all$bc
  pi_ci <- pi_ci_all$percentile
  A_ci_low <- A_ci_all$percentile$low
  A_ci_high <- A_ci_all$percentile$high
  mu_low <- mu_ci_all$percentile$low
  mu_high <- mu_ci_all$percentile$high
  mu_basic_low <- mu_ci_all$basic$low
  mu_basic_high <- mu_ci_all$basic$high
  mu_bc_low <- mu_ci_all$bc$low
  mu_bc_high <- mu_ci_all$bc$high

  list(
    B_requested = B,
    B_success = length(res_ok),
    success_rate = length(res_ok) / B,
    level = level,
    phi_boot = phi_boot,
    A_boot = A_boot,
    pi_boot = pi_boot,
    mu_boot = mu_boot,
    phi_ci = phi_ci,
    A_ci_low = A_ci_low,
    A_ci_high = A_ci_high,
    pi_ci = pi_ci,
    mu_low = mu_low,
    mu_high = mu_high,
    mu_basic_low = mu_basic_low,
    mu_basic_high = mu_basic_high,
    mu_bc_low = mu_bc_low,
    mu_bc_high = mu_bc_high,
    phi_ci_all = phi_ci_all,
    A_ci_all = A_ci_all,
    pi_ci_all = pi_ci_all,
    mu_ci_all = mu_ci_all,
    age = age_grid
  )
}
