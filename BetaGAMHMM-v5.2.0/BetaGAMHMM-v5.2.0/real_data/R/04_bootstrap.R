# Parametric bootstrap for the selected year-level Beta-GAM HMM.
#
# The reported intervals are conditional on selected K and lambda. The official
# interval convention is pre-specified in config/real_data_config.R. Percentile,
# basic, and bias-corrected intervals are all computed from the same costly refits
# so interval sensitivity never requires repeating the bootstrap.

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
                                          seed = 20262001L,
                                          interval_methods = list(phi = "bc", transition = "percentile",
                                                                  initial = "percentile", curve = "percentile"),
                                          retry_failed = TRUE,
                                          retry_n_start = 8L) {
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

  success_first <- !vapply(res, is.null, logical(1))
  retry_idx <- integer(0)
  retry_recovered <- integer(0)

  # Retry only failed datasets once. The simulated dataset and replicate seed
  # are unchanged; only the random-start budget is modestly increased.
  # Successful first-pass refits are never recomputed.
  if (isTRUE(retry_failed) && any(!success_first)) {
    retry_idx <- which(!success_first)
    retry_n_start <- max(as.integer(retry_n_start), as.integer(n_start) + 1L)
    message(
      "Retrying ", length(retry_idx),
      " failed conditional-bootstrap refit(s) once with ", retry_n_start,
      " starts; successful first-pass refits are retained."
    )
    retry_jobs <- lapply(retry_idx, function(b) list(
      fit = fit,
      b = b,
      seed = seed,
      n_start = retry_n_start,
      max_iter = max_iter,
      tol = tol,
      phi_lower = phi_lower,
      phi_upper = phi_upper,
      optimizer_maxit = optimizer_maxit,
      precision_saturation_fraction = precision_saturation_fraction
    ))
    retry_res <- if (n_cores > 1L) {
      parallel_lapply(
        retry_jobs,
        beta_bootstrap_worker,
        n_cores = n_cores,
        project_root = getwd(),
        seed = seed + 910000L
      )
    } else {
      lapply(retry_jobs, beta_bootstrap_worker)
    }
    for (j in seq_along(retry_idx)) {
      if (!is.null(retry_res[[j]])) {
        res[[retry_idx[j]]] <- retry_res[[j]]
        retry_recovered <- c(retry_recovered, retry_idx[j])
      }
    }
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
  phi_ci_all <- setNames(lapply(methods, function(m) {
    t(vapply(seq_len(K), function(k) bootstrap_interval(fit$phi[k], phi_boot[, k], m, level), numeric(2)))
  }), methods)
  pi_ci_all <- setNames(lapply(methods, function(m) {
    t(vapply(seq_len(K), function(k) bootstrap_interval(fit$pi[k], pi_boot[, k], m, level), numeric(2)))
  }), methods)
  A_ci_all <- setNames(lapply(methods, function(m) {
    lo <- hi <- matrix(NA_real_, K, K)
    for (i in seq_len(K)) for (j in seq_len(K)) {
      ci <- bootstrap_interval(fit$A[i, j], A_boot[, i, j], m, level)
      lo[i, j] <- ci[1L]; hi[i, j] <- ci[2L]
    }
    list(low = lo, high = hi)
  }), methods)
  mu_ci_all <- setNames(lapply(methods, function(m) {
    lo <- hi <- matrix(NA_real_, nrow = length(age_grid), ncol = K)
    for (j in seq_along(age_grid)) for (k in seq_len(K)) {
      ci <- bootstrap_interval(fit$mu_age[j, k], mu_boot[, j, k], m, level)
      lo[j, k] <- ci[1L]; hi[j, k] <- ci[2L]
    }
    list(low = lo, high = hi)
  }), methods)

  normalize_method <- function(x, default) {
    x <- tolower(as.character(x %||% default)[1L])
    if (!x %in% methods) stop("Unknown bootstrap interval method: ", x)
    x
  }
  interval_methods <- list(
    phi = normalize_method(interval_methods$phi, "bc"),
    transition = normalize_method(interval_methods$transition, "percentile"),
    initial = normalize_method(interval_methods$initial, "percentile"),
    curve = normalize_method(interval_methods$curve, "percentile")
  )
  phi_ci <- phi_ci_all[[interval_methods$phi]]
  pi_ci <- pi_ci_all[[interval_methods$initial]]
  A_ci_low <- A_ci_all[[interval_methods$transition]]$low
  A_ci_high <- A_ci_all[[interval_methods$transition]]$high
  mu_low <- mu_ci_all[[interval_methods$curve]]$low
  mu_high <- mu_ci_all[[interval_methods$curve]]$high

  list(
    B_requested = B,
    B_success_first_pass = sum(success_first),
    B_retry_attempted = length(retry_idx),
    B_retry_recovered = length(retry_recovered),
    B_success = length(res_ok),
    success_rate = length(res_ok) / B,
    retry_n_start = if (length(retry_idx)) retry_n_start else NA_integer_,
    success_first_pass = success_first,
    success_final = success,
    retry_indices = retry_idx,
    retry_recovered_indices = retry_recovered,
    level = level,
    interval_methods = interval_methods,
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
    phi_ci_all = phi_ci_all,
    A_ci_all = A_ci_all,
    pi_ci_all = pi_ci_all,
    mu_ci_all = mu_ci_all,
    age = age_grid
  )
}
