# Real-data robustness analyses used in the revision.

run_basis_selection_sensitivity <- function(panel, primary_grid, K_grid, lambda_grid,
                                            inner_knots_grid = c(4L, 6L, 8L),
                                            primary_inner_knots = 6L, degree = 3L,
                                            n_start = 6L, n_cores = 1L,
                                            max_iter = 300L, tol = 1e-6,
                                            phi_lower = 1, phi_upper = 2000,
                                            optimizer_maxit = 150L,
                                            precision_saturation_fraction = 0.98,
                                            seed = 1L) {
  rows <- list(); fits <- list(); grids <- list(); idx <- 0L
  primary_sel <- select_beta_model(primary_grid)
  for (m in as.integer(inner_knots_grid)) {
    if (m == as.integer(primary_inner_knots)) {
      grid <- primary_grid
    } else {
      bs <- build_age_basis(panel$ages, degree = degree, n_inner_knots = m)
      grid <- run_beta_grid(
        panel = panel, basis_spec = bs, K_grid = K_grid, lambda_grid = lambda_grid,
        n_start = n_start, n_cores = n_cores, max_iter = max_iter, tol = tol,
        phi_lower = phi_lower, phi_upper = phi_upper,
        optimizer_maxit = optimizer_maxit,
        precision_saturation_fraction = precision_saturation_fraction,
        seed = seed + 100000L * m, progress = TRUE
      )
    }
    sel <- select_beta_model(grid)
    grids[[as.character(m)]] <- grid
    fits[[as.character(m)]] <- sel$fit
    aligned <- if (sel$fit$K == primary_sel$fit$K) align_fit_to_reference(sel$fit, primary_sel$fit) else sel$fit
    idx <- idx + 1L
    rows[[idx]] <- data.frame(
      n_inner_knots = m, selected_K = sel$fit$K,
      selected_lambda = sel$fit$lambda[1L], BIC_EDF = sel$selected_row$BIC_EDF,
      min_state_share = sel$fit$occupancy$min_share,
      precision_flag = sel$fit$precision_diagnostic$flagged,
      same_K_as_M6 = sel$fit$K == primary_sel$fit$K,
      viterbi_ARI_vs_M6 = adjusted_rand_index(sel$fit$viterbi, primary_sel$fit$viterbi),
      exact_viterbi_agreement = if (sel$fit$K == primary_sel$fit$K) mean(aligned$viterbi == primary_sel$fit$viterbi) else NA_real_,
      curve_RMSE_vs_M6 = if (sel$fit$K == primary_sel$fit$K) sqrt(mean((aligned$mu_age - primary_sel$fit$mu_age)^2)) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  list(summary = do.call(rbind, rows), fits = fits, grids = grids)
}

summarize_state_persistence <- function(fit, years, boundary_tol = 1e-6) {
  K <- fit$K
  if (K == 1L) return(list(
    transition = data.frame(
      state = 1L, self_transition = 1,
      implied_mean_sojourn = NA_real_,
      raw_implied_mean_sojourn = Inf,
      boundary_transition_estimate = TRUE,
      terminally_right_censored = TRUE,
      stringsAsFactors = FALSE
    ),
    decoded_runs = data.frame(
      state = 1L, start_year = min(years), end_year = max(years),
      duration_years = length(years), right_censored_at_series_end = TRUE,
      stringsAsFactors = FALSE
    )
  ))

  akk <- diag(fit$A)
  raw_sojourn <- 1 / pmax(1 - akk, .Machine$double.eps)
  boundary <- akk >= (1 - boundary_tol)
  terminal_state <- tail(fit$viterbi, 1L)

  trans <- data.frame(
    state = seq_len(K),
    self_transition = akk,
    implied_mean_sojourn = ifelse(boundary, NA_real_, raw_sojourn),
    raw_implied_mean_sojourn = raw_sojourn,
    boundary_transition_estimate = boundary,
    terminally_right_censored = seq_len(K) == terminal_state,
    stringsAsFactors = FALSE
  )

  rr <- rle(fit$viterbi)
  ends <- cumsum(rr$lengths)
  starts <- c(1L, head(ends, -1L) + 1L)
  runs <- data.frame(
    state = rr$values,
    start_year = years[starts],
    end_year = years[ends],
    duration_years = rr$lengths,
    right_censored_at_series_end = ends == length(years),
    stringsAsFactors = FALSE
  )
  list(transition = trans, decoded_runs = runs)
}

phi_bound_refit_worker <- function(job) {
  fit <- tryCatch(
    fit_beta_year_hmm(
      panel = job$panel, basis_spec = job$basis_spec,
      K = job$K, lambda = job$lambda,
      n_start = job$n_start, n_cores = 1L,
      max_iter = job$max_iter, tol = job$tol,
      phi_lower = job$phi_lower, phi_upper = job$phi_upper,
      optimizer_maxit = job$optimizer_maxit,
      precision_saturation_fraction = job$precision_saturation_fraction,
      seed = job$seed,
      warm_start_fit = job$warm_start_fit
    ), error = function(e) NULL
  )
  if (is.null(fit)) return(data.frame(
    setting = job$setting, phi_lower = job$phi_lower, phi_upper = job$phi_upper,
    K = job$K, lambda = job$lambda, converged = FALSE, monotonicity_ok = FALSE,
    BIC_EDF = NA_real_, AIC_EDF = NA_real_, loglik = NA_real_, max_phi = NA_real_,
    precision_flag = NA, warm_start_used = !is.null(job$warm_start_fit),
    error = TRUE, stringsAsFactors = FALSE
  ))
  data.frame(
    setting = job$setting, phi_lower = job$phi_lower, phi_upper = job$phi_upper,
    K = fit$K, lambda = fit$lambda[1L], converged = isTRUE(fit$converged),
    monotonicity_ok = isTRUE(fit$monotonicity_ok),
    BIC_EDF = fit$criteria$bic_edf, AIC_EDF = fit$criteria$aic_edf,
    loglik = fit$loglik_obs, max_phi = max(fit$phi),
    precision_flag = isTRUE(fit$precision_diagnostic$flagged),
    warm_start_used = !is.null(job$warm_start_fit),
    error = FALSE, stringsAsFactors = FALSE
  )
}

run_phi_bound_sensitivity <- function(panel, basis_spec, primary_selection,
                                      bound_sets = list(relaxed = c(lower = 0.1, upper = 5000)),
                                      n_start = 4L, n_cores = 1L,
                                      max_iter = 300L, tol = 1e-6,
                                      optimizer_maxit = 150L,
                                      precision_saturation_fraction = 0.98,
                                      seed = 1L,
                                      primary_grid = NULL) {
  base <- primary_selection$best_by_K
  if (!nrow(base)) stop("Primary best-by-K table is empty.")

  find_primary_fit <- function(K, lambda) {
    if (is.null(primary_grid) || !length(primary_grid)) return(NULL)
    idx <- which(vapply(primary_grid, function(f) {
      !inherits(f, "grid_error") &&
        as.integer(f$K) == as.integer(K) &&
        abs(as.numeric(f$lambda[1L]) - as.numeric(lambda)) < 1e-12
    }, logical(1)))
    if (!length(idx)) return(NULL)
    primary_grid[[idx[1L]]]
  }

  primary_fits <- lapply(seq_len(nrow(base)), function(i) {
    find_primary_fit(base$K[i], base$lambda[i])
  })
  if (!is.null(primary_grid) && any(vapply(primary_fits, is.null, logical(1)))) {
    miss <- which(vapply(primary_fits, is.null, logical(1)))
    stop("Could not recover primary warm-start fit(s) for K=",
         paste(base$K[miss], collapse = ", "), ".")
  }

  rows <- list(data.frame(
    setting = "primary", phi_lower = NA_real_, phi_upper = NA_real_,
    K = base$K, lambda = base$lambda, converged = TRUE, monotonicity_ok = TRUE,
    BIC_EDF = base$BIC_EDF, AIC_EDF = base$AIC_EDF,
    loglik = base$loglik, max_phi = NA_real_,
    precision_flag = base$precision_flag, warm_start_used = FALSE,
    error = FALSE, stringsAsFactors = FALSE
  ))

  jobs <- list(); jj <- 0L
  for (setting in names(bound_sets)) {
    bounds <- bound_sets[[setting]]
    if (is.null(names(bounds)) || !all(c("lower", "upper") %in% names(bounds)))
      stop("Each phi-bound setting must contain named lower and upper values.")
    for (i in seq_len(nrow(base))) {
      jj <- jj + 1L
      warm <- primary_fits[[i]]
      if (!is.null(warm) &&
          (any(warm$phi < as.numeric(bounds[["lower"]])) ||
           any(warm$phi > as.numeric(bounds[["upper"]])))) {
        warm <- NULL
      }
      jobs[[jj]] <- list(
        setting = setting,
        phi_lower = as.numeric(bounds[["lower"]]), phi_upper = as.numeric(bounds[["upper"]]),
        panel = panel, basis_spec = basis_spec,
        K = as.integer(base$K[i]), lambda = as.numeric(base$lambda[i]),
        n_start = n_start, max_iter = max_iter, tol = tol,
        optimizer_maxit = optimizer_maxit,
        precision_saturation_fraction = precision_saturation_fraction,
        warm_start_fit = warm,
        seed = seed + 100000L * match(setting, names(bound_sets)) + 1000L * as.integer(base$K[i])
      )
    }
  }

  if (length(jobs)) {
    alt <- parallel_lapply(jobs, phi_bound_refit_worker, n_cores = n_cores,
                           project_root = getwd(), seed = seed + 900000L)
    rows <- c(rows, alt)
  }
  grid <- do.call(rbind, rows)

  summary <- do.call(rbind, lapply(unique(grid$setting), function(st) {
    d <- grid[grid$setting == st & !grid$error & grid$converged & grid$monotonicity_ok &
                is.finite(grid$BIC_EDF), , drop = FALSE]
    if (!nrow(d)) return(data.frame(
      setting = st, selected_K = NA_integer_, selected_lambda = NA_real_,
      min_BIC_EDF = NA_real_, same_K_as_primary = NA,
      all_primary_warm_starts_available = NA, stringsAsFactors = FALSE
    ))
    z <- d[order(d$BIC_EDF, d$K), , drop = FALSE][1L, ]
    data.frame(
      setting = st, selected_K = z$K, selected_lambda = z$lambda,
      min_BIC_EDF = z$BIC_EDF,
      same_K_as_primary = z$K == primary_selection$fit$K,
      all_primary_warm_starts_available =
        if (st == "primary") NA else all(d$warm_start_used),
      stringsAsFactors = FALSE
    )
  }))

  list(
    summary = summary, grid = grid,
    note = paste(
      "Alternative bounds refit the primary AIC-selected lambda within each K.",
      "For each K, the corresponding primary fit is supplied as a warm start and",
      "retained as a feasible candidate whenever it lies inside the relaxed bounds;",
      "the number of expensive starts is unchanged because the warm refinement",
      "replaces one random start. The full lambda grid is not repeated."
    )
  )
}

