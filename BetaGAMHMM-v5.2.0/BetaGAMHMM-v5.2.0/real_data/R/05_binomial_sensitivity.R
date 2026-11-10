# Denominator-aware binomial-kernel sensitivity model.
#
# The hidden-state process is identical to the primary Beta model: one latent
# state per calendar year. The SAME HMD estimated death counts are used.
#
# When Female/Male are integer counts, the emission below is the ordinary
# Binomial log-likelihood. HMD Complete Data Series may contain fractional
# estimated death counts after standard HMD adjustments. Rather than rounding
# those estimates, we retain them exactly and use the gamma-extended binomial
# log-kernel
#
#   lgamma(N+1)-lgamma(F+1)-lgamma(M+1) + F log(p) + M log(1-p),
#
# where N=F+M. The combinatorial term is data-only and therefore does not affect
# state posterior probabilities, parameter estimation, or within-family model
# comparisons. This analysis is explicitly a denominator-aware robustness
# sensitivity, not a claim that fractional HMD estimates are literal Binomial
# observations.

binomial_kernel_cell_loglik <- function(F, M, prob) {
  prob <- clip_unit(prob, 1e-12)
  N <- F + M
  lgamma(N + 1) - lgamma(F + 1) - lgamma(M + 1) + F * log(prob) + M * log1p(-prob)
}

binomial_log_emission_year <- function(F, M, B_age, beta_list) {
  Tn <- nrow(F)
  J <- ncol(F)
  K <- length(beta_list)
  mu_age <- vapply(seq_len(K), function(k) inv_logit(drop(B_age %*% beta_list[[k]])), numeric(J))
  if (K == 1L) mu_age <- matrix(mu_age, ncol = 1L)

  const <- lgamma(F + M + 1) - lgamma(F + 1) - lgamma(M + 1)
  out <- matrix(0, Tn, K)
  for (k in seq_len(K)) {
    mu <- clip_unit(mu_age[, k], 1e-12)
    ll <- const + sweep(F, 2L, log(mu), FUN = "*") + sweep(M, 2L, log1p(-mu), FUN = "*")
    out[, k] <- rowSums(ll)
  }
  list(log_emission = out, mu_age = mu_age)
}

initialize_binomial_year_hmm <- function(F, M, B_age, K, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Tn <- nrow(F)
  N <- F + M
  yearly_ratio <- rowSums(F) / rowSums(N)

  if (K == 1L) {
    z0 <- rep.int(1L, Tn)
  } else {
    s <- stats::qlogis(clip_unit(yearly_ratio, 1e-8))
    s <- s + stats::rnorm(Tn, 0, max(stats::sd(s), 1e-4) * 0.05)
    z0 <- as.integer(cut(
      rank(s, ties.method = "random"),
      breaks = seq(0, Tn, length.out = K + 1L),
      include.lowest = TRUE,
      labels = FALSE
    ))
  }

  Mcoef <- ncol(B_age)
  beta_list <- vector("list", K)
  for (k in seq_len(K)) {
    idx <- which(z0 == k)
    if (!length(idx)) idx <- seq_len(Tn)
    f_age <- colSums(F[idx, , drop = FALSE])
    n_age <- colSums(N[idx, , drop = FALSE])
    p_age <- (f_age + 0.5) / (n_age + 1)
    beta_list[[k]] <- tryCatch(
      as.numeric(qr.solve(B_age, stats::qlogis(clip_unit(p_age, 1e-8)), tol = 1e-10)),
      error = function(e) stats::rnorm(Mcoef, 0, 0.1)
    )
  }

  if (K == 1L) {
    A <- matrix(1, 1, 1)
    pi <- 1
  } else {
    A <- matrix(stats::runif(K * K, 0.01, 0.06), K, K)
    diag(A) <- diag(A) + stats::runif(K, 0.80, 0.95)
    A <- normalize_rows(A)
    pi <- stationary_distribution(A)
  }
  list(beta = beta_list, A = A, pi = pi)
}

binomial_state_objective <- function(beta, F, M, B_age, gamma_k, lambda, penalty) {
  mu <- clip_unit(inv_logit(drop(B_age %*% beta)), 1e-12)
  const <- lgamma(F + M + 1) - lgamma(F + 1) - lgamma(M + 1)
  ll <- const + sweep(F, 2L, log(mu), FUN = "*") + sweep(M, 2L, log1p(-mu), FUN = "*")
  ll_year <- rowSums(ll)
  -sum(gamma_k * ll_year) + lambda * drop(crossprod(beta, penalty %*% beta))
}

binomial_state_gradient <- function(beta, F, M, B_age, gamma_k, lambda, penalty) {
  mu <- clip_unit(inv_logit(drop(B_age %*% beta)), 1e-12)
  N <- F + M
  score_age <- colSums((F - sweep(N, 2L, mu, FUN = "*")) * gamma_k)
  -drop(crossprod(B_age, score_age)) + 2 * lambda * drop(penalty %*% beta)
}

binomial_penalized_observed <- function(F, M, B_age, beta_list, pi, A, lambda, penalty) {
  em <- binomial_log_emission_year(F, M, B_age, beta_list)
  fb <- forward_backward_hmm(em$log_emission, pi, A)
  pen <- sum(vapply(seq_along(beta_list), function(k) {
    lambda[k] * drop(crossprod(beta_list[[k]], penalty %*% beta_list[[k]]))
  }, numeric(1)))
  list(value = fb$loglik - pen, fb = fb, em = em, penalty = pen)
}

interpolate_binomial_parameters <- function(old, new, step) {
  K <- length(old$beta)
  list(
    beta = lapply(seq_len(K), function(k) old$beta[[k]] + step * (new$beta[[k]] - old$beta[[k]])),
    pi = normalize_vec(old$pi + step * (new$pi - old$pi)),
    A = normalize_rows(old$A + step * (new$A - old$A))
  )
}

fit_single_binomial_year_hmm <- function(panel,
                                         basis_spec,
                                         K,
                                         lambda,
                                         max_iter = 300L,
                                         tol = 1e-6,
                                         optimizer_maxit = 150L,
                                         monotone_tol = 1e-8,
                                         seed = NULL) {
  F <- panel$F
  M <- panel$M
  N <- panel$N
  B_age <- basis_spec$B
  penalty <- basis_spec$penalty
  lambda_vec <- if (length(lambda) == 1L) rep(as.numeric(lambda), K) else as.numeric(lambda)
  init <- initialize_binomial_year_hmm(F, M, B_age, K, seed = seed)
  beta_list <- init$beta
  pi <- init$pi
  A <- init$A

  current <- binomial_penalized_observed(F, M, B_age, beta_list, pi, A, lambda_vec, penalty)
  objective_prev <- current$value
  converged <- FALSE
  monotonicity_ok <- TRUE
  history <- vector("list", max_iter)

  for (iter in seq_len(max_iter)) {
    gamma <- current$fb$gamma
    xi <- current$fb$xi
    old <- list(beta = beta_list, pi = pi, A = A)

    pi_new <- normalize_vec(gamma[1L, ])
    A_new <- A
    if (panel$n_years >= 2L && K > 1L) {
      A_new <- normalize_rows(apply(xi[seq_len(panel$n_years - 1L), , , drop = FALSE], c(2, 3), sum) + 1e-12)
    }

    beta_new <- beta_list
    for (k in seq_len(K)) {
      old_q <- binomial_state_objective(beta_list[[k]], F, M, B_age, gamma[, k], lambda_vec[k], penalty)
      opt <- tryCatch(
        stats::optim(
          par = beta_list[[k]],
          fn = binomial_state_objective,
          gr = binomial_state_gradient,
          method = "BFGS",
          F = F,
          M = M,
          B_age = B_age,
          gamma_k = gamma[, k],
          lambda = lambda_vec[k],
          penalty = penalty,
          control = list(maxit = optimizer_maxit, reltol = 1e-9)
        ),
        error = function(e) NULL
      )
      if (!is.null(opt) && is.finite(opt$value) && opt$value <= old_q + 1e-7) beta_new[[k]] <- opt$par
    }

    proposal <- list(beta = beta_new, pi = pi_new, A = A_new)
    proposed <- binomial_penalized_observed(F, M, B_age, proposal$beta, proposal$pi, proposal$A, lambda_vec, penalty)
    accepted <- proposed$value >= objective_prev - monotone_tol
    step <- 1

    if (!accepted) {
      for (j in seq_len(8L)) {
        step <- 0.5^j
        damped <- interpolate_binomial_parameters(old, proposal, step)
        damped_eval <- binomial_penalized_observed(F, M, B_age, damped$beta, damped$pi, damped$A, lambda_vec, penalty)
        if (damped_eval$value >= objective_prev - monotone_tol) {
          proposal <- damped
          proposed <- damped_eval
          accepted <- TRUE
          break
        }
      }
    }
    if (!accepted) {
      monotonicity_ok <- FALSE
      break
    }

    beta_list <- proposal$beta
    pi <- proposal$pi
    A <- proposal$A
    delta <- proposed$value - objective_prev
    current <- proposed
    history[[iter]] <- data.frame(iter = iter, loglik = current$fb$loglik, objective = current$value, delta = delta, step = step)
    rel <- abs(delta) / (1 + abs(objective_prev))
    objective_prev <- current$value
    if (rel < tol) {
      converged <- TRUE
      break
    }
  }

  fit <- list(
    family = "binomial_kernel",
    K = K,
    lambda = lambda_vec,
    beta = beta_list,
    phi = NULL,
    pi = pi,
    A = A,
    gamma = current$fb$gamma,
    xi = current$fb$xi,
    alpha = current$fb$alpha,
    log_emission = current$em$log_emission,
    mu_age = current$em$mu_age,
    loglik_obs = current$fb$loglik,
    objective_pen = current$value,
    converged = converged,
    monotonicity_ok = monotonicity_ok,
    n_iter = length(Filter(Negate(is.null), history)),
    panel = panel,
    basis_spec = basis_spec,
    interpretation = "denominator-aware binomial-kernel sensitivity using original HMD estimated counts"
  )
  fit <- relabel_fit(fit, order_states_by_profile(fit))
  fit$viterbi <- viterbi_from_log_emission(fit$log_emission, fit$pi, fit$A)
  fit$occupancy <- state_occupancy(fit$gamma)
  fit
}

binomial_multistart_worker <- function(job) {
  tryCatch(
    do.call(fit_single_binomial_year_hmm, job),
    error = function(e) structure(
      list(error = conditionMessage(e)),
      class = "fit_error"
    )
  )
}

fit_binomial_year_hmm <- function(panel,
                                  basis_spec,
                                  K,
                                  lambda,
                                  n_start = 15L,
                                  n_cores = 1L,
                                  max_iter = 300L,
                                  tol = 1e-6,
                                  optimizer_maxit = 150L,
                                  seed = 20263001L) {
  started <- proc.time()[[3L]]
  seeds <- seed + seq_len(n_start)
  jobs <- lapply(seq_len(n_start), function(s) list(
    panel = panel,
    basis_spec = basis_spec,
    K = K,
    lambda = lambda,
    max_iter = max_iter,
    tol = tol,
    optimizer_maxit = optimizer_maxit,
    seed = seeds[s]
  ))

  fits <- if (n_cores > 1L) {
    parallel_lapply(
      jobs,
      binomial_multistart_worker,
      n_cores = n_cores,
      project_root = getwd(),
      seed = seed + 80000L
    )
  } else {
    lapply(jobs, binomial_multistart_worker)
  }

  good <- vapply(fits, function(x) !inherits(x, "fit_error"), logical(1))
  if (!any(good)) {
    worker_errors <- unique(vapply(fits, function(x) {
      if (inherits(x, "fit_error")) x$error else NA_character_
    }, character(1)))
    worker_errors <- worker_errors[!is.na(worker_errors) & nzchar(worker_errors)]
    detail <- if (length(worker_errors)) paste0(
      " Worker error(s): ", paste(worker_errors, collapse = " | ")
    ) else ""
    stop("All binomial-kernel multi-start fits failed for K=", K, ", lambda=", lambda, ".", detail)
  }

  fits <- fits[good]
  reliable <- vapply(fits, function(f) isTRUE(f$converged) && isTRUE(f$monotonicity_ok), logical(1))
  pool <- if (any(reliable)) which(reliable) else seq_along(fits)
  best <- fits[[pool[which.max(vapply(fits[pool], `[[`, numeric(1), "loglik_obs"))]]]

  best$edf_spline <- vapply(seq_len(K), function(k) {
    binomial_edf_state(
      B = basis_spec$B,
      gamma_k = best$gamma[, k],
      mu_age = best$mu_age[, k],
      N = panel$N,
      lambda = best$lambda[k],
      penalty = basis_spec$penalty
    )
  }, numeric(1))

  best$criteria <- information_criteria(
    loglik = best$loglik_obs,
    gamma = best$gamma,
    K = K,
    edf_spline = best$edf_spline,
    n_years = panel$n_years,
    n_cells = panel$n_cells,
    has_phi = FALSE
  )
  best$start_diagnostics <- data.frame(
    n_start = n_start,
    n_successful = length(fits),
    n_converged = sum(vapply(fits, function(f) isTRUE(f$converged), logical(1))),
    n_monotone = sum(vapply(fits, function(f) isTRUE(f$monotonicity_ok), logical(1))),
    best_loglik = max(vapply(fits, `[[`, numeric(1), "loglik_obs")),
    median_loglik = stats::median(vapply(fits, `[[`, numeric(1), "loglik_obs")),
    stringsAsFactors = FALSE
  )
  best$timing_seconds <- proc.time()[[3L]] - started
  best
}

run_binomial_grid <- function(panel,
                              basis_spec,
                              K_grid,
                              lambda_grid,
                              n_start = 15L,
                              n_cores = 1L,
                              max_iter = 300L,
                              tol = 1e-6,
                              optimizer_maxit = 150L,
                              seed = 20264001L,
                              progress = TRUE) {
  grid <- expand.grid(K = K_grid, lambda = lambda_grid, KEEP.OUT.ATTRS = FALSE)
  results <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    if (progress) message(sprintf("Binomial-kernel grid %d/%d: K=%d, lambda=%g", i, nrow(grid), grid$K[i], grid$lambda[i]))
    results[[i]] <- tryCatch(
      fit_binomial_year_hmm(
        panel = panel,
        basis_spec = basis_spec,
        K = grid$K[i],
        lambda = grid$lambda[i],
        n_start = n_start,
        n_cores = n_cores,
        max_iter = max_iter,
        tol = tol,
        optimizer_maxit = optimizer_maxit,
        seed = seed + 10000L * grid$K[i] + 100L * i
      ),
      error = function(e) structure(list(K = grid$K[i], lambda = grid$lambda[i], error = conditionMessage(e)), class = "grid_error")
    )
  }
  results
}

summarize_binomial_grid <- function(results) {
  rows <- lapply(results, function(fit) {
    if (inherits(fit, "grid_error")) return(data.frame(
      K = fit$K, lambda = fit$lambda, converged = FALSE, monotonicity_ok = FALSE,
      min_state_share = NA_real_, loglik_kernel = NA_real_, AIC_EDF = NA_real_,
      BIC_EDF = NA_real_, ICL_EDF = NA_real_, npar_edf = NA_real_,
      edf_spline = NA_real_, timing_seconds = NA_real_, error = fit$error,
      stringsAsFactors = FALSE
    ))
    cr <- fit$criteria
    data.frame(K = fit$K, lambda = fit$lambda[1L], converged = fit$converged,
      monotonicity_ok = fit$monotonicity_ok, min_state_share = fit$occupancy$min_share,
      loglik_kernel = fit$loglik_obs, AIC_EDF = cr$aic_edf, BIC_EDF = cr$bic_edf,
      ICL_EDF = cr$icl_edf, npar_edf = cr$npar_edf,
      edf_spline = sum(fit$edf_spline), timing_seconds = fit$timing_seconds,
      error = NA_character_, stringsAsFactors = FALSE)
  })
  ans <- do.call(rbind, rows); ans[order(ans$K, ans$lambda), , drop = FALSE]
}

select_binomial_model <- function(results) {
  tab <- summarize_binomial_grid(results)
  tab$admissible <- with(tab, is.na(error) & converged & monotonicity_ok & is.finite(AIC_EDF) & is.finite(BIC_EDF))
  best_by_K <- do.call(rbind, lapply(split(tab, tab$K), function(dd) {
    ok <- dd[dd$admissible, , drop = FALSE]; if (!nrow(ok)) return(NULL)
    ok[order(ok$AIC_EDF, ok$lambda), , drop = FALSE][1L, ]
  }))
  if (is.null(best_by_K) || !nrow(best_by_K)) stop("No numerically valid denominator-aware model found.")
  rownames(best_by_K) <- NULL
  sel <- best_by_K[order(best_by_K$BIC_EDF, best_by_K$K), , drop = FALSE][1L, ]
  idx <- which(vapply(results, function(f) !inherits(f, "grid_error") && f$K == sel$K &&
                       abs(f$lambda[1L] - sel$lambda) < 1e-12, logical(1)))[1L]
  list(fit = results[[idx]], grid_summary = tab, best_by_K = best_by_K,
       selected_row = sel, criterion_label = "BIC", criterion_column = "BIC_EDF")
}

compare_beta_binomial <- function(beta_selection, binomial_selection) {
  beta_fit <- beta_selection$fit; bin_fit <- binomial_selection$fit
  aligned <- if (beta_fit$K == bin_fit$K) align_fit_to_reference(bin_fit, beta_fit) else bin_fit
  data.frame(beta_K = beta_fit$K, beta_lambda = beta_fit$lambda[1L],
    denominator_K = bin_fit$K, denominator_lambda = bin_fit$lambda[1L],
    viterbi_ARI = adjusted_rand_index(beta_fit$viterbi, bin_fit$viterbi),
    exact_yearly_state_agreement = if (beta_fit$K == bin_fit$K) mean(beta_fit$viterbi == aligned$viterbi) else NA_real_,
    curve_RMSE_if_same_K = if (beta_fit$K == bin_fit$K) sqrt(mean((beta_fit$mu_age - aligned$mu_age)^2)) else NA_real_,
    stringsAsFactors = FALSE)
}


run_binomial_grid_with_boundary_extension <- function(panel,
                                                       basis_spec,
                                                       K_grid,
                                                       lambda_grid,
                                                       K_auto_max = max(K_grid),
                                                       K_auto_step = 2L,
                                                       n_start = 15L,
                                                       n_cores = 1L,
                                                       max_iter = 300L,
                                                       tol = 1e-6,
                                                       optimizer_maxit = 150L,
                                                       seed = 1L,
                                                       progress = TRUE) {
  tested <- sort(unique(as.integer(K_grid)))
  res <- run_binomial_grid(
    panel = panel, basis_spec = basis_spec,
    K_grid = tested, lambda_grid = lambda_grid,
    n_start = n_start, n_cores = n_cores,
    max_iter = max_iter, tol = tol,
    optimizer_maxit = optimizer_maxit,
    seed = seed, progress = progress
  )
  repeat {
    sel <- select_binomial_model(res)
    current_max <- max(tested)
    if (sel$fit$K < current_max || current_max >= K_auto_max) break
    next_max <- min(K_auto_max, current_max + as.integer(K_auto_step))
    new_K <- seq.int(current_max + 1L, next_max)
    message("      Binomial-kernel BIC minimum is at upper boundary K=", current_max,
            "; extending search to K=", next_max, ".")
    add <- run_binomial_grid(
      panel = panel, basis_spec = basis_spec,
      K_grid = new_K, lambda_grid = lambda_grid,
      n_start = n_start, n_cores = n_cores,
      max_iter = max_iter, tol = tol,
      optimizer_maxit = optimizer_maxit,
      seed = seed + 1000000L + 10000L * current_max,
      progress = progress
    )
    res <- append_grid_results(res, add)
    tested <- sort(unique(c(tested, new_K)))
  }
  sel <- select_binomial_model(res)
  attr(res, "K_tested") <- tested
  attr(res, "boundary_hit") <- identical(sel$fit$K, max(tested))
  res
}
