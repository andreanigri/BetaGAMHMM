# Core estimation routines for the Beta-GAM hidden Markov model.
#
# The implementation uses a scaled forward-backward recursion, analytical
# gradients for the Beta emission M-step, and an explicit monotonicity safeguard
# for the generalized EM iterations.  These choices substantially reduce run
# time while making the numerical behaviour easier to diagnose and reproduce.

compute_log_emission <- function(y, B, beta_list, phi) {
  K <- length(phi)
  T <- length(y)
  eta <- vapply(seq_len(K), function(k) drop(B %*% beta_list[[k]]), numeric(T))
  if (K == 1L) eta <- matrix(eta, ncol = 1L)
  mu <- inv_logit(eta)
  log_emission <- vapply(
    seq_len(K),
    function(k) log_beta_density(y, mu[, k], phi[k]),
    numeric(T)
  )
  if (K == 1L) log_emission <- matrix(log_emission, ncol = 1L)
  list(eta = eta, mu = mu, log_emission = log_emission)
}

# Scaled forward-backward algorithm.
#
# Emission log-densities are shifted by their row maxima before exponentiation.
# The resulting recursion is carried out with matrix-vector products and scaling
# constants, which is both stable and materially faster than repeated log-sum-exp
# calculations inside nested R loops.
forward_backward_hmm <- function(log_emission, pi, A) {
  T <- nrow(log_emission)
  K <- ncol(log_emission)

  row_shift <- apply(log_emission, 1L, max)
  emission <- exp(log_emission - row_shift)

  alpha <- matrix(0, nrow = T, ncol = K)
  beta <- matrix(1, nrow = T, ncol = K)
  scale <- numeric(T)

  a1 <- as.numeric(pi) * emission[1L, ]
  scale[1L] <- sum(a1)
  if (!is.finite(scale[1L]) || scale[1L] <= 0) {
    stop("Forward recursion failed at the first observation.")
  }
  alpha[1L, ] <- a1 / scale[1L]

  if (T >= 2L) {
    for (t in 2:T) {
      at <- drop(alpha[t - 1L, ] %*% A) * emission[t, ]
      scale[t] <- sum(at)
      if (!is.finite(scale[t]) || scale[t] <= 0) {
        stop("Forward recursion produced a non-positive scaling constant at t = ", t, ".")
      }
      alpha[t, ] <- at / scale[t]
    }

    for (t in (T - 1L):1L) {
      beta[t, ] <- drop(A %*% (emission[t + 1L, ] * beta[t + 1L, ])) / scale[t + 1L]
    }
  }

  gamma <- alpha * beta
  gamma <- gamma / rowSums(gamma)

  xi <- array(0, dim = c(max(T - 1L, 1L), K, K))
  if (T >= 2L) {
    for (t in seq_len(T - 1L)) {
      right <- emission[t + 1L, ] * beta[t + 1L, ]
      tmp <- A * tcrossprod(alpha[t, ], right)
      den <- sum(tmp)
      if (is.finite(den) && den > 0) tmp <- tmp / den
      xi[t, , ] <- tmp
    }
  }

  loglik <- sum(row_shift + log(scale))
  list(alpha = alpha, beta = beta, gamma = gamma, xi = xi, loglik = loglik, scale = scale)
}

initialize_parameters <- function(y, B, K, phi_lower = 1, phi_upper = 500, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  T <- length(y)
  M <- ncol(B)

  if (K == 1L) {
    z0 <- rep.int(1L, T)
  } else {
    # Balanced, mildly randomized initial partitions are more stable than purely
    # random labels and reduce the number of multi-start fits needed in practice.
    qy <- stats::qlogis(y)
    jitter_sd <- max(stats::sd(qy), 1e-3) * 0.05
    score <- qy + stats::rnorm(T, sd = jitter_sd)
    z0 <- as.integer(cut(rank(score, ties.method = "random"),
                         breaks = seq(0, T, length.out = K + 1L),
                         include.lowest = TRUE, labels = FALSE))
  }

  qy <- stats::qlogis(y)
  beta_list <- vector("list", K)
  phi <- numeric(K)

  for (k in seq_len(K)) {
    w <- 0.05 + as.numeric(z0 == k)
    BW <- B * sqrt(w)
    XtWX <- crossprod(BW)
    XtWy <- crossprod(BW, qy * sqrt(w))
    beta_list[[k]] <- tryCatch(
      as.numeric(solve(XtWX + 1e-6 * diag(M), XtWy)),
      error = function(e) stats::rnorm(M, 0, 0.20)
    )

    mu_k <- inv_logit(drop(B %*% beta_list[[k]]))
    s2 <- stats::weighted.mean((y - mu_k)^2, w = w)
    mbar <- stats::weighted.mean(mu_k, w = w)
    phi_mom <- mbar * (1 - mbar) / max(s2, 1e-6) - 1
    if (!is.finite(phi_mom)) phi_mom <- 20
    phi[k] <- min(max(phi_mom, phi_lower), phi_upper)
  }

  if (K == 1L) {
    A <- matrix(1, 1L, 1L)
    pi <- 1
  } else {
    A <- matrix(stats::runif(K * K, min = 0.01, max = 0.10), K, K)
    diag(A) <- diag(A) + stats::runif(K, 0.75, 0.95)
    A <- normalize_rows(A)
    pi <- stationary_distribution(A)
  }

  list(pi = pi, A = A, beta = beta_list, phi = phi)
}

state_penalized_objective <- function(par,
                                      y,
                                      B,
                                      gamma_k,
                                      lambda,
                                      penalty_matrix,
                                      phi_lower = 1,
                                      phi_upper = 500) {
  M <- ncol(B)
  beta_k <- par[seq_len(M)]
  phi_k <- exp(par[M + 1L])
  mu_k <- inv_logit(drop(B %*% beta_k))
  nll <- -sum(gamma_k * log_beta_density(y, mu_k, phi_k))
  penalty <- as.numeric(lambda) * drop(crossprod(beta_k, penalty_matrix %*% beta_k))
  nll + penalty
}

# Analytical gradient of the negative penalized expected log-likelihood.
# The final component is with respect to log(phi), matching the optimization
# parameterization and enforcing positivity without ad-hoc clipping.
state_penalized_gradient <- function(par,
                                     y,
                                     B,
                                     gamma_k,
                                     lambda,
                                     penalty_matrix,
                                     phi_lower = 1,
                                     phi_upper = 500) {
  M <- ncol(B)
  beta_k <- par[seq_len(M)]
  phi_k <- exp(par[M + 1L])
  eta <- drop(B %*% beta_k)
  mu <- clip_unit(inv_logit(eta), 1e-10)
  dmu_deta <- mu * (1 - mu)

  dloglik_dmu <- phi_k * (
    log(y) - log1p(-y) -
      digamma(mu * phi_k) + digamma((1 - mu) * phi_k)
  )
  score_eta <- gamma_k * dloglik_dmu * dmu_deta
  grad_beta <- -drop(crossprod(B, score_eta)) +
    2 * as.numeric(lambda) * drop(penalty_matrix %*% beta_k)

  dloglik_dphi <- digamma(phi_k) -
    mu * digamma(mu * phi_k) -
    (1 - mu) * digamma((1 - mu) * phi_k) +
    mu * log(y) + (1 - mu) * log1p(-y)
  grad_log_phi <- -sum(gamma_k * dloglik_dphi) * phi_k

  c(grad_beta, grad_log_phi)
}

state_penalty_components <- function(beta_list, penalty_matrix) {
  vapply(beta_list, function(b) drop(crossprod(b, penalty_matrix %*% b)), numeric(1))
}

penalized_observed_objective <- function(y, B, beta_list, phi, pi, A, lambda, penalty_matrix) {
  em <- compute_log_emission(y, B, beta_list, phi)
  fb <- forward_backward_hmm(em$log_emission, pi, A)
  lambda_vec <- if (length(lambda) == 1L) rep(as.numeric(lambda), length(phi)) else as.numeric(lambda)
  penalty <- sum(lambda_vec * state_penalty_components(beta_list, penalty_matrix))
  list(value = fb$loglik - penalty, fb = fb, em = em, penalty = penalty)
}

# Interpolate between two parameter sets when a numerical GEM step causes a tiny
# decrease in the penalized observed-data objective.  This is a numerical
# safeguard, not a separate estimation method.
interpolate_parameters <- function(old, new, step_size) {
  K <- length(old$phi)
  beta <- lapply(seq_len(K), function(k) old$beta[[k]] + step_size * (new$beta[[k]] - old$beta[[k]]))
  log_phi <- log(old$phi) + step_size * (log(new$phi) - log(old$phi))
  pi <- normalize_vec(old$pi + step_size * (new$pi - old$pi))
  A <- normalize_rows(old$A + step_size * (new$A - old$A))
  list(beta = beta, phi = exp(log_phi), pi = pi, A = A)
}

fit_single_start_beta_gam_hmm <- function(y,
                                          B,
                                          K,
                                          lambda,
                                          penalty_matrix,
                                          max_iter = 200L,
                                          tol = 1e-6,
                                          phi_lower = 1,
                                          phi_upper = 500,
                                          optimizer_maxit = 100L,
                                          mstep_tol = 1e-8,
                                          monotone_tol = 1e-8,
                                          boundary_handling = "error",
                                          verbose = FALSE,
                                          seed = NULL,
                                          initial_parameters = NULL) {
  y <- prepare_unit_interval_response(y, boundary_handling = boundary_handling)
  T <- length(y)
  M <- ncol(B)
  lambda_vec <- if (length(lambda) == 1L) rep(as.numeric(lambda), K) else as.numeric(lambda)
  if (length(lambda_vec) != K) stop("lambda must be scalar or length K.")

  init <- initial_parameters %||% initialize_parameters(
    y = y, B = B, K = K,
    phi_lower = phi_lower, phi_upper = phi_upper, seed = seed
  )
  pi_vec <- normalize_vec(init$pi)
  A_mat <- normalize_rows(init$A)
  beta_list <- init$beta
  phi_vec <- pmin(pmax(as.numeric(init$phi), phi_lower), phi_upper)

  current <- penalized_observed_objective(
    y, B, beta_list, phi_vec, pi_vec, A_mat, lambda_vec, penalty_matrix
  )
  objective_prev <- current$value
  converged <- FALSE
  monotonicity_ok <- TRUE
  stalled <- FALSE
  total_mstep_rejections <- 0L
  history_rows <- vector("list", max_iter)

  for (iter in seq_len(max_iter)) {
    gamma <- current$fb$gamma
    xi <- current$fb$xi
    old <- list(beta = beta_list, phi = phi_vec, pi = pi_vec, A = A_mat)

    pi_new <- normalize_vec(gamma[1L, ])
    A_new <- A_mat
    if (T >= 2L) {
      A_num <- apply(xi[seq_len(T - 1L), , , drop = FALSE], c(2, 3), sum) + 1e-12
      A_new <- normalize_rows(A_num)
    }

    beta_new <- beta_list
    phi_new <- phi_vec
    iter_rejections <- 0L

    for (k in seq_len(K)) {
      par0 <- c(beta_list[[k]], log(phi_vec[k]))
      old_q <- state_penalized_objective(
        par0, y, B, gamma[, k], lambda_vec[k], penalty_matrix,
        phi_lower, phi_upper
      )
      opt <- tryCatch(
        stats::optim(
          par = par0,
          fn = state_penalized_objective,
          gr = state_penalized_gradient,
          method = "L-BFGS-B",
          lower = c(rep(-Inf, M), log(phi_lower)),
          upper = c(rep(Inf, M), log(phi_upper)),
          y = y,
          B = B,
          gamma_k = gamma[, k],
          lambda = lambda_vec[k],
          penalty_matrix = penalty_matrix,
          phi_lower = phi_lower,
          phi_upper = phi_upper,
          control = list(maxit = optimizer_maxit, factr = 1e7)
        ),
        error = function(e) NULL
      )

      accept <- !is.null(opt) && is.finite(opt$value) && opt$value <= old_q + mstep_tol
      if (accept) {
        beta_new[[k]] <- opt$par[seq_len(M)]
        phi_new[k] <- exp(opt$par[M + 1L])
      } else {
        iter_rejections <- iter_rejections + 1L
      }
    }
    total_mstep_rejections <- total_mstep_rejections + iter_rejections

    proposal <- list(beta = beta_new, phi = phi_new, pi = pi_new, A = A_new)
    proposed_eval <- penalized_observed_objective(
      y, B, proposal$beta, proposal$phi, proposal$pi, proposal$A,
      lambda_vec, penalty_matrix
    )

    step_size <- 1
    accepted <- proposed_eval$value >= objective_prev - monotone_tol
    if (!accepted) {
      # Backtracking is rarely needed when the statewise Q updates are accepted,
      # but protects against accumulated numerical error in the generalized M-step.
      for (j in seq_len(8L)) {
        step_size <- 0.5^j
        damped <- interpolate_parameters(old, proposal, step_size)
        damped_eval <- penalized_observed_objective(
          y, B, damped$beta, damped$phi, damped$pi, damped$A,
          lambda_vec, penalty_matrix
        )
        if (damped_eval$value >= objective_prev - monotone_tol) {
          proposal <- damped
          proposed_eval <- damped_eval
          accepted <- TRUE
          break
        }
      }
    }

    if (!accepted) {
      monotonicity_ok <- FALSE
      stalled <- TRUE
      history_rows[[iter]] <- data.frame(
        iter = iter,
        loglik_obs = current$fb$loglik,
        objective_pen = objective_prev,
        delta_objective = 0,
        step_size = 0,
        mstep_rejections = iter_rejections,
        stringsAsFactors = FALSE
      )
      break
    }

    beta_list <- proposal$beta
    phi_vec <- proposal$phi
    pi_vec <- proposal$pi
    A_mat <- proposal$A
    delta_obj <- proposed_eval$value - objective_prev
    current <- proposed_eval

    history_rows[[iter]] <- data.frame(
      iter = iter,
      loglik_obs = current$fb$loglik,
      objective_pen = current$value,
      delta_objective = delta_obj,
      step_size = step_size,
      mstep_rejections = iter_rejections,
      stringsAsFactors = FALSE
    )

    if (verbose) {
      message(sprintf(
        "Iter %03d | logLik=%.4f | penalized=%.4f | delta=%.3e | step=%.3f",
        iter, current$fb$loglik, current$value, delta_obj, step_size
      ))
    }

    rel_change <- abs(delta_obj) / (1 + abs(objective_prev))
    objective_prev <- current$value
    if (rel_change < tol) {
      converged <- TRUE
      break
    }
  }

  history_rows <- Filter(Negate(is.null), history_rows)
  history <- if (length(history_rows)) do.call(rbind, history_rows) else data.frame()

  list(
    beta = beta_list,
    phi = phi_vec,
    pi = normalize_vec(pi_vec),
    A = normalize_rows(A_mat),
    gamma = current$fb$gamma,
    xi = current$fb$xi,
    eta = current$em$eta,
    mu = current$em$mu,
    log_emission = current$em$log_emission,
    loglik_obs = current$fb$loglik,
    objective_pen = current$value,
    converged = converged,
    monotonicity_ok = monotonicity_ok,
    stalled = stalled,
    mstep_rejections = total_mstep_rejections,
    n_iter = if (nrow(history) > 0L) max(history$iter) else 0L,
    history = history,
    K = K,
    lambda = lambda_vec,
    penalty_matrix = penalty_matrix,
    y = y,
    B = B,
    basis_spec = NULL
  )
}

format_lambda <- function(lambda, digits = 3L) {
  lam <- as.numeric(lambda)
  if (length(lam) == 1L) return(format(round(lam, digits), nsmall = digits, trim = TRUE))
  paste(format(round(lam, digits), nsmall = digits, trim = TRUE), collapse = " | ")
}

summarize_start_fits <- function(fits, diagnostic) {
  rows <- lapply(seq_along(fits), function(i) {
    fit_i <- fits[[i]]
    diag_i <- do.call(diagnose_phi_degeneracy, c(list(phi = sort(fit_i$phi)), diagnostic))
    data.frame(
      start_id = i,
      converged = isTRUE(fit_i$converged),
      monotonicity_ok = isTRUE(fit_i$monotonicity_ok),
      stalled = isTRUE(fit_i$stalled),
      loglik_obs = fit_i$loglik_obs,
      objective_pen = fit_i$objective_pen,
      n_iter = fit_i$n_iter,
      mstep_rejections = fit_i$mstep_rejections %||% 0L,
      phi_min = min(fit_i$phi),
      phi_max = max(fit_i$phi),
      phi_range = max(fit_i$phi) - min(fit_i$phi),
      lambda = format_lambda(fit_i$lambda),
      degenerate = isTRUE(diag_i$flagged),
      n_sat = diag_i$n_sat,
      delta_tail = diag_i$delta_tail,
      large_jump = diag_i$large_jump,
      phi_sorted = format_phi(sort(fit_i$phi)),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$loglik_gap_from_best <- max(out$loglik_obs, na.rm = TRUE) - out$loglik_obs
  out[order(!out$converged, !out$monotonicity_ok, -out$loglik_obs), , drop = FALSE]
}

fit_beta_gam_hmm <- function(y,
                             basis_spec,
                             K,
                             lambda,
                             n_start = 8L,
                             n_start_cores = 1L,
                             max_iter = 200L,
                             tol = 1e-6,
                             phi_lower = 1,
                             phi_upper = 500,
                             optimizer_maxit = 100L,
                             mstep_tol = 1e-8,
                             monotone_tol = 1e-8,
                             boundary_handling = "error",
                             store_all_starts = FALSE,
                             verbose = FALSE,
                             diagnostic = list(
                               phi_max = 500,
                               eps = 1e-6,
                               s_thresh = 2L,
                               delta_abs = 50,
                               delta_sum = 100
                             ),
                             seed = NULL) {
  started <- proc.time()[[3L]]
  y <- prepare_unit_interval_response(y, boundary_handling = boundary_handling)
  B <- basis_spec$B
  penalty_matrix <- basis_spec$penalty
  seeds <- if (is.null(seed)) sample.int(1e7, n_start) else seed + seq_len(n_start)

  fit_one_start <- function(s) {
    tryCatch(
      fit_single_start_beta_gam_hmm(
        y = y,
        B = B,
        K = K,
        lambda = lambda,
        penalty_matrix = penalty_matrix,
        max_iter = max_iter,
        tol = tol,
        phi_lower = phi_lower,
        phi_upper = phi_upper,
        optimizer_maxit = optimizer_maxit,
        mstep_tol = mstep_tol,
        monotone_tol = monotone_tol,
        boundary_handling = "error",
        verbose = verbose,
        seed = seeds[s]
      ),
      error = function(e) NULL
    )
  }

  fits <- parallel_lapply_crossplatform(
    X = seq_len(n_start),
    FUN = fit_one_start,
    n_cores = min(as.integer(n_start_cores %||% 1L), n_start),
    exports = c(
      "fit_one_start", "y", "B", "K", "lambda", "penalty_matrix", "max_iter", "tol",
      "phi_lower", "phi_upper", "optimizer_maxit", "mstep_tol", "monotone_tol", "verbose", "seeds"
    ),
    project_root = getwd(),
    seed = if (is.null(seed)) NULL else seed + 50000L,
    progress = FALSE
  )

  fits <- Filter(Negate(is.null), fits)
  if (length(fits) == 0L) stop("All multi-start fits failed for K = ", K, ".")

  start_summary <- summarize_start_fits(fits, diagnostic = diagnostic)
  logliks <- vapply(fits, function(f) f$loglik_obs, numeric(1))
  reliable <- vapply(fits, function(f) isTRUE(f$converged) && isTRUE(f$monotonicity_ok), logical(1))
  candidate_idx <- if (any(reliable)) which(reliable) else seq_along(fits)
  best_idx <- candidate_idx[which.max(logliks[candidate_idx])]
  best_fit <- fits[[best_idx]]

  ord <- order(best_fit$phi)
  best_fit <- relabel_fit(best_fit, ord)
  best_fit$basis_spec <- basis_spec
  best_fit$y <- y
  best_fit$diagnostic <- do.call(diagnose_phi_degeneracy, c(list(phi = best_fit$phi), diagnostic))
  best_fit$occupancy <- state_occupancy_stats(best_fit$gamma)
  best_fit$edf_spline <- compute_fit_edf(best_fit)
  best_fit$penalty_components <- state_penalty_components(best_fit$beta, best_fit$penalty_matrix)
  best_fit$criteria <- calculate_information_criteria(
    loglik_obs = best_fit$loglik_obs,
    gamma = best_fit$gamma,
    K = K,
    M = ncol(B),
    T = length(y),
    edf_spline = best_fit$edf_spline
  )
  best_fit$all_starts <- if (isTRUE(store_all_starts)) fits else NULL
  best_fit$best_start_index <- best_idx
  best_fit$start_summary <- start_summary
  best_fit$start_diagnostics <- list(
    n_start = n_start,
    n_successful_starts = nrow(start_summary),
    n_failed_starts = n_start - nrow(start_summary),
    n_converged = sum(start_summary$converged),
    n_monotone = sum(start_summary$monotonicity_ok),
    n_degenerate = sum(start_summary$degenerate),
    best_loglik = max(start_summary$loglik_obs),
    median_loglik = stats::median(start_summary$loglik_obs),
    loglik_gap_best_vs_median = max(start_summary$loglik_obs) - stats::median(start_summary$loglik_obs)
  )
  best_fit$timing_seconds <- proc.time()[[3L]] - started
  best_fit
}

viterbi_decode <- function(y, basis_spec, fit, boundary_handling = "error") {
  y <- prepare_unit_interval_response(y, boundary_handling = boundary_handling)
  B <- basis_spec$B
  em <- compute_log_emission(y, B, fit$beta, fit$phi)
  log_emission <- em$log_emission
  T <- nrow(log_emission)
  K <- ncol(log_emission)

  delta <- matrix(-Inf, nrow = T, ncol = K)
  psi <- matrix(NA_integer_, nrow = T, ncol = K)
  delta[1L, ] <- log(pmax(fit$pi, 1e-16)) + log_emission[1L, ]

  if (T >= 2L) {
    log_A <- log(pmax(fit$A, 1e-16))
    for (t in 2:T) {
      for (k in seq_len(K)) {
        vals <- delta[t - 1L, ] + log_A[, k]
        psi[t, k] <- which.max(vals)
        delta[t, k] <- max(vals) + log_emission[t, k]
      }
    }
  }

  path <- integer(T)
  path[T] <- which.max(delta[T, ])
  if (T >= 2L) {
    for (t in (T - 1L):1L) path[t] <- psi[t + 1L, path[t + 1L]]
  }
  path
}
