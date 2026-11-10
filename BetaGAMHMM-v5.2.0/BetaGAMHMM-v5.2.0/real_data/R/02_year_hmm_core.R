# Core estimation routines for the year-level vector-emission HMM.
#
# For calendar year t, the latent state S_t is shared by all ages. Conditional
# on S_t = k, age-specific proportions are modelled by a product of Beta
# densities with a smooth state-specific mean curve mu_k(age) and one precision
# phi_k. A Binomial analogue is provided later as a denominator-aware
# sensitivity model.

beta_log_emission_year <- function(Y, B_age, beta_list, phi) {
  Tn <- nrow(Y)
  J <- ncol(Y)
  K <- length(phi)
  mu_age <- vapply(seq_len(K), function(k) inv_logit(drop(B_age %*% beta_list[[k]])), numeric(J))
  if (K == 1L) mu_age <- matrix(mu_age, ncol = 1L)
  out <- matrix(0, Tn, K)
  for (k in seq_len(K)) {
    a <- mu_age[, k] * phi[k]
    b <- (1 - mu_age[, k]) * phi[k]
    ld <- vapply(seq_len(J), function(j) stats::dbeta(Y[, j], a[j], b[j], log = TRUE), numeric(Tn))
    if (J == 1L) ld <- matrix(ld, ncol = 1L)
    out[, k] <- rowSums(ld)
  }
  list(log_emission = out, mu_age = mu_age)
}

initialize_beta_year_hmm <- function(Y, B_age, K, phi_lower, phi_upper, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Tn <- nrow(Y)
  M <- ncol(B_age)
  yearly_score <- rowMeans(stats::qlogis(clip_unit(Y)))

  if (K == 1L) {
    z0 <- rep.int(1L, Tn)
  } else {
    score <- yearly_score + stats::rnorm(Tn, 0, max(stats::sd(yearly_score), 1e-4) * 0.05)
    z0 <- as.integer(cut(
      rank(score, ties.method = "random"),
      breaks = seq(0, Tn, length.out = K + 1L),
      include.lowest = TRUE,
      labels = FALSE
    ))
  }

  beta_list <- vector("list", K)
  phi <- numeric(K)
  for (k in seq_len(K)) {
    idx <- which(z0 == k)
    if (!length(idx)) idx <- seq_len(Tn)
    ybar <- colMeans(Y[idx, , drop = FALSE])
    eta_bar <- stats::qlogis(clip_unit(ybar, 1e-6))
    beta_list[[k]] <- tryCatch(
      as.numeric(qr.solve(B_age, eta_bar, tol = 1e-10)),
      error = function(e) stats::rnorm(M, 0, 0.1)
    )
    mu <- inv_logit(drop(B_age %*% beta_list[[k]]))
    resid2 <- as.vector(sweep(Y[idx, , drop = FALSE], 2L, mu, FUN = "-")^2)
    mu_rep <- rep(mu, each = length(idx))
    mbar <- mean(mu_rep)
    s2 <- mean(resid2)
    phi_mom <- mbar * (1 - mbar) / max(s2, 1e-8) - 1
    if (!is.finite(phi_mom)) phi_mom <- 50
    phi[k] <- min(max(phi_mom, phi_lower), 0.8 * phi_upper)
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
  list(beta = beta_list, phi = phi, A = A, pi = pi)
}

beta_state_objective <- function(par, Y, B_age, gamma_k, lambda, penalty, phi_lower, phi_upper) {
  M <- ncol(B_age)
  beta <- par[seq_len(M)]
  phi <- exp(par[M + 1L])
  mu <- clip_unit(inv_logit(drop(B_age %*% beta)))
  J <- ncol(Y)
  ll_year <- numeric(nrow(Y))
  for (j in seq_len(J)) {
    ll_year <- ll_year + stats::dbeta(Y[, j], mu[j] * phi, (1 - mu[j]) * phi, log = TRUE)
  }
  -sum(gamma_k * ll_year) + lambda * drop(crossprod(beta, penalty %*% beta))
}

beta_state_gradient <- function(par, Y, B_age, gamma_k, lambda, penalty, phi_lower, phi_upper) {
  M <- ncol(B_age)
  beta <- par[seq_len(M)]
  phi <- exp(par[M + 1L])
  eta <- drop(B_age %*% beta)
  mu <- clip_unit(inv_logit(eta))
  dmu <- mu * (1 - mu)
  J <- ncol(Y)

  score_eta_age <- numeric(J)
  score_phi <- 0
  for (j in seq_len(J)) {
    yj <- Y[, j]
    d_ll_dmu <- phi * (
      log(yj) - log1p(-yj) - digamma(mu[j] * phi) + digamma((1 - mu[j]) * phi)
    )
    score_eta_age[j] <- sum(gamma_k * d_ll_dmu) * dmu[j]

    d_ll_dphi <- digamma(phi) -
      mu[j] * digamma(mu[j] * phi) -
      (1 - mu[j]) * digamma((1 - mu[j]) * phi) +
      mu[j] * log(yj) + (1 - mu[j]) * log1p(-yj)
    score_phi <- score_phi + sum(gamma_k * d_ll_dphi)
  }

  grad_beta <- -drop(crossprod(B_age, score_eta_age)) + 2 * lambda * drop(penalty %*% beta)
  grad_log_phi <- -score_phi * phi
  c(grad_beta, grad_log_phi)
}

beta_penalized_observed <- function(Y, B_age, beta_list, phi, pi, A, lambda, penalty) {
  em <- beta_log_emission_year(Y, B_age, beta_list, phi)
  fb <- forward_backward_hmm(em$log_emission, pi, A)
  penalty_value <- sum(vapply(seq_along(beta_list), function(k) {
    lambda[k] * drop(crossprod(beta_list[[k]], penalty %*% beta_list[[k]]))
  }, numeric(1)))
  list(value = fb$loglik - penalty_value, fb = fb, em = em, penalty = penalty_value)
}

interpolate_beta_parameters <- function(old, new, step) {
  K <- length(old$phi)
  list(
    beta = lapply(seq_len(K), function(k) old$beta[[k]] + step * (new$beta[[k]] - old$beta[[k]])),
    phi = exp(log(old$phi) + step * (log(new$phi) - log(old$phi))),
    pi = normalize_vec(old$pi + step * (new$pi - old$pi)),
    A = normalize_rows(old$A + step * (new$A - old$A))
  )
}

fit_single_beta_year_hmm <- function(panel,
                                     basis_spec,
                                     K,
                                     lambda,
                                     max_iter = 300L,
                                     tol = 1e-6,
                                     phi_lower = 1,
                                     phi_upper = 2000,
                                     optimizer_maxit = 150L,
                                     monotone_tol = 1e-8,
                                     seed = NULL,
                                     init_fit = NULL) {
  Y <- panel$Y
  B_age <- basis_spec$B
  penalty <- basis_spec$penalty
  M <- ncol(B_age)
  lambda_vec <- if (length(lambda) == 1L) rep(as.numeric(lambda), K) else as.numeric(lambda)
  if (length(lambda_vec) != K) stop("lambda must be scalar or length K.")

  if (is.null(init_fit)) {
    init <- initialize_beta_year_hmm(Y, B_age, K, phi_lower, phi_upper, seed = seed)
    init_source <- "random"
  } else {
    if (is.null(init_fit$K) || as.integer(init_fit$K) != as.integer(K))
      stop("Warm-start fit has incompatible K.")
    if (length(init_fit$beta) != K || any(vapply(init_fit$beta, length, integer(1)) != M))
      stop("Warm-start fit has incompatible spline coefficient dimension.")
    if (length(init_fit$phi) != K || any(!is.finite(init_fit$phi)) ||
        any(init_fit$phi < phi_lower) || any(init_fit$phi > phi_upper))
      stop("Warm-start precision parameters fall outside the requested bounds.")
    if (is.null(init_fit$A) || is.null(dim(init_fit$A)) ||
        !identical(as.integer(dim(init_fit$A)), c(as.integer(K), as.integer(K))) ||
        is.null(init_fit$pi) || length(init_fit$pi) != K)
      stop("Warm-start transition/initial probabilities are incompatible.")
    init <- list(
      beta = lapply(init_fit$beta, as.numeric),
      phi = as.numeric(init_fit$phi),
      A = normalize_rows(init_fit$A),
      pi = normalize_vec(init_fit$pi)
    )
    init_source <- "warm_start"
  }
  beta_list <- init$beta
  phi <- init$phi
  pi <- init$pi
  A <- init$A

  current <- beta_penalized_observed(Y, B_age, beta_list, phi, pi, A, lambda_vec, penalty)
  objective_prev <- current$value
  converged <- FALSE
  monotonicity_ok <- TRUE
  history <- vector("list", max_iter)
  total_rejections <- 0L

  for (iter in seq_len(max_iter)) {
    gamma <- current$fb$gamma
    xi <- current$fb$xi
    old <- list(beta = beta_list, phi = phi, pi = pi, A = A)

    pi_new <- normalize_vec(gamma[1L, ])
    A_new <- A
    if (panel$n_years >= 2L && K > 1L) {
      A_new <- normalize_rows(apply(xi[seq_len(panel$n_years - 1L), , , drop = FALSE], c(2, 3), sum) + 1e-12)
    }

    beta_new <- beta_list
    phi_new <- phi
    iter_rejections <- 0L
    for (k in seq_len(K)) {
      par0 <- c(beta_list[[k]], log(phi[k]))
      old_q <- beta_state_objective(par0, Y, B_age, gamma[, k], lambda_vec[k], penalty, phi_lower, phi_upper)
      opt <- tryCatch(
        stats::optim(
          par = par0,
          fn = beta_state_objective,
          gr = beta_state_gradient,
          method = "L-BFGS-B",
          lower = c(rep(-Inf, M), log(phi_lower)),
          upper = c(rep(Inf, M), log(phi_upper)),
          Y = Y,
          B_age = B_age,
          gamma_k = gamma[, k],
          lambda = lambda_vec[k],
          penalty = penalty,
          phi_lower = phi_lower,
          phi_upper = phi_upper,
          control = list(maxit = optimizer_maxit, factr = 1e7)
        ),
        error = function(e) NULL
      )
      if (!is.null(opt) && is.finite(opt$value) && opt$value <= old_q + 1e-7) {
        beta_new[[k]] <- opt$par[seq_len(M)]
        phi_new[k] <- exp(opt$par[M + 1L])
      } else {
        iter_rejections <- iter_rejections + 1L
      }
    }
    total_rejections <- total_rejections + iter_rejections

    proposal <- list(beta = beta_new, phi = phi_new, pi = pi_new, A = A_new)
    proposed <- beta_penalized_observed(Y, B_age, proposal$beta, proposal$phi, proposal$pi, proposal$A, lambda_vec, penalty)
    accepted <- proposed$value >= objective_prev - monotone_tol
    step <- 1
    if (!accepted) {
      for (j in seq_len(8L)) {
        step <- 0.5^j
        damped <- interpolate_beta_parameters(old, proposal, step)
        damped_eval <- beta_penalized_observed(Y, B_age, damped$beta, damped$phi, damped$pi, damped$A, lambda_vec, penalty)
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
    phi <- proposal$phi
    pi <- proposal$pi
    A <- proposal$A
    delta <- proposed$value - objective_prev
    current <- proposed
    history[[iter]] <- data.frame(
      iter = iter,
      loglik = current$fb$loglik,
      penalized_objective = current$value,
      delta = delta,
      step = step,
      mstep_rejections = iter_rejections
    )
    rel <- abs(delta) / (1 + abs(objective_prev))
    objective_prev <- current$value
    if (rel < tol) {
      converged <- TRUE
      break
    }
  }

  mu_age <- current$em$mu_age
  fit <- list(
    family = "beta",
    K = K,
    lambda = lambda_vec,
    beta = beta_list,
    phi = phi,
    pi = pi,
    A = A,
    gamma = current$fb$gamma,
    xi = current$fb$xi,
    alpha = current$fb$alpha,
    log_emission = current$em$log_emission,
    mu_age = mu_age,
    loglik_obs = current$fb$loglik,
    objective_pen = current$value,
    converged = converged,
    monotonicity_ok = monotonicity_ok,
    n_iter = length(Filter(Negate(is.null), history)),
    mstep_rejections = total_rejections,
    initialization = init_source,
    history = Filter(Negate(is.null), history),
    panel = panel,
    basis_spec = basis_spec
  )

  ord <- order_states_by_profile(fit)
  fit <- relabel_fit(fit, ord)
  fit$viterbi <- viterbi_from_log_emission(fit$log_emission, fit$pi, fit$A)
  fit$occupancy <- state_occupancy(fit$gamma)
  fit$precision_diagnostic <- NULL
  fit
}

beta_multistart_worker <- function(job) {
  tryCatch(
    do.call(fit_single_beta_year_hmm, job),
    error = function(e) structure(
      list(error = conditionMessage(e)),
      class = "fit_error"
    )
  )
}

fit_beta_year_hmm <- function(panel,
                              basis_spec,
                              K,
                              lambda,
                              n_start = 15L,
                              n_cores = 1L,
                              max_iter = 300L,
                              tol = 1e-6,
                              phi_lower = 1,
                              phi_upper = 2000,
                              optimizer_maxit = 150L,
                              precision_saturation_fraction = 0.98,
                              seed = 1001L,
                              warm_start_fit = NULL) {
  started <- proc.time()[[3L]]
  n_start <- max(1L, as.integer(n_start))
  seeds <- seed + seq_len(n_start)

  if (!is.null(warm_start_fit)) {
    if (is.null(warm_start_fit$K) || as.integer(warm_start_fit$K) != as.integer(K))
      stop("warm_start_fit has incompatible K.")
    if (any(!is.finite(warm_start_fit$phi)) ||
        any(warm_start_fit$phi < phi_lower) || any(warm_start_fit$phi > phi_upper))
      stop("warm_start_fit is not feasible under the requested phi bounds.")
  }

  # Each PSOCK task is self-contained. When a feasible warm start is supplied,
  # it replaces the first random start, so the number of expensive fits does
  # not increase. The original primary fit is retained as a zero-cost candidate
  # below, preventing a poorer local optimum from displacing an already
  # feasible solution.
  jobs <- lapply(seq_len(n_start), function(s) list(
    panel = panel,
    basis_spec = basis_spec,
    K = K,
    lambda = lambda,
    max_iter = max_iter,
    tol = tol,
    phi_lower = phi_lower,
    phi_upper = phi_upper,
    optimizer_maxit = optimizer_maxit,
    seed = seeds[s],
    init_fit = if (!is.null(warm_start_fit) && s == 1L) warm_start_fit else NULL
  ))

  if (n_cores > 1L) {
    fits <- parallel_lapply(
      jobs,
      beta_multistart_worker,
      n_cores = n_cores,
      project_root = getwd(),
      seed = seed + 90000L
    )
  } else {
    fits <- lapply(jobs, beta_multistart_worker)
  }

  if (!is.null(warm_start_fit)) {
    warm_candidate <- warm_start_fit
    warm_candidate$converged <- TRUE
    warm_candidate$monotonicity_ok <- TRUE
    warm_candidate$initialization <- "primary_fit_candidate"
    fits <- c(list(warm_candidate), fits)
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
    stop("All multi-start fits failed for K=", K, ", lambda=", lambda, ".", detail)
  }

  fits_good <- fits[good]
  reliable <- vapply(fits_good, function(f) isTRUE(f$converged) && isTRUE(f$monotonicity_ok), logical(1))
  pool <- if (any(reliable)) which(reliable) else seq_along(fits_good)
  best <- fits_good[[pool[which.max(vapply(fits_good[pool], `[[`, numeric(1), "loglik_obs"))]]]

  best$precision_diagnostic <- precision_safeguard(best$phi, phi_upper, precision_saturation_fraction)
  best$edf_spline <- vapply(seq_len(K), function(k) {
    beta_edf_state(
      B = basis_spec$B,
      gamma_k = best$gamma[, k],
      mu_age = best$mu_age[, k],
      phi = best$phi[k],
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
    has_phi = TRUE
  )
  best$start_diagnostics <- data.frame(
    n_start = n_start,
    n_candidates = length(fits),
    warm_start_used = !is.null(warm_start_fit),
    n_successful = sum(good),
    n_converged = sum(vapply(fits_good, function(f) isTRUE(f$converged), logical(1))),
    n_monotone = sum(vapply(fits_good, function(f) isTRUE(f$monotonicity_ok), logical(1))),
    best_loglik = max(vapply(fits_good, `[[`, numeric(1), "loglik_obs")),
    median_loglik = stats::median(vapply(fits_good, `[[`, numeric(1), "loglik_obs")),
    stringsAsFactors = FALSE
  )
  best$timing_seconds <- proc.time()[[3L]] - started
  best
}
