# Scalar annual aggregate benchmark.
#
# This is a scientifically valid scalar proportion time series: for each year,
# female deaths are aggregated over the selected ages and divided by total
# deaths over the same ages.  It is fitted with an intercept-only Beta HMM
# (state-specific mean and precision).  The benchmark deliberately does not
# flatten age-year cells into a fictitious 2,255-step HMM sequence.

prepare_scalar_aggregate <- function(panel) {
  F_year <- rowSums(panel$F)
  N_year <- rowSums(panel$N)
  y <- F_year / N_year
  if (any(y <= 0 | y >= 1)) {
    stop("Scalar annual aggregate contains exact 0/1 values; Beta benchmark is undefined.")
  }
  list(
    years = panel$years,
    y = as.numeric(y),
    female = as.numeric(F_year),
    total = as.numeric(N_year),
    n_years = length(y)
  )
}

scalar_beta_log_emission <- function(y, mu, phi) {
  K <- length(mu)
  out <- vapply(seq_len(K), function(k) {
    stats::dbeta(y, mu[k] * phi[k], (1 - mu[k]) * phi[k], log = TRUE)
  }, numeric(length(y)))
  if (K == 1L) out <- matrix(out, ncol = 1L)
  out
}

initialize_scalar_beta_hmm <- function(y, K, phi_lower, phi_upper, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Tn <- length(y)
  score <- stats::qlogis(y)
  if (K == 1L) {
    z <- rep.int(1L, Tn)
  } else {
    score <- score + stats::rnorm(Tn, 0, max(stats::sd(score), 1e-5) * 0.05)
    z <- as.integer(cut(rank(score, ties.method = "random"),
                        breaks = seq(0, Tn, length.out = K + 1L),
                        include.lowest = TRUE, labels = FALSE))
  }
  mu <- phi <- numeric(K)
  for (k in seq_len(K)) {
    yy <- y[z == k]
    if (!length(yy)) yy <- y
    m <- min(max(mean(yy), 1e-5), 1 - 1e-5)
    v <- stats::var(yy)
    if (!is.finite(v) || v <= 1e-10) v <- max(m * (1 - m) / 100, 1e-8)
    ph <- m * (1 - m) / v - 1
    mu[k] <- m
    phi[k] <- min(max(ph, phi_lower), 0.8 * phi_upper)
  }
  if (K == 1L) {
    A <- matrix(1, 1, 1); pi <- 1
  } else {
    A <- matrix(stats::runif(K * K, 0.01, 0.06), K, K)
    diag(A) <- diag(A) + stats::runif(K, 0.80, 0.95)
    A <- normalize_rows(A)
    pi <- stationary_distribution(A)
  }
  list(mu = mu, phi = phi, A = A, pi = pi)
}

scalar_state_objective <- function(par, y, w, phi_lower, phi_upper) {
  mu <- stats::plogis(par[1L])
  phi <- exp(par[2L])
  -sum(w * stats::dbeta(y, mu * phi, (1 - mu) * phi, log = TRUE))
}

scalar_observed_loglik <- function(y, mu, phi, pi, A) {
  em <- scalar_beta_log_emission(y, mu, phi)
  fb <- forward_backward_hmm(em, pi, A)
  list(value = fb$loglik, fb = fb, log_emission = em)
}

fit_single_scalar_beta_hmm <- function(series, K, max_iter = 300L, tol = 1e-7,
                                       phi_lower = 1, phi_upper = 2000,
                                       optimizer_maxit = 150L,
                                       monotone_tol = 1e-8, seed = NULL) {
  y <- series$y
  init <- initialize_scalar_beta_hmm(y, K, phi_lower, phi_upper, seed)
  mu <- init$mu; phi <- init$phi; pi <- init$pi; A <- init$A
  current <- scalar_observed_loglik(y, mu, phi, pi, A)
  prev <- current$value
  converged <- FALSE
  monotonicity_ok <- TRUE
  history <- vector("list", max_iter)

  for (iter in seq_len(max_iter)) {
    gamma <- current$fb$gamma
    xi <- current$fb$xi
    old <- list(mu = mu, phi = phi, pi = pi, A = A)
    pi_new <- normalize_vec(gamma[1L, ])
    A_new <- A
    if (length(y) >= 2L && K > 1L) {
      A_new <- normalize_rows(apply(xi[seq_len(length(y) - 1L), , , drop = FALSE], c(2, 3), sum) + 1e-12)
    }
    mu_new <- mu; phi_new <- phi
    for (k in seq_len(K)) {
      p0 <- c(stats::qlogis(mu[k]), log(phi[k]))
      old_q <- scalar_state_objective(p0, y, gamma[, k], phi_lower, phi_upper)
      opt <- tryCatch(stats::optim(
        p0, scalar_state_objective, method = "L-BFGS-B",
        lower = c(-Inf, log(phi_lower)), upper = c(Inf, log(phi_upper)),
        y = y, w = gamma[, k], phi_lower = phi_lower, phi_upper = phi_upper,
        control = list(maxit = optimizer_maxit, factr = 1e7)
      ), error = function(e) NULL)
      if (!is.null(opt) && is.finite(opt$value) && opt$value <= old_q + 1e-7) {
        mu_new[k] <- stats::plogis(opt$par[1L])
        phi_new[k] <- exp(opt$par[2L])
      }
    }
    proposal <- scalar_observed_loglik(y, mu_new, phi_new, pi_new, A_new)
    accepted <- proposal$value >= prev - monotone_tol
    step <- 1
    if (!accepted) {
      for (j in seq_len(8L)) {
        step <- 0.5^j
        mu_d <- stats::plogis(stats::qlogis(old$mu) + step * (stats::qlogis(mu_new) - stats::qlogis(old$mu)))
        phi_d <- exp(log(old$phi) + step * (log(phi_new) - log(old$phi)))
        pi_d <- normalize_vec(old$pi + step * (pi_new - old$pi))
        A_d <- normalize_rows(old$A + step * (A_new - old$A))
        cand <- scalar_observed_loglik(y, mu_d, phi_d, pi_d, A_d)
        if (cand$value >= prev - monotone_tol) {
          mu_new <- mu_d; phi_new <- phi_d; pi_new <- pi_d; A_new <- A_d
          proposal <- cand; accepted <- TRUE; break
        }
      }
    }
    if (!accepted) { monotonicity_ok <- FALSE; break }
    delta <- proposal$value - prev
    mu <- mu_new; phi <- phi_new; pi <- pi_new; A <- A_new; current <- proposal
    history[[iter]] <- data.frame(iter = iter, loglik = current$value, delta = delta, step = step)
    rel <- abs(delta) / (1 + abs(prev)); prev <- current$value
    if (rel < tol) { converged <- TRUE; break }
  }

  ord <- order(mu)
  mu <- mu[ord]; phi <- phi[ord]; pi <- pi[ord]; A <- A[ord, ord, drop = FALSE]
  gamma <- current$fb$gamma[, ord, drop = FALSE]
  log_em <- current$log_emission[, ord, drop = FALSE]
  viterbi <- viterbi_from_log_emission(log_em, pi, A)
  occupancy <- state_occupancy(gamma)
  npar <- 2L * K + K * (K - 1L) + (K - 1L)
  entropy <- entropy_from_gamma(gamma)
  bic <- -2 * current$value + log(length(y)) * npar
  icl <- bic + 2 * entropy
  list(
    family = "scalar_beta", K = K, mu = mu, phi = phi, pi = pi, A = A,
    gamma = gamma, log_emission = log_em, viterbi = viterbi,
    loglik_obs = current$value, converged = converged,
    monotonicity_ok = monotonicity_ok, occupancy = occupancy,
    precision_diagnostic = precision_safeguard(phi, phi_upper, 0.98),
    criteria = list(bic = bic, icl = icl, aic = -2 * current$value + 2 * npar,
                    npar = npar, entropy = entropy),
    history = Filter(Negate(is.null), history), series = series
  )
}

scalar_multistart_worker <- function(job) {
  tryCatch(do.call(fit_single_scalar_beta_hmm, job),
           error = function(e) structure(list(error = conditionMessage(e)), class = "fit_error"))
}

fit_scalar_beta_hmm <- function(series, K, n_start = 12L, n_cores = 1L,
                                max_iter = 300L, tol = 1e-7,
                                phi_lower = 1, phi_upper = 2000,
                                optimizer_maxit = 150L, seed = 1L) {
  started <- proc.time()[[3L]]
  jobs <- lapply(seq_len(n_start), function(s) list(
    series = series, K = K, max_iter = max_iter, tol = tol,
    phi_lower = phi_lower, phi_upper = phi_upper,
    optimizer_maxit = optimizer_maxit, seed = seed + s
  ))
  fits <- if (n_cores > 1L) {
    parallel_lapply(jobs, scalar_multistart_worker, n_cores = n_cores,
                    project_root = getwd(), seed = seed + 90000L)
  } else lapply(jobs, scalar_multistart_worker)
  good <- vapply(fits, function(x) !inherits(x, "fit_error"), logical(1))
  if (!any(good)) stop("All scalar benchmark starts failed for K=", K, ".")
  ff <- fits[good]
  reliable <- vapply(ff, function(f) f$converged && f$monotonicity_ok && !f$precision_diagnostic$flagged, logical(1))
  pool <- if (any(reliable)) which(reliable) else seq_along(ff)
  best <- ff[[pool[which.max(vapply(ff[pool], `[[`, numeric(1), "loglik_obs"))]]]
  best$timing_seconds <- proc.time()[[3L]] - started
  best
}

run_scalar_beta_grid <- function(series, K_grid = 1:8, n_start = 12L,
                                 n_cores = 1L, max_iter = 300L, tol = 1e-7,
                                 phi_lower = 1, phi_upper = 2000,
                                 optimizer_maxit = 150L, seed = 1L,
                                 progress = TRUE) {
  out <- vector("list", length(K_grid))
  for (i in seq_along(K_grid)) {
    K <- K_grid[i]
    if (progress) message("Scalar annual benchmark: K=", K, " (", i, "/", length(K_grid), ")")
    out[[i]] <- tryCatch(fit_scalar_beta_hmm(
      series, K, n_start, n_cores, max_iter, tol, phi_lower, phi_upper,
      optimizer_maxit, seed + 10000L * K
    ), error = function(e) structure(list(K = K, error = conditionMessage(e)), class = "grid_error"))
  }
  names(out) <- paste0("K", K_grid)
  out
}

summarize_scalar_beta_grid <- function(grid) {
  do.call(rbind, lapply(grid, function(f) {
    if (inherits(f, "grid_error")) return(data.frame(
      K = f$K, converged = FALSE, monotonicity_ok = FALSE, precision_flag = NA,
      loglik = NA_real_, aic = NA_real_, bic = NA_real_, icl = NA_real_,
      npar = NA_real_, min_state_share = NA_real_, timing_seconds = NA_real_,
      error = f$error, stringsAsFactors = FALSE))
    data.frame(
      K = f$K, converged = f$converged, monotonicity_ok = f$monotonicity_ok,
      precision_flag = f$precision_diagnostic$flagged, loglik = f$loglik_obs,
      aic = f$criteria$aic, bic = f$criteria$bic, icl = f$criteria$icl,
      npar = f$criteria$npar, min_state_share = f$occupancy$min_share,
      timing_seconds = f$timing_seconds, error = NA_character_, stringsAsFactors = FALSE
    )
  }))
}

select_scalar_beta_model <- function(grid) {
  tab <- summarize_scalar_beta_grid(grid)
  ok <- tab[is.na(tab$error) & tab$converged & tab$monotonicity_ok & is.finite(tab$bic), , drop = FALSE]
  if (!nrow(ok)) stop("No admissible scalar annual Beta HMM fit.")
  row <- ok[order(ok$bic, ok$icl, ok$K), , drop = FALSE][1L, ]
  fit <- grid[[which(vapply(grid, function(f) !inherits(f, "grid_error") && f$K == row$K, logical(1)))[1L]]]
  list(fit = fit, grid_summary = tab, selected_row = row)
}

compute_scalar_beta_pseudo_residuals <- function(fit) {
  y <- fit$series$y
  Tn <- length(y); K <- fit$K
  pred_state <- matrix(NA_real_, Tn, K)
  pred_state[1L, ] <- fit$pi
  if (Tn >= 2L) {
    # Reconstruct filtered alpha from the fitted emission model.
    fb <- forward_backward_hmm(fit$log_emission, fit$pi, fit$A)
    for (t in 2:Tn) pred_state[t, ] <- drop(fb$alpha[t - 1L, ] %*% fit$A)
  }
  u <- numeric(Tn)
  for (t in seq_len(Tn)) {
    cdf_k <- vapply(seq_len(K), function(k) {
      stats::pbeta(y[t], fit$mu[k] * fit$phi[k], (1 - fit$mu[k]) * fit$phi[k])
    }, numeric(1))
    u[t] <- sum(pred_state[t, ] * cdf_k)
  }
  q <- stats::qnorm(clip_unit(u, 1e-10))
  lag_lb <- min(10L, floor(Tn / 3L))
  data.frame(
    mean = mean(q), sd = stats::sd(q),
    ljung_box_p = if (lag_lb >= 1L) stats::Box.test(q, lag = lag_lb, type = "Ljung-Box")$p.value else NA_real_,
    stringsAsFactors = FALSE
  )
}
