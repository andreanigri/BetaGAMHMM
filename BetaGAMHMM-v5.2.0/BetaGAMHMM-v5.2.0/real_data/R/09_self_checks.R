# Lightweight invariants executed at the start of RUN_ALL.R.

run_internal_checks <- function() {
  # The gamma-extended kernel must reproduce dbinom exactly for integer counts.
  F <- c(0, 2, 5, 11)
  M <- c(5, 3, 0, 9)
  p <- c(0.2, 0.4, 0.7, 0.55)
  lhs <- binomial_kernel_cell_loglik(F, M, p)
  rhs <- stats::dbinom(F, size = F + M, prob = p, log = TRUE)
  if (max(abs(lhs - rhs)) > 1e-10) {
    stop("Internal check failed: binomial kernel does not reproduce dbinom for integer counts.")
  }

  # Forward-backward posterior probabilities must be normalized.
  log_em <- matrix(c(-0.2, -1.1, -0.4, -0.8, -1.0, -0.3), nrow = 3L, byrow = TRUE)
  A <- matrix(c(0.9, 0.1, 0.2, 0.8), 2L, 2L, byrow = TRUE)
  fb <- forward_backward_hmm(log_em, pi = c(0.5, 0.5), A = A)
  if (max(abs(rowSums(fb$gamma) - 1)) > 1e-10) {
    stop("Internal check failed: forward-backward posterior probabilities are not normalized.")
  }

  # Relabelling must also transform an already decoded Viterbi sequence.
  toy <- list(
    beta = list(1, 2), phi = c(10, 20), lambda = c(1, 1),
    pi = c(0.4, 0.6), A = A, gamma = fb$gamma, alpha = fb$alpha,
    xi = fb$xi, mu_age = matrix(c(0.2, 0.4, 0.3, 0.5), 2, 2),
    log_emission = log_em, viterbi = c(1L, 2L, 1L)
  )
  toy2 <- relabel_fit(toy, c(2L, 1L))
  if (!identical(as.integer(toy2$viterbi), c(2L, 1L, 2L))) {
    stop("Internal check failed: Viterbi labels were not transformed with state relabelling.")
  }

  # Scalable assignment must recover a known minimum-cost relabelling.
  C <- matrix(c(9, 1, 8, 2, 7, 3, 6, 5, 0), 3L, 3L, byrow = TRUE)
  ass <- min_cost_assignment_dp(C)
  if (!identical(as.integer(ass), c(2L, 1L, 3L))) {
    stop("Internal check failed: dynamic-programming state assignment is incorrect.")
  }

  # The final BIC must be finite for a simple repeated-profile example.
  toy_gamma <- matrix(c(0.8, 0.2, 0.3, 0.7), 2L, 2L, byrow = TRUE)
  ic <- information_criteria(loglik = -10, gamma = toy_gamma, K = 2L,
                             edf_spline = c(2, 2), n_years = 55L,
                             n_cells = 2255L, has_phi = TRUE)
  if (!is.finite(ic$bic_edf) || !is.finite(ic$aic_edf) || !is.finite(ic$icl_edf)) {
    stop("Internal check failed: information criteria are not finite.")
  }

  invisible(TRUE)
}
