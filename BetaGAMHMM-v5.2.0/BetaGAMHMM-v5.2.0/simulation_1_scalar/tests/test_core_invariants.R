# Lightweight numerical checks for the core implementation.
# Loaded automatically by RUN_ALL.R.

assert_true <- function(x, message) {
  if (!isTRUE(x)) stop(message, call. = FALSE)
}

set.seed(101)

# 1. Strict boundary handling: a standard Beta likelihood must reject 0 and 1.
boundary_error <- tryCatch({
  prepare_unit_interval_response(c(0, 0.5, 1), boundary_handling = "error")
  FALSE
}, error = function(e) TRUE)
assert_true(boundary_error, "Boundary validation did not reject y = 0 or y = 1.")

# 2. Small end-to-end fit, including K = 1 support.
sim <- simulate_default_scenario(
  T = 120L, K = 2L, diagonal = 0.90,
  phi_true = c(12, 28), seed = 202
)

fit2 <- fit_beta_gam_hmm(
  sim$y, sim$basis_spec,
  K = 2L, lambda = 0.1,
  n_start = 2L, max_iter = 60L, optimizer_maxit = 40L,
  seed = 303
)
assert_true(is.finite(fit2$loglik_obs), "Two-state fit returned a non-finite log-likelihood.")
assert_true(all(is.finite(fit2$edf_spline)), "Beta-specific EDF contains non-finite values.")
assert_true(all(fit2$edf_spline >= 0), "Beta-specific EDF contains negative values.")
assert_true(all(fit2$edf_spline <= ncol(sim$basis_spec$B) + 1e-8), "EDF exceeds the spline basis dimension.")
assert_true(isTRUE(fit2$monotonicity_ok), "GEM monotonicity safeguard failed in the smoke test.")
expected_bic <- -2 * fit2$loglik_obs + fit2$criteria$npar_edf * log(length(sim$y))
assert_true(abs(fit2$criteria$bic_edf - expected_bic) < 1e-8, "BIC does not use the intended log(T) penalty.")
assert_true(fit2$criteria$icl_edf + 1e-10 >= fit2$criteria$bic_edf, "ICL entropy penalty has the wrong sign.")

fit1 <- fit_beta_gam_hmm(
  sim$y, sim$basis_spec,
  K = 1L, lambda = 0.1,
  n_start = 1L, max_iter = 50L, optimizer_maxit = 40L,
  seed = 304
)
assert_true(fit1$K == 1L, "K = 1 model did not fit correctly.")
assert_true(n_parameters_hmm_edf(1L, fit1$edf_spline) >= 1, "K = 1 parameter count is invalid.")

# 3. Transition matrices and posterior probabilities must be normalized.
assert_true(max(abs(rowSums(fit2$A) - 1)) < 1e-8, "Transition matrix rows do not sum to one.")
assert_true(max(abs(rowSums(fit2$gamma) - 1)) < 1e-8, "Posterior state probabilities do not sum to one.")

# 4. Model selection must consider K = 1 rather than forcing regime switching.
grid <- run_grid_search_hmm(
  sim$y, sim$basis_spec,
  K_grid = 1:3, lambda_grid = c(0.05, 0.1),
  n_start = 1L, max_iter = 45L, optimizer_maxit = 35L,
  n_cores = 1L, n_start_cores = 1L,
  progress = FALSE, seed = 404
)
summary_grid <- summarize_grid_results(grid)
assert_true(all(1:3 %in% summary_grid$K), "Grid search did not evaluate all requested K values.")

# 5. Additive basis supports more than one covariate without changing the fitter.
set.seed(505)
dat <- data.frame(x1 = runif(100), x2 = rnorm(100))
add_basis <- build_additive_basis(dat, n_inner_knots = 3L)
assert_true(nrow(add_basis$B) == 100L, "Additive basis has the wrong number of rows.")
assert_true(nrow(add_basis$penalty) == ncol(add_basis$B), "Additive penalty has the wrong dimension.")

# 6. Forecast pseudo-residual diagnostics must be finite.
pr <- compute_pseudo_residuals(sim$y, fit2, sim$basis_spec)
pr_diag <- pseudo_residual_diagnostics(pr, lag = 10L)
assert_true(all(is.finite(unlist(pr_diag))), "Pseudo-residual diagnostics contain non-finite values.")

cat("All core invariant checks passed.\n")
