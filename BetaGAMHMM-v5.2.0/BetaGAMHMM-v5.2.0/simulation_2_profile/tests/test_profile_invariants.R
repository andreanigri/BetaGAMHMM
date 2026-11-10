message("Running profile-model smoke/invariant checks...")
sc <- profile_scenarios("quick")$baseline_T55
sc$T <- 25L; sc$J <- 15L; sc$K <- 2L; sc$phi <- c(120, 220); sc$weights <- c(0.5, 0.5)
sim <- simulate_profile_scenario(sc, seed = 101)
bs <- build_age_basis(sim$panel$ages, degree = 3L, n_inner_knots = 4L)
fit <- fit_beta_year_hmm(sim$panel, bs, K = 2L, lambda = 1,
  n_start = 2L, n_cores = 1L, max_iter = 50L,
  phi_lower = 1, phi_upper = 800, optimizer_maxit = 50L, seed = 202)
stopifnot(is.finite(fit$loglik_obs))
stopifnot(all(abs(rowSums(fit$gamma) - 1) < 1e-8))
stopifnot(all(abs(rowSums(fit$A) - 1) < 1e-8))
stopifnot(is.finite(fit$criteria$bic_edf))
stopifnot(is.finite(fit$criteria$icl_edf))
stopifnot(abs(fit$criteria$bic_edf - (-2 * fit$loglik_obs + log(sim$panel$n_years) * fit$criteria$npar_edf)) < 1e-8)
stopifnot(fit$criteria$icl_edf + 1e-10 >= fit$criteria$bic_edf)
stopifnot(fit$criteria$markov_par == 2L * 1L + 1L)
stopifnot(!any(c("bic_edf_cells", "bic_edf_block") %in% names(fit$criteria)))
message("Profile-model smoke/invariant checks passed.")
