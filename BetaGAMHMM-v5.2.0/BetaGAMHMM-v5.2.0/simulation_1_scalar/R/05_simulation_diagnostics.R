# Performance diagnostics for simulated data sets with known generating values.

compare_fit_to_truth <- function(fit, sim, x_grid = NULL) {
  true_K <- sim$K_true %||% length(sim$phi_true)
  out <- list(
    selected_K = fit$K,
    true_K = true_K,
    phi_rmse = NA_real_,
    A_rmse = NA_real_,
    curve_rmse = NA_real_,
    pi_l1 = NA_real_,
    viterbi_accuracy = NA_real_
  )
  if (fit$K != true_K) return(out)

  fit_aligned <- align_fit_to_reference(
    fit,
    list(beta = sim$beta_true, phi = sim$phi_true, pi = sim$pi_true, A = sim$A_true),
    basis_spec = sim$basis_spec,
    x_grid = x_grid
  )

  if (identical(sim$basis_spec$type, "additive_bs")) {
    B_eval <- sim$basis_spec$B
  } else {
    x_grid <- x_grid %||% seq(min(sim$x), max(sim$x), length.out = 250L)
    B_eval <- predict_basis(x_grid, sim$basis_spec)
  }

  mu_true <- vapply(
    seq_len(true_K),
    function(k) inv_logit(drop(B_eval %*% sim$beta_true[[k]])),
    numeric(nrow(B_eval))
  )
  mu_hat <- vapply(
    seq_len(true_K),
    function(k) inv_logit(drop(B_eval %*% fit_aligned$beta[[k]])),
    numeric(nrow(B_eval))
  )
  if (true_K == 1L) {
    mu_true <- matrix(mu_true, ncol = 1L)
    mu_hat <- matrix(mu_hat, ncol = 1L)
  }

  z_hat <- viterbi_decode(sim$y, sim$basis_spec, fit_aligned)
  out$phi_rmse <- sqrt(mean((fit_aligned$phi - sim$phi_true)^2))
  out$A_rmse <- sqrt(mean((fit_aligned$A - sim$A_true)^2))
  out$curve_rmse <- sqrt(mean((mu_hat - mu_true)^2))
  out$pi_l1 <- sum(abs(fit_aligned$pi - sim$pi_true))
  out$viterbi_accuracy <- mean(z_hat == sim$z)
  out
}
