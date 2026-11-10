# Diagnostics for the selected year-level Beta model.

compute_beta_pseudo_residuals <- function(fit) {
  panel <- fit$panel
  Tn <- panel$n_years
  J <- panel$n_ages
  K <- fit$K
  pred_state <- matrix(NA_real_, Tn, K)
  pred_state[1L, ] <- fit$pi
  if (Tn >= 2L) {
    for (t in 2:Tn) pred_state[t, ] <- drop(fit$alpha[t - 1L, ] %*% fit$A)
  }

  u <- matrix(NA_real_, Tn, J)
  for (t in seq_len(Tn)) {
    for (j in seq_len(J)) {
      cdf_k <- vapply(seq_len(K), function(k) {
        stats::pbeta(
          panel$Y[t, j],
          shape1 = fit$mu_age[j, k] * fit$phi[k],
          shape2 = (1 - fit$mu_age[j, k]) * fit$phi[k]
        )
      }, numeric(1))
      u[t, j] <- sum(pred_state[t, ] * cdf_k)
    }
  }
  u <- clip_unit(u, 1e-10)
  q <- stats::qnorm(u)
  q_year <- rowMeans(q)
  acf_obj <- stats::acf(q_year, plot = FALSE, lag.max = min(20L, length(q_year) - 1L))
  max_abs_acf <- if (length(acf_obj$acf) > 1L) max(abs(acf_obj$acf[-1L])) else NA_real_
  lag_lb <- min(10L, floor(length(q_year) / 3L))
  lb_p <- if (lag_lb >= 1L) stats::Box.test(q_year, lag = lag_lb, type = "Ljung-Box")$p.value else NA_real_

  list(
    residual_matrix = q,
    year_mean = q_year,
    summary = data.frame(
      metric = c("mean_cell_residual", "sd_cell_residual", "mean_year_residual", "sd_year_residual", "max_abs_acf_year_mean", "ljung_box_p_year_mean"),
      value = c(mean(q), stats::sd(as.vector(q)), mean(q_year), stats::sd(q_year), max_abs_acf, lb_p)
    )
  )
}

save_beta_diagnostic_plots <- function(diag, years, output_dir) {
  safe_dir_create(output_dir)

  grDevices::pdf(file.path(output_dir, "pseudo_residual_qq.pdf"), width = 6.2, height = 5.2)
  stats::qqnorm(as.vector(diag$residual_matrix), main = "Beta one-step-ahead pseudo-residuals")
  stats::qqline(as.vector(diag$residual_matrix), lty = 2)
  grDevices::dev.off()

  grDevices::png(file.path(output_dir, "pseudo_residual_qq.png"), width = 1200, height = 950, res = 170)
  stats::qqnorm(as.vector(diag$residual_matrix), main = "Beta one-step-ahead pseudo-residuals")
  stats::qqline(as.vector(diag$residual_matrix), lty = 2)
  grDevices::dev.off()

  grDevices::pdf(file.path(output_dir, "year_mean_pseudo_residual_acf.pdf"), width = 6.2, height = 5.2)
  stats::acf(diag$year_mean, lag.max = min(20L, length(diag$year_mean) - 1L), main = "ACF of yearly mean pseudo-residuals")
  grDevices::dev.off()

  grDevices::png(file.path(output_dir, "year_mean_pseudo_residual_acf.png"), width = 1200, height = 950, res = 170)
  stats::acf(diag$year_mean, lag.max = min(20L, length(diag$year_mean) - 1L), main = "ACF of yearly mean pseudo-residuals")
  grDevices::dev.off()

  invisible(output_dir)
}

# Predictive standardized residuals for the denominator-aware binomial-kernel
# sensitivity. These residuals remain well-defined for HMD fractional estimated
# counts, unlike randomized discrete Binomial residuals that require integer
# support.
compute_binomial_kernel_residuals <- function(fit) {
  panel <- fit$panel
  F <- panel$F
  N <- panel$N
  Tn <- panel$n_years
  J <- panel$n_ages
  K <- fit$K

  pred_state <- matrix(NA_real_, Tn, K)
  pred_state[1L, ] <- fit$pi
  if (Tn >= 2L) {
    for (t in 2:Tn) pred_state[t, ] <- drop(fit$alpha[t - 1L, ] %*% fit$A)
  }

  residual <- matrix(NA_real_, Tn, J)
  for (t in seq_len(Tn)) {
    for (j in seq_len(J)) {
      p_k <- fit$mu_age[j, ]
      n <- N[t, j]
      mean_k <- n * p_k
      var_k <- n * p_k * (1 - p_k)
      pred_mean <- sum(pred_state[t, ] * mean_k)
      # Law of total variance under the predictive state mixture.
      pred_var <- sum(pred_state[t, ] * (var_k + mean_k^2)) - pred_mean^2
      residual[t, j] <- (F[t, j] - pred_mean) / sqrt(max(pred_var, 1e-10))
    }
  }

  year_mean <- rowMeans(residual, na.rm = TRUE)
  acf_obj <- stats::acf(year_mean, plot = FALSE, lag.max = min(20L, length(year_mean) - 1L), na.action = stats::na.pass)
  max_abs_acf <- if (length(acf_obj$acf) > 1L) max(abs(acf_obj$acf[-1L]), na.rm = TRUE) else NA_real_
  lag_lb <- min(10L, floor(length(year_mean) / 3L))
  lb_p <- if (lag_lb >= 1L) stats::Box.test(year_mean, lag = lag_lb, type = "Ljung-Box")$p.value else NA_real_

  list(
    residual_matrix = residual,
    year_mean = year_mean,
    summary = data.frame(
      metric = c(
        "mean_cell_standardized_residual", "sd_cell_standardized_residual",
        "mean_year_standardized_residual", "sd_year_standardized_residual",
        "max_abs_acf_year_mean", "ljung_box_p_year_mean"
      ),
      value = c(
        mean(residual, na.rm = TRUE), stats::sd(as.vector(residual), na.rm = TRUE),
        mean(year_mean, na.rm = TRUE), stats::sd(year_mean, na.rm = TRUE),
        max_abs_acf, lb_p
      ),
      stringsAsFactors = FALSE
    )
  )
}

save_binomial_kernel_diagnostic_plots <- function(diag, output_dir) {
  safe_dir_create(output_dir)
  vals <- as.vector(diag$residual_matrix)
  vals <- vals[is.finite(vals)]

  grDevices::pdf(file.path(output_dir, "binomial_kernel_standardized_residual_qq.pdf"), width = 6.2, height = 5.2)
  stats::qqnorm(vals, main = "Binomial-kernel predictive standardized residuals")
  stats::qqline(vals, lty = 2)
  grDevices::dev.off()

  grDevices::png(file.path(output_dir, "binomial_kernel_standardized_residual_qq.png"), width = 1200, height = 950, res = 170)
  stats::qqnorm(vals, main = "Binomial-kernel predictive standardized residuals")
  stats::qqline(vals, lty = 2)
  grDevices::dev.off()

  grDevices::pdf(file.path(output_dir, "binomial_kernel_year_mean_residual_acf.pdf"), width = 6.2, height = 5.2)
  stats::acf(diag$year_mean, lag.max = min(20L, length(diag$year_mean) - 1L), main = "ACF of binomial-kernel yearly mean residuals")
  grDevices::dev.off()

  grDevices::png(file.path(output_dir, "binomial_kernel_year_mean_residual_acf.png"), width = 1200, height = 950, res = 170)
  stats::acf(diag$year_mean, lag.max = min(20L, length(diag$year_mean) - 1L), main = "ACF of binomial-kernel yearly mean residuals")
  grDevices::dev.off()
  invisible(output_dir)
}
