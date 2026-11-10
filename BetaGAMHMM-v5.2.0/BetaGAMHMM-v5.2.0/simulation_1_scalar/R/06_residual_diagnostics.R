# Forecast pseudo-residual diagnostics used in the final revision study.
#
# These functions evaluate one-step-ahead forecast residual calibration and
# residual serial dependence for a fitted Beta-GAM-HMM.

compute_pseudo_residuals <- function(y, fit, basis_spec = fit$basis_spec) {
  y <- prepare_unit_interval_response(y, boundary_handling = "error")
  T <- length(y)
  K <- fit$K
  B <- basis_spec$B
  
  mu <- sapply(seq_len(K), function(k) inv_logit(drop(B %*% fit$beta[[k]])))
  if (is.null(dim(mu))) mu <- matrix(mu, nrow = T, ncol = K)
  
  log_A <- log(pmax(fit$A, 1e-16))
  log_pi <- log(pmax(fit$pi, 1e-16))
  
  # forecast probabilities and CDF
  u <- numeric(T)
  log_alpha <- numeric(K)
  
  # t = 1: forecast state probs are the stationary distribution
  forecast_prob <- fit$pi
  cdf_components <- sapply(seq_len(K), function(k) {
    stats::pbeta(y[1], mu[1, k] * fit$phi[k], (1 - mu[1, k]) * fit$phi[k])
  })
  u[1] <- sum(forecast_prob * cdf_components)
  
  # initialise log-alpha after observing y[1]
  log_emission_1 <- sapply(seq_len(K), function(k) {
    log_beta_density(y[1], mu[1, k], fit$phi[k])
  })
  log_alpha <- log_pi + log_emission_1
  
  for (t in 2:T) {
    # forecast: P(S_t = k | y_{1:t-1}) via prediction step
    log_forecast <- numeric(K)
    for (k in seq_len(K)) {
      log_forecast[k] <- log_sum_exp(log_alpha + log_A[, k])
    }
    # normalise to get forecast state probabilities
    log_forecast <- log_forecast - log_sum_exp(log_forecast)
    forecast_prob <- exp(log_forecast)
    
    # CDF of y_t under each state
    cdf_k <- sapply(seq_len(K), function(k) {
      stats::pbeta(y[t], mu[t, k] * fit$phi[k], (1 - mu[t, k]) * fit$phi[k])
    })
    u[t] <- sum(forecast_prob * cdf_k)
    
    # update: incorporate y_t
    log_emission_t <- sapply(seq_len(K), function(k) {
      log_beta_density(y[t], mu[t, k], fit$phi[k])
    })
    log_alpha <- log_forecast + log_emission_t
  }
  
  # clip to avoid Inf in qnorm
  u <- pmin(pmax(u, 1e-10), 1 - 1e-10)
  residuals <- stats::qnorm(u)
  
  list(
    u = u,
    residuals = residuals,
    T = T
  )
}

pseudo_residual_diagnostics <- function(pr, lag = NULL) {
  r <- as.numeric(pr$residuals)
  n <- length(r)
  lag <- lag %||% max(1L, min(20L, floor(sqrt(n))))
  lb <- stats::Box.test(r, lag = lag, type = "Ljung-Box")
  acf_obj <- stats::acf(r, lag.max = lag, plot = FALSE)$acf
  data.frame(
    n = n,
    mean = mean(r),
    sd = stats::sd(r),
    max_abs_acf = max(abs(as.numeric(acf_obj[-1L])), na.rm = TRUE),
    ljung_box_lag = lag,
    ljung_box_statistic = unname(lb$statistic),
    ljung_box_p_value = lb$p.value,
    stringsAsFactors = FALSE
  )
}

plot_pseudo_residuals <- function(pr) {
  check_required_packages("ggplot2")
  n <- pr$T
  res <- pr$residuals
  
  p_acf <- {
    acf_obj <- stats::acf(res, lag.max = min(30L, n - 1L), plot = FALSE)
    acf_df <- data.frame(
      lag = as.numeric(acf_obj$lag[-1]),
      acf = as.numeric(acf_obj$acf[-1])
    )
    ci <- stats::qnorm(0.975) / sqrt(n)
    ggplot2::ggplot(acf_df, ggplot2::aes(x = lag, y = acf)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey50") +
      ggplot2::geom_hline(yintercept = c(-ci, ci), linetype = 2, colour = "blue") +
      ggplot2::geom_segment(ggplot2::aes(xend = lag, yend = 0), linewidth = 0.7) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::labs(title = "ACF of pseudo-residuals", x = "Lag", y = "ACF")
  }
  
  p_qq <- {
    qq <- stats::qqnorm(res, plot.it = FALSE)
    qq_df <- data.frame(theoretical = qq$x, sample = qq$y)
    ggplot2::ggplot(qq_df, ggplot2::aes(x = theoretical, y = sample)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "red") +
      ggplot2::geom_point(size = 0.8, alpha = 0.5) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::labs(title = "QQ-plot of pseudo-residuals", x = "Theoretical quantiles", y = "Sample quantiles")
  }
  
  p_hist <- {
    ggplot2::ggplot(data.frame(r = res), ggplot2::aes(x = r)) +
      ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                              bins = 30, fill = "grey70", colour = "grey30") +
      ggplot2::stat_function(fun = stats::dnorm, colour = "red", linewidth = 0.8) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::labs(title = "Histogram of pseudo-residuals", x = "Pseudo-residual", y = "Density")
  }
  
  list(acf = p_acf, qq = p_qq, hist = p_hist)
}

save_plot <- function(plot_obj, file, width = 8, height = 5, dpi = 300) {
  check_required_packages("ggplot2")
  ggplot2::ggsave(
    filename = file, plot = plot_obj, width = width, height = height, dpi = dpi
  )
  invisible(file)
}
