# Plot/table helpers shared by pre-final and final reporting.

make_table4 <- function(fit, boot) {
  rows <- list(); r <- 1L
  for (k in seq_len(fit$K)) {
    rows[[r]] <- data.frame(
      parameter = paste0("phi", k), estimate = fit$phi[k],
      lower = boot$phi_ci[k, 1L], upper = boot$phi_ci[k, 2L],
      interval = "BC percentile", stringsAsFactors = FALSE
    ); r <- r + 1L
  }
  if (fit$K > 1L) {
    for (i in seq_len(fit$K)) for (j in seq_len(fit$K)) {
      rows[[r]] <- data.frame(
        parameter = paste0("A", i, j), estimate = fit$A[i, j],
        lower = boot$A_ci_low[i, j], upper = boot$A_ci_high[i, j],
        interval = "percentile", stringsAsFactors = FALSE
      ); r <- r + 1L
    }
  }
  do.call(rbind, rows)
}

plot_state_profiles <- function(fit, boot, panel, pdf_path, png_path) {
  ages <- panel$ages; K <- fit$K
  ncol_panel <- if (K <= 2L) K else 2L
  nrow_panel <- ceiling(K / ncol_panel)
  draw <- function() {
    old <- graphics::par(mfrow = c(nrow_panel, ncol_panel), mar = c(4.1, 4.3, 2.3, 1.0))
    on.exit(graphics::par(old), add = TRUE)
    for (k in seq_len(K)) {
      ylim <- range(c(panel$Y_raw, boot$mu_low[, k], boot$mu_high[, k]), finite = TRUE)
      ylim <- c(max(0, ylim[1L]), min(1, ylim[2L]))
      graphics::plot(ages, fit$mu_age[, k], type = "n", ylim = ylim,
                     xlab = "Age", ylab = "Female-to-total death ratio",
                     main = paste("State", k))
      # Reviewer-facing display: observed profiles are shown only in the
      # panel corresponding to their Viterbi-decoded year-level state.
      year_idx <- which(fit$viterbi == k)
      for (t in year_idx) {
        graphics::points(ages, panel$Y_raw[t, ], pch = 16, cex = 0.24,
                         col = grDevices::adjustcolor("grey35", alpha.f = 0.16))
      }
      graphics::polygon(c(ages, rev(ages)), c(boot$mu_low[, k], rev(boot$mu_high[, k])),
                        border = NA, col = grDevices::adjustcolor("grey55", alpha.f = 0.25))
      graphics::lines(ages, fit$mu_age[, k], lwd = 2.2)
    }
  }
  grDevices::pdf(pdf_path, width = 8.5, height = max(4.5, 3.7 * nrow_panel)); draw(); grDevices::dev.off()
  grDevices::png(png_path, width = 1600, height = max(900, 700 * nrow_panel), res = 180); draw(); grDevices::dev.off()
}

plot_viterbi_years <- function(fit, panel, pdf_path, png_path) {
  draw <- function() {
    oldpar <- graphics::par(mar = c(4.2, 6.2, 3.0, 1.0)); on.exit(graphics::par(oldpar), add = TRUE)
    graphics::plot(panel$years, fit$viterbi, type = "s", lwd = 2,
                   xlab = "Calendar year", ylab = "", yaxt = "n",
                   ylim = c(0.75, fit$K + 0.25), main = "Year-level latent regime sequence")
    graphics::axis(2, at = seq_len(fit$K), labels = paste("State", seq_len(fit$K)), las = 1)
    graphics::mtext("Viterbi-decoded state", side = 2, line = 4.7)
    graphics::points(panel$years, fit$viterbi, pch = 16, cex = 0.55)
  }
  grDevices::pdf(pdf_path, width = 8.2, height = 4.6); draw(); grDevices::dev.off()
  grDevices::png(png_path, width = 1550, height = 850, res = 180); draw(); grDevices::dev.off()
}

compare_beta_binomial_extended <- function(beta_selection, binomial_selection) {
  bf <- beta_selection$fit; nf <- binomial_selection$fit
  ari <- adjusted_rand_index(bf$viterbi, nf$viterbi)
  exact <- curve_rmse <- NA_real_
  if (bf$K == nf$K) {
    nf2 <- align_fit_to_reference(nf, bf)
    exact <- mean(bf$viterbi == nf2$viterbi)
    curve_rmse <- sqrt(mean((bf$mu_age - nf2$mu_age)^2))
  }
  data.frame(
    beta_K = bf$K, beta_lambda = bf$lambda[1L],
    denominator_K = nf$K, denominator_lambda = nf$lambda[1L],
    viterbi_ARI = ari, exact_yearly_state_agreement = exact,
    fitted_curve_RMSE = curve_rmse,
    same_HMD_cells = TRUE, same_HMD_death_estimates = TRUE,
    note = "Same HMD cells and death estimates; raw information-criterion values are not compared across emission families.",
    stringsAsFactors = FALSE
  )
}
