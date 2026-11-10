# Grid search and final state-number selection for the repeated-profile Beta model.
#
# Workflow fixed before the final real-data rerun:
#   1) for each K, select lambda by minimum EDF-adjusted AIC;
#   2) select K by minimum EDF-adjusted BIC = -2 log L + nu log(T).
# Numerical diagnostics are reported, but no occupancy threshold, +2 BIC rule,
# or multi-BIC voting rule determines K.

run_beta_grid <- function(panel, basis_spec, K_grid, lambda_grid,
                          n_start = 15L, n_cores = 1L,
                          max_iter = 300L, tol = 1e-6,
                          phi_lower = 1, phi_upper = 2000,
                          optimizer_maxit = 150L,
                          precision_saturation_fraction = 0.98,
                          seed = 20261001L, progress = TRUE) {
  grid <- expand.grid(K = K_grid, lambda = lambda_grid, KEEP.OUT.ATTRS = FALSE)
  results <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    if (progress) message(sprintf("Beta profile grid %d/%d: K=%d, lambda=%g",
                                  i, nrow(grid), grid$K[i], grid$lambda[i]))
    results[[i]] <- tryCatch(
      fit_beta_year_hmm(
        panel = panel, basis_spec = basis_spec, K = grid$K[i], lambda = grid$lambda[i],
        n_start = n_start, n_cores = n_cores, max_iter = max_iter, tol = tol,
        phi_lower = phi_lower, phi_upper = phi_upper,
        optimizer_maxit = optimizer_maxit,
        precision_saturation_fraction = precision_saturation_fraction,
        seed = seed + 10000L * grid$K[i] + 100L * i
      ),
      error = function(e) structure(list(K = grid$K[i], lambda = grid$lambda[i],
                                         error = conditionMessage(e)), class = "grid_error")
    )
  }
  names(results) <- paste0("K", grid$K, "_lam", grid$lambda)
  results
}

summarize_beta_grid <- function(results) {
  rows <- lapply(results, function(fit) {
    if (inherits(fit, "grid_error")) {
      return(data.frame(
        K = fit$K, lambda = fit$lambda, converged = FALSE,
        monotonicity_ok = FALSE, precision_flag = NA,
        min_state_share = NA_real_, loglik = NA_real_,
        AIC_EDF = NA_real_, BIC_EDF = NA_real_, ICL_EDF = NA_real_,
        npar_edf = NA_real_, edf_spline = NA_real_, phi = NA_character_,
        timing_seconds = NA_real_, error = fit$error, stringsAsFactors = FALSE
      ))
    }
    cr <- fit$criteria
    data.frame(
      K = fit$K, lambda = fit$lambda[1L], converged = fit$converged,
      monotonicity_ok = fit$monotonicity_ok,
      precision_flag = fit$precision_diagnostic$flagged,
      min_state_share = fit$occupancy$min_share,
      loglik = fit$loglik_obs,
      AIC_EDF = cr$aic_edf,
      BIC_EDF = cr$bic_edf,
      ICL_EDF = cr$icl_edf,
      npar_edf = cr$npar_edf,
      edf_spline = sum(fit$edf_spline),
      phi = paste(formatC(fit$phi, digits = 2, format = "f"), collapse = " | "),
      timing_seconds = fit$timing_seconds,
      error = NA_character_, stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(out$K, out$lambda), , drop = FALSE]
}

select_beta_model <- function(results) {
  tab <- summarize_beta_grid(results)
  # State selection is not filtered by posterior occupancy or heuristic
  # precision-jump thresholds. Convergence and GEM monotonicity are numerical
  # validity requirements; the precision-bound flag is retained for reporting.
  tab$admissible <- with(tab,
    is.na(error) & converged & monotonicity_ok & is.finite(AIC_EDF) & is.finite(BIC_EDF)
  )
  best_by_K <- do.call(rbind, lapply(split(tab, tab$K), function(dd) {
    ok <- dd[dd$admissible, , drop = FALSE]
    if (!nrow(ok)) return(NULL)
    ok[order(ok$AIC_EDF, ok$lambda), , drop = FALSE][1L, ]
  }))
  if (is.null(best_by_K) || !nrow(best_by_K)) stop("No numerically valid Beta profile fit remains.")
  rownames(best_by_K) <- NULL
  selected_row <- best_by_K[order(best_by_K$BIC_EDF, best_by_K$K), , drop = FALSE][1L, ]
  idx <- which(vapply(results, function(f) {
    !inherits(f, "grid_error") && f$K == selected_row$K &&
      abs(f$lambda[1L] - selected_row$lambda) < 1e-12
  }, logical(1)))
  if (!length(idx)) stop("Internal error: selected Beta profile fit not found.")
  list(fit = results[[idx[1L]]], grid_summary = tab,
       best_by_K = best_by_K, selected_row = selected_row,
       criterion_label = "BIC", criterion_column = "BIC_EDF")
}

append_grid_results <- function(x, y) {
  out <- c(x, y)
  key <- vapply(out, function(f) {
    if (inherits(f, "grid_error")) paste(f$K, f$lambda, sep = "|")
    else paste(f$K, f$lambda[1L], sep = "|")
  }, character(1))
  out[!duplicated(key)]
}
