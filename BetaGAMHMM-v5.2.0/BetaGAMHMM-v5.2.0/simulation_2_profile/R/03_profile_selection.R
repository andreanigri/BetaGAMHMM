# Grid search and final model-selection utilities for profile HMMs.
#
# The manuscript uses one state-number criterion only:
#   BIC_nu = -2 log L + nu log(T).
# ICL_nu is retained only as a classification-oriented secondary comparison in
# simulation. Lambda is selected by AIC_nu within each K for both summaries.
# Precision/occupancy diagnostics are reported; they do not select K.

run_profile_grid <- function(panel,
                             basis_spec,
                             K_grid,
                             lambda_grid,
                             n_start = 3L,
                             n_cores = 1L,
                             max_iter = 150L,
                             tol = 1e-6,
                             phi_lower = 1,
                             phi_upper = 800,
                             optimizer_maxit = 80L,
                             precision_saturation_fraction = 0.98,
                             seed = 1L,
                             progress = FALSE) {
  grid <- expand.grid(K = K_grid, lambda = lambda_grid, KEEP.OUT.ATTRS = FALSE)
  out <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    if (progress) message(sprintf("Profile grid %d/%d: K=%d, lambda=%g", i, nrow(grid), grid$K[i], grid$lambda[i]))
    out[[i]] <- tryCatch(
      fit_beta_year_hmm(
        panel = panel, basis_spec = basis_spec,
        K = grid$K[i], lambda = grid$lambda[i],
        n_start = n_start, n_cores = n_cores,
        max_iter = max_iter, tol = tol,
        phi_lower = phi_lower, phi_upper = phi_upper,
        optimizer_maxit = optimizer_maxit,
        precision_saturation_fraction = precision_saturation_fraction,
        seed = seed + 10000L * grid$K[i] + 100L * i
      ),
      error = function(e) structure(
        list(K = grid$K[i], lambda = grid$lambda[i], error = conditionMessage(e)),
        class = "profile_grid_error"
      )
    )
  }
  names(out) <- paste0("K", grid$K, "_lam", grid$lambda)
  out
}

summarize_profile_grid <- function(results) {
  rows <- lapply(results, function(f) {
    if (inherits(f, "profile_grid_error")) {
      return(data.frame(
        K = f$K, lambda = f$lambda, converged = FALSE, monotonicity_ok = FALSE,
        precision_flag = NA, min_state_share = NA_real_, loglik = NA_real_,
        aic_edf = NA_real_, bic_edf = NA_real_, icl_edf = NA_real_,
        npar_edf = NA_real_, markov_par = NA_real_, emission_edf = NA_real_,
        phi = NA_character_, timing_seconds = NA_real_, error = f$error,
        stringsAsFactors = FALSE
      ))
    }
    cr <- f$criteria
    data.frame(
      K = f$K, lambda = f$lambda[1L], converged = f$converged,
      monotonicity_ok = f$monotonicity_ok,
      precision_flag = f$precision_diagnostic$flagged,
      min_state_share = f$occupancy$min_share,
      loglik = f$loglik_obs,
      aic_edf = cr$aic_edf,
      bic_edf = cr$bic_edf,
      icl_edf = cr$icl_edf,
      npar_edf = cr$npar_edf,
      markov_par = cr$markov_par,
      emission_edf = cr$emission_edf,
      phi = paste(formatC(f$phi, digits = 2, format = "f"), collapse = " | "),
      timing_seconds = f$timing_seconds,
      error = NA_character_, stringsAsFactors = FALSE
    )
  })
  ans <- do.call(rbind, rows)
  ans[order(ans$K, ans$lambda), , drop = FALSE]
}

profile_criterion_map <- function() c(BIC = "bic_edf", ICL = "icl_edf")

select_profile_model <- function(results, criterion = c("bic_edf", "icl_edf")) {
  criterion <- match.arg(criterion)
  tab <- summarize_profile_grid(results)
  # A model is admissible when the numerical fit succeeded and the safeguarded
  # GEM monotonicity check passed. Precision/occupancy flags are diagnostics and
  # are deliberately NOT model-selection filters.
  tab$admissible <- with(tab, is.na(error) & is.finite(aic_edf) & converged & monotonicity_ok)

  best_by_K <- do.call(rbind, lapply(split(tab, tab$K), function(dd) {
    ok <- dd[dd$admissible, , drop = FALSE]
    if (!nrow(ok)) return(NULL)
    ok[order(ok$aic_edf, ok$lambda), , drop = FALSE][1L, ]
  }))
  if (is.null(best_by_K) || !nrow(best_by_K)) stop("No admissible profile fit.")
  rownames(best_by_K) <- NULL
  selrow <- best_by_K[order(best_by_K[[criterion]], best_by_K$K), , drop = FALSE][1L, ]
  idx <- which(vapply(results, function(f) {
    !inherits(f, "profile_grid_error") && f$K == selrow$K && abs(f$lambda[1L] - selrow$lambda) < 1e-12
  }, logical(1)))
  if (!length(idx)) stop("Selected profile fit not found.")
  list(fit = results[[idx[1L]]], selected_row = selrow, best_by_K = best_by_K, grid_summary = tab)
}
