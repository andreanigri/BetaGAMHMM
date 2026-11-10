# Model fitting over (K, lambda) grids and model-selection utilities.

run_grid_search_hmm <- function(y,
                                basis_spec,
                                K_grid,
                                lambda_grid,
                                n_start = 8L,
                                max_iter = 200L,
                                tol = 1e-6,
                                phi_lower = 1,
                                phi_upper = 500,
                                optimizer_maxit = 100L,
                                mstep_tol = 1e-8,
                                monotone_tol = 1e-8,
                                boundary_handling = "error",
                                store_all_starts = FALSE,
                                verbose = FALSE,
                                diagnostic = list(
                                  phi_max = 500,
                                  eps = 1e-6,
                                  s_thresh = 2L,
                                  delta_abs = 50,
                                  delta_sum = 100
                                ),
                                n_cores = 1L,
                                n_start_cores = 1L,
                                progress = interactive(),
                                seed = NULL) {
  grid <- expand.grid(K = as.integer(K_grid), lambda = as.numeric(lambda_grid), KEEP.OUT.ATTRS = FALSE)
  grid$job_id <- seq_len(nrow(grid))

  fit_one <- function(i) {
    K_i <- grid$K[i]
    lambda_i <- grid$lambda[i]
    seed_i <- model_candidate_seed(seed, K_i, lambda_i)
    tryCatch({
      fit <- fit_beta_gam_hmm(
        y = y,
        basis_spec = basis_spec,
        K = K_i,
        lambda = lambda_i,
        n_start = n_start,
        n_start_cores = n_start_cores,
        max_iter = max_iter,
        tol = tol,
        phi_lower = phi_lower,
        phi_upper = phi_upper,
        optimizer_maxit = optimizer_maxit,
        mstep_tol = mstep_tol,
        monotone_tol = monotone_tol,
        boundary_handling = boundary_handling,
        store_all_starts = store_all_starts,
        verbose = verbose,
        diagnostic = diagnostic,
        seed = seed_i
      )
      list(K = K_i, lambda = lambda_i, fit = fit, error = NA_character_)
    }, error = function(e) {
      list(K = K_i, lambda = lambda_i, fit = NULL, error = conditionMessage(e))
    })
  }

  idx <- seq_len(nrow(grid))
  parallel_lapply_crossplatform(
    X = idx,
    FUN = fit_one,
    n_cores = min(as.integer(n_cores %||% 1L), length(idx)),
    exports = c(
      "fit_one", "grid", "seed", "y", "basis_spec", "n_start", "max_iter", "tol",
      "phi_lower", "phi_upper", "optimizer_maxit", "mstep_tol", "monotone_tol",
      "boundary_handling", "store_all_starts", "verbose", "diagnostic", "n_start_cores"
    ),
    project_root = getwd(),
    seed = seed,
    progress = progress,
    progress_label = sprintf("Grid search: %d candidate (K, lambda) models", length(idx))
  )
}

summarize_grid_results <- function(results, diagnostic_override = NULL) {
  rows <- lapply(results, function(res) {
    if (is.null(res$fit)) {
      return(data.frame(
        K = res$K, lambda = res$lambda, lambda_fit = NA_character_,
        converged = FALSE, monotonicity_ok = FALSE,
        loglik_obs = NA_real_, objective_pen = NA_real_,
        aic = NA_real_, bic = NA_real_, icl = NA_real_,
        aic_edf = NA_real_, bic_edf = NA_real_, icl_edf = NA_real_,
        entropy = NA_real_, npar = NA_real_, npar_edf = NA_real_,
        edf_total = NA_real_, edf_by_state = NA_character_,
        degenerate = NA, n_sat = NA_integer_, delta_tail = NA_real_, large_jump = NA,
        min_state_share = NA_real_, min_state_mass = NA_real_, state_share = NA_character_,
        occupancy_flag = NA, phi_sorted = NA_character_,
        best_start_index = NA_integer_, n_successful_starts = NA_integer_, n_failed_starts = NA_integer_,
        n_converged_starts = NA_integer_, n_monotone_starts = NA_integer_, n_degenerate_starts = NA_integer_,
        start_loglik_gap = NA_real_, n_iter = NA_integer_,
        mstep_rejections = NA_integer_, timing_seconds = NA_real_,
        error = res$error, stringsAsFactors = FALSE
      ))
    }

    fit <- res$fit
    occ <- fit$occupancy
    edf <- fit$edf_spline
    diag_i <- if (is.null(diagnostic_override)) {
      fit$diagnostic
    } else {
      do.call(diagnose_phi_degeneracy, c(list(phi = fit$phi), diagnostic_override))
    }
    data.frame(
      K = res$K,
      lambda = res$lambda,
      lambda_fit = format_lambda(fit$lambda),
      converged = isTRUE(fit$converged),
      monotonicity_ok = isTRUE(fit$monotonicity_ok %||% TRUE),
      loglik_obs = fit$loglik_obs,
      objective_pen = fit$objective_pen,
      aic = fit$criteria$aic,
      bic = fit$criteria$bic,
      icl = fit$criteria$icl,
      aic_edf = fit$criteria$aic_edf,
      bic_edf = fit$criteria$bic_edf,
      icl_edf = fit$criteria$icl_edf,
      entropy = fit$criteria$entropy,
      npar = fit$criteria$npar,
      npar_edf = fit$criteria$npar_edf,
      edf_total = sum(edf),
      edf_by_state = format_phi(edf),
      degenerate = isTRUE(diag_i$flagged),
      n_sat = diag_i$n_sat,
      delta_tail = diag_i$delta_tail,
      large_jump = isTRUE(diag_i$large_jump),
      min_state_share = occ$min_share,
      min_state_mass = occ$min_mass,
      state_share = format_phi(occ$share),
      occupancy_flag = FALSE,
      phi_sorted = format_phi(diag_i$phi_sorted),
      best_start_index = fit$best_start_index,
      n_successful_starts = fit$start_diagnostics$n_successful_starts %||% NA_integer_,
      n_failed_starts = fit$start_diagnostics$n_failed_starts %||% NA_integer_,
      n_converged_starts = fit$start_diagnostics$n_converged,
      n_monotone_starts = fit$start_diagnostics$n_monotone %||% NA_integer_,
      n_degenerate_starts = fit$start_diagnostics$n_degenerate,
      start_loglik_gap = fit$start_diagnostics$loglik_gap_best_vs_median,
      n_iter = fit$n_iter,
      mstep_rejections = fit$mstep_rejections %||% NA_integer_,
      timing_seconds = fit$timing_seconds %||% NA_real_,
      error = NA_character_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(out$K, out$lambda), , drop = FALSE]
}

apply_selection_filters <- function(summary_df,
                                    exclude_degenerate = FALSE,
                                    require_converged = TRUE,
                                    require_monotone = TRUE,
                                    enforce_occupancy = FALSE,
                                    min_state_share = 0,
                                    min_state_mass = 0,
                                    max_start_loglik_gap = Inf) {
  out <- summary_df
  share_thr <- min_state_share
  mass_thr <- min_state_mass
  out$occupancy_flag <- with(
    out,
    ifelse(is.na(min_state_share) | is.na(min_state_mass), TRUE,
           (min_state_share < share_thr) | (min_state_mass < mass_thr))
  )
  keep <- !is.na(out$loglik_obs) & is.na(out$error)
  if (exclude_degenerate) keep <- keep & !out$degenerate
  if (require_converged) keep <- keep & out$converged
  if (require_monotone && "monotonicity_ok" %in% names(out)) keep <- keep & out$monotonicity_ok
  if (enforce_occupancy) keep <- keep & !out$occupancy_flag
  if (is.finite(max_start_loglik_gap)) keep <- keep & (out$start_loglik_gap <= max_start_loglik_gap)
  out$passes_filter <- keep
  out
}

rank_candidates <- function(df,
                            criterion_primary = "bic_edf",
                            criterion_secondary = "icl_edf",
                            equivalence_margin = 0) {
  if (nrow(df) == 0L) return(df)
  good <- is.finite(df[[criterion_primary]])
  df <- df[good, , drop = FALSE]
  if (nrow(df) == 0L) return(df)
  secondary <- if (criterion_secondary %in% names(df)) df[[criterion_secondary]] else rep(0, nrow(df))
  secondary[!is.finite(secondary)] <- Inf
  df <- df[order(df[[criterion_primary]], secondary, -df$n_converged_starts,
                 df$start_loglik_gap, df$lambda), , drop = FALSE]

  if (equivalence_margin > 0 && nrow(df) > 1L) {
    best_val <- df[[criterion_primary]][1L]
    near <- df[df[[criterion_primary]] <= best_val + equivalence_margin, , drop = FALSE]
    if (nrow(near) > 1L) {
      near_secondary <- near[[criterion_secondary]]
      near_secondary[!is.finite(near_secondary)] <- Inf
      near <- near[order(-near$K, near_secondary, -near$n_converged_starts,
                         near$start_loglik_gap, near$lambda), , drop = FALSE]
      rest_keys <- paste(df$K, df$lambda)
      near_keys <- paste(near$K, near$lambda)
      rest <- df[!rest_keys %in% near_keys, , drop = FALSE]
      df <- rbind(near, rest)
    }
  }
  df
}

select_best_model_hierarchical <- function(results,
                                           criterion_primary = "bic_edf",
                                           criterion_secondary = "icl_edf",
                                           lambda_criterion_primary = "aic_edf",
                                           lambda_criterion_secondary = "bic_edf",
                                           exclude_degenerate = FALSE,
                                           require_converged = TRUE,
                                           require_monotone = TRUE,
                                           enforce_occupancy = FALSE,
                                           min_state_share = 0,
                                           min_state_mass = 0,
                                           max_start_loglik_gap = Inf,
                                           equivalence_margin = 0,
                                           prefer_richer_K_within_margin = FALSE,
                                           diagnostic_override = NULL) {
  summary_df <- summarize_grid_results(results, diagnostic_override = diagnostic_override)
  filtered <- apply_selection_filters(
    summary_df,
    exclude_degenerate = exclude_degenerate,
    require_converged = require_converged,
    require_monotone = require_monotone,
    enforce_occupancy = enforce_occupancy,
    min_state_share = min_state_share,
    min_state_mass = min_state_mass,
    max_start_loglik_gap = max_start_loglik_gap
  )
  candidates <- filtered[filtered$passes_filter, , drop = FALSE]
  if (nrow(candidates) == 0L) stop("No valid candidate models after filtering.")

  by_k_list <- lapply(split(candidates, candidates$K), function(dd) {
    ranked <- rank_candidates(dd, lambda_criterion_primary, lambda_criterion_secondary)
    if (nrow(ranked) == 0L) NULL else ranked[1L, , drop = FALSE]
  })
  by_k <- do.call(rbind, Filter(Negate(is.null), by_k_list))
  if (is.null(by_k) || nrow(by_k) == 0L) stop("No finite information criterion values after lambda selection.")
  by_k <- by_k[order(by_k$K), , drop = FALSE]

  ranked <- rank_candidates(
    by_k,
    criterion_primary = criterion_primary,
    criterion_secondary = criterion_secondary,
    equivalence_margin = if (prefer_richer_K_within_margin) equivalence_margin else 0
  )
  if (nrow(ranked) == 0L) stop("No finite information criterion values for K selection.")
  winner <- ranked[1L, , drop = FALSE]

  idx <- which(vapply(results, function(r) {
    isTRUE(all.equal(as.numeric(r$K), as.numeric(winner$K[1L]))) &&
      isTRUE(all.equal(as.numeric(r$lambda), as.numeric(winner$lambda[1L]))) &&
      !is.null(r$fit)
  }, logical(1)))[1L]

  list(
    selected = results[[idx]],
    summary = filtered,
    candidate_by_K = by_k,
    ranked_candidates = ranked,
    criterion_primary = criterion_primary,
    criterion_secondary = criterion_secondary,
    lambda_criterion_primary = lambda_criterion_primary,
    lambda_criterion_secondary = lambda_criterion_secondary,
    winner_row = winner
  )
}

select_best_model <- function(results,
                              criterion = c("bic_edf", "icl_edf", "bic", "icl", "aic", "aic_edf"),
                              exclude_degenerate = FALSE,
                              require_converged = TRUE,
                              prefer_smallest_K = FALSE,
                              enforce_occupancy = FALSE,
                              min_state_share = 0,
                              min_state_mass = 0,
                              max_start_loglik_gap = Inf,
                              equivalence_margin = 0,
                              diagnostic_override = NULL) {
  criterion <- match.arg(criterion)
  secondary <- if (criterion == "bic_edf") "icl_edf" else if (criterion == "bic") "icl" else criterion
  out <- select_best_model_hierarchical(
    results,
    criterion_primary = criterion,
    criterion_secondary = secondary,
    lambda_criterion_primary = "aic_edf",
    lambda_criterion_secondary = criterion,
    exclude_degenerate = exclude_degenerate,
    require_converged = require_converged,
    enforce_occupancy = enforce_occupancy,
    min_state_share = min_state_share,
    min_state_mass = min_state_mass,
    max_start_loglik_gap = max_start_loglik_gap,
    equivalence_margin = if (prefer_smallest_K) 0 else equivalence_margin,
    prefer_richer_K_within_margin = !prefer_smallest_K,
    diagnostic_override = diagnostic_override
  )
  out$criterion <- criterion
  out
}

extract_k_diagnostic_table <- function(results, K_target) {
  summary_df <- summarize_grid_results(results)
  summary_df[summary_df$K == K_target, , drop = FALSE]
}

diagnostic_sensitivity_table <- function(results, diagnostic_rules) {
  rows <- lapply(seq_along(results), function(i) {
    rr <- results[[i]]
    if (is.null(rr$fit)) return(NULL)
    do.call(rbind, lapply(names(diagnostic_rules), function(rule_name) {
      rule <- diagnostic_rules[[rule_name]]
      diag_args <- rule$diagnostic %||% rule
      diag_i <- do.call(diagnose_phi_degeneracy, c(list(phi = rr$fit$phi), diag_args))
      data.frame(
        K = rr$K,
        lambda = rr$lambda,
        rule = rule_name,
        flagged = diag_i$flagged,
        n_sat = diag_i$n_sat,
        delta_tail = diag_i$delta_tail,
        large_jump = diag_i$large_jump,
        phi_sorted = format_phi(diag_i$phi_sorted),
        stringsAsFactors = FALSE
      )
    }))
  })
  merge_data_frames(rows)
}

# Apply several selection rules to the same fitted grid.  This makes sensitivity
# analysis inexpensive: changing diagnostic/occupancy thresholds does not require
# refitting the HMMs.
evaluate_selection_rules <- function(results, rules, K_true = NULL) {
  rows <- lapply(names(rules), function(rule_name) {
    r <- rules[[rule_name]]
    sel <- tryCatch(
      select_best_model_hierarchical(
        results,
        criterion_primary = r$criterion_primary %||% "bic_edf",
        criterion_secondary = r$criterion_secondary %||% "icl_edf",
        lambda_criterion_primary = r$lambda_criterion_primary %||% "aic_edf",
        lambda_criterion_secondary = r$lambda_criterion_secondary %||% "bic_edf",
        exclude_degenerate = r$exclude_degenerate %||% FALSE,
        require_converged = r$require_converged %||% TRUE,
        require_monotone = r$require_monotone %||% TRUE,
        enforce_occupancy = r$enforce_occupancy %||% FALSE,
        min_state_share = r$min_state_share %||% 0,
        min_state_mass = r$min_state_mass %||% 0,
        max_start_loglik_gap = r$max_start_loglik_gap %||% Inf,
        equivalence_margin = r$equivalence_margin %||% 0,
        prefer_richer_K_within_margin = r$prefer_richer_K_within_margin %||% FALSE,
        diagnostic_override = r$diagnostic %||% NULL
      ),
      error = function(e) NULL
    )
    if (is.null(sel)) {
      return(data.frame(
        rule = rule_name, selected_K = NA_integer_, selected_lambda = NA_real_,
        selected_true_K = NA, criterion_value = NA_real_, stringsAsFactors = FALSE
      ))
    }
    w <- sel$winner_row
    data.frame(
      rule = rule_name,
      selected_K = as.integer(w$K[1L]),
      selected_lambda = as.numeric(w$lambda[1L]),
      selected_true_K = if (is.null(K_true)) NA else as.integer(w$K[1L]) == K_true,
      criterion_value = as.numeric(w[[r$criterion_primary %||% "bic_edf"]][1L]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

selection_diagnostic_report <- function(results, K_true = NULL, selection_cfg = NULL) {
  summary_df <- summarize_grid_results(results)
  selection_cfg <- selection_cfg %||% list()
  filtered <- apply_selection_filters(
    summary_df,
    exclude_degenerate = selection_cfg$exclude_degenerate %||% FALSE,
    require_converged = selection_cfg$require_converged %||% TRUE,
    require_monotone = selection_cfg$require_monotone %||% TRUE,
    enforce_occupancy = selection_cfg$enforce_occupancy %||% FALSE,
    min_state_share = selection_cfg$min_state_share %||% 0,
    min_state_mass = selection_cfg$min_state_mass %||% 0,
    max_start_loglik_gap = selection_cfg$max_start_loglik_gap %||% Inf
  )
  filtered$converged_num <- as.integer(filtered$converged %in% TRUE)
  filtered$nondeg_num <- as.integer(!is.na(filtered$degenerate) & !filtered$degenerate)
  filtered$occupancy_ok_num <- as.integer(!is.na(filtered$occupancy_flag) & !filtered$occupancy_flag)
  filtered$passed_num <- as.integer(filtered$passes_filter %in% TRUE)
  filtered$valid_num <- as.integer(is.na(filtered$error))
  by_K <- aggregate(
    cbind(converged = filtered$converged_num,
          nondeg = filtered$nondeg_num,
          occupancy_ok = filtered$occupancy_ok_num,
          passed = filtered$passed_num,
          valid = filtered$valid_num) ~ K,
    data = filtered,
    FUN = sum
  )
  by_K$total <- as.numeric(table(filtered$K)[as.character(by_K$K)])
  if (!is.null(K_true)) by_K$is_true_K <- by_K$K == K_true
  list(summary = filtered, by_K = by_K)
}
