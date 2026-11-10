# Optional comparison with hmmTMB (Michelot, 2025).
#
# This file is not required by the core implementation.  It provides a transparent
# external benchmark when the CRAN package `hmmTMB` is installed.  The comparator
# fits a Beta HMM in which both shape parameters may vary smoothly with x, making
# it at least as flexible as the mean-only smooth specification used in the paper.

hmmtmb_available <- function() {
  requireNamespace("hmmTMB", quietly = TRUE)
}

best_label_accuracy <- function(reference, estimated, K) {
  perms <- all_permutations(seq_len(K))
  acc <- apply(perms, 1L, function(p) {
    mapped <- p[estimated]
    mean(mapped == reference)
  })
  max(acc)
}

fit_hmmtmb_beta_gam <- function(y,
                                x,
                                K,
                                basis_k = 8L,
                                silent = TRUE) {
  if (!hmmtmb_available()) {
    stop("Package 'hmmTMB' is not installed. Run install.packages('hmmTMB') to use this benchmark.")
  }
  y <- prepare_unit_interval_response(y, boundary_handling = "error")
  dat <- data.frame(ID = 1L, y = y, x = as.numeric(x))

  # Stable positive starting values spanning the observed response distribution.
  mu0 <- as.numeric(stats::quantile(y, probs = seq(0.2, 0.8, length.out = K), names = FALSE))
  phi0 <- rep(20, K)
  par0 <- list(y = list(
    shape1 = pmax(mu0 * phi0, 0.2),
    shape2 = pmax((1 - mu0) * phi0, 0.2)
  ))
  smooth_formula <- stats::as.formula(sprintf("~ s(x, k = %d, bs = 'cs')", as.integer(basis_k)))
  formulas <- list(y = list(shape1 = smooth_formula, shape2 = smooth_formula))

  started <- proc.time()[[3L]]
  hid <- hmmTMB::MarkovChain$new(
    data = dat,
    n_states = K,
    initial_state = "estimated"
  )
  obs <- hmmTMB::Observation$new(
    data = dat,
    dists = list(y = "beta"),
    n_states = K,
    par = par0,
    formulas = formulas
  )
  hmm <- hmmTMB::HMM$new(obs = obs, hid = hid)
  hmm$fit(silent = silent)

  list(
    model = hmm,
    viterbi = as.integer(hmm$viterbi()),
    aic_conditional = tryCatch(as.numeric(stats::AIC(hmm)), error = function(e) NA_real_),
    aic_marginal = tryCatch(as.numeric(hmm$AIC_marginal()), error = function(e) NA_real_),
    elapsed_seconds = proc.time()[[3L]] - started
  )
}

run_hmmtmb_comparison_study <- function(n_rep = 20L,
                                         T = 700L,
                                         K = 4L,
                                         basis_k = 8L,
                                         lambda_grid = c(0.05, 0.1, 0.5),
                                         n_start_ours = 2L,
                                         n_cores = 1L,
                                         seed = 7400L) {
  if (!hmmtmb_available()) {
    return(data.frame(
      status = "skipped",
      reason = "hmmTMB is not installed",
      stringsAsFactors = FALSE
    ))
  }

  run_one <- function(r) {
    sim <- simulate_default_scenario(
      T = T, K = K, diagonal = 0.90,
      seed = seed + r * 100L
    )

    # Proposed method: K is held fixed so that this comparison focuses on
    # estimation and state decoding rather than on state-number selection.
    t0 <- proc.time()[[3L]]
    grid <- tryCatch(
      run_grid_search_hmm(
        sim$y, sim$basis_spec,
        K_grid = K, lambda_grid = lambda_grid,
        n_start = n_start_ours, max_iter = 120L,
        optimizer_maxit = 70L, n_cores = 1L,
        n_start_cores = 1L, progress = FALSE,
        seed = seed + r * 100L + 1L
      ),
      error = function(e) NULL
    )
    ours_time <- proc.time()[[3L]] - t0

    ours_row <- data.frame(
      rep = r, method = "Beta-GAM-HMM",
      converged = FALSE, viterbi_accuracy = NA_real_,
      information_criterion = NA_real_, elapsed_seconds = ours_time,
      stringsAsFactors = FALSE
    )
    if (!is.null(grid)) {
      sel <- tryCatch(
        select_best_model_hierarchical(
          grid,
          criterion_primary = "bic_edf",
          criterion_secondary = "icl_edf",
          lambda_criterion_primary = "aic_edf",
          lambda_criterion_secondary = "bic_edf",
          exclude_degenerate = FALSE,
          require_converged = TRUE,
          require_monotone = TRUE,
          enforce_occupancy = FALSE,
          equivalence_margin = 0,
          prefer_richer_K_within_margin = FALSE
        ),
        error = function(e) NULL
      )
      if (!is.null(sel)) {
        fit <- sel$selected$fit
        zhat <- viterbi_decode(sim$y, sim$basis_spec, fit)
        ours_row$converged <- isTRUE(fit$converged)
        ours_row$viterbi_accuracy <- best_label_accuracy(sim$z, zhat, K)
        ours_row$information_criterion <- as.numeric(sel$winner_row$aic_edf[1L])
      }
    }

    cmp <- tryCatch(
      fit_hmmtmb_beta_gam(sim$y, sim$x, K = K, basis_k = basis_k),
      error = function(e) NULL
    )
    hmm_row <- data.frame(
      rep = r, method = "hmmTMB",
      converged = !is.null(cmp), viterbi_accuracy = NA_real_,
      information_criterion = NA_real_, elapsed_seconds = NA_real_,
      stringsAsFactors = FALSE
    )
    if (!is.null(cmp)) {
      hmm_row$viterbi_accuracy <- best_label_accuracy(sim$z, cmp$viterbi, K)
      hmm_row$information_criterion <- cmp$aic_conditional
      hmm_row$elapsed_seconds <- cmp$elapsed_seconds
    }

    rbind(ours_row, hmm_row)
  }

  raw <- parallel_lapply_crossplatform(
    seq_len(n_rep), run_one,
    n_cores = min(as.integer(n_cores %||% 1L), n_rep),
    exports = c(
      "run_one", "T", "K", "basis_k", "lambda_grid",
      "n_start_ours", "seed"
    ),
    project_root = getwd(), seed = seed, progress = TRUE,
    progress_label = sprintf("External-method comparison: %d replicates", n_rep)
  )
  do.call(rbind, raw)
}

summarize_hmmtmb_comparison <- function(x) {
  required <- c("rep", "method", "converged", "viterbi_accuracy", "elapsed_seconds")
  if (!all(required %in% names(x))) {
    stop("Comparison table does not contain the required columns.")
  }

  method_summary <- do.call(rbind, lapply(split(x, x$method), function(dd) {
    data.frame(
      method = unique(dd$method),
      n = nrow(dd),
      convergence_rate = mean(dd$converged, na.rm = TRUE),
      mean_viterbi_accuracy = mean(dd$viterbi_accuracy, na.rm = TRUE),
      median_viterbi_accuracy = stats::median(dd$viterbi_accuracy, na.rm = TRUE),
      sd_viterbi_accuracy = stats::sd(dd$viterbi_accuracy, na.rm = TRUE),
      mean_elapsed_seconds = mean(dd$elapsed_seconds, na.rm = TRUE),
      median_elapsed_seconds = stats::median(dd$elapsed_seconds, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))

  ours <- x[x$method == "Beta-GAM-HMM" & x$converged, , drop = FALSE]
  comparator <- x[x$method == "hmmTMB" & x$converged, , drop = FALSE]
  paired <- merge(
    ours[, c("rep", "viterbi_accuracy", "elapsed_seconds")],
    comparator[, c("rep", "viterbi_accuracy", "elapsed_seconds")],
    by = "rep", suffixes = c("_ours", "_hmmtmb")
  )

  if (nrow(paired) >= 2L) {
    accuracy_difference <- paired$viterbi_accuracy_ours - paired$viterbi_accuracy_hmmtmb
    tt <- stats::t.test(accuracy_difference, mu = 0)
    speed_ratio <- paired$elapsed_seconds_hmmtmb / paired$elapsed_seconds_ours

    paired_summary <- data.frame(
      n_paired = nrow(paired),
      mean_accuracy_ours = mean(paired$viterbi_accuracy_ours),
      mean_accuracy_hmmtmb = mean(paired$viterbi_accuracy_hmmtmb),
      mean_accuracy_difference = mean(accuracy_difference),
      accuracy_difference_ci_lower = unname(tt$conf.int[1L]),
      accuracy_difference_ci_upper = unname(tt$conf.int[2L]),
      paired_accuracy_p_value = tt$p.value,
      mean_speed_ratio_hmmtmb_over_ours = mean(speed_ratio),
      median_speed_ratio_hmmtmb_over_ours = stats::median(speed_ratio),
      stringsAsFactors = FALSE
    )
  } else {
    paired_summary <- data.frame(
      n_paired = nrow(paired),
      mean_accuracy_ours = NA_real_,
      mean_accuracy_hmmtmb = NA_real_,
      mean_accuracy_difference = NA_real_,
      accuracy_difference_ci_lower = NA_real_,
      accuracy_difference_ci_upper = NA_real_,
      paired_accuracy_p_value = NA_real_,
      mean_speed_ratio_hmmtmb_over_ours = NA_real_,
      median_speed_ratio_hmmtmb_over_ours = NA_real_,
      stringsAsFactors = FALSE
    )
  }

  list(method_summary = method_summary, paired_summary = paired_summary)
}
