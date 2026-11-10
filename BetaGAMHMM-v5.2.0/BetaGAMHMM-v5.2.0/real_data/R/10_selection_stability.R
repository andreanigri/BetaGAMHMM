# Separate model-selection bootstrap with criterion re-selection.
# K and lambda are reselected in each bootstrap sample. Parameter confidence
# intervals remain a separate conditional bootstrap with K and lambda fixed.

selection_stability_worker <- function(job) {
  sim <- simulate_beta_panel_from_fit(job$fit, seed = job$seed + job$b)
  grid <- tryCatch(
    run_beta_grid(
      panel = sim$panel, basis_spec = job$fit$basis_spec,
      K_grid = job$K_grid, lambda_grid = job$lambda_grid,
      n_start = job$n_start, n_cores = 1L,
      max_iter = job$max_iter, tol = job$tol,
      phi_lower = job$phi_lower, phi_upper = job$phi_upper,
      optimizer_maxit = job$optimizer_maxit,
      precision_saturation_fraction = job$precision_saturation_fraction,
      seed = job$seed + 100000L + 1000L * job$b, progress = FALSE
    ), error = function(e) NULL
  )
  if (is.null(grid)) return(data.frame(
    replicate_id = job$b, success = FALSE, selected_K = NA_integer_, selected_lambda = NA_real_,
    selected_value = NA_real_, selected_min_state_share = NA_real_, boundary_selection = NA,
    stringsAsFactors = FALSE
  ))
  sel <- tryCatch(select_beta_model(grid), error = function(e) NULL)
  if (is.null(sel)) return(data.frame(
    replicate_id = job$b, success = FALSE, selected_K = NA_integer_, selected_lambda = NA_real_,
    selected_value = NA_real_, selected_min_state_share = NA_real_, boundary_selection = NA,
    stringsAsFactors = FALSE
  ))
  data.frame(
    replicate_id = job$b, success = TRUE,
    selected_K = sel$fit$K, selected_lambda = sel$fit$lambda[1L],
    selected_value = sel$selected_row[[sel$criterion_column]][1L],
    selected_min_state_share = sel$fit$occupancy$min_share,
    boundary_selection = sel$fit$K == max(job$K_grid), stringsAsFactors = FALSE
  )
}

run_selection_stability_bootstrap <- function(fit,
                                              B = 100L,
                                              K_grid,
                                              lambda_grid,
                                              criterion,
                                              n_start = 5L,
                                              n_cores = 1L,
                                              max_iter = 300L,
                                              tol = 1e-6,
                                              phi_lower = 1,
                                              phi_upper = 2000,
                                              optimizer_maxit = 150L,
                                              precision_saturation_fraction = 0.98,
                                              seed = 20265001L) {
  jobs <- lapply(seq_len(B), function(b) list(
    fit = fit, b = b, seed = seed, criterion = criterion,
    K_grid = as.integer(K_grid), lambda_grid = lambda_grid,
    n_start = n_start, max_iter = max_iter, tol = tol,
    phi_lower = phi_lower, phi_upper = phi_upper,
    optimizer_maxit = optimizer_maxit,
    precision_saturation_fraction = precision_saturation_fraction
  ))
  res <- if (n_cores > 1L) {
    parallel_lapply(jobs, selection_stability_worker, n_cores = n_cores,
                    project_root = getwd(), seed = seed + 900000L)
  } else lapply(jobs, selection_stability_worker)

  details <- do.call(rbind, res)
  ok <- details[details$success, , drop = FALSE]
  if (!nrow(ok)) stop("All model-selection bootstrap replicates failed.")

  k_freq <- as.data.frame(table(factor(ok$selected_K, levels = K_grid)), stringsAsFactors = FALSE)
  names(k_freq) <- c("K", "frequency")
  k_freq$K <- as.integer(as.character(k_freq$K)); k_freq$proportion <- k_freq$frequency / nrow(ok)

  lam_levels <- sort(unique(lambda_grid))
  lam_freq <- as.data.frame(table(factor(ok$selected_lambda, levels = lam_levels)), stringsAsFactors = FALSE)
  names(lam_freq) <- c("lambda", "frequency")
  lam_freq$lambda <- as.numeric(as.character(lam_freq$lambda)); lam_freq$proportion <- lam_freq$frequency / nrow(ok)

  joint <- as.data.frame(table(ok$selected_K, ok$selected_lambda), stringsAsFactors = FALSE)
  names(joint) <- c("K", "lambda", "frequency")
  joint <- joint[joint$frequency > 0L, , drop = FALSE]
  joint$K <- as.integer(as.character(joint$K)); joint$lambda <- as.numeric(as.character(joint$lambda))
  joint$proportion <- joint$frequency / nrow(ok)

  summary <- data.frame(
    criterion = criterion, B_requested = B, B_success = nrow(ok), success_rate = nrow(ok) / B,
    reference_K = fit$K, reference_lambda = fit$lambda[1L],
    same_K_rate = mean(ok$selected_K == fit$K),
    same_K_lambda_rate = mean(ok$selected_K == fit$K & abs(ok$selected_lambda - fit$lambda[1L]) < 1e-12),
    boundary_selection_rate = mean(ok$boundary_selection), stringsAsFactors = FALSE
  )
  list(summary = summary, details = details, K_frequency = k_freq,
       lambda_frequency = lam_freq, joint_frequency = joint,
       K_grid = K_grid, lambda_grid = lambda_grid, criterion = criterion)
}
