# Computational benchmarking and reproducibility helpers.

system_information_table <- function() {
  info <- Sys.info()
  data.frame(
    item = c(
      "R_version", "platform", "system", "release", "machine",
      "logical_cores", "physical_cores"
    ),
    value = c(
      R.version.string,
      R.version$platform,
      unname(info[["sysname"]] %||% NA_character_),
      unname(info[["release"]] %||% NA_character_),
      unname(info[["machine"]] %||% NA_character_),
      as.character(parallel::detectCores(logical = TRUE)),
      as.character(parallel::detectCores(logical = FALSE))
    ),
    stringsAsFactors = FALSE
  )
}

recommended_core_count <- function(reserve = 1L, max_cores = 12L) {
  n <- parallel::detectCores(logical = FALSE)
  if (is.na(n)) n <- parallel::detectCores(logical = TRUE)
  if (is.na(n)) n <- 2L
  max(1L, min(as.integer(max_cores), as.integer(n) - as.integer(reserve)))
}

# Empirical scaling benchmark for one fixed (K, lambda) fit.  It is intentionally
# small and can be run on the same machine used for the paper to report wall-clock
# times and software/hardware information.
run_computational_scaling_benchmark <- function(T_grid = c(300L, 700L, 1200L),
                                                K_grid = c(2L, 4L, 5L),
                                                n_start = 3L,
                                                max_iter = 120L,
                                                lambda = 0.1,
                                                seed = 8500L) {
  design <- expand.grid(T = T_grid, K = K_grid, KEEP.OUT.ATTRS = FALSE)
  rows <- lapply(seq_len(nrow(design)), function(i) {
    T_i <- design$T[i]
    K_i <- design$K[i]
    phi_i <- seq(10, 10 + 8 * (K_i - 1L), length.out = K_i)
    sim <- simulate_default_scenario(
      T = T_i, K = K_i, diagonal = 0.90,
      phi_true = phi_i, seed = seed + i * 100L
    )
    started <- proc.time()[[3L]]
    fit <- tryCatch(
      fit_beta_gam_hmm(
        sim$y, sim$basis_spec,
        K = K_i, lambda = lambda,
        n_start = n_start,
        max_iter = max_iter,
        optimizer_maxit = 80L,
        seed = seed + i * 100L + 1L
      ),
      error = function(e) NULL
    )
    elapsed <- proc.time()[[3L]] - started
    data.frame(
      T = T_i,
      K = K_i,
      basis_dimension = ncol(sim$basis_spec$B),
      n_start = n_start,
      elapsed_seconds = elapsed,
      converged = if (is.null(fit)) FALSE else isTRUE(fit$converged),
      iterations = if (is.null(fit)) NA_integer_ else fit$n_iter,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# Compare serial and outer-parallel execution for a modest grid.  Nested
# parallelism is deliberately avoided because it usually oversubscribes CPUs and
# slows down HMM estimation.
benchmark_grid_parallelism <- function(n_cores = recommended_core_count(),
                                       T = 700L,
                                       K_grid = 1:4,
                                       lambda_grid = c(0.05, 0.1, 0.5),
                                       n_start = 2L,
                                       seed = 8600L) {
  sim <- simulate_default_scenario(T = T, K = 4L, seed = seed)
  run_mode <- function(cores) {
    started <- proc.time()[[3L]]
    out <- run_grid_search_hmm(
      sim$y, sim$basis_spec,
      K_grid = K_grid, lambda_grid = lambda_grid,
      n_start = n_start, max_iter = 100L, optimizer_maxit = 60L,
      n_cores = cores, n_start_cores = 1L,
      progress = FALSE, seed = seed + cores
    )
    data.frame(
      mode = if (cores == 1L) "serial" else "outer_parallel",
      n_cores = cores,
      elapsed_seconds = proc.time()[[3L]] - started,
      successful_models = sum(is.na(summarize_grid_results(out)$error)),
      stringsAsFactors = FALSE
    )
  }
  rbind(run_mode(1L), run_mode(max(1L, as.integer(n_cores))))
}

# Compact summaries for replicate-level simulation output.  Each requested
# metric is reported with the number of finite replicates, mean, standard
# deviation, median and interquartile range.
summarize_replicate_metrics <- function(df, metrics, group_cols = character(0)) {
  if (nrow(df) == 0L) return(data.frame())
  missing_metrics <- setdiff(metrics, names(df))
  if (length(missing_metrics)) stop("Missing metric columns: ", paste(missing_metrics, collapse = ", "))
  missing_groups <- setdiff(group_cols, names(df))
  if (length(missing_groups)) stop("Missing group columns: ", paste(missing_groups, collapse = ", "))

  if (length(group_cols)) {
    key <- interaction(df[group_cols], drop = TRUE, lex.order = TRUE)
    parts <- split(df, key)
  } else {
    parts <- list(all = df)
  }

  rows <- lapply(parts, function(dd) {
    group_values <- if (length(group_cols)) dd[1L, group_cols, drop = FALSE] else data.frame(dummy = 1)[, FALSE, drop = FALSE]
    metric_rows <- lapply(metrics, function(m) {
      x <- as.numeric(dd[[m]])
      x <- x[is.finite(x)]
      data.frame(
        metric = m,
        n = length(x),
        mean = if (length(x)) mean(x) else NA_real_,
        sd = if (length(x) > 1L) stats::sd(x) else NA_real_,
        median = if (length(x)) stats::median(x) else NA_real_,
        q25 = if (length(x)) stats::quantile(x, 0.25, names = FALSE) else NA_real_,
        q75 = if (length(x)) stats::quantile(x, 0.75, names = FALSE) else NA_real_,
        stringsAsFactors = FALSE
      )
    })
    out <- do.call(rbind, metric_rows)
    if (length(group_cols)) {
      group_values <- group_values[rep(1L, nrow(out)), , drop = FALSE]
      rownames(group_values) <- NULL
      out <- cbind(group_values, out)
    }
    out
  })
  do.call(rbind, rows)
}
