
check_required_packages <- function(pkgs = c("ggplot2", "pbapply")) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop(
      "Missing packages: ", paste(missing, collapse = ", "),
      ". Run scripts/install_packages.R first."
    )
  }
  invisible(TRUE)
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    suppressWarnings(dir.create(path, recursive = TRUE, showWarnings = FALSE))
  }
  if (!dir.exists(path)) {
    base <- tryCatch(normalizePath(getwd(), winslash = "/", mustWork = TRUE), error = function(e) getwd())
    target <- file.path(base, path)
    hint <- if (.Platform$OS.type == "windows")
      " On Windows this is commonly caused by an overly long path. Move/extract the project to a short location such as C:/BGHMM and rerun; checkpoints are relative and will be preserved when the whole project folder is moved."
    else " Check write permissions and the parent directory."
    stop("Could not create output directory: ", target, ".", hint, call. = FALSE)
  }
  invisible(path)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

# Derive a reproducible random seed from a model candidate rather than from its
# position in a grid.  This makes a given (K, lambda) fit use the same random
# starts in a full-grid search, an adaptive refinement, or a standalone fit.
# Stable candidate-level seeds are important when benchmarking search strategies:
# otherwise differences in selected models can be caused by different local
# optima rather than by the search strategy itself.
model_candidate_seed <- function(base_seed, K, lambda) {
  if (is.null(base_seed)) return(NULL)
  lambda_txt <- paste(formatC(as.numeric(lambda), digits = 12, format = "fg", flag = "#"), collapse = ",")
  key <- paste0("K=", as.integer(K), "|lambda=", lambda_txt)
  code <- utf8ToInt(key)
  # A lightweight deterministic checksum that does not require an extra package.
  checksum <- sum((seq_along(code) * code) %% 1000003L) %% 1000003L
  modulus <- .Machine$integer.max - 1
  out <- (as.double(base_seed) + 1000003 * as.integer(K) + checksum) %% modulus
  as.integer(out + 1)
}

clip_unit <- function(x, eps = 1e-8) {
  pmin(pmax(x, eps), 1 - eps)
}

# Stable inverse-logit implementation.
inv_logit <- function(eta) {
  stats::plogis(eta)
}

# Validate proportion responses before fitting.  The default policy is deliberately
# strict: a standard Beta likelihood is defined for observations in the open unit
# interval.  Boundary values should therefore be modelled explicitly rather than
# silently perturbed.  The optional `clip` mode is retained only for reproducibility
# of legacy analyses and emits a warning.
prepare_unit_interval_response <- function(y,
                                           boundary_handling = c("error", "clip"),
                                           eps = 1e-8) {
  boundary_handling <- match.arg(boundary_handling)
  y <- as.numeric(y)
  if (any(!is.finite(y))) stop("Response contains non-finite values.")
  at_boundary <- y <= 0 | y >= 1
  if (any(at_boundary)) {
    if (boundary_handling == "error") {
      stop(
        "The Beta likelihood requires y in (0, 1). Found ", sum(at_boundary),
        " boundary observation(s). Use a zero/one-inflated model, a count model ",
        "when denominators are available, or set boundary_handling = 'clip' only ",
        "for a documented sensitivity analysis."
      )
    }
    warning(
      "Clipping boundary observations for numerical compatibility. This changes ",
      "the likelihood and should be reported as a sensitivity analysis.",
      call. = FALSE
    )
    y <- clip_unit(y, eps = eps)
  }
  y
}

log_sum_exp <- function(v) {
  m <- max(v)
  m + log(sum(exp(v - m)))
}

normalize_vec <- function(x, eps = 1e-12) {
  x <- pmax(x, eps)
  x / sum(x)
}

normalize_rows <- function(A, eps = 1e-12) {
  A <- pmax(A, eps)
  A / rowSums(A)
}

log_beta_density <- function(y, mu, phi, eps = 1e-10) {
  if (any(y <= 0 | y >= 1 | !is.finite(y))) {
    stop("log_beta_density requires finite y strictly inside (0, 1).")
  }
  mu <- clip_unit(mu, eps)
  phi <- pmax(phi, eps)
  a <- mu * phi
  b <- (1 - mu) * phi
  lgamma(phi) - lgamma(a) - lgamma(b) +
    (a - 1) * log(y) + (b - 1) * log1p(-y)
}

pspline_penalty <- function(M, diff_order = 2L) {
  if (M <= diff_order) {
    return(diag(M))
  }
  D <- diff(diag(M), differences = diff_order)
  crossprod(D)
}

compute_boundary_knots <- function(x, boundary_padding_frac = 0, boundary_knots = NULL) {
  x <- as.numeric(x)
  if (!is.null(boundary_knots)) return(as.numeric(boundary_knots))
  rng <- range(x)
  pad <- boundary_padding_frac * diff(rng)
  c(rng[1] - pad, rng[2] + pad)
}

compute_inner_knots <- function(x,
                                n_inner_knots = 6L,
                                knots = NULL,
                                boundary_knots = range(x),
                                knot_quantile_range = c(0, 1)) {
  if (!is.null(knots)) return(unique(as.numeric(knots)))
  probs <- seq(knot_quantile_range[1], knot_quantile_range[2], length.out = n_inner_knots + 2L)
  probs <- probs[probs > 0 & probs < 1]
  knots <- if (length(probs) > 0L) stats::quantile(x, probs = probs, names = FALSE) else NULL
  knots <- unique(as.numeric(knots))
  knots <- knots[knots > boundary_knots[1] & knots < boundary_knots[2]]
  if (length(knots) == 0L) knots <- NULL
  knots
}

build_spline_basis <- function(x,
                               type = "bs",
                               degree = 3L,
                               n_inner_knots = 6L,
                               knots = NULL,
                               boundary_knots = NULL,
                               boundary_padding_frac = 0,
                               knot_quantile_range = c(0, 1)) {
  type <- "bs"
  x <- as.numeric(x)
  boundary_knots <- compute_boundary_knots(x, boundary_padding_frac = boundary_padding_frac, boundary_knots = boundary_knots)
  knots <- compute_inner_knots(x,
                               n_inner_knots = n_inner_knots,
                               knots = knots,
                               boundary_knots = boundary_knots,
                               knot_quantile_range = knot_quantile_range)

  B_obj <- splines::bs(
    x = x,
    knots = knots,
    degree = degree,
    intercept = TRUE,
    Boundary.knots = boundary_knots
  )

  B <- unclass(B_obj)
  list(
    x = x,
    B = B,
    type = "bs",
    degree = degree,
    knots = attr(B_obj, "knots"),
    boundary_knots = attr(B_obj, "Boundary.knots"),
    intercept = TRUE,
    boundary_padding_frac = boundary_padding_frac,
    knot_quantile_range = knot_quantile_range,
    penalty = pspline_penalty(ncol(B), diff_order = 2L)
  )
}

build_bspline_basis <- function(x,
                                degree = 3L,
                                n_inner_knots = 6L,
                                knots = NULL,
                                boundary_knots = NULL,
                                boundary_padding_frac = 0,
                                knot_quantile_range = c(0, 1),
                                type = "bs") {
  build_spline_basis(
    x = x,
    type = "bs",
    degree = degree,
    n_inner_knots = n_inner_knots,
    knots = knots,
    boundary_knots = boundary_knots,
    boundary_padding_frac = boundary_padding_frac,
    knot_quantile_range = knot_quantile_range
  )
}

block_diag_matrices <- function(mats) {
  mats <- Filter(Negate(is.null), mats)
  if (length(mats) == 0L) return(matrix(0, 0, 0))
  nr <- sum(vapply(mats, nrow, integer(1)))
  nc <- sum(vapply(mats, ncol, integer(1)))
  out <- matrix(0, nr, nc)
  r0 <- 0L
  c0 <- 0L
  for (m in mats) {
    rr <- seq.int(r0 + 1L, r0 + nrow(m))
    cc <- seq.int(c0 + 1L, c0 + ncol(m))
    out[rr, cc] <- m
    r0 <- r0 + nrow(m)
    c0 <- c0 + ncol(m)
  }
  out
}

# Construct an additive B-spline design for p >= 1 covariates.  Each covariate
# receives its own basis block and second-difference penalty.  The current model
# can use one common smoothing parameter or a state-specific smoothing parameter;
# the design itself is fully multivariable.
build_additive_basis <- function(data,
                                 covariates = names(data),
                                 degree = 3L,
                                 n_inner_knots = 6L,
                                 boundary_padding_frac = 0,
                                 knot_quantile_range = c(0, 1)) {
  data <- as.data.frame(data)
  if (length(covariates) < 1L) stop("At least one covariate is required.")
  if (!all(covariates %in% names(data))) stop("Unknown covariate name in `covariates`.")

  term_specs <- lapply(covariates, function(v) {
    x <- as.numeric(data[[v]])
    boundary_knots <- compute_boundary_knots(x, boundary_padding_frac = boundary_padding_frac)
    knots <- compute_inner_knots(
      x, n_inner_knots = n_inner_knots, boundary_knots = boundary_knots,
      knot_quantile_range = knot_quantile_range
    )
    B_obj <- splines::bs(
      x = x, knots = knots, degree = degree, intercept = FALSE,
      Boundary.knots = boundary_knots
    )
    B_raw <- unclass(B_obj)
    center <- colMeans(B_raw)
    B_centered <- sweep(B_raw, 2L, center, FUN = "-")
    list(
      B = B_centered,
      center = center,
      degree = degree,
      knots = attr(B_obj, "knots"),
      boundary_knots = attr(B_obj, "Boundary.knots"),
      intercept = FALSE,
      penalty = pspline_penalty(ncol(B_centered), diff_order = 2L)
    )
  })
  names(term_specs) <- covariates
  B_blocks <- lapply(term_specs, `[[`, "B")
  P_blocks <- lapply(term_specs, `[[`, "penalty")
  widths <- vapply(B_blocks, ncol, integer(1))
  ends <- 1L + cumsum(widths)
  starts <- c(2L, head(ends, -1L) + 1L)
  term_columns <- Map(seq.int, starts, ends)
  names(term_columns) <- covariates

  list(
    x = data,
    data = data,
    covariates = covariates,
    B = cbind(`(Intercept)` = 1, do.call(cbind, B_blocks)),
    type = "additive_bs",
    term_specs = term_specs,
    term_columns = term_columns,
    penalty = block_diag_matrices(c(list(matrix(0, 1L, 1L)), P_blocks)),
    degree = degree,
    n_inner_knots = n_inner_knots,
    boundary_padding_frac = boundary_padding_frac,
    knot_quantile_range = knot_quantile_range
  )
}

predict_basis <- function(x_new, basis_spec) {
  if (identical(basis_spec$type, "additive_bs")) {
    new_data <- as.data.frame(x_new)
    if (!all(basis_spec$covariates %in% names(new_data))) {
      stop("New data must contain: ", paste(basis_spec$covariates, collapse = ", "))
    }
    blocks <- lapply(basis_spec$covariates, function(v) {
      spec <- basis_spec$term_specs[[v]]
      B_new <- unclass(splines::bs(
        x = new_data[[v]],
        knots = spec$knots,
        degree = spec$degree,
        intercept = FALSE,
        Boundary.knots = spec$boundary_knots
      ))
      sweep(B_new, 2L, spec$center, FUN = "-")
    })
    return(cbind(`(Intercept)` = 1, do.call(cbind, blocks)))
  }

  unclass(
    splines::bs(
      x = x_new,
      knots = basis_spec$knots,
      degree = basis_spec$degree,
      intercept = TRUE,
      Boundary.knots = basis_spec$boundary_knots
    )
  )
}

predict_bspline_basis <- function(x_new, basis_spec) {
  predict_basis(x_new, basis_spec)
}

stationary_distribution <- function(A) {
  ev <- eigen(t(A))
  idx <- which.min(abs(ev$values - 1))
  vec <- Re(ev$vectors[, idx])
  vec <- vec / sum(vec)
  vec <- pmax(vec, 0)
  normalize_vec(vec)
}

flatten_transition <- function(A) {
  as.vector(t(A))
}

n_parameters_hmm <- function(K, M) {
  K * M + K + K * (K - 1L) + (K - 1L)
}

n_parameters_hmm_edf <- function(K, edf_spline) {
  sum(edf_spline) + K + K * (K - 1L) + (K - 1L)
}

entropy_from_gamma <- function(gamma, eps = 1e-12) {
  -sum(gamma * log(pmax(gamma, eps)))
}

state_occupancy_stats <- function(gamma) {
  T <- nrow(gamma)
  mass <- colSums(gamma)
  share <- mass / T
  list(
    mass = as.numeric(mass),
    share = as.numeric(share),
    min_mass = min(mass),
    min_share = min(share),
    max_share = max(share)
  )
}

diagnose_state_occupancy <- function(gamma,
                                     min_state_share = 0.05,
                                     min_state_mass = 50) {
  occ <- state_occupancy_stats(gamma)
  flagged <- any(occ$share < min_state_share) || any(occ$mass < min_state_mass)
  list(
    flagged = flagged,
    min_state_share = occ$min_share,
    min_state_mass = occ$min_mass,
    state_share = occ$share,
    state_mass = occ$mass,
    threshold_share = min_state_share,
    threshold_mass = min_state_mass
  )
}

effective_df_state <- function(B,
                               gamma_k,
                               mu_k,
                               phi_k,
                               lambda,
                               penalty_matrix,
                               ridge = 1e-8) {
  # Fisher information for the additive predictor in a Beta mean--precision
  # parameterization.  This replaces the earlier working-weight approximation
  # that omitted the precision and trigamma terms.
  mu_k <- clip_unit(as.numeric(mu_k), 1e-10)
  phi_k <- as.numeric(phi_k)[1]
  dmu_deta <- mu_k * (1 - mu_k)
  fisher_eta <- phi_k^2 *
    (trigamma(mu_k * phi_k) + trigamma((1 - mu_k) * phi_k)) *
    dmu_deta^2
  w <- pmax(as.numeric(gamma_k) * fisher_eta, 1e-12)

  BW <- B * sqrt(w)
  XtWX <- crossprod(BW)
  # The penalized objective contains lambda * beta' P beta; hence its Hessian
  # contribution is 2 * lambda * P.
  S <- 2 * as.numeric(lambda)[1] * penalty_matrix + ridge * diag(ncol(B))
  A <- XtWX + S
  H <- tryCatch(
    solve(A, XtWX),
    error = function(e) qr.solve(A, XtWX, tol = 1e-10)
  )
  edf <- sum(diag(H))
  as.numeric(min(max(edf, 0), ncol(B)))
}

compute_fit_edf <- function(fit) {
  required <- c("B", "gamma", "mu", "phi", "penalty_matrix")
  if (any(vapply(required, function(nm) is.null(fit[[nm]]), logical(1)))) {
    return(rep(NA_real_, fit$K))
  }
  vapply(seq_len(fit$K), function(k) {
    effective_df_state(
      B = fit$B,
      gamma_k = fit$gamma[, k],
      mu_k = fit$mu[, k],
      phi_k = fit$phi[k],
      lambda = if (length(fit$lambda) == 1L) fit$lambda else fit$lambda[k],
      penalty_matrix = fit$penalty_matrix
    )
  }, numeric(1))
}

calculate_information_criteria <- function(loglik_obs, gamma, K, M, T, edf_spline = NULL) {
  npar_raw <- n_parameters_hmm(K, M)
  entropy <- entropy_from_gamma(gamma)
  aic_raw <- -2 * loglik_obs + 2 * npar_raw
  bic_raw <- -2 * loglik_obs + npar_raw * log(T)
  icl_raw <- bic_raw + 2 * entropy

  if (is.null(edf_spline) || any(!is.finite(edf_spline))) {
    npar_edf <- NA_real_
    aic_edf <- NA_real_
    bic_edf <- NA_real_
    icl_edf <- NA_real_
  } else {
    npar_edf <- n_parameters_hmm_edf(K, edf_spline)
    aic_edf <- -2 * loglik_obs + 2 * npar_edf
    bic_edf <- -2 * loglik_obs + npar_edf * log(T)
    icl_edf <- bic_edf + 2 * entropy
  }

  list(
    aic = aic_raw,
    bic = bic_raw,
    icl = icl_raw,
    aic_edf = aic_edf,
    bic_edf = bic_edf,
    icl_edf = icl_edf,
    entropy = entropy,
    npar = npar_raw,
    npar_edf = npar_edf,
    edf_spline = edf_spline
  )
}

diagnose_phi_degeneracy <- function(phi,
                                    phi_max = 500,
                                    eps = 1e-6,
                                    s_thresh = 2L,
                                    delta_abs = 50,
                                    delta_sum = 100) {
  phi <- sort(as.numeric(phi))
  K <- length(phi)
  n_sat <- sum(phi >= (phi_max - eps))
  deltas <- if (K >= 2L) diff(phi) else numeric(0)
  tail_idx <- if (length(deltas) >= 3L) seq.int(length(deltas) - 2L, length(deltas)) else seq_along(deltas)
  delta_tail <- if (length(tail_idx) > 0L) sum(deltas[tail_idx]) else 0
  large_jump <- if (length(tail_idx) > 0L) any(deltas[tail_idx] > delta_abs) else FALSE
  flagged <- (n_sat >= s_thresh) || large_jump || (delta_tail > delta_sum)
  list(
    flagged = flagged,
    phi_sorted = phi,
    n_sat = n_sat,
    deltas = deltas,
    delta_tail = delta_tail,
    large_jump = large_jump
  )
}

all_permutations <- function(v) {
  v <- as.integer(v)
  if (length(v) == 1L) {
    return(matrix(v, nrow = 1L))
  }
  out <- vector("list", length(v))
  for (i in seq_along(v)) {
    rest <- v[-i]
    p_rest <- all_permutations(rest)
    out[[i]] <- cbind(v[i], p_rest)
  }
  do.call(rbind, out)
}

relabel_fit <- function(fit, perm) {
  perm <- as.integer(perm)
  fit$beta <- fit$beta[perm]
  fit$phi <- fit$phi[perm]
  fit$pi <- fit$pi[perm]
  fit$A <- fit$A[perm, perm, drop = FALSE]
  if (!is.null(fit$gamma)) fit$gamma <- fit$gamma[, perm, drop = FALSE]
  if (!is.null(fit$xi)) fit$xi <- fit$xi[, perm, perm, drop = FALSE]
  if (!is.null(fit$z_viterbi)) fit$z_viterbi <- match(fit$z_viterbi, perm)
  fit
}

state_mean_curves <- function(fit, basis_spec, x_grid) {
  B_grid <- predict_bspline_basis(x_grid, basis_spec)
  K <- length(fit$beta)
  mu <- sapply(seq_len(K), function(k) inv_logit(drop(B_grid %*% fit$beta[[k]])))
  list(x = x_grid, B_grid = B_grid, mu = mu)
}

align_fit_to_reference <- function(fit, reference_fit, basis_spec, x_grid = NULL) {
  K_fit <- length(fit$phi)
  K_ref <- length(reference_fit$phi)
  if (K_fit != K_ref) return(fit)

  # When no prediction grid is supplied, align states on the observed design
  # matrix.  This also works for additive models with more than one covariate.
  if (is.null(x_grid)) {
    B_eval <- basis_spec$B
    curves_fit <- sapply(seq_len(K_fit), function(k) inv_logit(drop(B_eval %*% fit$beta[[k]])))
    curves_ref <- sapply(seq_len(K_ref), function(k) inv_logit(drop(B_eval %*% reference_fit$beta[[k]])))
  } else {
    curves_fit <- state_mean_curves(fit, basis_spec, x_grid)$mu
    curves_ref <- state_mean_curves(reference_fit, basis_spec, x_grid)$mu
  }

  perms <- all_permutations(seq_len(K_fit))
  costs <- apply(perms, 1L, function(perm) {
    perm <- as.integer(perm)
    curve_cost <- sum((curves_fit[, perm, drop = FALSE] - curves_ref)^2)
    phi_cost <- sum((fit$phi[perm] - reference_fit$phi)^2) / max(1, K_fit)
    pi_cost <- sum((fit$pi[perm] - reference_fit$pi)^2)
    curve_cost + 0.01 * phi_cost + 0.01 * pi_cost
  })
  best_perm <- perms[which.min(costs), ]
  relabel_fit(fit, best_perm)
}

safe_quantile <- function(x, probs = c(0.025, 0.975)) {
  stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE)
}

project_r_files <- function() {
  # Discover numbered project modules automatically.  Keeping this dynamic avoids
  # silently omitting newly added modules (for example, targeted follow-up code).
  files <- list.files(
    "R",
    pattern = "^[0-9]{2}_[A-Za-z0-9_]+\\.R$",
    full.names = TRUE
  )
  sort(files)
}

setup_parallel_workers <- function(cl, project_root = getwd()) {
  project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  parallel::clusterExport(cl, varlist = c("project_root"), envir = environment())
  parallel::clusterEvalQ(cl, {
    setwd(project_root)

    # Load utilities first, then discover and load every remaining numbered
    # module.  This mirrors the main R session and is important on Windows,
    # where PSOCK workers start from clean R sessions.
    source("R/00_utils.R", local = .GlobalEnv)
    module_files <- setdiff(project_r_files(), "R/00_utils.R")
    for (f in module_files) source(f, local = .GlobalEnv)
    NULL
  })
  invisible(cl)
}

parallel_lapply_crossplatform <- function(X,
                                          FUN,
                                          n_cores = 1L,
                                          exports = NULL,
                                          project_root = getwd(),
                                          seed = NULL,
                                          use_fork_on_unix = TRUE,
                                          progress = interactive(),
                                          progress_label = NULL) {
  n_cores <- as.integer(n_cores %||% 1L)
  use_pb <- isTRUE(progress) && requireNamespace("pbapply", quietly = TRUE)

  if (!is.null(progress_label) && nzchar(progress_label)) {
    message(progress_label)
  }

  if (n_cores <= 1L || length(X) <= 1L) {
    if (use_pb) return(pbapply::pblapply(X, FUN))
    return(lapply(X, FUN))
  }

  if (.Platform$OS.type != "windows" && isTRUE(use_fork_on_unix)) {
    if (use_pb) return(pbapply::pblapply(X, FUN, cl = n_cores))
    return(parallel::mclapply(X, FUN, mc.cores = n_cores))
  }

  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  setup_parallel_workers(cl, project_root = project_root)
  if (!is.null(seed)) parallel::clusterSetRNGStream(cl, iseed = as.integer(seed)[1])
  if (!is.null(exports) && length(exports) > 0L) {
    parallel::clusterExport(cl, varlist = unique(exports), envir = parent.frame())
  }
  if (use_pb) return(pbapply::pblapply(X, FUN, cl = cl))
  parallel::parLapply(cl, X, FUN)
}

format_phi <- function(phi, digits = 3L) {
  paste(format(round(phi, digits), nsmall = digits, trim = TRUE), collapse = " | ")
}

save_rds_safe <- function(object, file) {
  ensure_dir(dirname(file))
  saveRDS(object, file = file)
  invisible(file)
}

write_csv_safe <- function(df, file, row.names = FALSE) {
  ensure_dir(dirname(file))
  utils::write.csv(df, file = file, row.names = row.names)
  invisible(file)
}

chunk_indices <- function(n, chunk_size) {
  idx <- seq_len(n)
  split(idx, ceiling(idx / chunk_size))
}

merge_data_frames <- function(lst) {
  lst <- Filter(Negate(is.null), lst)
  if (length(lst) == 0L) return(data.frame())

  all_names <- unique(unlist(lapply(lst, names)))
  lst_aligned <- lapply(lst, function(df) {
    miss <- setdiff(all_names, names(df))
    if (length(miss) > 0L) {
      for (nm in miss) df[[nm]] <- NA
    }
    df[, all_names, drop = FALSE]
  })

  do.call(rbind, lst_aligned)
}

rule_grid_from_cfg <- function(cfg) {
  cfg$diagnostic_sensitivity %||% list(
    baseline = cfg$diagnostic,
    permissive = modifyList(cfg$diagnostic, list(s_thresh = cfg$diagnostic$s_thresh + 1L, delta_abs = cfg$diagnostic$delta_abs + 25, delta_sum = cfg$diagnostic$delta_sum + 50)),
    strict = modifyList(cfg$diagnostic, list(s_thresh = max(1L, cfg$diagnostic$s_thresh - 1L), delta_abs = max(10, cfg$diagnostic$delta_abs - 15), delta_sum = max(25, cfg$diagnostic$delta_sum - 30)))
  )
}
