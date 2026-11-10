# Utility functions for the year-level Beta-GAM hidden Markov model.
#
# The real-data application treats calendar year as the HMM time index.
# Ages within a year form a vector-valued emission and share the same latent
# state. This file contains numerical helpers, spline construction, state
# relabelling, and information-criterion calculations used throughout the
# analysis.

`%||%` <- function(x, y) if (is.null(x)) y else x

inv_logit <- function(x) stats::plogis(x)

clip_unit <- function(x, eps = 1e-10) {
  pmin(pmax(x, eps), 1 - eps)
}

normalize_vec <- function(x) {
  x <- pmax(as.numeric(x), 1e-15)
  x / sum(x)
}

normalize_rows <- function(A) {
  A <- pmax(A, 1e-15)
  A / rowSums(A)
}

stationary_distribution <- function(A) {
  if (nrow(A) == 1L) return(1)
  ev <- eigen(t(A))
  idx <- which.min(abs(ev$values - 1))
  v <- Re(ev$vectors[, idx])
  v <- pmax(v, 0)
  if (sum(v) <= 0) v <- rep(1, length(v))
  normalize_vec(v)
}

pspline_penalty <- function(M, diff_order = 2L) {
  D <- diff(diag(M), differences = diff_order)
  crossprod(D)
}

build_age_basis <- function(age,
                            degree = 3L,
                            n_inner_knots = 6L,
                            knots = NULL,
                            boundary_knots = range(age)) {
  age <- as.numeric(age)
  if (is.null(knots)) {
    probs <- seq(0, 1, length.out = n_inner_knots + 2L)
    probs <- probs[probs > 0 & probs < 1]
    knots <- unique(as.numeric(stats::quantile(age, probs = probs, names = FALSE)))
    knots <- knots[knots > boundary_knots[1] & knots < boundary_knots[2]]
    if (!length(knots)) knots <- NULL
  }
  B_obj <- splines::bs(
    x = age,
    knots = knots,
    degree = degree,
    intercept = TRUE,
    Boundary.knots = boundary_knots
  )
  B <- unclass(B_obj)
  list(
    age = age,
    B = B,
    degree = degree,
    knots = attr(B_obj, "knots"),
    boundary_knots = attr(B_obj, "Boundary.knots"),
    penalty = pspline_penalty(ncol(B), diff_order = 2L)
  )
}

predict_age_basis <- function(age_new, basis_spec) {
  unclass(splines::bs(
    x = age_new,
    knots = basis_spec$knots,
    degree = basis_spec$degree,
    intercept = TRUE,
    Boundary.knots = basis_spec$boundary_knots
  ))
}

entropy_from_gamma <- function(gamma, eps = 1e-12) {
  -sum(gamma * log(pmax(gamma, eps)))
}

all_permutations <- function(v) {
  v <- as.integer(v)
  if (length(v) == 1L) return(matrix(v, nrow = 1L))
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
  K <- length(perm)
  fit$beta <- fit$beta[perm]
  if (!is.null(fit$phi)) fit$phi <- fit$phi[perm]
  if (!is.null(fit$lambda) && length(fit$lambda) == K) fit$lambda <- fit$lambda[perm]
  if (!is.null(fit$edf_spline) && length(fit$edf_spline) == K) fit$edf_spline <- fit$edf_spline[perm]
  fit$pi <- fit$pi[perm]
  fit$A <- fit$A[perm, perm, drop = FALSE]
  if (!is.null(fit$gamma)) fit$gamma <- fit$gamma[, perm, drop = FALSE]
  if (!is.null(fit$alpha)) fit$alpha <- fit$alpha[, perm, drop = FALSE]
  if (!is.null(fit$xi)) fit$xi <- fit$xi[, perm, perm, drop = FALSE]
  if (!is.null(fit$mu_age)) fit$mu_age <- fit$mu_age[, perm, drop = FALSE]
  if (!is.null(fit$log_emission)) fit$log_emission <- fit$log_emission[, perm, drop = FALSE]
  # perm gives old-state indices in the new order. Convert any already decoded
  # old labels to their corresponding new labels.
  if (!is.null(fit$viterbi)) fit$viterbi <- match(as.integer(fit$viterbi), perm)
  fit
}

order_states_by_profile <- function(fit) {
  profile_mean <- colMeans(fit$mu_age)
  order(profile_mean, seq_along(profile_mean))
}

forward_backward_hmm <- function(log_emission, pi, A) {
  Tn <- nrow(log_emission)
  K <- ncol(log_emission)
  row_shift <- apply(log_emission, 1L, max)
  E <- exp(log_emission - row_shift)

  alpha <- matrix(0, Tn, K)
  beta <- matrix(1, Tn, K)
  scale <- numeric(Tn)

  a1 <- as.numeric(pi) * E[1L, ]
  scale[1L] <- sum(a1)
  if (!is.finite(scale[1L]) || scale[1L] <= 0) stop("Forward recursion failed at year 1.")
  alpha[1L, ] <- a1 / scale[1L]

  if (Tn >= 2L) {
    for (t in 2:Tn) {
      at <- drop(alpha[t - 1L, ] %*% A) * E[t, ]
      scale[t] <- sum(at)
      if (!is.finite(scale[t]) || scale[t] <= 0) stop("Forward recursion failed at year index ", t, ".")
      alpha[t, ] <- at / scale[t]
    }
    for (t in (Tn - 1L):1L) {
      beta[t, ] <- drop(A %*% (E[t + 1L, ] * beta[t + 1L, ])) / scale[t + 1L]
    }
  }

  gamma <- alpha * beta
  gamma <- gamma / rowSums(gamma)

  xi <- array(0, dim = c(max(Tn - 1L, 1L), K, K))
  if (Tn >= 2L) {
    for (t in seq_len(Tn - 1L)) {
      right <- E[t + 1L, ] * beta[t + 1L, ]
      tmp <- A * tcrossprod(alpha[t, ], right)
      den <- sum(tmp)
      if (is.finite(den) && den > 0) tmp <- tmp / den
      xi[t, , ] <- tmp
    }
  }

  list(
    alpha = alpha,
    beta = beta,
    gamma = gamma,
    xi = xi,
    loglik = sum(row_shift + log(scale)),
    scale = scale
  )
}

viterbi_from_log_emission <- function(log_emission, pi, A) {
  Tn <- nrow(log_emission)
  K <- ncol(log_emission)
  delta <- matrix(-Inf, Tn, K)
  psi <- matrix(NA_integer_, Tn, K)
  delta[1L, ] <- log(pmax(pi, 1e-16)) + log_emission[1L, ]
  log_A <- log(pmax(A, 1e-16))

  if (Tn >= 2L) {
    for (t in 2:Tn) {
      for (k in seq_len(K)) {
        vals <- delta[t - 1L, ] + log_A[, k]
        psi[t, k] <- which.max(vals)
        delta[t, k] <- max(vals) + log_emission[t, k]
      }
    }
  }

  path <- integer(Tn)
  path[Tn] <- which.max(delta[Tn, ])
  if (Tn >= 2L) {
    for (t in (Tn - 1L):1L) path[t] <- psi[t + 1L, path[t + 1L]]
  }
  path
}

beta_edf_state <- function(B, gamma_k, mu_age, phi, lambda, penalty, ridge = 1e-8) {
  mu_age <- clip_unit(mu_age)
  dmu <- mu_age * (1 - mu_age)
  fisher_eta <- phi^2 *
    (trigamma(mu_age * phi) + trigamma((1 - mu_age) * phi)) *
    dmu^2
  effective_year_mass <- sum(gamma_k)
  w <- pmax(effective_year_mass * fisher_eta, 1e-12)
  BW <- B * sqrt(w)
  XtWX <- crossprod(BW)
  Hpen <- 2 * lambda * penalty + ridge * diag(ncol(B))
  H <- tryCatch(
    solve(XtWX + Hpen, XtWX),
    error = function(e) qr.solve(XtWX + Hpen, XtWX, tol = 1e-10)
  )
  as.numeric(min(max(sum(diag(H)), 0), ncol(B)))
}

binomial_edf_state <- function(B, gamma_k, mu_age, N, lambda, penalty, ridge = 1e-8) {
  # N is a T x J denominator matrix. Posterior state weights are at the year level.
  weighted_N_by_age <- colSums(N * gamma_k)
  w <- pmax(weighted_N_by_age * mu_age * (1 - mu_age), 1e-12)
  BW <- B * sqrt(w)
  XtWX <- crossprod(BW)
  Hpen <- 2 * lambda * penalty + ridge * diag(ncol(B))
  H <- tryCatch(
    solve(XtWX + Hpen, XtWX),
    error = function(e) qr.solve(XtWX + Hpen, XtWX, tol = 1e-10)
  )
  as.numeric(min(max(sum(diag(H)), 0), ncol(B)))
}

information_criteria <- function(loglik,
                                 gamma,
                                 K,
                                 edf_spline,
                                 n_years,
                                 n_cells = NULL,
                                 has_phi = TRUE) {
  # Final information criteria. The HMM sample size is the number of time
  # points T. This retains the original manuscript's BIC structure while using
  # the corrected Beta EDF for the spline component.
  trans_par <- K * (K - 1L)
  init_par <- K - 1L
  emission_extra <- if (has_phi) K else 0L
  npar_edf <- trans_par + init_par + sum(edf_spline) + emission_extra
  entropy <- entropy_from_gamma(gamma)
  aic <- -2 * loglik + 2 * npar_edf
  bic <- -2 * loglik + log(n_years) * npar_edf
  list(
    aic_edf = aic,
    bic_edf = bic,
    icl_edf = bic + 2 * entropy,
    entropy = entropy,
    npar_edf = npar_edf,
    edf_spline = edf_spline
  )
}


precision_safeguard <- function(phi, phi_upper, saturation_fraction = 0.98) {
  if (is.null(phi)) return(list(flagged = FALSE, max_phi = NA_real_, threshold = NA_real_))
  threshold <- saturation_fraction * phi_upper
  list(
    flagged = any(phi >= threshold),
    max_phi = max(phi),
    threshold = threshold
  )
}

state_occupancy <- function(gamma) {
  mass <- colSums(gamma)
  list(mass = mass, share = mass / nrow(gamma), min_share = min(mass / nrow(gamma)))
}

adjusted_rand_index <- function(x, y) {
  x <- as.factor(x)
  y <- as.factor(y)
  tab <- table(x, y)
  comb2 <- function(z) z * (z - 1) / 2
  a <- sum(comb2(rowSums(tab)))
  b <- sum(comb2(colSums(tab)))
  ab <- sum(comb2(tab))
  n2 <- comb2(sum(tab))
  if (n2 <= 0) return(1)
  expected <- a * b / n2
  max_index <- 0.5 * (a + b)
  denom <- max_index - expected
  if (abs(denom) < 1e-12) return(1)
  (ab - expected) / denom
}

bc_interval <- function(theta_hat, boot, level = 0.95) {
  boot <- boot[is.finite(boot)]
  if (length(boot) < 10L) return(c(NA_real_, NA_real_))
  alpha <- 1 - level
  p_less <- mean(boot < theta_hat)
  p_less <- min(max(p_less, 1 / (2 * length(boot))), 1 - 1 / (2 * length(boot)))
  z0 <- stats::qnorm(p_less)
  zlo <- stats::qnorm(alpha / 2)
  zhi <- stats::qnorm(1 - alpha / 2)
  probs <- stats::pnorm(c(2 * z0 + zlo, 2 * z0 + zhi))
  probs <- pmin(pmax(probs, 0), 1)
  as.numeric(stats::quantile(boot, probs = probs, na.rm = TRUE, names = FALSE, type = 8))
}

percentile_interval <- function(boot, level = 0.95) {
  alpha <- 1 - level
  as.numeric(stats::quantile(boot, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE, names = FALSE, type = 8))
}

basic_interval <- function(theta_hat, boot, level = 0.95) {
  boot <- boot[is.finite(boot)]
  if (length(boot) < 2L) return(c(NA_real_, NA_real_))
  alpha <- 1 - level
  q <- as.numeric(stats::quantile(boot, probs = c(alpha / 2, 1 - alpha / 2),
                                  na.rm = TRUE, names = FALSE, type = 8))
  c(2 * theta_hat - q[2L], 2 * theta_hat - q[1L])
}

bootstrap_interval <- function(theta_hat, boot, method = c("percentile", "basic", "bc"), level = 0.95) {
  method <- match.arg(method)
  if (method == "percentile") return(percentile_interval(boot, level = level))
  if (method == "basic") return(basic_interval(theta_hat, boot, level = level))
  bc_interval(theta_hat, boot, level = level)
}

write_latex_table <- function(df, path, caption = NULL, label = NULL, digits = 3L) {
  fmt <- function(x) {
    if (is.numeric(x)) return(ifelse(is.na(x), "--", formatC(x, digits = digits, format = "f")))
    x <- ifelse(is.na(x), "--", as.character(x))
    x <- gsub("_", "\\\\_", x, fixed = TRUE)
    x
  }
  dat <- as.data.frame(lapply(df, fmt), stringsAsFactors = FALSE)
  align <- paste0("l", paste(rep("r", max(ncol(dat) - 1L, 0L)), collapse = ""))
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines("\\begin{table}[!ht]", con)
  writeLines("\\centering", con)
  if (!is.null(caption)) writeLines(paste0("\\caption{", caption, "}"), con)
  if (!is.null(label)) writeLines(paste0("\\label{", label, "}"), con)
  writeLines(paste0("\\begin{tabular}{", align, "}"), con)
  writeLines("\\hline", con)
  writeLines(paste(paste(names(dat), collapse = " & "), "\\\\"), con)
  writeLines("\\hline", con)
  for (i in seq_len(nrow(dat))) {
    writeLines(paste(paste(dat[i, ], collapse = " & "), "\\\\"), con)
  }
  writeLines("\\hline", con)
  writeLines("\\end{tabular}", con)
  writeLines("\\end{table}", con)
  invisible(path)
}

safe_dir_create <- function(path) {
  if (!dir.exists(path)) suppressWarnings(dir.create(path, recursive = TRUE, showWarnings = FALSE))
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

parallel_lapply <- function(X, FUN, n_cores = 1L, project_root = getwd(), seed = 1L) {
  # Windows-safe PSOCK helper.  All model workers passed here are top-level
  # functions and each element of X contains the complete job specification.
  # This avoids relying on serialization of nested closures/environments,
  # which is fragile on PSOCK workers (especially under RStudio/Windows).
  n_cores <- max(1L, min(as.integer(n_cores), length(X)))
  if (n_cores <= 1L) return(lapply(X, FUN))

  project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  r_files <- sort(list.files(file.path(project_root, "R"), pattern = "\\.R$", full.names = TRUE))
  if (!length(r_files)) stop("No R modules found under: ", file.path(project_root, "R"))

  cl <- parallel::makePSOCKcluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterSetRNGStream(cl, iseed = seed)
  parallel::clusterExport(cl, "r_files", envir = environment())
  parallel::clusterEvalQ(cl, {
    for (f in r_files) source(f, local = .GlobalEnv)
    NULL
  })

  # parLapplyLB is deliberately used instead of pbapply on PSOCK clusters.
  # The grid itself prints progress, and load balancing is useful because
  # multi-start fits can have heterogeneous run times.
  parallel::parLapplyLB(cl, X, FUN)
}


# -----------------------------------------------------------------------------
# Scalable state alignment for bootstrap refits
# -----------------------------------------------------------------------------
# Dynamic-programming solution of the square linear assignment problem.
# Complexity is O(K 2^K), which is trivial for the K ranges used here and avoids
# factorial enumeration of all label permutations when K > 5.
min_cost_assignment_dp <- function(cost) {
  cost <- as.matrix(cost)
  K <- nrow(cost)
  if (K != ncol(cost)) stop("Assignment cost matrix must be square.")
  if (K == 1L) return(1L)
  if (K > 20L) stop("DP assignment is intended for K <= 20.")

  nmask <- bitwShiftL(1L, K)
  dp <- rep(Inf, nmask)
  prev_mask <- rep(NA_integer_, nmask)
  prev_col <- rep(NA_integer_, nmask)
  dp[1L] <- 0 # mask 0 is stored at index 1

  popcount <- function(mask) {
    sum(as.integer(intToBits(mask))[seq_len(K)])
  }

  for (mask in 0:(nmask - 1L)) {
    idx <- mask + 1L
    if (!is.finite(dp[idx])) next
    row <- popcount(mask) + 1L
    if (row > K) next
    for (col in seq_len(K)) {
      bit <- bitwShiftL(1L, col - 1L)
      if (bitwAnd(mask, bit) != 0L) next
      new_mask <- bitwOr(mask, bit)
      new_idx <- new_mask + 1L
      value <- dp[idx] + cost[row, col]
      if (value < dp[new_idx]) {
        dp[new_idx] <- value
        prev_mask[new_idx] <- mask
        prev_col[new_idx] <- col
      }
    }
  }

  perm <- integer(K)
  mask <- nmask - 1L
  for (row in K:1L) {
    idx <- mask + 1L
    col <- prev_col[idx]
    if (is.na(col)) stop("Assignment backtracking failed.")
    perm[row] <- col
    mask <- prev_mask[idx]
  }
  perm
}

align_fit_to_reference <- function(fit, reference_fit) {
  K <- reference_fit$K
  if (fit$K != K) return(fit)
  if (K == 1L) return(fit)

  ref_mu <- reference_fit$mu_age
  ref_phi <- reference_fit$phi
  cost <- matrix(0, K, K)
  for (new_state in seq_len(K)) {
    for (old_state in seq_len(K)) {
      curve_cost <- sum((fit$mu_age[, old_state] - ref_mu[, new_state])^2)
      phi_cost <- 0
      if (!is.null(ref_phi) && !is.null(fit$phi)) {
        scale_phi <- max(stats::sd(ref_phi), 1)
        phi_cost <- ((fit$phi[old_state] - ref_phi[new_state]) / scale_phi)^2
      }
      cost[new_state, old_state] <- curve_cost + 0.01 * phi_cost
    }
  }
  perm <- min_cost_assignment_dp(cost)
  relabel_fit(fit, perm)
}
