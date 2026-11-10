# Simulation utilities used in the main paper and revision studies.

default_eta_functions <- function(K, scale = 1) {
  funs <- list(
    function(x) scale * (-0.9 + 0.8 * sin(pi * x / 2)),
    function(x) scale * ( 0.5 + 0.6 * x - 0.35 * x^2),
    function(x) scale * (-0.15 + 1.05 * exp(-x^2) - 0.20 * x),
    function(x) scale * ( 0.55 * cos(pi * x / 2) + 0.18 * x),
    function(x) scale * (-0.4 - 0.65 * sin(pi * x) + 0.10 * x^2)
  )
  if (K > length(funs)) {
    extra <- lapply(seq_len(K - length(funs)), function(j) {
      force(j)
      function(x) scale * (0.2 * cos(j * pi * x / 2) - 0.1 * j * x)
    })
    funs <- c(funs, extra)
  }
  funs[seq_len(K)]
}

project_eta_functions_to_basis <- function(x, basis_spec, eta_funs) {
  B <- basis_spec$B
  lapply(eta_funs, function(f) qr.solve(B, f(x), tol = 1e-10))
}

make_transition_matrix <- function(K, diagonal = 0.90) {
  if (K == 1L) return(matrix(1, 1L, 1L))
  off <- (1 - diagonal) / (K - 1L)
  A <- matrix(off, nrow = K, ncol = K)
  diag(A) <- diagonal
  A
}

# An asymmetric transition matrix is useful for testing rare or short-lived
# regimes, providing a deliberately imbalanced stress scenario.
make_imbalanced_transition_matrix <- function(K = 4L) {
  if (K != 4L) stop("The default imbalanced transition design is defined for K = 4.")
  normalize_rows(matrix(c(
    0.965, 0.020, 0.010, 0.005,
    0.090, 0.865, 0.030, 0.015,
    0.100, 0.060, 0.820, 0.020,
    0.140, 0.060, 0.030, 0.770
  ), nrow = 4L, byrow = TRUE))
}

simulate_from_parameters <- function(x,
                                     basis_spec,
                                     beta_list,
                                     phi,
                                     A,
                                     pi = NULL,
                                     seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  B <- basis_spec$B
  T <- nrow(B)
  K <- length(phi)
  pi <- pi %||% stationary_distribution(A)

  z <- integer(T)
  z[1L] <- sample.int(K, 1L, prob = pi)
  if (T >= 2L) {
    for (t in 2:T) z[t] <- sample.int(K, 1L, prob = A[z[t - 1L], ])
  }

  eta_by_state <- vapply(seq_len(K), function(k) drop(B %*% beta_list[[k]]), numeric(T))
  if (K == 1L) eta_by_state <- matrix(eta_by_state, ncol = 1L)
  mu_by_state <- inv_logit(eta_by_state)
  mu <- mu_by_state[cbind(seq_len(T), z)]
  y <- stats::rbeta(T, mu * phi[z], (1 - mu) * phi[z])

  list(
    x = x,
    y = y,
    z = z,
    B = B,
    basis_spec = basis_spec,
    eta_by_state = eta_by_state,
    mu_by_state = mu_by_state,
    mu = mu,
    beta_true = beta_list,
    phi_true = phi,
    A_true = A,
    pi_true = pi,
    K_true = K
  )
}

simulate_default_scenario <- function(T = 1200L,
                                      K = 4L,
                                      degree = 3L,
                                      n_inner_knots = 6L,
                                      diagonal = 0.92,
                                      phi_true = c(10, 18, 25, 34),
                                      seed = 1L,
                                      basis_type = "bs",
                                      boundary_padding_frac = 0,
                                      knot_quantile_range = c(0, 1),
                                      eta_scale = 1,
                                      A_true = NULL) {
  set.seed(seed)
  x <- stats::runif(T, -2, 2)
  basis_spec <- build_bspline_basis(
    x = x,
    degree = degree,
    n_inner_knots = n_inner_knots,
    type = basis_type,
    boundary_padding_frac = boundary_padding_frac,
    knot_quantile_range = knot_quantile_range
  )
  eta_funs <- default_eta_functions(K, scale = eta_scale)
  beta_true <- project_eta_functions_to_basis(x, basis_spec, eta_funs)
  if (length(phi_true) != K) {
    if (length(phi_true) > K) {
      phi_true <- phi_true[seq_len(K)]
    } else {
      phi_true <- c(phi_true, rep(max(phi_true), K - length(phi_true)))
    }
  }
  user_supplied_A <- !is.null(A_true)
  if (is.null(A_true)) A_true <- make_transition_matrix(K, diagonal = diagonal)
  pi_true <- stationary_distribution(A_true)
  sim <- simulate_from_parameters(
    x = x,
    basis_spec = basis_spec,
    beta_list = beta_true,
    phi = phi_true,
    A = A_true,
    pi = pi_true,
    seed = seed + 1L
  )
  sim$scenario <- list(
    T = T, K = K, diagonal = diagonal, phi_true = phi_true,
    eta_scale = eta_scale, asymmetric_A = user_supplied_A
  )
  sim
}

# Two-covariate design used to verify that the additive implementation scales
# beyond the single-covariate examples in the original manuscript.
simulate_multicovariate_scenario <- function(T = 800L,
                                             K = 3L,
                                             degree = 3L,
                                             n_inner_knots = 5L,
                                             diagonal = 0.90,
                                             phi_true = c(12, 22, 36),
                                             seed = 101L,
                                             boundary_padding_frac = 0.05) {
  set.seed(seed)
  dat <- data.frame(
    x1 = stats::runif(T, -2, 2),
    x2 = stats::runif(T, -1, 1)
  )
  basis_spec <- build_additive_basis(
    dat,
    covariates = c("x1", "x2"),
    degree = degree,
    n_inner_knots = n_inner_knots,
    boundary_padding_frac = boundary_padding_frac
  )

  eta_true <- cbind(
    -0.6 + 0.75 * sin(pi * dat$x1 / 2) + 0.35 * dat$x2,
     0.2 + 0.45 * dat$x1 - 0.20 * dat$x1^2 - 0.50 * dat$x2^2,
     0.5 * cos(pi * dat$x1 / 2) + 0.55 * dat$x2
  )
  if (K != 3L) stop("The default two-covariate scenario is defined for K = 3.")
  beta_true <- lapply(seq_len(K), function(k) qr.solve(basis_spec$B, eta_true[, k], tol = 1e-10))
  A_true <- make_transition_matrix(K, diagonal = diagonal)
  pi_true <- stationary_distribution(A_true)
  sim <- simulate_from_parameters(
    x = dat,
    basis_spec = basis_spec,
    beta_list = beta_true,
    phi = phi_true,
    A = A_true,
    pi = pi_true,
    seed = seed + 1L
  )
  sim$scenario <- list(T = T, K = K, p = 2L, diagonal = diagonal, phi_true = phi_true)
  sim
}

# Scenarios used in the repeated K-selection study.  They deliberately include
# a shorter series and more overlap rather than only a long, symmetric, well-separated benchmark.
revision_simulation_scenarios <- function() {
  list(
    baseline = list(
      T = 1200L, K_true = 4L, diagonal = 0.92,
      phi_true = c(10, 18, 28, 40), eta_scale = 1.00, A_true = NULL
    ),
    short_series = list(
      T = 300L, K_true = 4L, diagonal = 0.85,
      phi_true = c(10, 18, 28, 40), eta_scale = 0.85, A_true = NULL
    ),
    overlapping_states = list(
      T = 700L, K_true = 4L, diagonal = 0.88,
      phi_true = c(14, 18, 22, 27), eta_scale = 0.65, A_true = NULL
    )
  )
}
