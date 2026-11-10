# Simulation designs for repeated proportion profiles.
#
# At HMM time t we observe J proportions sharing one latent state. The profile
# coordinate is x_j (age in the mortality-style design). The truth is analytical,
# rather than generated from the exact spline basis used for estimation.

make_profile_transition <- function(K, refresh_persistence = 0.70, weights = NULL) {
  # Sticky-refresh chain. With probability p the process stays in its current
  # state; otherwise it refreshes from `weights`. Therefore the actual diagonal
  # is p + (1-p) * weight_k, while `weights` is the stationary distribution.
  if (is.null(weights)) weights <- rep(1 / K, K)
  weights <- normalize_vec(weights)
  p <- as.numeric(refresh_persistence)
  if (!is.finite(p) || p < 0 || p >= 1) stop("refresh_persistence must lie in [0,1).")
  A <- matrix(rep((1 - p) * weights, each = K), K, K, byrow = FALSE)
  diag(A) <- diag(A) + p
  A <- normalize_rows(A)
  list(A = A, pi = weights, stationary = weights,
       refresh_persistence = p, diagonal = diag(A))
}

profile_true_curves <- function(age, K, separation = 1) {
  age <- as.numeric(age)
  z <- (age - mean(range(age))) / diff(range(age))
  base_eta <- stats::qlogis(0.50) -
    1.05 * exp(-((age - 25) / 8.5)^2) +
    0.16 * exp(-((age - 2) / 3.5)^2) +
    0.08 * z

  u <- if (K == 1L) 0 else seq(-1, 1, length.out = K)
  mu <- vapply(seq_len(K), function(k) {
    eta <- base_eta + separation * (
      0.32 * u[k] +
      0.28 * u[k] * exp(-((age - 24) / 9)^2) +
      0.07 * sin((age + 3 * k) / 8) +
      0.06 * u[k] * z
    )
    clip_unit(inv_logit(eta), 1e-6)
  }, numeric(length(age)))
  if (K == 1L) mu <- matrix(mu, ncol = 1L)
  mu
}

profile_scenarios <- function(mode = c("paper", "quick")) {
  mode <- match.arg(mode)
  nrep <- if (mode == "paper") 200L else 2L

  # The final profile study keeps K=4 and J=41 fixed. It varies temporal
  # information (T) and includes one intentionally adverse separation/occupancy
  # scenario. This directly addresses the reviewers' concerns about the original
  # favourable, very-long-series simulation design.
  list(
    baseline_T55 = list(
      name = "baseline_T55", K = 4L, T = 55L, J = 41L,
      refresh_persistence = 0.70, weights = rep(0.25, 4), separation = 1.00,
      phi = c(120, 170, 230, 300), n_rep = nrep
    ),
    baseline_T250 = list(
      name = "baseline_T250", K = 4L, T = 250L, J = 41L,
      refresh_persistence = 0.70, weights = rep(0.25, 4), separation = 1.00,
      phi = c(120, 170, 230, 300), n_rep = nrep
    ),
    overlap_rare_T55 = list(
      name = "overlap_rare_T55", K = 4L, T = 55L, J = 41L,
      # 8% stationary mass remains genuinely rare, but the refresh rate is high
      # enough that the state is usually observed in T=55. This tests rare
      # occupancy rather than recovery of a latent component never seen in data.
      refresh_persistence = 0.40, weights = c(0.34, 0.31, 0.27, 0.08), separation = 0.45,
      phi = c(120, 150, 190, 240), n_rep = nrep
    )
  )
}

simulate_profile_scenario <- function(scenario, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  K <- scenario$K
  Tn <- scenario$T
  J <- scenario$J
  ages <- seq(0, 40, length.out = J)
  trans <- make_profile_transition(K, scenario$refresh_persistence, scenario$weights)
  mu_true <- profile_true_curves(ages, K, scenario$separation)
  phi_true <- scenario$phi
  if (length(phi_true) != K) stop("phi length must equal K.")

  z <- integer(Tn)
  z[1L] <- sample.int(K, 1L, prob = trans$pi)
  if (Tn >= 2L) {
    for (t in 2:Tn) z[t] <- sample.int(K, 1L, prob = trans$A[z[t - 1L], ])
  }

  Y <- matrix(NA_real_, Tn, J)
  for (t in seq_len(Tn)) {
    k <- z[t]
    Y[t, ] <- stats::rbeta(J, mu_true[, k] * phi_true[k], (1 - mu_true[, k]) * phi_true[k])
  }

  panel <- list(
    Y = Y, Y_raw = Y,
    years = seq_len(Tn), ages = ages,
    n_years = Tn, n_ages = J, n_cells = Tn * J,
    boundary_cells = matrix(integer(0), ncol = 2L),
    application_label = scenario$name, country_code = "SIM"
  )

  state_counts <- tabulate(z, nbins = K)
  list(
    panel = panel, z = z, mu_true = mu_true, phi_true = phi_true,
    A_true = trans$A, pi_true = trans$pi, scenario = scenario,
    transition_diagonal_true = trans$diagonal,
    true_state_counts = state_counts,
    true_state_shares = state_counts / Tn,
    n_states_realized = sum(state_counts > 0L),
    all_states_visited = all(state_counts > 0L),
    min_true_state_count = min(state_counts),
    min_true_state_share = min(state_counts / Tn)
  )
}

profile_scenario_design_summary <- function(scenarios) {
  rows <- lapply(scenarios, function(sc) {
    tr <- make_profile_transition(sc$K, sc$refresh_persistence, sc$weights)
    no_visit <- (1 - tr$pi) * (1 - (1 - tr$refresh_persistence) * tr$pi)^(sc$T - 1L)
    data.frame(
      scenario = sc$name,
      T = sc$T,
      K = sc$K,
      state = seq_len(sc$K),
      stationary_share = tr$pi,
      expected_state_count = sc$T * tr$pi,
      transition_diagonal = diag(tr$A),
      implied_mean_sojourn = 1 / pmax(1 - diag(tr$A), 1e-12),
      theoretical_no_visit_probability = no_visit,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
