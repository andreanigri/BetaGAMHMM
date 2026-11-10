# Frozen configuration for the final scalar simulation study.

make_study_config <- function(mode = Sys.getenv("BGHMM_SCALAR_MODE", unset = "paper")) {
  mode <- match.arg(tolower(mode), c("paper", "quick"))
  quick <- identical(mode, "quick")
  diagnostic <- list(
    phi_max = 500, eps = 1e-6, s_thresh = 2L, delta_abs = 50, delta_sum = 100
  )
  override <- suppressWarnings(as.integer(Sys.getenv("BGHMM_N_CORES", unset = "")))
  ncores <- if (is.finite(override) && override >= 1L) override else min(12L, recommended_core_count())

  cfg <- list(
    mode = mode,
    seed = 20260921L, # frozen final seed
    output_root = file.path("outputs", paste0("s1_", mode)),
    parallel = list(n_cores = ncores, n_cores_fit_starts = 1L),
    basis = list(
      type = "bs", degree = 3L, n_inner_knots = 6L,
      boundary_padding_frac = 0.05, knot_quantile_range = c(0, 1)
    ),
    fit = list(
      K_grid = 1:5,
      lambda_grid = c(0.01, 0.05, 0.1, 0.5, 1),
      n_start = if (quick) 2L else 3L,
      max_iter = if (quick) 80L else 150L,
      tol = 1e-6, optimizer_maxit = if (quick) 60L else 80L,
      phi_lower = 1, phi_upper = 500
    ),
    diagnostic = diagnostic,
    scenarios = revision_simulation_scenarios(),
    primary_rule = "BIC",
    selection_study = list(
      n_rep = if (quick) 2L else 200L,
      chunk_size = if (quick) 2L else max(5L, ncores),
      scenarios = c("baseline", "short_series", "overlapping_states")
    ),
    basis_sensitivity = list(
      n_inner_knots = c(4L, 6L, 8L), n_start = if (quick) 2L else 5L
    ),
    multicovariate_followup = list(
      n_rep = if (quick) 2L else 30L, T = if (quick) 250L else 800L,
      n_start = if (quick) 2L else 5L
    ),
    smoothing_followup = list(
      T = if (quick) 300L else 1200L, K = 4L,
      n_start = if (quick) 2L else 4L, n_cycles = if (quick) 1L else 2L
    ),
    coverage_followup = list(
      n_rep = if (quick) 2L else 30L,
      B_boot = if (quick) 12L else 100L,
      T = if (quick) 250L else 700L, K = 4L, lambda = 0.1,
      n_start_fit = if (quick) 2L else 4L, n_start_boot = 2L
    ),
    hmmtmb_benchmark = list(
      n_rep = if (quick) 2L else 20L,
      T = if (quick) 250L else 700L, K = 4L, basis_k = 8L,
      lambda_grid = c(0.05, 0.1, 0.5, 1), n_start_ours = if (quick) 2L else 3L
    ),
    computational_benchmark = list(
      n_start = if (quick) 2L else 3L, max_iter = if (quick) 80L else 120L
    ),
    paper_reporting = list(
      # 100 fixed-K replications reproduce the original conditional-estimation table;
      # 200 replications are reserved for the reviewer-critical K-selection study.
      n_rep = if (quick) 2L else 100L,
      chunk_size = if (quick) 2L else max(5L, ncores),
      n_start = if (quick) 2L else 3L,
      max_iter = if (quick) 80L else 150L,
      x_grid_n = 250L,
      include_observed_points = TRUE
    )
  )
  cfg$selection_rules <- revision_selection_rules(cfg)
  cfg
}
