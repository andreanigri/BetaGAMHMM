# Frozen configuration for the final repeated-profile simulation study.
#
# Modes:
#   paper (default): manuscript-quality Monte Carlo study
#   quick:           smoke test only; never use quick outputs in the manuscript

profile_auto_cores <- function(max_workers = 8L) {
  override <- suppressWarnings(as.integer(Sys.getenv("BGHMM_N_CORES", unset = "")))
  if (is.finite(override) && override >= 1L) return(as.integer(override))
  detected <- parallel::detectCores(logical = FALSE)
  if (!is.finite(detected)) detected <- parallel::detectCores(logical = TRUE)
  if (!is.finite(detected)) detected <- 2L
  max(1L, min(as.integer(max_workers), as.integer(detected) - 1L))
}

make_profile_config <- function(mode = Sys.getenv("BGHMM_PROFILE_MODE", unset = "paper")) {
  mode <- match.arg(tolower(mode), c("paper", "quick"))
  quick <- identical(mode, "quick")
  list(
    mode = mode,
    seed = 20260920L, # frozen final seed
    output_root = file.path("outputs", paste0("s2_", mode)),
    n_cores = profile_auto_cores(12L),
    scenarios = profile_scenarios(mode),
    basis = list(degree = 3L, n_inner_knots = 6L),
    fit = list(
      K_grid = 1:8,
      lambda_grid = c(0.1, 0.5, 1, 5, 10),
      n_start = if (quick) 2L else 5L,
      max_iter = if (quick) 80L else 180L,
      tol = 1e-6,
      optimizer_maxit = if (quick) 60L else 100L,
      phi_lower = 1,
      phi_upper = 800,
      precision_saturation_fraction = 0.98
    ),
    # Paper-mode checkpoint chunks contain at least one full worker wave so
    # an i9-class CPU is not artificially limited to five concurrent replicates.
    selection = list(chunk_size = if (quick) 2L else max(5L, profile_auto_cores(12L))),
    coverage = list(
      n_outer = if (quick) 2L else 30L,
      B_boot = if (quick) 12L else 100L,
      n_start_fit = if (quick) 2L else 4L,
      n_start_boot = 2L,
      level = 0.95
    ),
    basis_sensitivity = list(
      n_rep = if (quick) 2L else 30L,
      n_inner_knots = c(4L, 6L, 8L),
      n_start = if (quick) 2L else 3L
    ),
    benchmark = list(
      n_repeat = if (quick) 1L else 3L,
      n_start = if (quick) 2L else 3L
    )
  )
}
