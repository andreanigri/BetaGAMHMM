# Final real-data configuration for the revised Beta-GAM HMM paper.
#
# The state-number criterion is the same BIC structure used in the original
# manuscript, now evaluated with the corrected Beta EDF and the correct HMM
# time unit T:
#   BIC_nu = -2 log L + nu log(T).
# Lambda is selected by minimum AIC_nu within each K; K is then selected by the
# minimum BIC_nu. No BIC tolerance rule and no diagnostic occupancy filter are
# used for state-number selection.

.env_chr <- function(name, default) {
  x <- Sys.getenv(name, unset = "")
  if (nzchar(x)) x else default
}
.env_int <- function(name, default) {
  x <- Sys.getenv(name, unset = "")
  if (!nzchar(x)) return(as.integer(default))
  out <- suppressWarnings(as.integer(x))
  if (!is.finite(out)) stop("Environment variable ", name, " must be an integer.")
  out
}
.env_bool <- function(name, default = FALSE) {
  x <- tolower(Sys.getenv(name, unset = if (default) "true" else "false"))
  x %in% c("1", "true", "yes", "y")
}

DATA_FILE <- .env_chr("BGHMM_DATA_FILE", file.path("data", "HMD_Deaths.zip"))
HMD_EXPECTED_LAST_MODIFIED <- .env_chr("BGHMM_HMD_EXPECTED_LAST_MODIFIED", "27 Aug 2026")
HMD_EXPECTED_DEATHS_MD5 <- .env_chr("BGHMM_HMD_EXPECTED_DEATHS_MD5", "08be9302bb6355fda6d9d89553507092")
HMD_ALLOW_VERSION_MISMATCH <- .env_bool("BGHMM_ALLOW_HMD_VERSION_MISMATCH", FALSE)
COUNTRY_CODE <- toupper(.env_chr("BGHMM_COUNTRY", "RUS"))
.default_labels <- c(RUS = "Russia", ITA = "Italy", USA = "United States", GBR_NP = "United Kingdom",
                     FRA = "France", DEUTNP = "Germany", ESP = "Spain", SWE = "Sweden", JPN = "Japan")
.default_label <- if (COUNTRY_CODE %in% names(.default_labels)) unname(.default_labels[COUNTRY_CODE]) else COUNTRY_CODE
APPLICATION_LABEL <- .env_chr("BGHMM_LABEL", .default_label)

YEAR_RANGE <- c(.env_int("BGHMM_START_YEAR", 1960L), .env_int("BGHMM_END_YEAR", 2014L))
AGE_RANGE <- c(.env_int("BGHMM_MIN_AGE", 0L), .env_int("BGHMM_MAX_AGE", 40L))
if (YEAR_RANGE[1L] > YEAR_RANGE[2L]) stop("Start year must not exceed end year.")
if (AGE_RANGE[1L] > AGE_RANGE[2L]) stop("Minimum age must not exceed maximum age.")

REAL_MODE <- tolower(.env_chr("BGHMM_REAL_MODE", "paper"))
if (!REAL_MODE %in% c("paper", "quick")) stop("BGHMM_REAL_MODE must be paper or quick.")

# Final, pre-specified state-selection workflow.
K_GRID <- 1:8
FINAL_K_AUTO_MAX <- 10L
LAMBDA_GRID <- c(0.1, 0.5, 1, 5, 10)

SPLINE_DEGREE <- 3L
N_INNER_KNOTS <- 6L
PHI_LOWER <- 1
PHI_UPPER <- 2000
PRECISION_SATURATION_FRACTION <- 0.98
BETA_BOUNDARY_ACTION <- .env_chr("BGHMM_BETA_BOUNDARY_ACTION", "error")
if (!BETA_BOUNDARY_ACTION %in% c("error", "continuity")) stop("BGHMM_BETA_BOUNDARY_ACTION must be 'error' or 'continuity'.")

PROFILE_N_START <- .env_int("BGHMM_PROFILE_N_START", 12L)
MAX_ITER <- 300L
FIT_TOL <- 1e-6
OPTIMIZER_MAXIT <- 150L
SCALAR_N_START <- .env_int("BGHMM_SCALAR_N_START", 15L)

RUN_BASIS_SENSITIVITY <- TRUE
BASIS_INNER_KNOTS_GRID <- c(4L, 6L, 8L)
BASIS_SENSITIVITY_N_START <- .env_int("BGHMM_BASIS_N_START", 6L)
RUN_DENOMINATOR_SENSITIVITY <- TRUE
DENOMINATOR_N_START <- .env_int("BGHMM_DENOM_N_START", 10L)

# Targeted precision-bound sensitivity. It reuses the primary AIC-selected
# lambda within each K, avoiding a second full K x lambda grid.
RUN_PHI_BOUND_SENSITIVITY <- .env_bool("BGHMM_RUN_PHI_BOUND_SENSITIVITY", TRUE)
PHI_BOUND_SENSITIVITY <- list(relaxed = c(lower = 0.1, upper = 5000))
PHI_BOUND_N_START <- .env_int("BGHMM_PHI_BOUND_N_START", 4L)

# Expensive inference is optional in the generic country run. Russia wrappers
# can enable it after model selection; the Italy robustness analysis defaults
# to selection/diagnostics only.
RUN_FINAL_BOOTSTRAP <- .env_bool("BGHMM_RUN_FINAL_BOOTSTRAP", FALSE)
BOOTSTRAP_B <- .env_int("BGHMM_BOOTSTRAP_B", 200L)
BOOTSTRAP_N_START <- .env_int("BGHMM_BOOTSTRAP_NSTART", 6L)
# Lightweight numerical safeguard: only failed conditional-bootstrap refits are
# retried once, with a modestly larger multistart budget. Successful refits are
# never repeated.
BOOTSTRAP_RETRY_FAILED <- .env_bool("BGHMM_BOOTSTRAP_RETRY_FAILED", TRUE)
BOOTSTRAP_RETRY_N_START <- .env_int("BGHMM_BOOTSTRAP_RETRY_NSTART", 8L)
BOOTSTRAP_LEVEL <- 0.95
# Frozen manuscript convention, specified before reading Simulation-2 coverage.
# All three formulas are retained from the same bootstrap refits for sensitivity.
BOOTSTRAP_INTERVAL_METHODS <- list(
  phi = tolower(.env_chr("BGHMM_CI_PHI", "bc")),
  transition = tolower(.env_chr("BGHMM_CI_TRANSITION", "percentile")),
  initial = tolower(.env_chr("BGHMM_CI_INITIAL", "percentile")),
  curve = tolower(.env_chr("BGHMM_CI_CURVE", "percentile"))
)
.valid_ci_methods <- c("percentile", "basic", "bc")
if (any(!unlist(BOOTSTRAP_INTERVAL_METHODS, use.names = FALSE) %in% .valid_ci_methods))
  stop("Bootstrap CI methods must be one of: percentile, basic, bc.")
SELECTION_BOOTSTRAP_B <- .env_int("BGHMM_SELECTION_BOOT_B", 100L)
SELECTION_BOOTSTRAP_N_START <- .env_int("BGHMM_SELECTION_BOOT_NSTART", 5L)

if (REAL_MODE == "quick") {
  K_GRID <- 1:3
  FINAL_K_AUTO_MAX <- 3L
  LAMBDA_GRID <- c(0.5, 1)
  PROFILE_N_START <- 2L
  SCALAR_N_START <- 2L
  BASIS_INNER_KNOTS_GRID <- c(4L, 6L)
  BASIS_SENSITIVITY_N_START <- 2L
  DENOMINATOR_N_START <- 2L
  PHI_BOUND_N_START <- 2L
  if (RUN_FINAL_BOOTSTRAP) {
    BOOTSTRAP_B <- 5L
    BOOTSTRAP_N_START <- 2L
    BOOTSTRAP_RETRY_N_START <- 3L
    SELECTION_BOOTSTRAP_B <- 5L
    SELECTION_BOOTSTRAP_N_START <- 2L
  }
}

SEED <- 20260916L + ifelse(COUNTRY_CODE == "ITA", 100000L, 0L)
OUTPUT_ROOT <- file.path("outputs", sprintf("%s_%s_%d_%d_a%d_%d",
  ifelse(REAL_MODE == "quick", "q", "p"), COUNTRY_CODE,
  YEAR_RANGE[1L], YEAR_RANGE[2L], AGE_RANGE[1L], AGE_RANGE[2L]))
