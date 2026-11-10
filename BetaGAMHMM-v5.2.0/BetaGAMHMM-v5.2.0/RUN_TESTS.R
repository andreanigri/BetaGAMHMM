# Lightweight invariant tests for the two model cores.
.run_tests <- function() {
  this_ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  this_file <- if (!is.null(this_ofile) && length(this_ofile) && nzchar(this_ofile[1L])) {
    tryCatch(normalizePath(this_ofile[1L], winslash = "/", mustWork = TRUE), error = function(e) NA_character_)
  } else NA_character_
  root <- if (length(this_file) == 1L && !is.na(this_file)) dirname(this_file) else getwd()
  oldwd <- getwd(); on.exit(setwd(oldwd), add = TRUE)

  cat("Running scalar invariant tests...\n")
  setwd(file.path(root, "simulation_1_scalar"))
  scalar_env <- new.env(parent = globalenv())
  for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f, local = scalar_env)
  source(file.path("tests", "test_core_invariants.R"), local = scalar_env)

  cat("Running profile invariant tests...\n")
  setwd(file.path(root, "simulation_2_profile"))
  profile_env <- new.env(parent = globalenv())
  for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f, local = profile_env)
  source(file.path("tests", "test_profile_invariants.R"), local = profile_env)

  cat("PASS: invariant tests completed.\n")
  invisible(TRUE)
}
.run_tests()
rm(.run_tests)
