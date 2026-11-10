# Deliberately remove generated outputs from all subprojects.
if (!identical(Sys.getenv("BGHMM_CONFIRM_RESET", unset = ""), "YES")) {
  stop('Refusing to delete outputs. First run: Sys.setenv(BGHMM_CONFIRM_RESET = "YES")', call. = FALSE)
}
.this_ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
.this_file <- if (!is.null(.this_ofile) && length(.this_ofile) && nzchar(.this_ofile[1L])) {
  tryCatch(normalizePath(.this_ofile[1L], winslash = "/", mustWork = TRUE), error = function(e) NA_character_)
} else NA_character_
.root <- if (length(.this_file) == 1L && !is.na(.this_file)) dirname(.this_file) else getwd()
for (project in c("simulation_1_scalar", "simulation_2_profile", "real_data")) {
  out <- file.path(.root, project, "outputs")
  if (dir.exists(out)) unlink(out, recursive = TRUE, force = TRUE)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  file.create(file.path(out, ".gitkeep"))
}
message("All generated outputs were removed. The repository is ready for a from-scratch run.")
rm(.this_ofile, .this_file, .root)
