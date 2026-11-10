.this_ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
.this_file <- if (!is.null(.this_ofile) && length(.this_ofile) && nzchar(.this_ofile[1L])) {
  tryCatch(normalizePath(.this_ofile[1L], winslash = "/", mustWork = TRUE), error = function(e) NA_character_)
} else NA_character_
.this_dir <- if (length(.this_file) == 1L && !is.na(.this_file)) dirname(.this_file) else getwd()
source(file.path(.this_dir, "LAUNCHER_UTILS.R"), local = .GlobalEnv)
.bghmm_run_final("real_data", "RUN_REAL_DATA_FINAL.R", "Real data", hardware = TRUE, root = .this_dir)
rm(.this_ofile, .this_file, .this_dir)
