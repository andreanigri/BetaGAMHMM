# Complete from-scratch manuscript workflow. On a fresh repository clone this
# executes every official analysis in the intended order.
.this_ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
.this_file <- if (!is.null(.this_ofile) && length(.this_ofile) && nzchar(.this_ofile[1L])) {
  tryCatch(normalizePath(.this_ofile[1L], winslash = "/", mustWork = TRUE), error = function(e) NA_character_)
} else NA_character_
.root <- if (length(.this_file) == 1L && !is.na(.this_file)) dirname(.this_file) else getwd()

entries <- c(
  "RUN_FINAL_SCALAR.R",
  "RUN_FINAL_PROFILE.R",
  "RUN_FINAL_REAL_DATA.R",
  "RUN_FINAL_INFERENCE_REAL_DATA.R"
)
for (entry in entries) {
  cat("\n\n############################################################\n")
  cat("### ", entry, "\n", sep = "")
  cat("############################################################\n")
  source(file.path(.root, entry), local = .GlobalEnv, chdir = FALSE)
}
cat("\n\n############################################################\n")
cat("### redraw_profile_figures.R (published state-profile figures)\n")
cat("############################################################\n")
.old_wd <- setwd(.root)
source(file.path(.root, "redraw_profile_figures.R"), local = new.env(), chdir = FALSE)
setwd(.old_wd)
rm(.old_wd)
cat("\nPASS: complete BetaGAMHMM final workflow finished.\n")
rm(.this_ofile, .this_file, .root, entries, entry)
