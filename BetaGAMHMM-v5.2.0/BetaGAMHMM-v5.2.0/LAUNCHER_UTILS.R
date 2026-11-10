# Root-level launcher utilities for the Beta-GAM HMM public reproducibility repository.
# These functions let users run any subproject without changing the R working directory manually.

.bghmm_norm <- function(x, mustWork = TRUE) {
  normalizePath(x, winslash = "/", mustWork = mustWork)
}

.bghmm_find_root <- function() {
  candidates <- character()

  # Prefer paths of files currently being sourced. This makes launchers self-locating.
  if (sys.nframe() > 0L) {
    for (i in rev(seq_len(sys.nframe()))) {
      of <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
      if (!is.null(of) && length(of) && nzchar(of[1L])) {
        p <- tryCatch(.bghmm_norm(of[1L], mustWork = TRUE), error = function(e) NA_character_)
        if (!is.na(p)) candidates <- c(candidates, dirname(p))
      }
    }
  }

  # Also search from the current working directory and its parents.
  cur <- tryCatch(.bghmm_norm(getwd(), mustWork = TRUE), error = function(e) getwd())
  repeat {
    candidates <- c(candidates, cur)
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }

  candidates <- unique(candidates)
  for (d in candidates) {
    marker <- file.path(d, "BetaGAMHMM.Rproj")
    ok_dirs <- all(dir.exists(file.path(d, c("simulation_1_scalar", "simulation_2_profile", "real_data"))))
    if (file.exists(marker) && ok_dirs) return(.bghmm_norm(d, mustWork = TRUE))
  }

  stop(
    paste0(
      "Cannot locate the BetaGAMHMM bundle root.\n",
      "Source this launcher from inside the extracted repository, or open BetaGAMHMM.Rproj first."
    ),
    call. = FALSE
  )
}

.bghmm_run_project <- function(project, scripts, label = project, root = NULL) {
  if (is.null(root)) {
    root <- .bghmm_find_root()
  } else {
    root <- .bghmm_norm(root, mustWork = TRUE)
    marker <- file.path(root, "BetaGAMHMM.Rproj")
    if (!file.exists(marker)) stop("Invalid bundle root: ", root, call. = FALSE)
  }
  project_dir <- file.path(root, project)
  if (!dir.exists(project_dir)) stop("Missing project directory: ", project_dir, call. = FALSE)

  oldwd <- getwd()
  on.exit(setwd(oldwd), add = TRUE)
  setwd(project_dir)

  cat("\n============================================================\n")
  cat("BetaGAMHMM: ", label, "\n", sep = "")
  cat("Bundle root: ", root, "\n", sep = "")
  cat("Project dir: ", .bghmm_norm(project_dir, mustWork = TRUE), "\n", sep = "")
  cat("============================================================\n")

  for (s in scripts) {
    if (!file.exists(s)) stop("Required entry-point not found: ", file.path(project_dir, s), call. = FALSE)
    cat("\n>>> source(\"", s, "\")\n", sep = "")
    source(s, local = .GlobalEnv, chdir = FALSE)
  }

  cat("\nPASS: ", label, " launcher completed.\n", sep = "")
  invisible(TRUE)
}

.bghmm_run_quick <- function(project, quick_script, label, root = NULL) {
  .bghmm_run_project(
    project = project,
    scripts = c("CHECK_CODE.R", "CHECK_REPRODUCIBILITY.R", quick_script),
    label = paste0(label, " — checks + quick run"),
    root = root
  )
}

.bghmm_run_final <- function(project, final_script, label, hardware = FALSE, root = NULL) {
  scripts <- c("CHECK_CODE.R")
  if (hardware) scripts <- c(scripts, "CHECK_HARDWARE.R")
  scripts <- c(scripts, "CHECK_REPRODUCIBILITY.R", final_script)
  .bghmm_run_project(project, scripts, paste0(label, " — preflight + final run"), root = root)
}
