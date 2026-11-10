# Optional reference-environment installer for the manuscript software versions.
# R 4.4.2 is recommended. On Windows, archived source packages may require Rtools.

options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

reference <- c(
  ggplot2 = "4.0.2",
  pbapply = "1.7-4",
  mgcv = "1.9-3",
  TMB = "1.9.17",
  tmbstan = "1.0.91",
  hmmTMB = "1.1.2"
)

for (pkg in names(reference)) {
  target <- unname(reference[[pkg]])
  current <- if (requireNamespace(pkg, quietly = TRUE)) as.character(utils::packageVersion(pkg)) else NA_character_
  if (identical(current, target)) {
    message(pkg, " ", target, ": already installed")
  } else {
    message("Installing ", pkg, " ", target, if (!is.na(current)) paste0(" (current: ", current, ")") else "")
    remotes::install_version(pkg, version = target, dependencies = TRUE, upgrade = "never")
  }
}

message("Reference package installation/check complete.")
