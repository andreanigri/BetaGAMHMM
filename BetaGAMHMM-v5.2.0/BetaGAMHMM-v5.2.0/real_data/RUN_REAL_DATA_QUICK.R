# Quick smoke analyses for the two 1960-2014 country applications.
run_one_quick <- function(country, label) {
  Sys.setenv(
    BGHMM_COUNTRY = country,
    BGHMM_LABEL = label,
    BGHMM_START_YEAR = "1960",
    BGHMM_END_YEAR = "2014",
    BGHMM_REAL_MODE = "quick"
  )
  source("RUN_COUNTRY.R", local = .GlobalEnv)
}
run_one_quick("RUS", "Russia")
run_one_quick("ITA", "Italy")
rm(run_one_quick)
