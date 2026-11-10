# Run expensive inference after the three country analyses are complete.
# The manuscript interval convention is already pre-specified in config; Simulation 2
# validates its coverage rather than selecting it post hoc.
run_inf <- function(country,label,start,end) {
  Sys.setenv(BGHMM_COUNTRY=country,BGHMM_LABEL=label,BGHMM_START_YEAR=as.character(start),BGHMM_END_YEAR=as.character(end),BGHMM_REAL_MODE="paper")
  source("RUN_FINAL_INFERENCE.R", local=.GlobalEnv)
}
run_inf("RUS","Russia",1960,2014)
run_inf("ITA","Italy",1960,2014)
run_inf("ITA","Italy extended",1950,2019)
