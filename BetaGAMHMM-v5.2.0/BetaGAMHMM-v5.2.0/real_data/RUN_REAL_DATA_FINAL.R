# Final three real-data analyses, run sequentially with identical protocol.
run_one <- function(country,label,start,end) {
  Sys.setenv(BGHMM_COUNTRY=country,BGHMM_LABEL=label,BGHMM_START_YEAR=as.character(start),BGHMM_END_YEAR=as.character(end),BGHMM_REAL_MODE="paper")
  source("RUN_COUNTRY.R", local=.GlobalEnv)
}
run_one("RUS","Russia",1960,2014)
run_one("ITA","Italy",1960,2014)
run_one("ITA","Italy extended",1950,2019)
