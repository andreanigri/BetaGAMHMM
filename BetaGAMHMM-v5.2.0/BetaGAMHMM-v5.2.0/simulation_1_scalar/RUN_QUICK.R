old_mode <- Sys.getenv("BGHMM_SCALAR_MODE", unset="")
Sys.setenv(BGHMM_SCALAR_MODE="quick")
tryCatch(source("RUN_FINAL.R", local=.GlobalEnv), finally={if(nzchar(old_mode)) Sys.setenv(BGHMM_SCALAR_MODE=old_mode) else Sys.unsetenv("BGHMM_SCALAR_MODE")})
rm(old_mode)
