old_mode <- Sys.getenv("BGHMM_PROFILE_MODE", unset="")
Sys.setenv(BGHMM_PROFILE_MODE="quick")
tryCatch(source("RUN_FINAL.R", local=.GlobalEnv), finally={if(nzchar(old_mode)) Sys.setenv(BGHMM_PROFILE_MODE=old_mode) else Sys.unsetenv("BGHMM_PROFILE_MODE")})
rm(old_mode)
