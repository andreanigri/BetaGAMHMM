# Expensive final inference for the already-selected country model.
# Uses the single pre-specified BIC rule from RUN_COUNTRY.R.
# Run after the corresponding final country analysis has completed.

options(stringsAsFactors = FALSE)
# Keep bootstrap parallelism at the outer-refit level; nested BLAS/OpenMP
# threading usually oversubscribes an i9 and slows the run.
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
           VECLIB_MAXIMUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
if (!file.exists(file.path("config","real_data_config.R")) || !dir.exists("R")) stop("Run from project root.")
PROJECT_ROOT <- normalizePath(getwd(),winslash="/",mustWork=TRUE)
source(file.path("config","real_data_config.R"),local=.GlobalEnv)
for (f in sort(list.files("R",pattern="\\.R$",full.names=TRUE))) source(f,local=.GlobalEnv)
if (REAL_MODE != "paper") stop("RUN_FINAL_INFERENCE.R requires paper mode.")
obj_path <- file.path(OUTPUT_ROOT,"analysis_object.rds")
if (!file.exists(obj_path)) stop("Run the country analysis first: analysis_object.rds not found.")
obj <- readRDS(obj_path)
if (is.null(obj$phi_bound_sensitivity) || is.null(obj$phi_bound_sensitivity$summary)) {
  stop(
    "The country analysis object is incomplete. Run the final real-data selection workflow before final inference.",
    call. = FALSE
  )
}
sel <- obj$profile_selection; panel <- obj$panel
hw <- choose_auto_workers(); N_CORES <- hw$selected_workers
if (N_CORES > 1L && !test_parallel_backend(N_CORES,PROJECT_ROOT,SEED+21L)) N_CORES <- 1L
FINAL_DIR <- safe_dir_create(file.path(OUTPUT_ROOT,"final_inference")); CHECKPOINT_DIR <- safe_dir_create(file.path(FINAL_DIR,"checkpoints"))
FORCE_RERUN <- identical(tolower(Sys.getenv("FORCE_RERUN","false")),"true")
load_or_run_final <- function(path,label,code) { if(file.exists(path)&&!FORCE_RERUN){message(label,": loading checkpoint");return(readRDS(path))}; message(label); x<-force(code); saveRDS(x,path); x }
message("Final inference: ",APPLICATION_LABEL,"; selected K=",sel$fit$K,", lambda=",sel$fit$lambda[1L])
selection_stability <- load_or_run_final(file.path(CHECKPOINT_DIR,"selection_stability.rds"),
  paste0("Selection-stability bootstrap B=",SELECTION_BOOTSTRAP_B),
  run_selection_stability_bootstrap(sel$fit,SELECTION_BOOTSTRAP_B,sort(unique(summarize_beta_grid(obj$profile_grid)$K)),
    LAMBDA_GRID,"BIC",SELECTION_BOOTSTRAP_N_START,N_CORES,MAX_ITER,FIT_TOL,PHI_LOWER,PHI_UPPER,
    OPTIMIZER_MAXIT,PRECISION_SATURATION_FRACTION,SEED+7000L))
conditional_bootstrap <- load_or_run_final(file.path(CHECKPOINT_DIR,"conditional_parametric_bootstrap.rds"),
  paste0("Conditional parametric bootstrap B=",BOOTSTRAP_B),
  run_beta_parametric_bootstrap(sel$fit,BOOTSTRAP_B,BOOTSTRAP_N_START,N_CORES,MAX_ITER,FIT_TOL,
    PHI_LOWER,PHI_UPPER,OPTIMIZER_MAXIT,PRECISION_SATURATION_FRACTION,BOOTSTRAP_LEVEL,SEED+8000L,
    interval_methods=BOOTSTRAP_INTERVAL_METHODS,
    retry_failed=BOOTSTRAP_RETRY_FAILED,retry_n_start=BOOTSTRAP_RETRY_N_START))
# Save raw bootstrap object so alternative curve interval formulas can be produced
# from the same refits without repeating the expensive bootstrap.
main_tables <- safe_dir_create(file.path(FINAL_DIR,"paper","main","tables"))
main_figs <- safe_dir_create(file.path(FINAL_DIR,"paper","main","figures"))
supp_tables <- safe_dir_create(file.path(FINAL_DIR,"paper","supp","tables"))
supp_figs <- safe_dir_create(file.path(FINAL_DIR,"paper","supp","figures"))
utils::write.csv(sel$best_by_K[,c("K","lambda","AIC_EDF","BIC_EDF","ICL_EDF","min_state_share","precision_flag")],file.path(main_tables,"Table_application_model_selection.csv"),row.names=FALSE)
utils::write.csv(selection_stability$summary,file.path(supp_tables,"selection_stability_summary.csv"),row.names=FALSE)
utils::write.csv(selection_stability$details,file.path(supp_tables,"selection_stability_replicates.csv"),row.names=FALSE)
utils::write.csv(selection_stability$K_frequency,file.path(supp_tables,"selection_stability_K_frequency.csv"),row.names=FALSE)
utils::write.csv(data.frame(
  B_requested=conditional_bootstrap$B_requested,
  B_success_first_pass=conditional_bootstrap$B_success_first_pass,
  B_retry_attempted=conditional_bootstrap$B_retry_attempted,
  B_retry_recovered=conditional_bootstrap$B_retry_recovered,
  B_success=conditional_bootstrap$B_success,
  success_rate=conditional_bootstrap$success_rate,
  retry_n_start=conditional_bootstrap$retry_n_start
),file.path(supp_tables,"conditional_bootstrap_success.csv"),row.names=FALSE)
utils::write.csv(data.frame(
  replicate=seq_len(conditional_bootstrap$B_requested),
  first_pass_success=conditional_bootstrap$success_first_pass,
  retried=seq_len(conditional_bootstrap$B_requested) %in% conditional_bootstrap$retry_indices,
  final_success=conditional_bootstrap$success_final
),file.path(supp_tables,"conditional_bootstrap_replicate_status.csv"),row.names=FALSE)
utils::write.csv(data.frame(parameter=names(conditional_bootstrap$interval_methods),
  method=unlist(conditional_bootstrap$interval_methods,use.names=FALSE)),
  file.path(supp_tables,"conditional_bootstrap_interval_methods.csv"),row.names=FALSE)
# Paper-style outputs use the pre-specified interval convention. Alternative
# formulas are retained from the same refits, so sensitivity needs no rerun.
tpar <- make_table4(sel$fit,conditional_bootstrap)
utils::write.csv(tpar,file.path(main_tables,"Table_parameter_estimates_CI.csv"),row.names=FALSE)
plot_state_profiles(sel$fit,conditional_bootstrap,panel,file.path(main_figs,"Figure_state_age_profiles.pdf"),file.path(main_figs,"Figure_state_age_profiles.png"))
plot_viterbi_years(sel$fit,panel,file.path(main_figs,"Figure_viterbi_by_year.pdf"),file.path(main_figs,"Figure_viterbi_by_year.png"))
saveRDS(list(selection=sel,selection_stability=selection_stability,conditional_bootstrap=conditional_bootstrap),file.path(FINAL_DIR,"final_inference_object.rds"))
writeLines(capture.output(sessionInfo()),file.path(FINAL_DIR,"sessionInfo.txt"))
message("Final inference complete. Raw bootstrap draws saved for interval recalibration without refitting.")
message("Outputs: ",normalizePath(FINAL_DIR,winslash="/",mustWork=FALSE))
