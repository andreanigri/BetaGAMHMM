# Reporting for the final single-BIC real-data pipeline.

write_simple_latex_table <- function(df, path, digits = 4L) {
  x <- df
  for (j in seq_along(x)) if (is.numeric(x[[j]])) x[[j]] <- formatC(x[[j]], digits = digits, format = "fg", flag = "#")
  esc <- function(z) gsub("_", "\\\\_", as.character(z), fixed = TRUE)
  hdr <- paste(esc(names(x)), collapse = " & ")
  rows <- apply(x, 1L, function(r) paste(esc(r), collapse = " & "))
  lines <- c(paste0("\\\\begin{tabular}{", paste(rep("l", ncol(x)), collapse = ""), "}"),
             "\\\\toprule", paste0(hdr, " \\\\\\"), "\\\\midrule",
             paste0(rows, " \\\\\\"), "\\\\bottomrule", "\\\\end{tabular}")
  writeLines(lines, path); invisible(path)
}

safe_file_tag <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

plot_profile_BIC_path <- function(best_by_k, output_dir) {
  safe_dir_create(output_dir)
  draw <- function() {
    d <- best_by_k[order(best_by_k$K), , drop=FALSE]
    y <- d$BIC_EDF - min(d$BIC_EDF, na.rm=TRUE)
    plot(d$K, y, type="b", pch=16, xlab="Number of states K", ylab="Delta BIC", xaxt="n")
    axis(1, at=d$K); abline(h=0, lty=3)
  }
  grDevices::pdf(file.path(output_dir,"profile_delta_BIC_by_K.pdf"),6.7,5.2); draw(); grDevices::dev.off()
  grDevices::png(file.path(output_dir,"profile_delta_BIC_by_K.png"),1250,950,res=170); draw(); grDevices::dev.off()
}

plot_profile_curves_selected <- function(fit, panel, output_dir) {
  safe_dir_create(output_dir)

  draw <- function() {
    K <- fit$K
    nr <- ceiling(sqrt(K))
    nc <- ceiling(K / nr)
    op <- par(mfrow = c(nr, nc), mar = c(3.8, 4.2, 2.2, 1.0), oma = c(0, 0, 2.2, 0))
    on.exit(par(op), add = TRUE)

    for (k in seq_len(K)) {
      idx <- which(fit$viterbi == k)
      vals <- c(fit$mu_age[, k], if (length(idx)) as.numeric(panel$Y[idx, , drop = FALSE]) else numeric(0))
      ylim <- range(vals, finite = TRUE)
      if (!all(is.finite(ylim)) || diff(ylim) <= 0) ylim <- c(0, 1)

      plot(panel$ages, fit$mu_age[, k], type = "n",
           xlab = "Age", ylab = "Female share of deaths",
           ylim = ylim, main = paste0("State ", k, " (n=", length(idx), ")"))
      if (length(idx)) {
        yy <- t(panel$Y[idx, , drop = FALSE])
        if (is.vector(yy)) yy <- matrix(yy, ncol = 1L)
        matlines(panel$ages, yy, lty = 1, lwd = 0.8,
                 col = grDevices::adjustcolor("grey40", alpha.f = 0.28))
      }
      lines(panel$ages, fit$mu_age[, k], lwd = 2.4, col = "black")
    }
    mtext(paste0("Observed annual age profiles by Viterbi state; fitted Beta-GAM curves (K=", K, ")"),
          outer = TRUE, cex = 1.05)
  }

  grDevices::pdf(file.path(output_dir, "selected_profile_curves.pdf"), 10, 7.5)
  draw()
  grDevices::dev.off()

  grDevices::png(file.path(output_dir, "selected_profile_curves.png"), 1800, 1350, res = 170)
  draw()
  grDevices::dev.off()
}

plot_scalar_aggregate <- function(scalar_selection, output_dir) {
  safe_dir_create(output_dir); f <- scalar_selection$fit; s <- f$series
  draw <- function() plot(s$years,s$y,type="o",pch=16,xlab="Year",ylab="Aggregated female share",
                          main=paste0("Scalar annual Beta-HMM (K=",f$K,")"))
  grDevices::pdf(file.path(output_dir,"scalar_annual_series.pdf"),7.2,4.8); draw(); grDevices::dev.off()
  grDevices::png(file.path(output_dir,"scalar_annual_series.png"),1350,900,res=170); draw(); grDevices::dev.off()
}

plot_state_sequence_comparison_single <- function(panel, scalar_selection, profile_selection, output_dir) {
  safe_dir_create(output_dir)
  states <- rbind(scalar_selection$fit$viterbi, profile_selection$fit$viterbi)
  draw <- function() {
    op <- par(mar=c(4.5,7.5,2,1)); on.exit(par(op))
    image(x=panel$years,y=1:2,z=t(states),axes=FALSE,xlab="Year",ylab="",main="Decoded annual state sequences")
    axis(1); axis(2,at=1:2,labels=c("Scalar aggregate","Profile"),las=1); box()
  }
  grDevices::pdf(file.path(output_dir,"scalar_vs_profile_state_sequences.pdf"),9,4.2); draw(); grDevices::dev.off()
  grDevices::png(file.path(output_dir,"scalar_vs_profile_state_sequences.png"),1600,800,res=170); draw(); grDevices::dev.off()
}

write_final_selection_outputs <- function(panel, scalar_selection, profile_grid, profile_selection,
                                          basis_sensitivity, beta_diagnostics, beta_persistence,
                                          binomial_grid=NULL, binomial_selection=NULL,
                                          binomial_diagnostics=NULL, phi_bound_sensitivity=NULL, output_root) {
  main_tables <- safe_dir_create(file.path(output_root,"main_text","tables"))
  main_figs <- safe_dir_create(file.path(output_root,"main_text","figures"))
  supp_tables <- safe_dir_create(file.path(output_root,"supplement","tables"))
  supp_figs <- safe_dir_create(file.path(output_root,"supplement","figures"))

  tsel <- profile_selection$best_by_K[,c("K","lambda","AIC_EDF","BIC_EDF","ICL_EDF","min_state_share","precision_flag"),drop=FALSE]
  names(tsel)[2] <- "lambda_star"
  utils::write.csv(tsel,file.path(main_tables,"Table_application_model_selection.csv"),row.names=FALSE)
  write_simple_latex_table(tsel,file.path(main_tables,"Table_application_model_selection.tex"))
  utils::write.csv(data.frame(selected_K=profile_selection$fit$K,selected_lambda=profile_selection$fit$lambda[1L],
    BIC_EDF=profile_selection$selected_row$BIC_EDF,loglik=profile_selection$fit$loglik_obs,
    min_state_share=profile_selection$fit$occupancy$min_share,
    precision_flag=profile_selection$fit$precision_diagnostic$flagged),
    file.path(main_tables,"selected_profile_model.csv"),row.names=FALSE)
  utils::write.csv(summarize_beta_grid(profile_grid),file.path(supp_tables,"profile_full_grid.csv"),row.names=FALSE)
  utils::write.csv(scalar_selection$grid_summary,file.path(supp_tables,"scalar_annual_model_selection.csv"),row.names=FALSE)
  utils::write.csv(compute_scalar_beta_pseudo_residuals(scalar_selection$fit),file.path(supp_tables,"scalar_pseudo_residual_summary.csv"),row.names=FALSE)
  utils::write.csv(beta_diagnostics$summary,file.path(supp_tables,"beta_pseudo_residual_summary.csv"),row.names=FALSE)
  save_beta_diagnostic_plots(beta_diagnostics,panel$years,supp_figs)
  utils::write.csv(beta_persistence$transition,file.path(supp_tables,"state_persistence_summary.csv"),row.names=FALSE)
  utils::write.csv(beta_persistence$decoded_runs,file.path(supp_tables,"decoded_state_runs.csv"),row.names=FALSE)
  utils::write.csv(data.frame(year=panel$years,scalar_state=scalar_selection$fit$viterbi,
                              profile_state=profile_selection$fit$viterbi),
                   file.path(main_tables,"annual_decoded_states_scalar_profile.csv"),row.names=FALSE)
  if (!is.null(basis_sensitivity)) utils::write.csv(basis_sensitivity$summary,file.path(supp_tables,"basis_dimension_sensitivity.csv"),row.names=FALSE)
  if (!is.null(phi_bound_sensitivity)) {
    utils::write.csv(phi_bound_sensitivity$summary,file.path(supp_tables,"phi_bound_sensitivity_summary.csv"),row.names=FALSE)
    utils::write.csv(phi_bound_sensitivity$grid,file.path(supp_tables,"phi_bound_sensitivity_grid.csv"),row.names=FALSE)
    writeLines(phi_bound_sensitivity$note,file.path(supp_tables,"phi_bound_sensitivity_note.txt"))
  }
  if (!is.null(binomial_selection)) {
    utils::write.csv(summarize_binomial_grid(binomial_grid),file.path(supp_tables,"denominator_aware_full_grid.csv"),row.names=FALSE)
    comp <- compare_beta_binomial(profile_selection,binomial_selection)
    utils::write.csv(comp,file.path(supp_tables,"beta_vs_denominator_aware_summary.csv"),row.names=FALSE)
    if (!is.null(binomial_diagnostics)) {
      utils::write.csv(binomial_diagnostics$summary,file.path(supp_tables,"denominator_aware_residual_summary.csv"),row.names=FALSE)
      save_binomial_kernel_diagnostic_plots(binomial_diagnostics,supp_figs)
    }
  }
  plot_profile_BIC_path(profile_selection$best_by_K,main_figs)
  plot_profile_curves_selected(profile_selection$fit,panel,main_figs)
  plot_scalar_aggregate(scalar_selection,main_figs)
  plot_state_sequence_comparison_single(panel,scalar_selection,profile_selection,main_figs)
  writeLines(c(
    "FINAL STATE-SELECTION RULE",
    "For each K, lambda is chosen by minimum EDF-adjusted AIC.",
    "K is chosen by minimum EDF-adjusted BIC = -2 log L + nu log(T).",
    "T is the number of HMM time points (calendar years in this application).",
    "No +2 BIC tolerance rule, occupancy threshold, or multi-BIC voting rule is used.",
    "Precision-bound and occupancy diagnostics are reported but do not determine K.",
    "The scalar benchmark is a valid annual aggregate series; the obsolete 2,255-cell flattening is not used.",
    "The denominator-aware analysis is a same-data sensitivity analysis; its BIC values are not compared numerically with Beta BIC values.",
    "Precision-bound sensitivity refits the primary AIC-selected lambda within each K under relaxed phi bounds using the corresponding primary fit as a warm-start/feasible candidate; it is a targeted numerical-safeguard check.",
    "The selected-profile figure shows only the observed annual age profiles assigned to each Viterbi state, overlaid with that state's fitted curve."
  ),file.path(output_root,"ANALYSIS_SCOPE.txt"))
  invisible(output_root)
}
