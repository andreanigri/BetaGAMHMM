# Redraw the real-data state-profile figures (manuscript Figs 3 and 5, Supplement Fig. S11)
# from the SAVED final-inference objects. No model is refitted and no bootstrap is rerun.
#
# Changes with respect to plot_state_profiles() in real_data/R/07_reporting.R:
#   * the bootstrap band is drawn BEFORE the observations, so points are no longer hidden;
#   * observed annual profiles are clearly visible (larger, darker points);
#   * larger axis labels / tick labels (reviewer request + colleague comment on Fig. 5);
#   * y-axis label "Female share of deaths" (consistent with "proportion" in the text).
#
# Run from the root of BetaGAMHMM-v5.2.0 (the folder containing real_data/).

plot_state_profiles_v2 <- function(fit, boot, panel, pdf_path, png_path) {
  ages <- panel$ages
  K <- fit$K
  ncol_panel <- if (K <= 2L) K else 2L
  nrow_panel <- ceiling(K / ncol_panel)
  # Common vertical scale within a country/window (as in the submitted figures)
  ylim <- range(c(panel$Y_raw, boot$mu_low, boot$mu_high), finite = TRUE)
  ylim <- c(max(0, ylim[1L]), min(1, ylim[2L]))
  draw <- function() {
    old <- graphics::par(mfrow = c(nrow_panel, ncol_panel),
                         mar = c(4.4, 4.8, 2.4, 1.0),
                         cex.lab = 1.35, cex.axis = 1.15, cex.main = 1.35)
    on.exit(graphics::par(old), add = TRUE)
    for (k in seq_len(K)) {
      graphics::plot(ages, fit$mu_age[, k], type = "n", ylim = ylim,
                     xlab = "Age", ylab = "Female share of deaths",
                     main = paste("State", k))
      # 1) pointwise 95% percentile bootstrap band
      graphics::polygon(c(ages, rev(ages)), c(boot$mu_low[, k], rev(boot$mu_high[, k])),
                        border = NA, col = grDevices::adjustcolor("grey70", alpha.f = 0.45))
      # 2) observed annual profiles of the years decoded to state k
      for (t in which(fit$viterbi == k)) {
        graphics::points(ages, panel$Y_raw[t, ], pch = 16, cex = 0.45,
                         col = grDevices::adjustcolor("grey25", alpha.f = 0.55))
      }
      # 3) fitted state-specific mean curve
      graphics::lines(ages, fit$mu_age[, k], lwd = 2.4)
    }
  }
  grDevices::pdf(pdf_path, width = 8.5, height = max(4.5, 3.7 * nrow_panel)); draw(); grDevices::dev.off()
  grDevices::png(png_path, width = 1600, height = max(900, 700 * nrow_panel), res = 180); draw(); grDevices::dev.off()
  invisible(TRUE)
}

cases <- c(
  ITA55 = "real_data/outputs/p_ITA_1960_2014_a0_40",
  RUS55 = "real_data/outputs/p_RUS_1960_2014_a0_40",
  ITA70 = "real_data/outputs/p_ITA_1950_2019_a0_40"
)
out_dir <- "redrawn_figures"
dir.create(out_dir, showWarnings = FALSE)

# Locate the saved final-inference object written by RUN_FINAL_INFERENCE_REAL_DATA.R
# (folder "final_inference"; any other "*final_inference" folder is accepted as a fallback).
find_inference_object <- function(case_dir) {
  candidates <- c(file.path(case_dir, "final_inference"),
                  list.dirs(case_dir, recursive = FALSE)[grepl("final_inference$", list.dirs(case_dir, recursive = FALSE))])
  files <- file.path(unique(candidates), "final_inference_object.rds")
  files <- files[file.exists(files)]
  if (!length(files)) stop("No final_inference_object.rds found under ", case_dir,
                           ". Run RUN_FINAL_INFERENCE_REAL_DATA.R first.")
  files[1L]
}

for (nm in names(cases)) {
  obj <- readRDS(find_inference_object(cases[[nm]]))
  fit <- obj$selection$fit
  boot <- obj$conditional_bootstrap
  panel <- fit$panel
  plot_state_profiles_v2(
    fit, boot, panel,
    pdf_path = file.path(out_dir, paste0("Figure_state_age_profiles_", nm, ".pdf")),
    png_path = file.path(out_dir, paste0("Figure_state_age_profiles_", nm, ".png"))
  )
  message("Redrawn: ", nm)
}
