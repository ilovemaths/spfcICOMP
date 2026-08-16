# Small integrity run before the full factorial experiment.
library(spfcICOMP)

dir.create("thesis_analysis/outputs_v9006", recursive = TRUE, showWarnings = FALSE)
dir.create("thesis_analysis/tables_v9006", recursive = TRUE, showWarnings = FALSE)
dir.create("thesis_analysis/figures_v9006", recursive = TRUE, showWarnings = FALSE)

set.seed(20260808)

design <- create_spfc_design_grid(
  response_type = "continuous",
  n_values = 50,
  p_values = c(20, 100),
  d_values = c(1, 2),
  s_values = 5,
  rho_x_values = c(0.2, 0.8),
  snr_values = c(1, 3),
  cov_methods = c("oas", "mec"),
  variable_methods = "adaptive_weighted_l1",
  nrep = 2
)

res <- run_spfc_design(
  design = design,
  d_grid = 1:5,
  criteria = c("AIC", "BIC", "CAIC", "ICOMP_IFIM", "CICOMP"),
  selection_rule = "c1f",
  nslices = 5,
  poly_degree = 2,
  save_every = 10,
  save_path = "thesis_analysis/outputs_v9006/smoke_partial.rds",
  verbose = TRUE
)

saveRDS(res, "thesis_analysis/outputs_v9006/smoke_results.rds")
write.csv(res, "thesis_analysis/outputs_v9006/smoke_results.csv", row.names = FALSE)

sm <- summarise_spfc_simulation(res)
write.csv(sm, "thesis_analysis/tables_v9006/smoke_summary.csv", row.names = FALSE)
# Print a focused diagnostic table instead of flooding the console with every
# available summary column. Four smoke replications per summary cell would
# indicate unintended pooling across SNR levels; the expected value is 2.
diag_cols <- intersect(
  c("cov_method", "criterion", "n", "p", "true_d", "rho_x", "snr", "nrep",
    "dimension_recovery_rate", "mean_selected_d", "mean_precision",
    "mean_recall", "mean_f1_variable", "mean_n_selected",
    "mean_c1f", "mean_c1f_hd_factor", "mean_c1f_complexity_multiplier",
    "mean_c1f_global_penalty", "mean_rmse"),
  names(sm)
)
print(sm[, diag_cols, drop = FALSE])

cat("\nSmoke-test structural checks:\n")
cat("  Expected nrep per summary cell = 2; observed range = ",
    paste(range(sm$nrep, na.rm = TRUE), collapse = " to "), "\n", sep = "")
cat("  Mean selected variables by p:\n")
print(aggregate(mean_n_selected ~ p, data = sm, FUN = mean))
cat("  Mean precision/recall/F1 by p:\n")
print(aggregate(cbind(mean_precision, mean_recall, mean_f1_variable) ~ p,
                data = sm, FUN = mean))

cat("\nICOMP-HD smoke diagnostics (descriptive only; no outcome-based tuning):\n")
cat("  The calibration is frozen before the definitive Monte Carlo study.\n")
cat("  These summaries are used only to detect coding/pathological scaling errors.\n")

cat("\nv0.0.0.9006 validity checks:\n")
stopifnot(all(res$prediction_sample == "independent_test"))
stopifnot(all(res$n_prediction == res$n))
cat("  Prediction sample: independent test data only.\n")
cat("  Test-size range: ", paste(range(res$n_prediction), collapse = " to "), "\n", sep = "")
