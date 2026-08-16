# v0.0.0.9006 simulation-validity check.
# This does not tune the method and does not replace the definitive study.
# It checks the corrected DGP, fixed degree-2 response basis, and independent
# test prediction path before launching the 200-replicate factorial experiment.

library(spfcICOMP)

dir.create("thesis_analysis/outputs_v9006", recursive = TRUE, showWarnings = FALSE)
dir.create("thesis_analysis/tables_v9006", recursive = TRUE, showWarnings = FALSE)

cat("spfcICOMP version: ", as.character(packageVersion("spfcICOMP")), "\n", sep = "")

# Direct algebraic checks on one rank-two population draw.
sim2 <- simulate_spfc_continuous(
  n = 200, n_test = 200, p = 20, d = 2, s = 5,
  rho_x = 0.2, snr = 3, seed = 20260809
)

cat("Rank of true inverse mean: ", qr(sim2$inverse_mean)$rank, "\n", sep = "")
cat("Requested SNR: 3; constructed SNR: ", signif(sim2$snr_check, 8), "\n", sep = "")
cat("True reduction-subspace check: ",
    signif(subspace_distance(sim2$B_true, solve(sim2$Sigma, sim2$Gamma_true)), 8),
    "\n", sep = "")
cat("Training/test matrices identical? ", identical(sim2$X, sim2$X_test), "\n", sep = "")

stopifnot(qr(sim2$inverse_mean)$rank == 2L)
stopifnot(abs(sim2$snr_check - 3) < 1e-10)
stopifnot(subspace_distance(sim2$B_true, solve(sim2$Sigma, sim2$Gamma_true)) < 1e-7)
stopifnot(!identical(sim2$X, sim2$X_test))

# Small descriptive validation grid. No performance threshold is used to alter
# or accept the method: the purpose is to expose gross implementation failures.
validity_design <- create_spfc_design_grid(
  response_type = "continuous",
  n_values = 100,
  p_values = 20,
  d_values = c(1, 2),
  s_values = 5,
  rho_x_values = 0.2,
  snr_values = 3,
  cov_methods = c("oas", "mec"),
  variable_methods = "adaptive_weighted_l1",
  nrep = 10
)

validity_results <- run_spfc_design(
  design = validity_design,
  base_seed = 20260808L,
  d_grid = 1:5,
  criteria = c("AIC", "BIC", "CAIC", "ICOMP_IFIM", "CICOMP"),
  selection_rule = "c1f",
  nslices = 5,
  poly_degree = 2,
  verbose = TRUE
)

stopifnot(all(validity_results$prediction_sample == "independent_test"))
stopifnot(all(validity_results$n_prediction == validity_results$n))

saveRDS(validity_results, "thesis_analysis/outputs_v9006/validity_results.rds")
write.csv(validity_results,
          "thesis_analysis/outputs_v9006/validity_results.csv",
          row.names = FALSE)

validity_summary <- aggregate(
  cbind(dimension_correct, rmse, mae, precision, recall, f1_variable) ~
    true_d + cov_method + criterion,
  data = validity_results,
  FUN = function(x) mean(x, na.rm = TRUE)
)

write.csv(validity_summary,
          "thesis_analysis/tables_v9006/validity_summary.csv",
          row.names = FALSE)
print(validity_summary)

cat("\nValidity check completed. No tuning decision is made from these values.\n")
