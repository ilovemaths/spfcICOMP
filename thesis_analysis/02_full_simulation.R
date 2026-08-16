# Full thesis simulation: continuous-response SPFC-ICOMP.
# Base factorial combinations:
# n={50,100}; p={20,100,500}; d0={1,2}; rho_x={0.2,0.8};
# SNR={1,3}; covariance={OAS,MEC}; 200 Monte Carlo replicates.
# This is 96 base configurations and 19,200 design rows.
# Each design row scores d=1,...,5 once and reports the selected dimension
# independently for AIC, BIC, CAIC, ICOMP(IFIM), and CICOMP.

library(spfcICOMP)

DESIGN_SEED <- 20260808L
NREP <- 200L
ACTIVE_S <- 5L

full_design <- create_spfc_design_grid(
  response_type = "continuous",
  n_values = c(50, 100),
  p_values = c(20, 100, 500),
  d_values = c(1, 2),
  s_values = ACTIVE_S,
  rho_x_values = c(0.2, 0.8),
  snr_values = c(1, 3),
  cov_methods = c("oas", "mec"),
  variable_methods = "adaptive_weighted_l1",
  nrep = NREP
)

write.csv(full_design, "thesis_analysis/outputs_v9006/full_design.csv", row.names = FALSE)

# Run in p-blocks to make restart/checkpoint management safer on Windows.
p_blocks <- split(full_design, full_design$p)
block_results <- vector("list", length(p_blocks))
names(block_results) <- names(p_blocks)

for (nm in names(p_blocks)) {
  message("Starting p = ", nm)
  block_results[[nm]] <- run_spfc_design(
    design = p_blocks[[nm]],
    base_seed = DESIGN_SEED,
    d_grid = 1:5,
    criteria = c("AIC", "BIC", "CAIC", "ICOMP_IFIM", "CICOMP"),
    selection_rule = "c1f",
    nslices = 5,
    poly_degree = 2,
    save_every = 100,
    save_path = paste0("thesis_analysis/outputs_v9006/partial_p", nm, ".rds"),
    verbose = TRUE
  )
  saveRDS(block_results[[nm]], paste0("thesis_analysis/outputs_v9006/results_p", nm, ".rds"))
}

results <- do.call(rbind, block_results)
rownames(results) <- NULL
saveRDS(results, "thesis_analysis/outputs_v9006/full_simulation_results.rds")
write.csv(results, "thesis_analysis/outputs_v9006/full_simulation_results.csv", row.names = FALSE)

summary_results <- summarise_spfc_simulation(results)
write.csv(summary_results, "thesis_analysis/tables_v9006/full_simulation_summary.csv", row.names = FALSE)
