# 07_real_data_audit_and_aggregate.R
# Audits and aggregates the definitive v0.0.0.9008 empirical analyses.
# No SPFC models are refitted.

EXPECTED_VERSION <- "0.0.0.9008"

output_dir <- "thesis_analysis/real_outputs_v9008"
table_dir  <- "thesis_analysis/real_tables_v9008"

required <- c(
  file.path(table_dir, "riboflavin_data_audit.csv"),
  file.path(table_dir, "riboflavin_full_data_summary.csv"),
  file.path(table_dir, "riboflavin_cv_summary.csv"),
  file.path(table_dir, "riboflavin_selected_variables.csv"),
  file.path(table_dir, "golub_data_audit.csv"),
  file.path(table_dir, "golub_full_data_summary.csv"),
  file.path(table_dir, "golub_cv_summary.csv"),
  file.path(table_dir, "golub_selected_variables.csv"),
  file.path(output_dir, "riboflavin_cv_predictions.csv"),
  file.path(output_dir, "golub_cv_predictions.csv")
)

if (!all(file.exists(required))) {
  missing <- required[!file.exists(required)]
  stop("Missing required empirical output(s):\n", paste(missing, collapse = "\n"))
}

rib_audit <- read.csv(file.path(table_dir, "riboflavin_data_audit.csv"), check.names = FALSE)
rib_full  <- read.csv(file.path(table_dir, "riboflavin_full_data_summary.csv"), check.names = FALSE)
rib_cv    <- read.csv(file.path(table_dir, "riboflavin_cv_summary.csv"), check.names = FALSE)
rib_vars  <- read.csv(file.path(table_dir, "riboflavin_selected_variables.csv"), check.names = FALSE)

gol_audit <- read.csv(file.path(table_dir, "golub_data_audit.csv"), check.names = FALSE)
gol_full  <- read.csv(file.path(table_dir, "golub_full_data_summary.csv"), check.names = FALSE)
gol_cv    <- read.csv(file.path(table_dir, "golub_cv_summary.csv"), check.names = FALSE)
gol_vars  <- read.csv(file.path(table_dir, "golub_selected_variables.csv"), check.names = FALSE)

rib_pred <- read.csv(file.path(output_dir, "riboflavin_cv_predictions.csv"), check.names = FALSE)
gol_pred <- read.csv(file.path(output_dir, "golub_cv_predictions.csv"), check.names = FALSE)

criteria <- c("AIC","BIC","CAIC","ICOMP_IFIM","CICOMP")
cov_methods <- c("oas","mec")

# ------------------------------------------------------------
# Structural checks
# ------------------------------------------------------------

stopifnot(
  nrow(rib_audit) == 1L,
  nrow(gol_audit) == 1L,
  setequal(unique(rib_full$criterion), criteria),
  setequal(unique(rib_full$cov_method), cov_methods),
  nrow(rib_full) == length(criteria) * length(cov_methods),
  setequal(unique(gol_full$cov_method), cov_methods),
  nrow(gol_full) == length(cov_methods),
  all(gol_full$fixed_d == 1L),
  all(rib_full$selected_d %in% 1:5),
  all(rib_full$n_selected >= 0),
  all(gol_full$n_selected >= 0),
  all(is.finite(rib_cv$mean_rmse)),
  all(is.finite(rib_cv$mean_mae)),
  all(gol_cv$mean_accuracy >= 0 & gol_cv$mean_accuracy <= 1),
  all(gol_cv$mean_balanced_accuracy >= 0 & gol_cv$mean_balanced_accuracy <= 1),
  all(gol_cv$mean_auc >= 0 & gol_cv$mean_auc <= 1)
)

# Expected prediction count:
# Riboflavin: n * repeats * 2 covariance * 5 criteria
expected_rib_predictions <- rib_audit$n[1] * 5L * 2L * 5L
stopifnot(nrow(rib_pred) == expected_rib_predictions)

# Golub: n * repeats * 2 covariance
expected_gol_predictions <- gol_audit$n[1] * 20L * 2L
stopifnot(nrow(gol_pred) == expected_gol_predictions)

# ------------------------------------------------------------
# Feature-selection overlap between OAS and MEC
# ------------------------------------------------------------

overlap_table <- function(vars, dataset, has_criterion = TRUE) {
  if (!nrow(vars)) return(data.frame())

  if (has_criterion) {
    out <- lapply(criteria, function(cr) {
      a <- unique(vars$variable_name[vars$cov_method == "oas" & vars$criterion == cr])
      b <- unique(vars$variable_name[vars$cov_method == "mec" & vars$criterion == cr])
      inter <- intersect(a, b)
      union_set <- union(a, b)
      data.frame(
        dataset = dataset,
        criterion = cr,
        oas_selected = length(a),
        mec_selected = length(b),
        common_selected = length(inter),
        jaccard = if (length(union_set)) length(inter) / length(union_set) else NA_real_
      )
    })
    do.call(rbind, out)
  } else {
    a <- unique(vars$variable_name[vars$cov_method == "oas"])
    b <- unique(vars$variable_name[vars$cov_method == "mec"])
    inter <- intersect(a, b)
    union_set <- union(a, b)
    data.frame(
      dataset = dataset,
      criterion = "d=1",
      oas_selected = length(a),
      mec_selected = length(b),
      common_selected = length(inter),
      jaccard = if (length(union_set)) length(inter) / length(union_set) else NA_real_
    )
  }
}

selection_overlap <- rbind(
  overlap_table(rib_vars, "Riboflavin", TRUE),
  overlap_table(gol_vars, "Golub leukaemia", FALSE)
)

write.csv(
  selection_overlap,
  file.path(table_dir, "real_data_selection_overlap.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Compact thesis tables
# ------------------------------------------------------------

rib_table <- merge(
  rib_full,
  rib_cv,
  by = c("cov_method", "criterion"),
  all.x = TRUE,
  sort = FALSE
)

rib_table <- rib_table[
  order(rib_table$cov_method, match(rib_table$criterion, criteria)),
]

write.csv(
  rib_table,
  file.path(table_dir, "table_riboflavin_definitive.csv"),
  row.names = FALSE
)

gol_table <- merge(
  gol_full,
  gol_cv,
  by = "cov_method",
  all.x = TRUE,
  sort = FALSE
)

write.csv(
  gol_table,
  file.path(table_dir, "table_golub_definitive.csv"),
  row.names = FALSE
)

dataset_summary <- data.frame(
  dataset = c("Riboflavin", "Golub leukaemia"),
  response = c("Continuous", "Binary"),
  n = c(rib_audit$n[1], gol_audit$n[1]),
  p = c(rib_audit$p[1], gol_audit$p[1]),
  p_over_n = c(rib_audit$p_over_n[1], gol_audit$p_over_n[1]),
  covariance_methods = "OAS, MEC",
  dimension_strategy = c(
    "AIC/BIC/CAIC/ICOMP_IFIM/CICOMP over d=1,...,5",
    "d=1 fixed (validated binary pathway)"
  ),
  validation = c(
    "Repeated 5-fold CV, 5 repeats",
    "Repeated stratified 5-fold CV, 20 repeats"
  )
)

write.csv(
  dataset_summary,
  file.path(table_dir, "table_real_data_design_summary.csv"),
  row.names = FALSE
)

manifest <- data.frame(
  item = c(
    "package_version",
    "riboflavin_prediction_protocol",
    "golub_prediction_protocol",
    "riboflavin_dimension_grid",
    "golub_dimension",
    "feature_selection",
    "covariance_methods",
    "riboflavin_prediction_rows",
    "golub_prediction_rows"
  ),
  value = c(
    EXPECTED_VERSION,
    "5-fold CV x 5 repeats; dimension selected inside training folds",
    "stratified 5-fold CV x 20 repeats; ridge logistic",
    "1,2,3,4,5",
    "1 (binary inverse-basis rank)",
    "C1F ICOMP-HD adaptive weighted L1",
    "OAS, MEC",
    nrow(rib_pred),
    nrow(gol_pred)
  )
)

write.csv(
  manifest,
  file.path(table_dir, "real_data_analysis_manifest.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    riboflavin = list(audit = rib_audit, full = rib_full, cv = rib_cv),
    golub = list(audit = gol_audit, full = gol_full, cv = gol_cv),
    selection_overlap = selection_overlap,
    manifest = manifest
  ),
  file.path(output_dir, "real_data_postprocessed.rds")
)

cat("\n============================================================\n")
cat("REAL-DATA AUDIT COMPLETE\n")
cat("============================================================\n")
print(dataset_summary)
cat("Riboflavin prediction rows:", nrow(rib_pred), "\n")
cat("Golub prediction rows:", nrow(gol_pred), "\n")
cat("Outputs written to", table_dir, "and", output_dir, "\n")
