# 08_real_data_figures.R
# Publication/thesis figures for definitive v0.0.0.9008 empirical analyses.
# No SPFC models are refitted.

stopifnot(requireNamespace("ggplot2", quietly = TRUE))
library(ggplot2)

table_dir  <- "thesis_analysis/real_tables_v9008"
figure_dir <- "thesis_analysis/real_figures_v9008"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

required <- c(
  file.path(table_dir, "table_riboflavin_definitive.csv"),
  file.path(table_dir, "table_golub_definitive.csv"),
  file.path(table_dir, "riboflavin_dimension_scores_oas.csv"),
  file.path(table_dir, "riboflavin_dimension_scores_mec.csv")
)
stopifnot(all(file.exists(required)))

rib <- read.csv(file.path(table_dir, "table_riboflavin_definitive.csv"), check.names = FALSE)
gol <- read.csv(file.path(table_dir, "table_golub_definitive.csv"), check.names = FALSE)
scores_oas <- read.csv(file.path(table_dir, "riboflavin_dimension_scores_oas.csv"), check.names = FALSE)
scores_mec <- read.csv(file.path(table_dir, "riboflavin_dimension_scores_mec.csv"), check.names = FALSE)

criterion_levels <- c("AIC","BIC","CAIC","ICOMP_IFIM","CICOMP")
rib$criterion <- factor(rib$criterion, levels = criterion_levels)
rib$cov_method <- toupper(rib$cov_method)
gol$cov_method <- toupper(gol$cov_method)

save_plot <- function(p, nm, w = 8, h = 5) {
  ggsave(
    filename = file.path(figure_dir, paste0(nm, ".pdf")),
    plot = p, width = w, height = h, units = "in"
  )
  ggsave(
    filename = file.path(figure_dir, paste0(nm, ".png")),
    plot = p, width = w, height = h, units = "in", dpi = 300
  )
}

# 1. Selected structural dimension, Riboflavin
p1 <- ggplot(
  rib,
  aes(
    x = criterion,
    y = selected_d,
    shape = cov_method,
    group = cov_method
  )
) +
  geom_point(size = 3) +
  geom_line() +
  scale_y_continuous(breaks = 1:5, limits = c(1,5)) +
  labs(
    x = "Information criterion",
    y = "Selected structural dimension",
    shape = "Covariance",
    title = "Riboflavin: selected structural dimension"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(p1, "riboflavin_selected_dimension")

# 2. Repeated-CV RMSE with ±1 SD, Riboflavin
p2 <- ggplot(
  rib,
  aes(
    x = criterion,
    y = mean_rmse,
    shape = cov_method,
    group = cov_method
  )
) +
  geom_errorbar(
    aes(ymin = pmax(mean_rmse - sd_rmse, 0), ymax = mean_rmse + sd_rmse),
    width = 0.15
  ) +
  geom_point(size = 3) +
  geom_line() +
  labs(
    x = "Information criterion",
    y = "Repeated-CV RMSE",
    shape = "Covariance",
    title = "Riboflavin: out-of-sample predictive performance"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(p2, "riboflavin_cv_rmse")

# 3. C1F-selected cardinality, Riboflavin
p3 <- ggplot(
  rib,
  aes(
    x = criterion,
    y = n_selected,
    shape = cov_method,
    group = cov_method
  )
) +
  geom_point(size = 3) +
  geom_line() +
  labs(
    x = "Information criterion",
    y = "Number of selected genes",
    shape = "Covariance",
    title = "Riboflavin: C1F-selected model cardinality"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(p3, "riboflavin_selected_genes")

# 4. Criterion curves over d, Riboflavin.
# Each criterion is standardised within covariance method to make unlike
# information-criterion scales visually comparable without changing minima.
score_long <- function(sc, cov_label) {
  crits <- criterion_levels
  out <- do.call(
    rbind,
    lapply(crits, function(cr) {
      z <- sc[[cr]]
      sdz <- stats::sd(z)
      z_std <- if (is.finite(sdz) && sdz > 0) (z - mean(z)) / sdz else z - mean(z)
      data.frame(
        cov_method = cov_label,
        d = sc$d,
        criterion = cr,
        standardised_score = z_std
      )
    })
  )
  out
}

curve_data <- rbind(
  score_long(scores_oas, "OAS"),
  score_long(scores_mec, "MEC")
)
curve_data$criterion <- factor(curve_data$criterion, levels = criterion_levels)

p4 <- ggplot(
  curve_data,
  aes(
    x = d,
    y = standardised_score,
    linetype = criterion,
    group = criterion
  )
) +
  geom_line() +
  geom_point() +
  facet_wrap(~ cov_method) +
  scale_x_continuous(breaks = 1:5) +
  labs(
    x = "Candidate structural dimension (d)",
    y = "Standardised criterion score",
    linetype = "Criterion",
    title = "Riboflavin: information-criterion profiles across structural dimension"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

save_plot(p4, "riboflavin_criterion_curves", 9, 5.5)

# 5. Golub classification metrics, mean ±1 SD.
metric_rows <- list(
  c("Accuracy", "mean_accuracy", "sd_accuracy"),
  c("Balanced accuracy", "mean_balanced_accuracy", "sd_balanced_accuracy"),
  c("Sensitivity", "mean_sensitivity", "sd_sensitivity"),
  c("Specificity", "mean_specificity", "sd_specificity"),
  c("F1", "mean_f1", "sd_f1"),
  c("AUC", "mean_auc", "sd_auc")
)

gol_long <- do.call(
  rbind,
  lapply(metric_rows, function(z) {
    data.frame(
      cov_method = gol$cov_method,
      metric = z[1],
      mean = gol[[z[2]]],
      sd = gol[[z[3]]]
    )
  })
)

p5 <- ggplot(
  gol_long,
  aes(
    x = metric,
    y = mean,
    shape = cov_method,
    group = cov_method
  )
) +
  geom_errorbar(
    aes(ymin = pmax(mean - sd, 0), ymax = pmin(mean + sd, 1)),
    width = 0.15
  ) +
  geom_point(size = 3) +
  geom_line() +
  coord_cartesian(ylim = c(0,1)) +
  labs(
    x = "Performance metric",
    y = "Repeated-CV performance",
    shape = "Covariance",
    title = "Golub leukaemia: classification performance"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

save_plot(p5, "golub_classification_performance", 9, 5)

# 6. Golub selected model cardinality
p6 <- ggplot(
  gol,
  aes(
    x = cov_method,
    y = n_selected,
    group = 1
  )
) +
  geom_point(size = 3) +
  geom_line() +
  labs(
    x = "Covariance estimator",
    y = "Number of selected genes",
    title = "Golub leukaemia: C1F-selected model cardinality"
  ) +
  theme_bw()

save_plot(p6, "golub_selected_genes", 7, 5)

# 7. Data-adaptive covariance shrinkage
rho_rib <- unique(rib[c("cov_method","shrinkage_rho")])
rho_rib$dataset <- "Riboflavin"
rho_gol <- unique(gol[c("cov_method","shrinkage_rho")])
rho_gol$dataset <- "Golub leukaemia"
rho <- rbind(rho_rib, rho_gol)

p7 <- ggplot(
  rho,
  aes(
    x = cov_method,
    y = shrinkage_rho,
    shape = dataset
  )
) +
  geom_point(size = 3, position = position_dodge(width = 0.25)) +
  coord_cartesian(ylim = c(0,1)) +
  labs(
    x = "Covariance estimator",
    y = "Data-adaptive shrinkage intensity",
    shape = "Dataset",
    title = "Empirical data-adaptive covariance shrinkage"
  ) +
  theme_bw()

save_plot(p7, "real_data_shrinkage_rho", 8, 5)

cat("Real-data figures generated in ", figure_dir, ".\n", sep = "")
