# 03_audit_and_aggregate.R
# Post-processing only: does NOT refit SPFC models.

dir.create("thesis_analysis/tables_v9006", recursive = TRUE, showWarnings = FALSE)
dir.create("thesis_analysis/outputs_v9006", recursive = TRUE, showWarnings = FALSE)

result_file <- "thesis_analysis/outputs_v9006/full_simulation_results.rds"
design_file <- "thesis_analysis/outputs_v9006/full_design.csv"
stopifnot(file.exists(result_file), file.exists(design_file))
res <- readRDS(result_file)
design <- read.csv(design_file, check.names = FALSE)

criteria_expected <- c("AIC", "BIC", "CAIC", "ICOMP_IFIM", "CICOMP")
cov_expected <- c("oas", "mec")

cat("\n============================================================\n")
cat("SPFC-ICOMP DEFINITIVE SIMULATION AUDIT\n")
cat("============================================================\n")

stopifnot(is.data.frame(res), is.data.frame(design))
stopifnot(nrow(design) == 19200L)
stopifnot(nrow(res) == 96000L)
stopifnot(setequal(unique(res$criterion), criteria_expected))
stopifnot(setequal(unique(res$cov_method), cov_expected))
stopifnot(setequal(sort(unique(res$n)), c(50L,100L)))
stopifnot(setequal(sort(unique(res$p)), c(20L,100L,500L)))
stopifnot(setequal(sort(unique(res$true_d)), c(1L,2L)))
stopifnot(setequal(sort(unique(res$rho_x)), c(0.2,0.8)))
stopifnot(setequal(sort(unique(res$snr)), c(1,3)))
stopifnot(all(res$s == 5L))
stopifnot(all(res$ytype == "continuous"))
stopifnot("prediction_sample" %in% names(res))
stopifnot("n_prediction" %in% names(res))
stopifnot(all(res$prediction_sample == "independent_test"))
stopifnot(all(res$n_prediction == res$n))

criterion_count <- aggregate(criterion ~ design_id, res, function(x) length(unique(x)))
stopifnot(all(criterion_count$criterion == 5L))
stopifnot(!anyDuplicated(paste(res$design_id, res$criterion, sep="::")))

base_rep <- unique(res[c("cov_method","n","p","true_d","rho_x","snr","replicate")])
base_counts <- aggregate(replicate ~ cov_method+n+p+true_d+rho_x+snr,
                         base_rep, function(x) length(unique(x)))
stopifnot(nrow(base_counts) == 96L, all(base_counts$replicate == 200L))

crit_counts <- aggregate(replicate ~ cov_method+criterion+n+p+true_d+rho_x+snr,
                         res, function(x) length(unique(x)))
stopifnot(nrow(crit_counts) == 480L, all(crit_counts$replicate == 200L))

cat("Design rows:", nrow(design), "\n")
cat("Criterion-level rows:", nrow(res), "\n")
cat("Base configurations:", nrow(base_counts), "\n")
cat("Criterion-specific cells:", nrow(crit_counts), "\n")

metrics <- c("fitted_d","dimension_correct","shrinkage_rho","c1f_l1",
             "c1f_global_penalty","runtime_sec","precision","recall",
             "f1_variable","n_selected","rmse","mae")
missing_tbl <- data.frame(
  variable = metrics,
  n_missing = vapply(metrics, function(v) sum(is.na(res[[v]])), integer(1)),
  pct_missing = vapply(metrics, function(v) 100*mean(is.na(res[[v]])), numeric(1))
)
write.csv(missing_tbl, "thesis_analysis/tables_v9006/audit_missingness.csv", row.names=FALSE)

if (any(!is.na(res$subspace_distance) & res$dimension_correct == 0))
  warning("Subspace distance found where fitted dimension is incorrect.")

stopifnot(all(res$fitted_d %in% 1:5))
stopifnot(all(res$dimension_correct %in% c(0,1)))
stopifnot(all(res$n_selected >= 0 & res$n_selected <= res$p, na.rm=TRUE))
stopifnot(all(res$precision >= 0 & res$precision <= 1, na.rm=TRUE))
stopifnot(all(res$recall >= 0 & res$recall <= 1, na.rm=TRUE))
stopifnot(all(res$f1_variable >= 0 & res$f1_variable <= 1, na.rm=TRUE))
stopifnot(all(res$rmse >= 0, na.rm=TRUE), all(res$mae >= 0, na.rm=TRUE))

safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm=TRUE)
safe_sd <- function(x) if (sum(!is.na(x)) < 2L) NA_real_ else sd(x, na.rm=TRUE)

summarise_group <- function(df) data.frame(
  nrep = length(unique(df$replicate)),
  dimension_recovery_rate = mean(df$dimension_correct, na.rm=TRUE),
  mean_selected_d = mean(df$fitted_d, na.rm=TRUE),
  sd_selected_d = sd(df$fitted_d, na.rm=TRUE),
  mean_rmse = safe_mean(df$rmse), sd_rmse = safe_sd(df$rmse),
  mean_mae = safe_mean(df$mae), sd_mae = safe_sd(df$mae),
  mean_subspace_distance = safe_mean(df$subspace_distance),
  sd_subspace_distance = safe_sd(df$subspace_distance),
  mean_precision = safe_mean(df$precision), sd_precision = safe_sd(df$precision),
  mean_recall = safe_mean(df$recall), sd_recall = safe_sd(df$recall),
  mean_f1_variable = safe_mean(df$f1_variable), sd_f1_variable = safe_sd(df$f1_variable),
  mean_n_selected = safe_mean(df$n_selected), sd_n_selected = safe_sd(df$n_selected),
  mean_shrinkage_rho = safe_mean(df$shrinkage_rho),
  mean_c1f = safe_mean(df$c1f_l1),
  mean_c1f_global_penalty = safe_mean(df$c1f_global_penalty),
  mean_runtime_sec = safe_mean(df$runtime_sec)
)

key <- interaction(res$cov_method,res$criterion,res$n,res$p,res$true_d,res$rho_x,res$snr,
                   drop=TRUE, lex.order=TRUE)
sp <- split(res,key)
full_summary <- do.call(rbind,lapply(sp,function(df)cbind(
  df[1,c("cov_method","criterion","n","p","true_d","rho_x","snr")],
  summarise_group(df)
)))
rownames(full_summary) <- NULL
full_summary <- full_summary[order(full_summary$p,full_summary$n,full_summary$true_d,
                                   full_summary$rho_x,full_summary$snr,
                                   full_summary$cov_method,
                                   match(full_summary$criterion,criteria_expected)),]
stopifnot(nrow(full_summary)==480L, all(full_summary$nrep==200L))
write.csv(full_summary,"thesis_analysis/tables_v9006/simulation_full_factorial_summary.csv",row.names=FALSE)

write.csv(full_summary[c("n","p","true_d","rho_x","snr","cov_method","criterion",
                         "dimension_recovery_rate","mean_selected_d","sd_selected_d")],
          "thesis_analysis/tables_v9006/table_dimension_recovery.csv",row.names=FALSE)

write.csv(full_summary[c("n","p","true_d","rho_x","snr","cov_method","criterion",
                         "mean_rmse","sd_rmse","mean_mae","sd_mae",
                         "mean_subspace_distance","sd_subspace_distance")],
          "thesis_analysis/tables_v9006/table_prediction_subspace.csv",row.names=FALSE)

write.csv(full_summary[c("n","p","true_d","rho_x","snr","cov_method","criterion",
                         "mean_precision","sd_precision","mean_recall","sd_recall",
                         "mean_f1_variable","sd_f1_variable","mean_n_selected","sd_n_selected")],
          "thesis_analysis/tables_v9006/table_variable_selection.csv",row.names=FALSE)

overall_key <- interaction(res$cov_method,res$criterion,drop=TRUE,lex.order=TRUE)
overall <- do.call(rbind,lapply(split(res,overall_key),function(df)data.frame(
  cov_method=df$cov_method[1], criterion=df$criterion[1], n_records=nrow(df),
  dimension_recovery_rate=mean(df$dimension_correct,na.rm=TRUE),
  mean_selected_d=mean(df$fitted_d,na.rm=TRUE),
  mean_rmse=safe_mean(df$rmse), mean_mae=safe_mean(df$mae),
  mean_subspace_distance=safe_mean(df$subspace_distance),
  mean_precision=safe_mean(df$precision), mean_recall=safe_mean(df$recall),
  mean_f1_variable=safe_mean(df$f1_variable),
  mean_n_selected=safe_mean(df$n_selected), mean_runtime_sec=safe_mean(df$runtime_sec)
)))
rownames(overall)<-NULL
overall <- overall[order(overall$cov_method,match(overall$criterion,criteria_expected)),]
write.csv(overall,"thesis_analysis/tables_v9006/table_overall_method_summary.csv",row.names=FALSE)

scenario_cols <- c("n","p","true_d","rho_x","snr","cov_method")
sk <- interaction(full_summary[scenario_cols],drop=TRUE,lex.order=TRUE)
winlist <- list(); ii <- 0L
for(df in split(full_summary,sk)){
  best_high <- function(v){z<-df[[v]]; if(all(is.na(z))) return(NA_character_); paste(df$criterion[z==max(z,na.rm=TRUE)],collapse=";")}
  best_low  <- function(v){z<-df[[v]]; if(all(is.na(z))) return(NA_character_); paste(df$criterion[z==min(z,na.rm=TRUE)],collapse=";")}
  ii<-ii+1L
  winlist[[ii]] <- cbind(df[1,scenario_cols],
    best_dimension=best_high("dimension_recovery_rate"),
    best_rmse=best_low("mean_rmse"), best_mae=best_low("mean_mae"),
    best_subspace=best_low("mean_subspace_distance"),
    best_precision=best_high("mean_precision"), best_recall=best_high("mean_recall"),
    best_f1=best_high("mean_f1_variable"))
}
winners <- do.call(rbind,winlist); rownames(winners)<-NULL
write.csv(winners,"thesis_analysis/tables_v9006/table_scenario_winners.csv",row.names=FALSE)

count_wins <- function(column){
  x<-winners[[column]]; out<-setNames(integer(length(criteria_expected)),criteria_expected)
  for(cr in criteria_expected) out[cr] <- sum(grepl(paste0("(^|;)",cr,"($|;)"),x))
  data.frame(criterion=names(out),wins=as.integer(out))
}
metrics_win <- c(dimension_recovery="best_dimension",rmse="best_rmse",mae="best_mae",
                 subspace="best_subspace",precision="best_precision",recall="best_recall",f1="best_f1")
wcount <- do.call(rbind,lapply(names(metrics_win),function(m)cbind(metric=m,count_wins(metrics_win[[m]]))))
rownames(wcount)<-NULL
write.csv(wcount,"thesis_analysis/tables_v9006/table_criterion_win_counts.csv",row.names=FALSE)

factor_summary <- function(fac){
  grp <- c(fac,"cov_method","criterion")
  kk <- interaction(res[grp],drop=TRUE,lex.order=TRUE)
  ans <- do.call(rbind,lapply(split(res,kk),function(df)cbind(df[1,grp],data.frame(
    dimension_recovery_rate=mean(df$dimension_correct,na.rm=TRUE),
    mean_rmse=safe_mean(df$rmse), mean_subspace_distance=safe_mean(df$subspace_distance),
    mean_precision=safe_mean(df$precision), mean_recall=safe_mean(df$recall),
    mean_f1_variable=safe_mean(df$f1_variable), mean_n_selected=safe_mean(df$n_selected)))))
  rownames(ans)<-NULL; ans
}
for(fac in c("n","p","true_d","rho_x","snr"))
  write.csv(factor_summary(fac),paste0("thesis_analysis/tables_v9006/marginal_by_",fac,".csv"),row.names=FALSE)

manifest <- data.frame(
  item=c("package_version_expected","design_rows","result_rows","base_configurations",
         "criterion_specific_cells","replicates_per_cell","criteria","covariance_methods",
         "candidate_dimensions","true_active_predictors"),
  value=c("0.0.0.9006",nrow(design),nrow(res),nrow(base_counts),nrow(crit_counts),
          paste(range(crit_counts$replicate),collapse="-"),paste(criteria_expected,collapse=", "),
          paste(cov_expected,collapse=", "),"1,2,3,4,5","5")
)
write.csv(manifest,"thesis_analysis/tables_v9006/simulation_audit_manifest.csv",row.names=FALSE)
saveRDS(list(results=res,full_summary=full_summary,overall=overall,winners=winners,
             winner_counts=wcount,audit_missingness=missing_tbl,manifest=manifest),
        "thesis_analysis/outputs_v9006/postprocessed_simulation_results.rds")

cat("Audit complete. Full factorial summary rows:",nrow(full_summary),"\n")
cat("Outputs written to thesis_analysis/tables and thesis_analysis/outputs.\n")
