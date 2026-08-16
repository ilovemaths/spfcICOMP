#' Safe Mean
#'
#' Computes the mean while safely handling all-NA vectors.
#'
#' @param x Numeric vector.
#'
#' @return Numeric mean or NA.
#'
#' @noRd
safe_mean <- function(x){

  if(all(is.na(x))){
    return(NA_real_)
  }

  mean(x, na.rm = TRUE)
}


#' Safe Standard Deviation
#'
#' Computes the standard deviation while safely handling all-NA vectors.
#'
#' @param x Numeric vector.
#'
#' @return Numeric standard deviation or NA.
#'
#' @noRd
safe_sd <- function(x){

  if(all(is.na(x))){
    return(NA_real_)
  }

  stats::sd(x, na.rm = TRUE)
}


#' Summarise SPFC Simulation Results
#'
#' Aggregates Monte Carlo simulation results by covariance estimator,
#' variable-selection method, and design factors.
#'
#' @param results Data frame returned by `run_spfc_simulation()`.
#'
#' @return A summary data frame.
#'
#' @examples
#' \donttest{
#' results <- run_spfc_simulation(
#'   response_type = 'continuous', nrep = 1, n = 30, p = 8,
#'   d = 1, s = 3, rho_x = 0.3, snr = 2,
#'   cov_methods = c('oas', 'sre'),
#'   variable_methods = 'adaptive_weighted_l1', seed = 123
#' )
#' summarise_spfc_simulation(results)
#' }
#' @export
summarise_spfc_simulation <- function(results){

  if(!is.data.frame(results)){
    stop("results must be a data frame.")
  }

  required_cols <- c(
    "ytype",
    "cov_method",
    "variable_method",
    "n",
    "p",
    "true_d",
    "fitted_d",
    "s",
    "rho_x",
    "runtime_sec",
    "subspace_distance",
    "precision",
    "recall",
    "f1_variable",
    "n_selected"
  )

  missing_cols <- setdiff(required_cols, names(results))

  if(length(missing_cols) > 0){
    stop(
      "results is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  grouping_cols <- c(
    "ytype",
    "cov_method",
    "variable_method",
    "n",
    "p",
    "true_d",
    "s",
    "rho_x"
  )

  if("criterion" %in% names(results)){
    grouping_cols <- append(grouping_cols, "criterion", after = 2)
  } else {
    grouping_cols <- c(grouping_cols, "fitted_d")
  }

  if("selection_rule" %in% names(results)){
    grouping_cols <- append(grouping_cols, "selection_rule", after = 3)
  }

  # Preserve all factorial response-strength factors. Continuous simulations
  # must never pool distinct SNR settings, and categorical simulations must
  # never pool distinct signal-strength settings.
  if("snr" %in% names(results) && any(!is.na(results$snr))){
    grouping_cols <- c(grouping_cols, "snr")
  }
  if("signal_strength" %in% names(results) && any(!is.na(results$signal_strength))){
    grouping_cols <- c(grouping_cols, "signal_strength")
  }

  split_factor <- interaction(
    results[, grouping_cols],
    drop = TRUE,
    lex.order = TRUE
  )

  groups <- split(results, split_factor)

  summary_list <- lapply(groups, function(df){

    base <- df[1, grouping_cols, drop = FALSE]

    out <- data.frame(
      base,
      nrep = nrow(df),
      dimension_recovery_rate = if("dimension_correct" %in% names(df)) safe_mean(df$dimension_correct) else NA_real_,
      mean_selected_d = if("fitted_d" %in% names(df)) safe_mean(df$fitted_d) else NA_real_,
      sd_selected_d = if("fitted_d" %in% names(df)) safe_sd(df$fitted_d) else NA_real_,

      mean_runtime_sec = safe_mean(df$runtime_sec),
      sd_runtime_sec = safe_sd(df$runtime_sec),

      mean_subspace_distance = safe_mean(df$subspace_distance),
      sd_subspace_distance = safe_sd(df$subspace_distance),

      mean_precision = safe_mean(df$precision),
      sd_precision = safe_sd(df$precision),

      mean_recall = safe_mean(df$recall),
      sd_recall = safe_sd(df$recall),

      mean_f1_variable = safe_mean(df$f1_variable),
      sd_f1_variable = safe_sd(df$f1_variable),

      mean_n_selected = safe_mean(df$n_selected),
      sd_n_selected = safe_sd(df$n_selected),

      mean_c1f = if("c1f_l1" %in% names(df)) safe_mean(df$c1f_l1) else NA_real_,
      mean_c1f_complexity_fraction = if("c1f_complexity_fraction" %in% names(df)) safe_mean(df$c1f_complexity_fraction) else NA_real_,
      mean_c1f_loading_scale = if("c1f_loading_scale" %in% names(df)) safe_mean(df$c1f_loading_scale) else NA_real_,
      mean_c1f_hd_factor = if("c1f_hd_factor" %in% names(df)) safe_mean(df$c1f_hd_factor) else NA_real_,
      mean_c1f_complexity_multiplier = if("c1f_complexity_multiplier" %in% names(df)) safe_mean(df$c1f_complexity_multiplier) else NA_real_,
      mean_c1f_global_penalty = if("c1f_global_penalty" %in% names(df)) safe_mean(df$c1f_global_penalty) else NA_real_
    )

    if("rmse" %in% names(df)){
      out$mean_rmse <- safe_mean(df$rmse)
      out$sd_rmse <- safe_sd(df$rmse)
    }

    if("mae" %in% names(df)){
      out$mean_mae <- safe_mean(df$mae)
      out$sd_mae <- safe_sd(df$mae)
    }

    if("accuracy" %in% names(df)){
      out$mean_accuracy <- safe_mean(df$accuracy)
      out$sd_accuracy <- safe_sd(df$accuracy)
    }

    if("sensitivity" %in% names(df)){
      out$mean_sensitivity <- safe_mean(df$sensitivity)
      out$sd_sensitivity <- safe_sd(df$sensitivity)
    }

    if("specificity" %in% names(df)){
      out$mean_specificity <- safe_mean(df$specificity)
      out$sd_specificity <- safe_sd(df$specificity)
    }

    if("f1_classification" %in% names(df)){
      out$mean_f1_classification <- safe_mean(df$f1_classification)
      out$sd_f1_classification <- safe_sd(df$f1_classification)
    }

    return(out)
  })

  summary <- do.call(rbind, summary_list)

  rownames(summary) <- NULL

  return(summary)
}


#' Rank Covariance Methods from Simulation Summaries
#'
#' Ranks covariance estimators according to a chosen performance metric.
#'
#' @param summary Data frame returned by `summarise_spfc_simulation()`.
#' @param metric Metric column used for ranking.
#' @param smaller_is_better Logical. If `TRUE`, lower metric values are better.
#'
#' @return Ranked data frame.
#'
#' @examples
#' summary <- data.frame(
#'   cov_method = c('oas', 'sre'),
#'   mean_rmse = c(0.50, 0.60)
#' )
#' rank_covariance_methods(summary, metric = 'mean_rmse')
#' @export
rank_covariance_methods <- function(
    summary,
    metric,
    smaller_is_better = TRUE
){

  if(!is.data.frame(summary)){
    stop("summary must be a data frame.")
  }

  if(!(metric %in% names(summary))){
    stop("metric not found in summary.")
  }

  metric_values <- summary[[metric]]

  if(smaller_is_better){
    ord <- order(metric_values, decreasing = FALSE, na.last = TRUE)
  } else {
    ord <- order(metric_values, decreasing = TRUE, na.last = TRUE)
  }

  ranked <- summary[ord, , drop = FALSE]

  ranked$rank <- seq_len(nrow(ranked))

  rownames(ranked) <- NULL

  return(ranked)
}


#' Summarise Results by Covariance Method
#'
#' Produces a compact method-level summary across all design settings.
#'
#' @param results Data frame returned by `run_spfc_simulation()`.
#'
#' @return A compact summary data frame.
#'
#' @noRd
summarise_by_method <- function(results){

  if(!is.data.frame(results)){
    stop("results must be a data frame.")
  }

  if(!("cov_method" %in% names(results))){
    stop("results must contain cov_method.")
  }

  groups <- split(results, results$cov_method)

  out_list <- lapply(groups, function(df){

    data.frame(
      cov_method = df$cov_method[1],
      nrep = nrow(df),

      mean_runtime_sec = safe_mean(df$runtime_sec),
      mean_subspace_distance = safe_mean(df$subspace_distance),

      mean_rmse = if("rmse" %in% names(df)) safe_mean(df$rmse) else NA_real_,
      mean_accuracy = if("accuracy" %in% names(df)) safe_mean(df$accuracy) else NA_real_,

      mean_precision = safe_mean(df$precision),
      mean_recall = safe_mean(df$recall),
      mean_f1_variable = safe_mean(df$f1_variable),

      mean_n_selected = safe_mean(df$n_selected)
    )
  })

  out <- do.call(rbind, out_list)

  rownames(out) <- NULL

  return(out)
}
