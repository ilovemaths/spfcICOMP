#' Benchmark SPFC Across Covariance Estimators
#'
#' Fits SPFC models using several covariance estimators and compares their
#' downstream predictive performance, runtime, shrinkage behaviour, and
#' variable-selection output.
#'
#' @param X Numeric predictor matrix.
#' @param y Response vector.
#' @param d Structural dimension.
#' @param ytype Character. One of `"auto"`, `"continuous"`, or `"categorical"`.
#' @param methods Character vector of covariance estimators to compare.
#' @param variable_method Variable-selection method.
#' @param quantile_cut Quantile threshold for adaptive variable selection.
#' @param validation Character. Validation strategy. One of `"resubstitution"`,
#' `"holdout"`, or `"cv"`.
#' @param test_fraction Proportion of observations allocated to the test set
#' when `validation = "holdout"`.
#' @param folds Number of folds used when `validation = "cv"`.
#' @param validation_seed Random seed used for train-test splitting or
#' cross-validation fold creation.
#' @param verbose Logical. If `TRUE`, print progress messages.
#' @param ... Additional arguments passed to `spfc()`.
#'
#' @return A list of class `spfc_benchmark`.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' X <- matrix(rnorm(180), nrow = 30, ncol = 6)
#' y <- X[, 1] - 0.7 * X[, 2] + rnorm(30, sd = 0.5)
#' bench <- benchmark_spfc(
#'   X, y, d = 1, ytype = 'continuous',
#'   methods = c('oas', 'sre'),
#'   variable_method = 'adaptive_weighted_l1',
#'   quantile_cut = 0.75, validation = 'resubstitution',
#'   verbose = FALSE
#' )
#' bench$summary
#' }
#' @export
benchmark_spfc <- function(
    X,
    y,
    d,
    ytype = c("auto", "continuous", "categorical"),
    methods = c("mec", "oas", "sre", "sde", "cse"),
    variable_method = "adaptive_weighted_l1",
    quantile_cut = 0.75,
    validation = c("resubstitution", "holdout", "cv"),
    test_fraction = 0.30,
    folds = 5,
    validation_seed = 123,
    verbose = TRUE,
    ...
){

  ytype <- match.arg(ytype)
  validation <- match.arg(validation)

  X <- as.matrix(X)

  detected_ytype <- detect_y_type(
    y = y,
    ytype = ytype
  )

  methods <- match.arg(
    methods,
    choices = c("mec", "oas", "sre", "sde", "cse"),
    several.ok = TRUE
  )

  results <- vector(
    mode = "list",
    length = length(methods)
  )

  summary_rows <- vector(
    mode = "list",
    length = length(methods)
  )

  for(i in seq_along(methods)){

    method_i <- methods[i]

    if(isTRUE(verbose)){
      message(
        "Running SPFC benchmark with covariance method: ",
        method_i
      )
    }

    t_start <- Sys.time()

    fit_i <- spfc(
      X = X,
      y = y,
      d = d,
      ytype = detected_ytype,
      cov_method = method_i,
      variable_method = variable_method,
      quantile_cut = quantile_cut,
      validation = validation,
      test_fraction = test_fraction,
      folds = folds,
      validation_seed = validation_seed,
      verbose = FALSE,
      ...
    )

    t_end <- Sys.time()

    runtime_sec <- as.numeric(
      difftime(
        t_end,
        t_start,
        units = "secs"
      )
    )

    ev <- fit_i$evaluation

    n_selected <- sum(
      fit_i$variable_selection$selected
    )

    if(detected_ytype == "continuous"){

      if(validation == "cv"){

        summary_rows[[i]] <- data.frame(
          cov_method = method_i,
          validation = validation,
          ytype = detected_ytype,
          d = fit_i$d,
          rho = fit_i$fit$rho,
          nslices_used = fit_i$fit$nslices_used,
          reduced_model = if(!is.null(fit_i$reduced_model)){
            fit_i$reduced_model$method
          } else {
            NA_character_
          },
          runtime_sec = runtime_sec,
          rmse = NA_real_,
          mae = NA_real_,
          mean_rmse = ev$mean_rmse,
          sd_rmse = ev$sd_rmse,
          mean_mae = ev$mean_mae,
          sd_mae = ev$sd_mae,
          accuracy = NA_real_,
          sensitivity = NA_real_,
          specificity = NA_real_,
          precision = NA_real_,
          f1 = NA_real_,
          balanced_accuracy = NA_real_,
          mean_accuracy = NA_real_,
          mean_sensitivity = NA_real_,
          mean_specificity = NA_real_,
          mean_precision = NA_real_,
          mean_f1 = NA_real_,
          mean_balanced_accuracy = NA_real_,
          n_selected = n_selected,
          stringsAsFactors = FALSE
        )

      } else {

        summary_rows[[i]] <- data.frame(
          cov_method = method_i,
          validation = validation,
          ytype = detected_ytype,
          d = fit_i$d,
          rho = fit_i$fit$rho,
          nslices_used = fit_i$fit$nslices_used,
          reduced_model = if(!is.null(fit_i$reduced_model)){
            fit_i$reduced_model$method
          } else {
            NA_character_
          },
          runtime_sec = runtime_sec,
          rmse = ev$rmse,
          mae = ev$mae,
          mean_rmse = NA_real_,
          sd_rmse = NA_real_,
          mean_mae = NA_real_,
          sd_mae = NA_real_,
          accuracy = NA_real_,
          sensitivity = NA_real_,
          specificity = NA_real_,
          precision = NA_real_,
          f1 = NA_real_,
          balanced_accuracy = NA_real_,
          mean_accuracy = NA_real_,
          mean_sensitivity = NA_real_,
          mean_specificity = NA_real_,
          mean_precision = NA_real_,
          mean_f1 = NA_real_,
          mean_balanced_accuracy = NA_real_,
          n_selected = n_selected,
          stringsAsFactors = FALSE
        )
      }

    } else {

      if(validation == "cv"){

        summary_rows[[i]] <- data.frame(
          cov_method = method_i,
          validation = validation,
          ytype = detected_ytype,
          d = fit_i$d,
          rho = fit_i$fit$rho,
          nslices_used = fit_i$fit$nslices_used,
          reduced_model = if(!is.null(fit_i$reduced_model)){
            fit_i$reduced_model$method
          } else {
            NA_character_
          },
          runtime_sec = runtime_sec,
          rmse = NA_real_,
          mae = NA_real_,
          mean_rmse = NA_real_,
          sd_rmse = NA_real_,
          mean_mae = NA_real_,
          sd_mae = NA_real_,
          accuracy = NA_real_,
          sensitivity = NA_real_,
          specificity = NA_real_,
          precision = NA_real_,
          f1 = NA_real_,
          balanced_accuracy = NA_real_,
          mean_accuracy = ev$mean_accuracy,
          mean_sensitivity = ev$mean_sensitivity,
          mean_specificity = ev$mean_specificity,
          mean_precision = ev$mean_precision,
          mean_f1 = ev$mean_f1,
          mean_balanced_accuracy = ev$mean_balanced_accuracy,
          n_selected = n_selected,
          stringsAsFactors = FALSE
        )

      } else {

        summary_rows[[i]] <- data.frame(
          cov_method = method_i,
          validation = validation,
          ytype = detected_ytype,
          d = fit_i$d,
          rho = fit_i$fit$rho,
          nslices_used = fit_i$fit$nslices_used,
          reduced_model = if(!is.null(fit_i$reduced_model)){
            fit_i$reduced_model$method
          } else {
            NA_character_
          },
          runtime_sec = runtime_sec,
          rmse = NA_real_,
          mae = NA_real_,
          mean_rmse = NA_real_,
          sd_rmse = NA_real_,
          mean_mae = NA_real_,
          sd_mae = NA_real_,
          accuracy = ev$accuracy,
          sensitivity = ev$sensitivity,
          specificity = ev$specificity,
          precision = ev$precision,
          f1 = ev$f1,
          balanced_accuracy = ev$balanced_accuracy,
          mean_accuracy = NA_real_,
          mean_sensitivity = NA_real_,
          mean_specificity = NA_real_,
          mean_precision = NA_real_,
          mean_f1 = NA_real_,
          mean_balanced_accuracy = NA_real_,
          n_selected = n_selected,
          stringsAsFactors = FALSE
        )
      }
    }

    results[[i]] <- fit_i
  }

  names(results) <- methods

  summary_table <- do.call(
    rbind,
    summary_rows
  )

  out <- list(
    call = match.call(),
    ytype = detected_ytype,
    d = d,
    methods = methods,
    variable_method = variable_method,
    quantile_cut = quantile_cut,
    validation = validation,
    test_fraction = test_fraction,
    folds = folds,
    validation_seed = validation_seed,
    fits = results,
    summary = summary_table
  )

  class(out) <- c(
    "spfc_benchmark",
    "list"
  )

  out
}


# Build a compact benchmark table containing only fields that apply to the
# response type and validation strategy. The complete rectangular schema
# remains available in `x$summary` for backwards compatibility.
#
# @param object Object returned by `benchmark_spfc()`.
# @return A data frame without structurally inapplicable columns.
# @noRd
benchmark_display_table <- function(object){

  tab <- object$summary

  metric_columns <- if(object$ytype == "continuous"){
    if(object$validation == "cv"){
      c("mean_rmse", "sd_rmse", "mean_mae", "sd_mae")
    } else {
      c("rmse", "mae")
    }
  } else {
    if(object$validation == "cv"){
      c(
        "mean_accuracy",
        "mean_sensitivity",
        "mean_specificity",
        "mean_precision",
        "mean_f1",
        "mean_balanced_accuracy"
      )
    } else {
      c(
        "accuracy",
        "sensitivity",
        "specificity",
        "precision",
        "f1",
        "balanced_accuracy"
      )
    }
  }

  display_columns <- c(
    "cov_method",
    "d",
    "reduced_model",
    "runtime_sec",
    metric_columns,
    "n_selected"
  )

  display_columns <- intersect(
    display_columns,
    names(tab)
  )

  tab[
    ,
    display_columns,
    drop = FALSE
  ]
}


# Collect method-specific covariance diagnostics in long form. Rows are
# created only where a diagnostic is defined, avoiding misleading missing
# values for estimators that do not use that quantity.
#
# @param object Object returned by `benchmark_spfc()`.
# @return A data frame with covariance method, diagnostic and value.
# @noRd
benchmark_covariance_details <- function(object){

  tab <- object$summary
  rows <- list()
  row_index <- 0L

  if("rho" %in% names(tab)){
    keep_rho <- is.finite(tab$rho)

    if(any(keep_rho)){
      row_index <- row_index + 1L
      rows[[row_index]] <- data.frame(
        cov_method = tab$cov_method[keep_rho],
        diagnostic = "rho",
        value = tab$rho[keep_rho],
        stringsAsFactors = FALSE
      )
    }
  }

  if("nslices_used" %in% names(tab)){
    keep_slices <- !is.na(tab$nslices_used)

    if(any(keep_slices)){
      row_index <- row_index + 1L
      rows[[row_index]] <- data.frame(
        cov_method = tab$cov_method[keep_slices],
        diagnostic = "nslices_used",
        value = as.numeric(tab$nslices_used[keep_slices]),
        stringsAsFactors = FALSE
      )
    }
  }

  if(length(rows) == 0L){
    return(data.frame(
      cov_method = character(),
      diagnostic = character(),
      value = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  details <- do.call(
    rbind,
    rows
  )

  rownames(details) <- NULL
  details
}


#' Print SPFC Benchmark
#'
#' @param x Object returned by `benchmark_spfc()`.
#' @param ... Additional arguments.
#'
#' @return Invisibly returns `x`.
#'
#' @export
print.spfc_benchmark <- function(x, ...){

  cat("\nSPFC Benchmark\n")
  cat("==============\n\n")

  cat("Response type:       ", x$ytype, "\n", sep = "")
  cat("Structural dimension:", x$d, "\n")
  cat("Variable method:     ", x$variable_method, "\n", sep = "")
  cat("Validation:          ", x$validation, "\n", sep = "")
  cat("Methods compared:    ", paste(x$methods, collapse = ", "), "\n", sep = "")

  cat("\nPerformance summary:\n")
  print(benchmark_display_table(x))

  covariance_details <- benchmark_covariance_details(x)

  if(nrow(covariance_details) > 0L){
    cat("\nApplicable covariance diagnostics:\n")
    print(covariance_details)
  }

  invisible(x)
}

#' Summarise SPFC Benchmark
#'
#' @param object Object returned by `benchmark_spfc()`.
#' @param metric Optional selection metric.
#' @param smaller_is_better Optional logical indicating whether smaller values
#' are better.
#' @param ... Additional arguments.
#'
#' @return A list of class `summary.spfc_benchmark`.
#'
#' @export
summary.spfc_benchmark <- function(
    object,
    metric = NULL,
    smaller_is_better = NULL,
    ...
){

  tab <- object$summary

  if(is.null(metric)){

    if(object$ytype == "continuous"){

      if(object$validation == "cv"){
        metric <- "mean_rmse"
      } else {
        metric <- "rmse"
      }

      smaller_is_better <- TRUE

    } else {

      if(object$validation == "cv"){
        metric <- "mean_accuracy"
      } else {
        metric <- "accuracy"
      }

      smaller_is_better <- FALSE
    }
  }

  if(!metric %in% names(tab)){
    stop("metric not found in benchmark summary.")
  }

  metric_values <- tab[[metric]]

  if(all(is.na(metric_values))){
    stop("metric contains only missing values.")
  }

  if(isTRUE(smaller_is_better)){

    best_idx <- which.min(metric_values)

  } else {

    best_idx <- which.max(metric_values)
  }

  out <- list(
    ytype = object$ytype,
    d = object$d,
    methods = object$methods,
    variable_method = object$variable_method,
    quantile_cut = object$quantile_cut,
    validation = object$validation,
    metric = metric,
    smaller_is_better = smaller_is_better,
    best_method = tab$cov_method[best_idx],
    best_row = tab[best_idx, , drop = FALSE],
    table = tab,
    display_table = benchmark_display_table(object),
    covariance_details = benchmark_covariance_details(object)
  )

  class(out) <- c(
    "summary.spfc_benchmark",
    "list"
  )

  out
}

#' Print Summary of SPFC Benchmark
#'
#' @param x Object returned by `summary.spfc_benchmark()`.
#' @param ... Additional arguments.
#'
#' @return Invisibly returns `x`.
#'
#' @export
print.summary.spfc_benchmark <- function(x, ...){

  cat("\nSummary of SPFC Benchmark\n")
  cat("=========================\n\n")

  cat("Response type:       ", x$ytype, "\n", sep = "")
  cat("Structural dimension:", x$d, "\n")
  cat("Variable method:     ", x$variable_method, "\n", sep = "")
  cat("Validation:          ", x$validation, "\n", sep = "")
  cat("Selection metric:    ", x$metric, "\n", sep = "")
  cat("Smaller is better:   ", x$smaller_is_better, "\n", sep = "")
  cat("Best method:         ", x$best_method, "\n", sep = "")

  cat("\nBenchmark table:\n")
  print(x$display_table)

  if(nrow(x$covariance_details) > 0L){
    cat("\nApplicable covariance diagnostics:\n")
    print(x$covariance_details)
  }

  invisible(x)
}
