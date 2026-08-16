#' Fit a Complete SPFC Workflow
#'
#' Fits a complete Shrinkage Principal Fitted Components workflow, including
#' optional structural dimension selection, SPFC fitting, reduced model fitting,
#' evaluation, and variable selection.
#'
#' @param X Numeric predictor matrix.
#' @param y Response vector.
#' @param d Optional structural dimension. If `NULL`, dimension selection is performed.
#' @param d_grid Candidate structural dimensions used when `d = NULL`.
#' @param ytype Character. One of `"auto"`, `"continuous"`, or `"categorical"`.
#' @param cov_method Covariance estimator. One of `"mec"`, `"oas"`, `"sre"`, `"sde"`, or `"cse"`.
#' @param selector Criterion for selecting structural dimension when `d = NULL`.
#' @param nslices Number of slices for continuous responses.
#' @param poly_degree Polynomial degree for continuous response basis.
#' @param variable_method Variable selection method.
#' @param selection_rule Character screening rule: `"c1f"`, `"quantile"`, or `"fixed"`.
#' @param threshold Optional finite numeric screening threshold. It is required for `selection_rule = "fixed"`; for `"quantile"`, `NULL` uses `quantile_cut`; it is ignored by the C1F-calibrated rule.
#' @param quantile_cut Quantile threshold for variable selection.
#' @param centre_x Logical. If `TRUE`, centre predictor columns.
#' @param scale_x Logical. If `TRUE`, scale predictor columns.
#' @param verbose Logical. If `TRUE`, print progress messages.
#' @param validation Character. Validation strategy. One of `"resubstitution"`,
#' `"holdout"`, or `"cv"`.
#' @param test_fraction Numeric value between 0 and 1 giving the test-set
#' proportion when `validation = "holdout"`.
#' @param folds Integer number of folds to use when `validation = "cv"`.
#' @param validation_seed Integer seed used for train-test splitting or
#' cross-validation fold creation.
#' @param ... Additional arguments passed to `spfc_fit()`.
#'
#' @return A list of class `spfc` containing the complete workflow.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' X <- matrix(rnorm(180), nrow = 30, ncol = 6)
#' y <- X[, 1] - 0.7 * X[, 2] + rnorm(30, sd = 0.5)
#' obj <- spfc(
#'   X, y, d = 1, ytype = 'continuous', cov_method = 'oas',
#'   poly_degree = 2, selection_rule = 'quantile',
#'   quantile_cut = 0.75, validation = 'resubstitution',
#'   verbose = FALSE
#' )
#' obj
#' }
#' @export
spfc <- function(
    X,
    y,
    d = NULL,
    d_grid = 1:3,
    ytype = c("auto", "continuous", "categorical"),
    cov_method = c("mec", "oas", "sre", "sde", "cse"),
    selector = c(
      "AIC",
      "BIC",
      "CAIC",
      "ICOMP_IFIM",
      "ICOMP_MISSPEC",
      "CICOMP"
    ),
    nslices = 5,
    poly_degree = 1,
    variable_method = c("adaptive_weighted_l1", "c1f_extension"),
    selection_rule = c("c1f", "quantile", "fixed"),
    threshold = NULL,
    quantile_cut = 0.75,
    validation = c("resubstitution", "holdout", "cv"),
    test_fraction = 0.30,
    folds = 5,
    validation_seed = 123,
    centre_x = TRUE,
    scale_x = FALSE,
    verbose = TRUE,
    ...
){

  ytype <- match.arg(ytype)
  cov_method <- match.arg(cov_method)
  selector <- match.arg(selector)
  variable_method <- match.arg(variable_method)
  selection_rule <- match.arg(selection_rule)
  validation <- match.arg(validation)

  detected_ytype <- detect_y_type(
    y = y,
    ytype = ytype
  )

  dimension_selection <- NULL

  if(is.null(d)){

    if(isTRUE(verbose)){
      message("Selecting structural dimension using ", selector, ".")
    }

    dimension_selection <- spfc_select_dimension(
      X = X,
      y = y,
      d_grid = d_grid,
      cov_method = cov_method,
      ytype = detected_ytype,
      nslices = nslices,
      poly_degree = poly_degree,
      verbose = verbose
    )

    selected_row <- dimension_selection$selected[
      dimension_selection$selected$criterion == selector,
      ,
      drop = FALSE
    ]

    if(nrow(selected_row) != 1){
      stop("Could not identify selected dimension for criterion: ", selector)
    }

    d <- selected_row$selected_d[1]

    if(isTRUE(verbose)){
      message("Selected structural dimension: ", d)
    }
  }

  if(isTRUE(verbose)){
    message("Fitting SPFC model using covariance method: ", cov_method)
  }

  fit <- spfc_fit(
    X = X,
    y = y,
    d = d,
    ytype = detected_ytype,
    cov_method = cov_method,
    nslices = nslices,
    poly_degree = poly_degree,
    centre_x = centre_x,
    scale_x = scale_x,
    ...
  )

  validation_results <- NULL

  if(validation == "resubstitution"){

    if(isTRUE(verbose)){
      message("Fitting reduced downstream model using resubstitution.")
    }

    reduced_model <- fit_reduced_model(
      Z = fit$Z,
      y = y,
      ytype = detected_ytype
    )

    evaluation <- evaluate_reduced_model(
      object = reduced_model,
      Z = fit$Z,
      y = y
    )

    if(isTRUE(verbose)){
      message("Running variable selection.")
    }

    variable_selection <- spfc_select_variables(
      fit = fit,
      method = variable_method,
      selection_rule = selection_rule,
      reduced_model = reduced_model,
      threshold = threshold,
      quantile_cut = quantile_cut
    )

  } else if(validation == "holdout"){

    if(isTRUE(verbose)){
      message("Running holdout validation.")
    }

    split <- create_train_test_split(
      y = y,
      test_fraction = test_fraction,
      seed = validation_seed,
      stratify = detected_ytype == "categorical"
    )

    validation_results <- evaluate_spfc_train_test(
      X = X,
      y = y,
      train_idx = split$train,
      test_idx = split$test,
      d = d,
      ytype = detected_ytype,
      cov_method = cov_method,
      nslices = nslices,
      poly_degree = poly_degree,
      variable_method = variable_method,
      selection_rule = selection_rule,
      threshold = threshold,
      quantile_cut = quantile_cut,
      centre_x = centre_x,
      scale_x = scale_x,
      ...
    )

    fit <- validation_results$fit
    reduced_model <- validation_results$reduced_model
    evaluation <- validation_results$evaluation
    variable_selection <- validation_results$variable_selection

  } else if(validation == "cv"){

    if(isTRUE(verbose)){
      message("Running ", folds, "-fold cross-validation.")
    }

    validation_results <- evaluate_spfc_cv(
      X = X,
      y = y,
      d = d,
      ytype = detected_ytype,
      cov_method = cov_method,
      folds = folds,
      seed = validation_seed,
      nslices = nslices,
      poly_degree = poly_degree,
      variable_method = variable_method,
      selection_rule = selection_rule,
      threshold = threshold,
      quantile_cut = quantile_cut,
      centre_x = centre_x,
      scale_x = scale_x,
      ...
    )

    fit <- validation_results$fit
    reduced_model <- validation_results$reduced_model
    evaluation <- validation_results$aggregate
    variable_selection <- validation_results$variable_selection
  }

  out <- list(
    call = match.call(),
    ytype = detected_ytype,
    cov_method = cov_method,
    selector = selector,
    d = d,
    dimension_selection = dimension_selection,
    fit = fit,
    reduced_model = reduced_model,
    evaluation = evaluation,
    cv_confusion_matrix = if(validation == "cv"){
      validation_results$confusion_matrix
    } else {
      NULL
    },
    variable_selection = variable_selection,
    variable_method = variable_method,
    selection_rule = selection_rule,
    validation = validation,
    validation_results = validation_results,
    test_fraction = test_fraction,
    folds = folds,
    validation_seed = validation_seed,
    quantile_cut = quantile_cut
  )

  class(out) <- c("spfc", "list")

  return(out)
}


#' Print Complete SPFC Workflow
#'
#' @param x Object returned by `spfc()`.
#' @param ... Additional arguments.
#'
#' @return Invisibly returns `x`.
#'
#' @export
print.spfc <- function(x, ...){

  cat("\nComplete SPFC Workflow\n")
  cat("======================\n\n")

  cat("Response type:       ", x$ytype, "\n", sep = "")
  cat("Covariance method:   ", x$cov_method, "\n", sep = "")
  cat("Structural dimension:", x$d, "\n")
  cat("Selector:            ", x$selector, "\n", sep = "")
  cat("Variable method:     ", x$variable_method, "\n", sep = "")
  cat("Validation:          ", x$validation, "\n", sep = "")

  if(!is.null(x$fit$rho)){
    cat("Shrinkage rho:       ", round(x$fit$rho, 6), "\n", sep = "")
  }

  cat("\nReduced model method:\n")

  if(!is.null(x$reduced_model)){
    print(x$reduced_model$method)
  } else {
    print(NA_character_)
  }



  cat("\nNumber of selected variables:\n")
  print(sum(x$variable_selection$selected))

  invisible(x)
}


#' Summarise Complete SPFC Workflow
#'
#' @param object Object returned by `spfc()`.
#' @param ... Additional arguments.
#'
#' @return A list of class `summary.spfc`.
#'
#' @export
summary.spfc <- function(object, ...){

  eval_summary <- NULL

  if(object$ytype == "continuous"){

    if(object$validation == "cv"){

      eval_summary <- data.frame(
        mean_rmse = object$evaluation$mean_rmse,
        sd_rmse = object$evaluation$sd_rmse,
        mean_mae = object$evaluation$mean_mae,
        sd_mae = object$evaluation$sd_mae
      )

    } else {

      eval_summary <- data.frame(
        rmse = object$evaluation$rmse,
        mae = object$evaluation$mae
      )
    }

  } else {

    if(object$validation == "cv"){

      eval_summary <- data.frame(
        mean_accuracy = object$evaluation$mean_accuracy,
        sd_accuracy = object$evaluation$sd_accuracy,
        mean_sensitivity = object$evaluation$mean_sensitivity,
        sd_sensitivity = object$evaluation$sd_sensitivity,
        mean_specificity = object$evaluation$mean_specificity,
        sd_specificity = object$evaluation$sd_specificity,
        mean_precision = object$evaluation$mean_precision,
        sd_precision = object$evaluation$sd_precision,
        mean_f1 = object$evaluation$mean_f1,
        sd_f1 = object$evaluation$sd_f1,
        mean_balanced_accuracy = object$evaluation$mean_balanced_accuracy,
        sd_balanced_accuracy = object$evaluation$sd_balanced_accuracy
      )

    } else {

      eval_summary <- data.frame(
        accuracy = object$evaluation$accuracy,
        sensitivity = object$evaluation$sensitivity,
        specificity = object$evaluation$specificity,
        precision = object$evaluation$precision,
        f1 = object$evaluation$f1,
        balanced_accuracy = object$evaluation$balanced_accuracy
      )
    }
  }

  out <- list(
    ytype = object$ytype,
    cov_method = object$cov_method,
    selector = object$selector,
    d = object$d,
    rho = object$fit$rho,
    reduced_model_method = if(!is.null(object$reduced_model)){
      object$reduced_model$method
    } else {
      NA_character_
    },
    evaluation = eval_summary,
    n_selected = sum(object$variable_selection$selected),
    selected_variables = object$variable_selection$variable[
      object$variable_selection$selected
    ],
    dimension_selection = object$dimension_selection,
    validation = object$validation
  )

  class(out) <- c("summary.spfc", "list")

  return(out)
}

#' Print Summary of Complete SPFC Workflow
#'
#' @param x Object returned by `summary.spfc()`.
#' @param ... Additional arguments.
#'
#' @return Invisibly returns `x`.
#'
#' @export
print.summary.spfc <- function(x, ...){

  cat("\nSummary of Complete SPFC Workflow\n")
  cat("=================================\n\n")

  cat("Response type:       ", x$ytype, "\n", sep = "")
  cat("Covariance method:   ", x$cov_method, "\n", sep = "")
  cat("Structural dimension:", x$d, "\n")
  cat("Selector:            ", x$selector, "\n", sep = "")

  if(!is.null(x$rho)){
    cat("Shrinkage rho:       ", round(x$rho, 6), "\n", sep = "")
  }

  cat("\nReduced model method:\n")
  print(x$reduced_model_method)

  cat("\nEvaluation:\n")
  print(x$evaluation)

  cat("\nNumber of selected variables:\n")
  print(x$n_selected)

  cat("\nSelected variables:\n")

  if (length(x$selected_variables) == 0) {

    cat("None\n")

  } else {

    max_print <- 20

    print(
      utils::head(
        x$selected_variables,
        max_print
      )
    )

    if (length(x$selected_variables) > max_print) {
      cat(
        "\nOutput truncated. Use summary(object)$selected_variables or object$variable_selection for the complete list.\n"
      )
    }
  }

  invisible(x)
}


#' Extract SPFC Fit from Complete Workflow
#'
#' @param object Object returned by `spfc()`.
#' @param ... Additional arguments.
#'
#' @return Estimated SPFC direction matrix from the fitted SPFC component.
#'
#' @export
coef.spfc <- function(object, ...){

  stats::coef(object$fit)
}


#' Extract Reduced Scores from Complete Workflow
#'
#' @param object Object returned by `spfc()`.
#' @param ... Additional arguments.
#'
#' @return Reduced SPFC score matrix.
#'
#' @export
fitted.spfc <- function(object, ...){

  stats::fitted(object$fit)
}


#' Predict SPFC Scores from Complete Workflow
#'
#' @param object Object returned by `spfc()`.
#' @param newdata Numeric matrix or data frame of new predictor observations.
#' @param ... Additional arguments.
#'
#' @return Numeric matrix of predicted SPFC scores.
#'
#' @export
predict.spfc <- function(object, newdata, ...){

  stats::predict(
    object$fit,
    newdata = newdata,
    ...
  )
}
