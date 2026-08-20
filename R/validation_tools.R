#' Determine Whether a Response Should Be Stratified
#'
#' @param y Response vector.
#'
#' @return Logical scalar.
#'
#' @noRd
is_stratifiable_response <- function(y){

  if(is.factor(y) || is.character(y) || is.logical(y)){
    return(TRUE)
  }

  if(is.numeric(y)){

    observed <- y[!is.na(y)]

    return(length(unique(observed)) == 2L)
  }

  FALSE
}

#' Create Train-test Split
#'
#' @param y Response vector.
#' @param test_fraction Proportion assigned to test set.
#' @param seed Random seed.
#' @param stratify Logical. If TRUE, stratify factor, character, logical, or numeric two-class responses.
#'
#' @return A list with train and test indices.
#'
#' @examples
#' y <- seq_len(20)
#' create_train_test_split(
#'   y, test_fraction = 0.25, seed = 123, stratify = FALSE
#' )
#' @export
create_train_test_split <- function(
    y,
    test_fraction = 0.30,
    seed = 123,
    stratify = TRUE
){

  if(test_fraction <= 0 || test_fraction >= 1){
    stop("test_fraction must lie between 0 and 1.")
  }

  set.seed(seed)

  n <- length(y)

  if(isTRUE(stratify) && is_stratifiable_response(y)){

    y_fac <- as.factor(y)

    test_idx <- unlist(
      lapply(
        split(seq_len(n), y_fac),
        function(idx){

          n_test <- max(
            1,
            floor(length(idx) * test_fraction)
          )

          sample(
            idx,
            size = n_test
          )
        }
      ),
      use.names = FALSE
    )

  } else {

    n_test <- floor(n * test_fraction)

    test_idx <- sample(
      seq_len(n),
      size = n_test
    )
  }

  train_idx <- setdiff(
    seq_len(n),
    test_idx
  )

  list(
    train = train_idx,
    test = test_idx
  )
}


#' Create K-fold Cross-validation Folds
#'
#' @param y Response vector.
#' @param folds Number of folds.
#' @param seed Random seed.
#' @param stratify Logical. If TRUE, stratify factor, character, logical, or numeric two-class responses.
#'
#' @return Integer vector of fold labels.
#'
#' @examples
#' y <- seq_len(12)
#' create_cv_folds(y, folds = 3, seed = 123, stratify = FALSE)
#' @export
create_cv_folds <- function(
    y,
    folds = 5,
    seed = 123,
    stratify = TRUE
){

  if(folds < 2){
    stop("folds must be at least 2.")
  }

  n <- length(y)

  if(folds > n){
    stop("folds cannot exceed length(y).")
  }

  set.seed(seed)

  fold_id <- integer(n)

  if(isTRUE(stratify) && is_stratifiable_response(y)){

    y_fac <- as.factor(y)

    for(cl in levels(y_fac)){

      idx <- which(y_fac == cl)
      idx <- sample(idx)

      fold_id[idx] <- rep(
        seq_len(folds),
        length.out = length(idx)
      )
    }

  } else {

    idx <- sample(seq_len(n))

    fold_id[idx] <- rep(
      seq_len(folds),
      length.out = n
    )
  }

  fold_id
}


#' Evaluate SPFC with Train-test Validation
#'
#' @param X Predictor matrix.
#' @param y Response vector.
#' @param train_idx Training indices.
#' @param test_idx Test indices.
#' @param d Structural dimension.
#' @param ytype Response type.
#' @param cov_method Covariance estimator.
#' @param nslices Number of response slices.
#' @param poly_degree Polynomial degree.
#' @param variable_method Variable selection method.
#' @param selection_rule Character screening rule: `"c1f"`, `"quantile"`, or `"fixed"`.
#' @param threshold Optional finite numeric screening threshold. It is required for `selection_rule = "fixed"`; for `"quantile"`, `NULL` uses `quantile_cut`; it is ignored by the C1F-calibrated rule.
#' @param quantile_cut Variable selection threshold.
#' @param centre_x Logical.
#' @param scale_x Logical.
#' @param ... Additional arguments passed to spfc_fit().
#'
#' @return A list containing fit, reduced model, predictions, and test metrics.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' X <- matrix(rnorm(180), nrow = 30, ncol = 6)
#' y <- X[, 1] - 0.7 * X[, 2] + rnorm(30, sd = 0.5)
#' split <- create_train_test_split(
#'   y, test_fraction = 0.25, seed = 123, stratify = FALSE
#' )
#' ans <- evaluate_spfc_train_test(
#'   X, y, train_idx = split$train, test_idx = split$test,
#'   d = 1, ytype = 'continuous', cov_method = 'oas',
#'   poly_degree = 2, selection_rule = 'quantile',
#'   quantile_cut = 0.75
#' )
#' ans$evaluation
#' }
#' @export
evaluate_spfc_train_test <- function(
    X,
    y,
    train_idx,
    test_idx,
    d,
    ytype = c("auto", "continuous", "categorical"),
    cov_method = c("mec", "oas", "sre", "sde", "cse"),
    nslices = 5,
    poly_degree = 1,
    variable_method = c("adaptive_weighted_l1", "c1f_extension"),
    selection_rule = c("c1f", "quantile", "fixed"),
    threshold = NULL,
    quantile_cut = 0.75,
    centre_x = TRUE,
    scale_x = FALSE,
    ...
){

  ytype <- match.arg(ytype)
  cov_method <- match.arg(cov_method)
  variable_method <- match.arg(variable_method)
  selection_rule <- match.arg(selection_rule)

  X <- as.matrix(X)

  detected_ytype <- detect_y_type(
    y = y,
    ytype = ytype
  )

  X_train <- X[train_idx, , drop = FALSE]
  X_test <- X[test_idx, , drop = FALSE]

  y_train <- y[train_idx]
  y_test <- y[test_idx]

  fit <- spfc_fit(
    X = X_train,
    y = y_train,
    d = d,
    ytype = detected_ytype,
    cov_method = cov_method,
    nslices = nslices,
    poly_degree = poly_degree,
    centre_x = centre_x,
    scale_x = scale_x,
    ...
  )

  Z_train <- fit$Z

  Z_test <- stats::predict(
    fit,
    newdata = X_test
  )

  reduced_model <- fit_reduced_model(
    Z = Z_train,
    y = y_train,
    ytype = detected_ytype
  )

  evaluation <- evaluate_reduced_model(
    object = reduced_model,
    Z = Z_test,
    y = y_test
  )

  variable_selection <- spfc_select_variables(
    fit = fit,
    method = variable_method,
    selection_rule = selection_rule,
    reduced_model = reduced_model,
    threshold = threshold,
    quantile_cut = quantile_cut
  )

  list(
    validation = "holdout",
    train_idx = train_idx,
    test_idx = test_idx,
    fit = fit,
    reduced_model = reduced_model,
    evaluation = evaluation,
    variable_selection = variable_selection,
    y_train = y_train,
    y_test = y_test,
    Z_train = Z_train,
    Z_test = Z_test
  )
}


#' Evaluate SPFC with K-fold Cross-validation
#'
#' @param X Predictor matrix.
#' @param y Response vector.
#' @param d Structural dimension.
#' @param ytype Response type.
#' @param cov_method Covariance estimator.
#' @param folds Number of folds.
#' @param seed Random seed.
#' @param nslices Number of response slices.
#' @param poly_degree Polynomial degree.
#' @param variable_method Variable selection method.
#' @param selection_rule Character screening rule: `"c1f"`, `"quantile"`, or `"fixed"`.
#' @param threshold Optional finite numeric screening threshold. It is required for `selection_rule = "fixed"`; for `"quantile"`, `NULL` uses `quantile_cut`; it is ignored by the C1F-calibrated rule.
#' @param quantile_cut Variable selection threshold.
#' @param centre_x Logical.
#' @param scale_x Logical.
#' @param ... Additional arguments passed to spfc_fit().
#'
#' @return A list containing fold-level and aggregate validation results.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' X <- matrix(rnorm(180), nrow = 30, ncol = 6)
#' y <- X[, 1] - 0.7 * X[, 2] + rnorm(30, sd = 0.5)
#' cv <- evaluate_spfc_cv(
#'   X, y, d = 1, ytype = 'continuous', cov_method = 'oas',
#'   folds = 3, seed = 123, poly_degree = 2,
#'   selection_rule = 'quantile', quantile_cut = 0.75
#' )
#' cv$aggregate
#' }
#' @export
evaluate_spfc_cv <- function(
    X,
    y,
    d,
    ytype = c("auto", "continuous", "categorical"),
    cov_method = c("mec", "oas", "sre", "sde", "cse"),
    folds = 5,
    seed = 123,
    nslices = 5,
    poly_degree = 1,
    variable_method = c("adaptive_weighted_l1", "c1f_extension"),
    selection_rule = c("c1f", "quantile", "fixed"),
    threshold = NULL,
    quantile_cut = 0.75,
    centre_x = TRUE,
    scale_x = FALSE,
    ...
){

  ytype <- match.arg(ytype)
  cov_method <- match.arg(cov_method)
  variable_method <- match.arg(variable_method)
  selection_rule <- match.arg(selection_rule)

  X <- as.matrix(X)

  detected_ytype <- detect_y_type(
    y = y,
    ytype = ytype
  )

  fold_id <- create_cv_folds(
    y = y,
    folds = folds,
    seed = seed,
    stratify = detected_ytype == "categorical"
  )

  fold_results <- vector(
    mode = "list",
    length = folds
  )

  fold_summary <- vector(
    mode = "list",
    length = folds
  )

  for(k in seq_len(folds)){

    test_idx <- which(fold_id == k)
    train_idx <- which(fold_id != k)

    fold_results[[k]] <- evaluate_spfc_train_test(
      X = X,
      y = y,
      train_idx = train_idx,
      test_idx = test_idx,
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

    ev <- fold_results[[k]]$evaluation

    if(detected_ytype == "continuous"){

      fold_summary[[k]] <- data.frame(
        fold = k,
        rmse = ev$rmse,
        mae = ev$mae,
        n_selected = sum(fold_results[[k]]$variable_selection$selected)
      )

    } else {

      fold_summary[[k]] <- data.frame(
        fold = k,
        accuracy = ev$accuracy,
        sensitivity = ev$sensitivity,
        specificity = ev$specificity,
        precision = ev$precision,
        f1 = ev$f1,
        balanced_accuracy = ev$balanced_accuracy,
        n_selected = sum(fold_results[[k]]$variable_selection$selected)
      )
    }
  }

  fold_summary <- do.call(
    rbind,
    fold_summary
  )

  if(detected_ytype == "continuous"){

    aggregate <- data.frame(
      mean_rmse = mean(fold_summary$rmse, na.rm = TRUE),
      sd_rmse = stats::sd(fold_summary$rmse, na.rm = TRUE),
      mean_mae = mean(fold_summary$mae, na.rm = TRUE),
      sd_mae = stats::sd(fold_summary$mae, na.rm = TRUE),
      mean_n_selected = mean(fold_summary$n_selected, na.rm = TRUE)
    )

  } else {

    aggregate <- data.frame(
      mean_accuracy = mean(fold_summary$accuracy, na.rm = TRUE),
      sd_accuracy = stats::sd(fold_summary$accuracy, na.rm = TRUE),
      mean_sensitivity = mean(fold_summary$sensitivity, na.rm = TRUE),
      sd_sensitivity = stats::sd(fold_summary$sensitivity, na.rm = TRUE),
      mean_specificity = mean(fold_summary$specificity, na.rm = TRUE),
      sd_specificity = stats::sd(fold_summary$specificity, na.rm = TRUE),
      mean_precision = mean(fold_summary$precision, na.rm = TRUE),
      sd_precision = stats::sd(fold_summary$precision, na.rm = TRUE),
      mean_f1 = mean(fold_summary$f1, na.rm = TRUE),
      sd_f1 = stats::sd(fold_summary$f1, na.rm = TRUE),
      mean_balanced_accuracy = mean(fold_summary$balanced_accuracy, na.rm = TRUE),
      sd_balanced_accuracy = stats::sd(fold_summary$balanced_accuracy, na.rm = TRUE),
      mean_n_selected = mean(fold_summary$n_selected, na.rm = TRUE)
    )
  }

  full_fit <- spfc_fit(
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

  full_reduced_model <- fit_reduced_model(
    Z = full_fit$Z,
    y = y,
    ytype = detected_ytype
  )

  full_variable_selection <- spfc_select_variables(
    fit = full_fit,
    method = variable_method,
    selection_rule = selection_rule,
    reduced_model = full_reduced_model,
    threshold = threshold,
    quantile_cut = quantile_cut
  )

  overall_confusion <- NULL

  if(detected_ytype == "categorical"){

    cms <- lapply(
      fold_results,
      function(x) x$evaluation$confusion_matrix
    )

    cms <- cms[
      !vapply(cms, is.null, logical(1))
    ]

    if(length(cms) > 0){

      overall_confusion <- Reduce(
        "+",
        cms
      )
    }
  }

  list(
    validation = "cv",
    folds = folds,
    fold_id = fold_id,
    fold_results = fold_results,
    fold_summary = fold_summary,
    aggregate = aggregate,
    confusion_matrix = overall_confusion,
    fit = full_fit,
    reduced_model = full_reduced_model,
    variable_selection = full_variable_selection
  )
}
