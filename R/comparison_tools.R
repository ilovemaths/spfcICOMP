#' Compare Covariance Estimators in SPFC
#'
#' Fits shrinkage principal fitted components using several covariance
#' estimators and evaluates the reduced model fitted on the resulting SPFC
#' scores.
#'
#' @param X Numeric predictor matrix.
#' @param y Response vector.
#' @param d Number of reduction directions.
#' @param methods Character vector of covariance estimators. Available options
#' are `"mec"`, `"oas"`, `"sre"`, `"sde"`, and `"cse"`.
#' @param ytype Character. One of `"auto"`, `"continuous"`, or `"categorical"`.
#' @param nslices Number of slices for continuous responses.
#' @param poly_degree Polynomial degree for continuous response basis.
#' @param gamma Ridge constant for SRE.
#' @param rho Convex shrinkage intensity for CSE.
#' @param centre_x Logical. If `TRUE`, centre predictor columns.
#' @param scale_x Logical. If `TRUE`, scale predictor columns.
#' @param classifier Character. One of `"auto"`, `"glm"`, or `"ridge_logistic"`.
#' @param eps Positive numerical floor for eigenvalues.
#' @param verbose Logical. If `TRUE`, print progress.
#'
#' @return A list containing a summary table and full fitted objects.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(180), nrow = 30, ncol = 6)
#' y <- X[, 1] - 0.7 * X[, 2] + rnorm(30, sd = 0.5)
#' cmp <- compare_covariance_methods(
#'   X, y, d = 1, methods = c('oas', 'sre'),
#'   ytype = 'continuous', poly_degree = 2, verbose = FALSE
#' )
#' cmp$summary
#' @export
compare_covariance_methods <- function(
    X,
    y,
    d,
    methods = c("mec", "oas", "sre", "sde", "cse"),
    ytype = c("auto", "continuous", "categorical"),
    nslices = 5,
    poly_degree = 1,
    gamma = 0.1,
    rho = 0.5,
    centre_x = TRUE,
    scale_x = FALSE,
    classifier = c("auto", "glm", "ridge_logistic"),
    eps = 1e-8,
    verbose = TRUE
){

  ytype <- match.arg(ytype)
  classifier <- match.arg(classifier)

  methods <- match.arg(
    methods,
    choices = c("mec", "oas", "sre", "sde", "cse"),
    several.ok = TRUE
  )

  detected_ytype <- detect_y_type(
    y = y,
    ytype = ytype
  )

  full_results <- vector("list", length(methods))
  names(full_results) <- methods

  summary_rows <- vector("list", length(methods))

  for(i in seq_along(methods)){

    method_i <- methods[i]

    if(verbose){
      message("Running SPFC with covariance method: ", method_i)
    }

    t1 <- Sys.time()

    fit_i <- spfc_fit(
      X = X,
      y = y,
      d = d,
      fy = NULL,
      ytype = detected_ytype,
      cov_method = method_i,
      nslices = nslices,
      poly_degree = poly_degree,
      gamma = gamma,
      rho = rho,
      centre_x = centre_x,
      scale_x = scale_x,
      eps = eps
    )

    t2 <- Sys.time()

    runtime_sec <- as.numeric(
      difftime(t2, t1, units = "secs")
    )

    reduced_i <- fit_reduced_model(
      Z = fit_i$Z,
      y = y,
      ytype = detected_ytype,
      classifier = classifier
    )

    eval_i <- evaluate_reduced_model(
      object = reduced_i,
      Z = fit_i$Z,
      y = y
    )

    if(detected_ytype == "continuous"){

      summary_rows[[i]] <- data.frame(
        cov_method = method_i,
        ytype = detected_ytype,
        d = d,
        rho = fit_i$rho,
        nslices_used = fit_i$nslices_used,
        runtime_sec = runtime_sec,
        runtime_min = runtime_sec / 60,
        reduced_model = reduced_i$method,
        rmse = eval_i$rmse,
        mae = eval_i$mae,
        accuracy = NA_real_,
        sensitivity = NA_real_,
        specificity = NA_real_,
        f1 = NA_real_,
        balanced_accuracy = NA_real_
      )

    } else {

      summary_rows[[i]] <- data.frame(
        cov_method = method_i,
        ytype = detected_ytype,
        d = d,
        rho = fit_i$rho,
        nslices_used = fit_i$nslices_used,
        runtime_sec = runtime_sec,
        runtime_min = runtime_sec / 60,
        reduced_model = reduced_i$method,
        rmse = NA_real_,
        mae = NA_real_,
        accuracy = eval_i$accuracy,
        sensitivity = eval_i$sensitivity,
        specificity = eval_i$specificity,
        f1 = eval_i$f1,
        balanced_accuracy = eval_i$balanced_accuracy
      )
    }

    full_results[[method_i]] <- list(
      spfc_fit = fit_i,
      reduced_model = reduced_i,
      evaluation = eval_i
    )
  }

  summary_table <- do.call(
    rbind,
    summary_rows
  )

  rownames(summary_table) <- NULL

  return(list(
    summary = summary_table,
    results = full_results
  ))
}
