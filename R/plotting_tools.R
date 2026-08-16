#' Plot SPFC Simulation Metric by Covariance Method
#'
#' Produces a base R boxplot for a selected simulation metric across covariance
#' estimators.
#'
#' @param results Data frame returned by `run_spfc_design()` or
#' `run_spfc_simulation()`.
#' @param metric Character. Name of the metric column to plot.
#' @param main Optional plot title.
#' @param ylab Optional y-axis label.
#'
#' @return Invisibly returns the formula used for plotting.
#'
#' @examples
#' results <- data.frame(
#'   cov_method = rep(c('oas', 'sre'), each = 2),
#'   rmse = c(0.50, 0.55, 0.62, 0.59),
#'   runtime_sec = c(0.01, 0.012, 0.008, 0.009),
#'   subspace_distance = c(0.20, 0.22, 0.28, 0.25)
#' )
#' plot_metric_by_covariance(results, metric = 'rmse')
#' @export
plot_metric_by_covariance <- function(
    results,
    metric,
    main = NULL,
    ylab = NULL
){

  if(!is.data.frame(results)){
    stop("results must be a data frame.")
  }

  if(!(metric %in% names(results))){
    stop("metric not found in results.")
  }

  if(!("cov_method" %in% names(results))){
    stop("results must contain cov_method.")
  }

  if(is.null(main)){
    main <- paste(metric, "by covariance estimator")
  }

  if(is.null(ylab)){
    ylab <- metric
  }

  form <- stats::as.formula(
    paste(metric, "~ cov_method")
  )

  graphics::boxplot(
    form,
    data = results,
    main = main,
    xlab = "Covariance estimator",
    ylab = ylab
  )

  invisible(form)
}


#' Plot RMSE by Covariance Method
#'
#' @param results Simulation results data frame.
#'
#' @return Invisibly returns the plotting formula.
#'
#' @examples
#' results <- data.frame(
#'   cov_method = rep(c('oas', 'sre'), each = 2),
#'   rmse = c(0.50, 0.55, 0.62, 0.59),
#'   runtime_sec = c(0.01, 0.012, 0.008, 0.009),
#'   subspace_distance = c(0.20, 0.22, 0.28, 0.25)
#' )
#' plot_rmse_by_covariance(results)
#' @export
plot_rmse_by_covariance <- function(results){

  plot_metric_by_covariance(
    results = results,
    metric = "rmse",
    main = "RMSE by covariance estimator",
    ylab = "RMSE"
  )
}


#' Plot Accuracy by Covariance Method
#'
#' @param results Simulation results data frame.
#'
#' @return Invisibly returns the plotting formula.
#'
#' @noRd
plot_accuracy_by_covariance <- function(results){

  plot_metric_by_covariance(
    results = results,
    metric = "accuracy",
    main = "Classification accuracy by covariance estimator",
    ylab = "Accuracy"
  )
}


#' Plot Subspace Distance by Covariance Method
#'
#' @param results Simulation results data frame.
#'
#' @return Invisibly returns the plotting formula.
#'
#' @examples
#' results <- data.frame(
#'   cov_method = rep(c('oas', 'sre'), each = 2),
#'   rmse = c(0.50, 0.55, 0.62, 0.59),
#'   runtime_sec = c(0.01, 0.012, 0.008, 0.009),
#'   subspace_distance = c(0.20, 0.22, 0.28, 0.25)
#' )
#' plot_subspace_distance_by_covariance(results)
#' @export
plot_subspace_distance_by_covariance <- function(results){

  plot_metric_by_covariance(
    results = results,
    metric = "subspace_distance",
    main = "Subspace distance by covariance estimator",
    ylab = "Subspace distance"
  )
}


#' Plot Runtime by Covariance Method
#'
#' @param results Simulation results data frame.
#'
#' @return Invisibly returns the plotting formula.
#'
#' @examples
#' results <- data.frame(
#'   cov_method = rep(c('oas', 'sre'), each = 2),
#'   rmse = c(0.50, 0.55, 0.62, 0.59),
#'   runtime_sec = c(0.01, 0.012, 0.008, 0.009),
#'   subspace_distance = c(0.20, 0.22, 0.28, 0.25)
#' )
#' plot_runtime_by_covariance(results)
#' @export
plot_runtime_by_covariance <- function(results){

  plot_metric_by_covariance(
    results = results,
    metric = "runtime_sec",
    main = "Runtime by covariance estimator",
    ylab = "Runtime in seconds"
  )
}


#' Plot Variable-Selection F1 by Covariance Method
#'
#' @param results Simulation results data frame.
#'
#' @return Invisibly returns the plotting formula.
#'
#' @noRd
plot_variable_f1_by_covariance <- function(results){

  plot_metric_by_covariance(
    results = results,
    metric = "f1_variable",
    main = "Variable-selection F1 by covariance estimator",
    ylab = "Variable-selection F1"
  )
}


#' Plot SPFC Reduced Scores
#'
#' Produces a reduced-score plot from an object returned by `spfc_fit()`.
#'
#' @param fit Object returned by `spfc_fit()`.
#' @param y Optional response vector used for colouring or grouping.
#' @param dims Integer vector of length 1 or 2 indicating SPFC score dimensions.
#' @param main Optional plot title.
#'
#' @return Invisibly returns the plotted data.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(180), nrow = 30, ncol = 6)
#' y <- X[, 1] - 0.7 * X[, 2] + rnorm(30, sd = 0.5)
#' fit <- spfc_fit(X, y, d = 1, ytype = 'continuous',
#'                 cov_method = 'oas', poly_degree = 2)
#' plot_spfc_scores(fit, y = y, dims = 1)
#' @export
plot_spfc_scores <- function(
    fit,
    y = NULL,
    dims = c(1, 2),
    main = NULL
){

  if(is.null(fit$Z)){
    stop("fit must contain Z.")
  }

  Z <- as.matrix(fit$Z)

  if(length(dims) == 1 || ncol(Z) == 1){

    dim1 <- dims[1]

    if(dim1 > ncol(Z)){
      stop("Requested dimension exceeds ncol(fit$Z).")
    }

    z1 <- Z[, dim1]

    if(is.null(main)){
      main <- paste("SPFC score", dim1)
    }

    if(is.null(y)){
      graphics::plot(
        z1,
        pch = 19,
        xlab = "Observation",
        ylab = paste0("SPFC", dim1),
        main = main
      )
    } else {
      graphics::boxplot(
        z1 ~ as.factor(y),
        xlab = "Response group",
        ylab = paste0("SPFC", dim1),
        main = main
      )
    }

    out <- data.frame(
      score = z1,
      y = if(is.null(y)) NA else y
    )

    return(invisible(out))
  }

  if(length(dims) >= 2){

    dim1 <- dims[1]
    dim2 <- dims[2]

    if(any(c(dim1, dim2) > ncol(Z))){
      stop("Requested dimensions exceed ncol(fit$Z).")
    }

    if(is.null(main)){
      main <- paste("SPFC score plot:", dim1, "vs", dim2)
    }

    if(is.null(y)){
      graphics::plot(
        Z[, dim1],
        Z[, dim2],
        pch = 19,
        xlab = paste0("SPFC", dim1),
        ylab = paste0("SPFC", dim2),
        main = main
      )
    } else {
      graphics::plot(
        Z[, dim1],
        Z[, dim2],
        pch = 19,
        col = as.numeric(as.factor(y)),
        xlab = paste0("SPFC", dim1),
        ylab = paste0("SPFC", dim2),
        main = main
      )

      graphics::legend(
        "topright",
        legend = levels(as.factor(y)),
        col = seq_along(levels(as.factor(y))),
        pch = 19,
        bty = "n"
      )
    }

    out <- data.frame(
      score1 = Z[, dim1],
      score2 = Z[, dim2],
      y = if(is.null(y)) NA else y
    )

    return(invisible(out))
  }
}


#' Plot Dimension-Selection Criteria
#'
#' Plots information criteria returned by `spfc_select_dimension()`.
#'
#' @param dsel Object returned by `spfc_select_dimension()`.
#' @param criterion Character. Criterion to plot.
#'
#' @return Invisibly returns the plotted data.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(180), nrow = 30, ncol = 6)
#' y <- X[, 1] - 0.7 * X[, 2] + rnorm(30, sd = 0.5)
#' dsel <- spfc_select_dimension(
#'   X, y, d_grid = 1:2, cov_method = 'oas',
#'   ytype = 'continuous', poly_degree = 2, verbose = FALSE
#' )
#' plot_dimension_selection(dsel, criterion = 'ICOMP_IFIM')
#' @export
plot_dimension_selection <- function(
    dsel,
    criterion = c(
      "AIC",
      "BIC",
      "CAIC",
      "ICOMP_IFIM",
      "ICOMP_MISSPEC",
      "CICOMP"
    )
){

  criterion <- match.arg(criterion)

  if(is.null(dsel$criteria)){
    stop("dsel must contain a criteria table.")
  }

  tab <- dsel$criteria

  if(!(criterion %in% names(tab))){
    stop("criterion not found in dsel$criteria.")
  }

  graphics::plot(
    tab$d,
    tab[[criterion]],
    type = "b",
    pch = 19,
    xlab = "Structural dimension d",
    ylab = criterion,
    main = paste("Dimension selection by", criterion)
  )

  invisible(tab[, c("d", criterion)])
}


#' Plot Shrinkage Intensities
#'
#' Plots shrinkage intensities across covariance estimators from simulation
#' results.
#'
#' @param results Simulation results data frame.
#'
#' @return Invisibly returns the plotting formula.
#'
#' @noRd
plot_shrinkage_rho_by_covariance <- function(results){

  plot_metric_by_covariance(
    results = results,
    metric = "shrinkage_rho",
    main = "Shrinkage intensity by covariance estimator",
    ylab = "Shrinkage intensity"
  )
}


#' Plot Simulation Summary as Barplot
#'
#' Produces a barplot of a mean performance metric from a simulation summary.
#'
#' @param summary Data frame returned by `summarise_spfc_simulation()`.
#' @param metric Character. Name of summary metric column.
#' @param main Optional plot title.
#' @param ylab Optional y-axis label.
#'
#' @return Invisibly returns the plotted data.
#'
#' @noRd
plot_summary_metric <- function(
    summary,
    metric,
    main = NULL,
    ylab = NULL
){

  if(!is.data.frame(summary)){
    stop("summary must be a data frame.")
  }

  if(!(metric %in% names(summary))){
    stop("metric not found in summary.")
  }

  labels <- paste(
    summary$cov_method,
    summary$ytype,
    sep = "_"
  )

  if(is.null(main)){
    main <- paste(metric, "summary")
  }

  if(is.null(ylab)){
    ylab <- metric
  }

  graphics::barplot(
    height = summary[[metric]],
    names.arg = labels,
    las = 2,
    main = main,
    ylab = ylab
  )

  invisible(summary[, c("cov_method", "ytype", metric)])
}
