#' Print an SPFC Fit
#'
#' @param x Object returned by `spfc_fit()`.
#' @param ... Additional arguments.
#'
#' @return Invisibly returns `x`.
#'
#' @export
print.spfc_fit <- function(x, ...){

  cat("\nShrinkage Principal Fitted Components Fit\n")
  cat("=========================================\n\n")

  cat("Response type:       ", x$ytype, "\n", sep = "")
  cat("Covariance method:   ", x$cov_method, "\n", sep = "")
  cat("Structural dimension:", x$d, "\n")

  if(!is.null(x$rho)){
    cat("Shrinkage rho:       ", round(x$rho, 6), "\n", sep = "")
  }

  if(!is.null(x$nslices_used)){
    cat("Slices used:         ", x$nslices_used, "\n", sep = "")
  }

  if(!is.null(x$eigvals)){
    eig_print <- x$eigvals
    eig_print[abs(eig_print) < 1e-10] <- 0
    eig_print <- pmax(eig_print, 0)

    cat("\nLeading eigenvalues:\n")
    print(round(utils::head(eig_print, 6), 6))
  }

  cat("\nReduced score matrix dimension:\n")
  print(dim(x$Z))

  invisible(x)
}


#' Summarise an SPFC Fit
#'
#' @param object Object returned by `spfc_fit()`.
#' @param ... Additional arguments.
#'
#' @return A list containing key SPFC summary components.
#'
#' @export
summary.spfc_fit <- function(object, ...){

  eigvals <- object$eigvals

  eig_summary <- NULL

  if(!is.null(eigvals)){

    eig_use <- eigvals[
      is.finite(eigvals)
    ]

    eig_use[abs(eig_use) < 1e-10] <- 0
    eig_use <- pmax(eig_use, 0)

    eig_sum <- sum(eig_use)

    if(eig_sum > 0){
      proportion <- eig_use / eig_sum
      cumulative_proportion <- cumsum(eig_use) / eig_sum
    } else {
      proportion <- rep(NA_real_, length(eig_use))
      cumulative_proportion <- rep(NA_real_, length(eig_use))
    }

    eig_summary <- data.frame(
      component = seq_along(eig_use),
      eigenvalue = eig_use,
      proportion = proportion,
      cumulative_proportion = cumulative_proportion
    )
  }

  out <- list(
    ytype = object$ytype,
    cov_method = object$cov_method,
    d = object$d,
    rho = object$rho,
    nslices_used = object$nslices_used,
    score_dimension = dim(object$Z),
    eigen_summary = eig_summary,
    preprocessing = object$preprocessing
  )

  class(out) <- c("summary.spfc_fit", "list")

  return(out)
}


#' Print Summary of an SPFC Fit
#'
#' @param x Object returned by `summary.spfc_fit()`.
#' @param ... Additional arguments.
#'
#' @return Invisibly returns `x`.
#'
#' @export
print.summary.spfc_fit <- function(x, ...){

  cat("\nSummary of Shrinkage Principal Fitted Components Fit\n")
  cat("====================================================\n\n")

  cat("Response type:       ", x$ytype, "\n", sep = "")
  cat("Covariance method:   ", x$cov_method, "\n", sep = "")
  cat("Structural dimension:", x$d, "\n")

  if(!is.null(x$rho)){
    cat("Shrinkage rho:       ", round(x$rho, 6), "\n", sep = "")
  }

  if(!is.null(x$nslices_used)){
    cat("Slices used:         ", x$nslices_used, "\n", sep = "")
  }

  if(!is.null(x$preprocessing)){
    cat("\nPreprocessing:\n")
    cat("Centred:             ", x$preprocessing$centred, "\n", sep = "")
    cat("Scaled:              ", x$preprocessing$scaled, "\n", sep = "")
  }

  cat("\nScore matrix dimension:\n")
  print(x$score_dimension)

  if(!is.null(x$eigen_summary)){
    cat("\nEigenvalue summary:\n")
    print(utils::head(x$eigen_summary, 10))
  }

  invisible(x)
}


#' Extract SPFC Coefficients
#'
#' Extracts the estimated SPFC direction matrix.
#'
#' @param object Object returned by `spfc_fit()`.
#' @param ... Additional arguments.
#'
#' @return Numeric matrix of estimated SPFC directions.
#'
#' @export
coef.spfc_fit <- function(object, ...){

  if(is.null(object$V)){
    stop("The object does not contain V, the SPFC direction matrix.")
  }

  return(object$V)
}


#' Fitted SPFC Scores
#'
#' Extracts the reduced SPFC score matrix.
#'
#' @param object Object returned by `spfc_fit()`.
#' @param ... Additional arguments.
#'
#' @return Numeric matrix of reduced SPFC scores.
#'
#' @export
fitted.spfc_fit <- function(object, ...){

  if(is.null(object$Z)){
    stop("The object does not contain Z, the SPFC score matrix.")
  }

  return(object$Z)
}


#' Predict SPFC Scores for New Data
#'
#' Projects new predictor observations into the estimated SPFC reduced space.
#'
#' @param object Object returned by `spfc_fit()`.
#' @param newdata Numeric matrix or data frame of new predictor observations.
#' @param ... Additional arguments.
#'
#' @return Numeric matrix of predicted SPFC scores.
#'
#' @export
predict.spfc_fit <- function(object, newdata, ...){

  if(missing(newdata)){
    stop("newdata must be supplied.")
  }

  if(is.null(object$V)){
    stop("The fitted object does not contain V, the SPFC direction matrix.")
  }

  Xnew <- as.matrix(newdata)

  if(ncol(Xnew) != nrow(object$V)){
    stop(
      "newdata has ",
      ncol(Xnew),
      " columns, but the fitted SPFC object expects ",
      nrow(object$V),
      " columns."
    )
  }

  if(is.null(object$preprocessing)){
    warning(
      "The fitted object does not contain preprocessing metadata. ",
      "Proceeding without centring or scaling newdata."
    )

    Znew <- Xnew %*% object$V

    colnames(Znew) <- paste0(
      "SPFC",
      seq_len(ncol(Znew))
    )

    return(Znew)
  }

  if(isTRUE(object$preprocessing$centred)){
    Xnew <- sweep(
      Xnew,
      2,
      object$preprocessing$means,
      FUN = "-"
    )
  }

  if(isTRUE(object$preprocessing$scaled)){
    Xnew <- sweep(
      Xnew,
      2,
      object$preprocessing$sds,
      FUN = "/"
    )
  }

  Znew <- Xnew %*% object$V

  colnames(Znew) <- paste0(
    "SPFC",
    seq_len(ncol(Znew))
  )

  return(Znew)
}
