#' Fit Ridge Logistic Regression
#'
#' Fits a ridge-penalised logistic regression model for binary classification.
#' This is useful when ordinary logistic regression suffers from complete or
#' quasi-complete separation.
#'
#' @param Z Numeric reduced predictor matrix.
#' @param y Binary response vector.
#' @param lambda Ridge penalty parameter.
#'
#' @return A list of class `ridge_logistic_fit`.
#'
#' @noRd
fit_ridge_logistic <- function(
    Z,
    y,
    lambda = 1
){

  Z <- as.matrix(Z)
  y <- as.factor(y)

  if (nlevels(y) != 2) {
    stop("fit_ridge_logistic() currently supports binary classification only.")
  }

  y_num <- as.numeric(y == levels(y)[2])

  X_design <- cbind(
    intercept = 1,
    Z
  )

  p <- ncol(X_design)

  objective <- function(beta){

    eta <- as.vector(
      X_design %*% beta
    )

    eta <- pmin(
      pmax(eta, -30),
      30
    )

    prob <- 1 / (1 + exp(-eta))

    nll <- -sum(
      y_num * log(prob + 1e-12) +
        (1 - y_num) * log(1 - prob + 1e-12)
    )

    penalty <- lambda * sum(
      beta[-1]^2
    )

    nll + penalty
  }

  start <- rep(
    0,
    p
  )

  opt <- stats::optim(
    par = start,
    fn = objective,
    method = "BFGS",
    hessian = TRUE,
    control = list(
      maxit = 1000
    )
  )

  # The inverse penalised observed Hessian provides a stable covariance-type
  # approximation for the ridge-logistic estimator. It is retained explicitly
  # because the C1F variable-selection route requires the covariance structure
  # of the fitted reduced-model parameters.
  hessian <- (opt$hessian + t(opt$hessian)) / 2
  hessian <- make_positive_definite(hessian, eps = 1e-10)
  vcov_penalised <- solve(hessian)
  vcov_penalised <- make_positive_definite(vcov_penalised, eps = 1e-10)

  out <- list(
    coefficients = opt$par,
    lambda = lambda,
    levels = levels(y),
    convergence = opt$convergence,
    value = opt$value,
    hessian = hessian,
    vcov = vcov_penalised,
    Z_names = colnames(Z)
  )

  class(out) <- c(
    "ridge_logistic_fit",
    "list"
  )

  out
}


#' Predict from Ridge Logistic Regression
#'
#' @param object Object returned by `fit_ridge_logistic()`.
#' @param newdata Numeric reduced predictor matrix.
#' @param type Character. One of `"class"` or `"response"`.
#' @param ... Additional arguments.
#'
#' @return Predicted classes or probabilities.
#'
#' @export
predict.ridge_logistic_fit <- function(
    object,
    newdata,
    type = c("class", "response"),
    ...
){

  type <- match.arg(type)

  Z <- as.matrix(newdata)

  X_design <- cbind(
    intercept = 1,
    Z
  )

  eta <- as.vector(
    X_design %*% object$coefficients
  )

  eta <- pmin(
    pmax(eta, -30),
    30
  )

  prob <- 1 / (1 + exp(-eta))

  if (type == "response") {
    return(prob)
  }

  pred <- ifelse(
    prob >= 0.5,
    object$levels[2],
    object$levels[1]
  )

  factor(
    pred,
    levels = object$levels
  )
}
