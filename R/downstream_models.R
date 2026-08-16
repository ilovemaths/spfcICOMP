#' Fit Reduced Downstream Model
#'
#' Fits a downstream predictive model using the reduced SPFC score matrix.
#'
#' @param Z Numeric reduced predictor matrix.
#' @param y Response vector.
#' @param ytype Character. One of `"auto"`, `"continuous"`, or `"categorical"`.
#' @param classifier Character. Classification model. One of `"auto"`, `"glm"`,
#' or `"ridge_logistic"`.
#' @param standardize Logical. If `TRUE`, standardise predictors internally in
#' ridge logistic regression.
#'
#' @return A list containing the fitted downstream model and metadata.
#'
#' @examples
#' set.seed(1)
#' Z <- matrix(rnorm(40), nrow = 20, ncol = 2)
#' y <- Z[, 1] - 0.5 * Z[, 2] + rnorm(20, sd = 0.4)
#' fit_reduced_model(Z, y, ytype = 'continuous')
#' @export
fit_reduced_model <- function(
    Z,
    y,
    ytype = c("auto", "continuous", "categorical"),
    classifier = c("auto", "glm", "ridge_logistic"),
    standardize = TRUE
){

  ytype <- match.arg(ytype)
  classifier <- match.arg(classifier)

  detected_ytype <- detect_y_type(
    y = y,
    ytype = ytype
  )

  Z <- as.matrix(Z)

  if(!is.numeric(Z)){
    stop("Z must be numeric.")
  }

  if(length(y) != nrow(Z)){
    stop("length(y) must equal nrow(Z).")
  }

  if(is.null(colnames(Z))){
    colnames(Z) <- paste0("Z", seq_len(ncol(Z)))
  }

  Zdf <- data.frame(
    y = y,
    Z
  )

  if(detected_ytype == "continuous"){

    fit <- stats::lm(
      y ~ .,
      data = Zdf
    )

    return(
      list(
        model = fit,
        ytype = "continuous",
        method = "lm",
        fitted_values = stats::fitted(fit),
        residuals = stats::residuals(fit)
      )
    )
  }

  yf <- as.factor(y)
  Zdf$y <- yf

  n_classes <- length(levels(yf))

  if(n_classes != 2){
    stop(
      "fit_reduced_model() currently supports binary categorical responses only. ",
      "Multiclass support will be added later."
    )
  }

  if(classifier == "auto"){

    classifier <- "ridge_logistic"

  }

  if(classifier == "glm"){

    fit <- stats::glm(
      y ~ .,
      data = Zdf,
      family = stats::binomial()
    )

    return(
      list(
        model = fit,
        ytype = "categorical",
        method = ifelse(ncol(Z) == 1, "glm_single_dimension", "glm"),
        fitted_probabilities = stats::fitted(fit),
        positive = levels(yf)[2],
        levels = levels(yf)
      )
    )
  }

  if(classifier == "ridge_logistic"){

    model <- fit_ridge_logistic(
      Z = Z,
      y = yf,
      lambda = 1
    )

    return(
      list(
        model = model,
        ytype = "categorical",
        method = "ridge_logistic",
        lambda = 1,
        positive = levels(yf)[2],
        levels = levels(yf)
      )
    )
  }

  stop("Unsupported classifier.")
}
