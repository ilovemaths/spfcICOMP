#' Compute Bozdogan's C1 Complexity
#'
#' Computes the maximal covariance complexity measure C1 from a positive
#' definite covariance-type matrix.
#'
#' @param Sigma Numeric positive definite covariance-type matrix.
#' @param eps Small positive numerical floor.
#'
#' @return Numeric C1 value.
#'
#' @noRd
compute_c1 <- function(Sigma, eps = 1e-10){

  Sigma <- make_positive_definite(Sigma, eps = eps)

  evals <- eigen(
    Sigma,
    symmetric = TRUE,
    only.values = TRUE
  )$values

  evals <- evals[is.finite(evals) & evals > eps]

  s <- length(evals)

  if(s == 0){
    stop("No positive eigenvalues available for C1 computation.")
  }

  lambda_bar <- mean(evals)

  c1 <- (s / 2) * log(lambda_bar) -
    0.5 * sum(log(evals))

  return(as.numeric(c1))
}


#' Compute Bozdogan's Scale-Invariant C1F Complexity
#'
#' Computes the scale-invariant C1F covariance complexity measure from a
#' positive definite covariance-type matrix.
#'
#' @param Sigma Numeric positive definite covariance-type matrix.
#' @param eps Small positive numerical floor.
#'
#' @return Numeric C1F value.
#'
#' @noRd
compute_c1f <- function(Sigma, eps = 1e-10){

  Sigma <- make_positive_definite(Sigma, eps = eps)

  evals <- eigen(
    Sigma,
    symmetric = TRUE,
    only.values = TRUE
  )$values

  evals <- evals[is.finite(evals) & evals > eps]

  if(length(evals) == 0){
    stop("No positive eigenvalues available for C1F computation.")
  }

  lambda_bar <- mean(evals)

  c1f <- sum((evals - lambda_bar)^2) /
    (4 * lambda_bar^2)

  return(as.numeric(c1f))
}


#' Extract Model Log-Likelihood
#'
#' @param model A fitted model object.
#'
#' @return Numeric log-likelihood value.
#'
#' @noRd
extract_loglik <- function(model){

  ll <- stats::logLik(model)

  return(as.numeric(ll))
}


#' Extract Number of Effective Parameters
#'
#' @param model A fitted model object.
#'
#' @return Integer number of estimated parameters.
#'
#' @noRd
extract_npar <- function(model){

  ll <- stats::logLik(model)

  k <- attr(ll, "df")

  if(is.null(k)){
    k <- length(stats::coef(model))
  }

  return(as.integer(k))
}


#' Compute AIC
#'
#' @param loglik Numeric maximised log-likelihood.
#' @param k Number of estimated parameters.
#'
#' @return Numeric AIC value.
#'
#' @noRd
compute_aic <- function(loglik, k){
  return(as.numeric(-2 * loglik + 2 * k))
}


#' Compute BIC
#'
#' @param loglik Numeric maximised log-likelihood.
#' @param k Number of estimated parameters.
#' @param n Sample size.
#'
#' @return Numeric BIC value.
#'
#' @noRd
compute_bic <- function(loglik, k, n){
  return(as.numeric(-2 * loglik + k * log(n)))
}


#' Compute CAIC
#'
#' @param loglik Numeric maximised log-likelihood.
#' @param k Number of estimated parameters.
#' @param n Sample size.
#'
#' @return Numeric CAIC value.
#'
#' @noRd
compute_caic <- function(loglik, k, n){
  return(as.numeric(-2 * loglik + k * (log(n) + 1)))
}


#' Compute ICOMP Based on the Inverse Fisher Information Matrix
#'
#' Computes ICOMP(IFIM) using the covariance matrix of the estimated model
#' parameters as a practical inverse Fisher information matrix.
#'
#' @param loglik Numeric maximised log-likelihood.
#' @param ifim Numeric inverse Fisher information matrix.
#' @param eps Small positive numerical floor.
#'
#' @return Numeric ICOMP(IFIM) value.
#'
#' @noRd
compute_icomp_ifim <- function(loglik, ifim, eps = 1e-10){

  ifim <- make_positive_definite(ifim, eps = eps)

  return(as.numeric(
    -2 * loglik + 2 * compute_c1(ifim, eps = eps)
  ))
}


#' Compute Misspecification-Resistant ICOMP
#'
#' Computes ICOMP(Misspec) from a sandwich covariance matrix.
#'
#' @param loglik Numeric maximised log-likelihood.
#' @param sandwich_cov Numeric sandwich covariance matrix.
#' @param eps Small positive numerical floor.
#'
#' @return Numeric ICOMP(Misspec) value.
#'
#' @noRd
compute_icomp_misspec <- function(loglik, sandwich_cov, eps = 1e-10){

  sandwich_cov <- make_positive_definite(sandwich_cov, eps = eps)

  return(as.numeric(
    -2 * loglik + 2 * compute_c1(sandwich_cov, eps = eps)
  ))
}


#' Compute Consistent ICOMP
#'
#' Computes CICOMP using C1F covariance complexity.
#'
#' @param loglik Numeric maximised log-likelihood.
#' @param k Number of estimated parameters.
#' @param n Sample size.
#' @param complexity_matrix Numeric covariance-type matrix used for C1F.
#' @param eps Small positive numerical floor.
#'
#' @return Numeric CICOMP value.
#'
#' @noRd
compute_cicomp <- function(
    loglik,
    k,
    n,
    complexity_matrix,
    eps = 1e-10
){

  complexity_matrix <- make_positive_definite(
    complexity_matrix,
    eps = eps
  )

  return(as.numeric(
    -2 * loglik +
      k +
      k * log(n) +
      2 * compute_c1f(complexity_matrix, eps = eps)
  ))
}


#' Extract IFIM Approximation from a Reduced Model
#'
#' Extracts a practical inverse Fisher information matrix approximation from
#' fitted reduced models. For `lm` and ordinary `glm`, this uses `vcov()`.
#' For the package's ridge-logistic fit, it uses the inverse penalised observed
#' Hessian stored at estimation time.
#'
#' @param reduced_model Object returned by `fit_reduced_model()`.
#' @param eps Small positive numerical floor.
#'
#' @return Numeric covariance-type matrix.
#'
#' @noRd
extract_ifim <- function(reduced_model, eps = 1e-10){

  if(reduced_model$method %in% c("lm", "glm", "glm_single_dimension")){

    V <- stats::vcov(reduced_model$model)

    V <- make_positive_definite(V, eps = eps)

    return(V)
  }

  if(identical(reduced_model$method, "ridge_logistic")){

    model <- reduced_model$model

    if(!is.null(model$vcov)){
      V <- as.matrix(model$vcov)
      V <- make_positive_definite(V, eps = eps)
      return(V)
    }

    if(!is.null(model$hessian)){
      H <- make_positive_definite(as.matrix(model$hessian), eps = eps)
      V <- solve(H)
      V <- make_positive_definite(V, eps = eps)
      return(V)
    }

    stop(
      "ridge_logistic reduced model does not contain a stored Hessian or ",
      "covariance approximation."
    )
  }

  stop(
    "IFIM extraction is not implemented for reduced-model method: ",
    reduced_model$method
  )
}


#' Compute Sandwich Covariance for Ordinary GLM or LM
#'
#' Computes a simple empirical sandwich covariance matrix using casewise
#' estimating scores. This is implemented for ordinary `lm` and `glm` objects.
#'
#' @param model A fitted `lm` or `glm` object.
#' @param eps Small positive numerical floor.
#'
#' @return Numeric sandwich covariance matrix.
#'
#' @noRd
compute_sandwich_cov <- function(model, eps = 1e-10){

  X <- stats::model.matrix(model)

  n <- nrow(X)

  if(inherits(model, "lm") && !inherits(model, "glm")){

    residuals <- stats::residuals(model)

    sigma2 <- sum(residuals^2) / n

    bread <- sigma2 * solve(crossprod(X))

    scores <- sweep(
      X,
      MARGIN = 1,
      STATS = residuals,
      FUN = "*"
    )

    meat <- crossprod(scores) / n

    sandwich <- bread %*% meat %*% bread

    sandwich <- make_positive_definite(sandwich, eps = eps)

    return(sandwich)
  }

  if(inherits(model, "glm")){

    family_name <- model$family$family

    if(family_name != "binomial"){
      stop("Sandwich covariance is currently implemented for binomial glm only.")
    }

    y <- model$y
    mu <- stats::fitted(model)

    scores <- sweep(
      X,
      MARGIN = 1,
      STATS = y - mu,
      FUN = "*"
    )

    W <- as.numeric(mu * (1 - mu))
    W <- pmax(W, eps)

    Fmat <- crossprod(X, X * W) / n

    Rmat <- crossprod(scores) / n

    Fmat <- make_positive_definite(Fmat, eps = eps)

    sandwich <- solve(Fmat) %*% Rmat %*% solve(Fmat) / n

    sandwich <- make_positive_definite(sandwich, eps = eps)

    return(sandwich)
  }

  stop("compute_sandwich_cov() supports ordinary lm and glm objects only.")
}


#' Score Information Criteria for a Reduced Model
#'
#' Computes AIC, BIC, CAIC, ICOMP(IFIM), ICOMP(Misspec), and CICOMP for a
#' model fitted on reduced SPFC scores.
#'
#' @param reduced_model Object returned by `fit_reduced_model()`.
#' @param n Sample size.
#' @param eps Small positive numerical floor.
#'
#' @return A one-row data frame of information criteria.
#'
#' @examples
#' set.seed(1)
#' Z <- matrix(rnorm(40), nrow = 20, ncol = 2)
#' y <- Z[, 1] - 0.5 * Z[, 2] + rnorm(20, sd = 0.4)
#' mod <- fit_reduced_model(Z, y, ytype = 'continuous')
#' score_information_criteria(mod, n = length(y))
#' @export
score_information_criteria <- function(
    reduced_model,
    n,
    eps = 1e-10
){

  if(reduced_model$method == "ridge_logistic"){
    stop(
      "Information criteria based on IFIM and sandwich covariance are not ",
      "currently computed for ridge_logistic objects. Use an ordinary glm ",
      "companion model on the reduced scores for dimension selection."
      )
  }

  model <- reduced_model$model

  loglik <- extract_loglik(model)
  k <- extract_npar(model)

  ifim <- extract_ifim(
    reduced_model = reduced_model,
    eps = eps
  )

  sandwich_cov <- compute_sandwich_cov(
    model = model,
    eps = eps
  )

  out <- data.frame(
    loglik = loglik,
    npar = k,
    AIC = compute_aic(loglik = loglik, k = k),
    BIC = compute_bic(loglik = loglik, k = k, n = n),
    CAIC = compute_caic(loglik = loglik, k = k, n = n),
    ICOMP_IFIM = compute_icomp_ifim(
      loglik = loglik,
      ifim = ifim,
      eps = eps
    ),
    ICOMP_MISSPEC = compute_icomp_misspec(
      loglik = loglik,
      sandwich_cov = sandwich_cov,
      eps = eps
    ),
    CICOMP = compute_cicomp(
      loglik = loglik,
      k = k,
      n = n,
      complexity_matrix = ifim,
      eps = eps
    )
  )

  return(out)
}


#' Select Structural Dimension for SPFC
#'
#' Fits SPFC models across candidate dimensions and scores classical and
#' information-complexity-based criteria.
#'
#' @param X Numeric predictor matrix.
#' @param y Response vector.
#' @param d_grid Integer vector of candidate structural dimensions.
#' @param cov_method Covariance estimator. One of `"mec"`, `"oas"`, `"sre"`,
#' `"sde"`, or `"cse"`.
#' @param ytype Character. One of `"auto"`, `"continuous"`, or `"categorical"`.
#' @param nslices Number of slices for MEC with continuous responses.
#' @param poly_degree Polynomial degree for continuous response basis.
#' @param gamma Ridge constant for SRE.
#' @param rho Convex shrinkage intensity for CSE.
#' @param centre_x Logical. If `TRUE`, centre predictor columns.
#' @param scale_x Logical. If `TRUE`, scale predictor columns.
#' @param classifier Character. One of `"auto"` or `"glm"` for categorical
#' responses when information criteria are required. Dimension selection uses
#' an ordinary logistic companion model for categorical responses because the
#' current information criteria require log-likelihood and covariance estimates.
#' @param eps Small positive numerical floor.
#' @param verbose Logical. If `TRUE`, print progress.
#'
#' @return A list containing the criterion table and selected dimensions by
#' criterion.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(180), nrow = 30, ncol = 6)
#' y <- X[, 1] - 0.7 * X[, 2] + rnorm(30, sd = 0.5)
#' dsel <- spfc_select_dimension(
#'   X, y, d_grid = 1:2, cov_method = 'oas',
#'   ytype = 'continuous', poly_degree = 2, verbose = FALSE
#' )
#' dsel$selected
#' @export
spfc_select_dimension <- function(
    X,
    y,
    d_grid = 1:3,
    cov_method = c("mec", "oas", "sre", "sde", "cse"),
    ytype = c("auto", "continuous", "categorical"),
    nslices = 5,
    poly_degree = 1,
    gamma = 0.1,
    rho = 0.5,
    centre_x = TRUE,
    scale_x = FALSE,
    classifier = c("auto", "glm"),
    eps = 1e-8,
    verbose = TRUE
){

  cov_method <- match.arg(cov_method)
  ytype <- match.arg(ytype)
  classifier <- match.arg(classifier)

  detected_ytype <- detect_y_type(
    y = y,
    ytype = ytype
  )

  if(detected_ytype == "categorical"){
    classifier <- "glm"
  }

  rows <- vector("list", length(d_grid))
  fits <- vector("list", length(d_grid))
  names(fits) <- paste0("d", d_grid)

  n <- nrow(as.matrix(X))

  for(i in seq_along(d_grid)){

    d_i <- d_grid[i]

    if(verbose){
      message("Scoring d = ", d_i, " using covariance method: ", cov_method)
    }

    fit_i <- spfc_fit(
      X = X,
      y = y,
      d = d_i,
      fy = NULL,
      ytype = detected_ytype,
      cov_method = cov_method,
      nslices = nslices,
      poly_degree = poly_degree,
      gamma = gamma,
      rho = rho,
      centre_x = centre_x,
      scale_x = scale_x,
      eps = eps
    )

    reduced_i <- suppressWarnings(
      fit_reduced_model(
        Z = fit_i$Z,
        y = y,
        ytype = detected_ytype,
        classifier = classifier
      )
    )

    crit_i <- score_information_criteria(
      reduced_model = reduced_i,
      n = n,
      eps = eps
    )

    rows[[i]] <- data.frame(
      cov_method = cov_method,
      ytype = detected_ytype,
      d = d_i,
      rho = fit_i$rho,
      nslices_used = fit_i$nslices_used,
      reduced_model = reduced_i$method,
      crit_i
    )

    fits[[i]] <- list(
      spfc_fit = fit_i,
      reduced_model = reduced_i,
      criteria = crit_i
    )
  }

  criteria_table <- do.call(rbind, rows)
  rownames(criteria_table) <- NULL

  criterion_names <- c(
    "AIC",
    "BIC",
    "CAIC",
    "ICOMP_IFIM",
    "ICOMP_MISSPEC",
    "CICOMP"
  )

  selected <- data.frame(
    criterion = criterion_names,
    selected_d = NA_integer_,
    minimum_value = NA_real_
  )

  for(j in seq_along(criterion_names)){

    crit <- criterion_names[j]

    idx <- which.min(criteria_table[[crit]])

    selected$selected_d[j] <- criteria_table$d[idx]
    selected$minimum_value[j] <- criteria_table[[crit]][idx]
  }

  return(list(
    criteria = criteria_table,
    selected = selected,
    fits = fits
  ))
}
