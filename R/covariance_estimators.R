#' Oracle Approximating Shrinkage-Type Covariance Estimator
#'
#' Computes an OAS-type stabilised covariance estimator with an additional
#' convex shrinkage step.
#'
#' @param x Numeric predictor matrix.
#'
#' @return A list with covariance matrix `covx`, shrinkage intensity `rho`,
#' and method label.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(120), nrow = 20, ncol = 6)
#' est <- oas(X)
#' c(method = est$method, rho = est$rho)
#' @export
oas <- function(x){

  X <- standardise_X(x, centre = TRUE, scale = FALSE)

  n <- nrow(X)
  p <- ncol(X)

  mle_cov <- safe_cov(X, centre = FALSE)

  eig <- eigen(mle_cov, symmetric = TRUE)

  vals <- eig$values
  vecs <- eig$vectors

  lambdastar <- diag(
    pmax(vals, mean(vals)),
    nrow = p
  )

  sta <- vecs %*% lambdastar %*% t(vecs)
  sta <- (sta + t(sta)) / 2

  tr_sta <- sum(diag(sta))
  tr_sta2 <- sum(diag(sta %*% sta))

  denom <- (p + 1 - (2 / p)) * (tr_sta2 + (tr_sta^2 / p))

  if(abs(denom) < .Machine$double.eps){
    rho_oas <- 0
  } else {
    rho_oas <- ((1 - 2 / p) * tr_sta2 + tr_sta^2) / denom
    rho_oas <- min(max(rho_oas, 0), 1)
  }

  target <- tr_sta * diag(p) / p

  oas_sta <- (1 - rho_oas) * sta + rho_oas * target
  oas_sta <- (oas_sta + t(oas_sta)) / 2

  lambda_target <- mean(diag(oas_sta)) * diag(p)

  beta <- (sum(diag(oas_sta)))^2 / sum(diag(oas_sta %*% oas_sta))

  if(abs(p - beta) < .Machine$double.eps){
    alpha <- 0
  } else {
    alpha <- 2 * (p * (1 + beta) - 2) / (p - beta)
  }

  mm <- 0.5 * alpha
  rho <- n / (n + mm)
  rho <- min(max(rho, 0), 1)

  covx <- rho * oas_sta + (1 - rho) * lambda_target
  covx <- make_positive_definite(covx)

  return(list(
    covx = covx,
    rho = rho,
    method = "oas",
    rho_oas = rho_oas
  ))
}


#' Stipulated Ridge Covariance Estimator
#'
#' Computes a ridge-stabilised covariance estimator by adding a positive
#' constant to the diagonal of the sample covariance matrix.
#'
#' @param x Numeric predictor matrix.
#' @param gamma Ridge constant.
#'
#' @return A list with covariance matrix `covx`, ridge intensity `rho`,
#' and method label.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(120), nrow = 20, ncol = 6)
#' est <- sre(X, gamma = 0.1)
#' c(method = est$method, rho = est$rho)
#' @export
sre <- function(x, gamma = 0.1){

  X <- standardise_X(x, centre = TRUE, scale = FALSE)

  p <- ncol(X)

  S <- safe_cov(X, centre = FALSE)

  covx <- S + gamma * diag(p)
  covx <- make_positive_definite(covx)

  return(list(
    covx = covx,
    rho = gamma,
    method = "sre"
  ))
}


#' Stipulated Diagonal Covariance Estimator
#'
#' Computes a diagonal covariance estimator by retaining only the diagonal
#' elements of the sample covariance matrix.
#'
#' @param x Numeric predictor matrix.
#'
#' @return A list with covariance matrix `covx`, shrinkage intensity `rho`,
#' and method label.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(120), nrow = 20, ncol = 6)
#' est <- sde(X)
#' c(method = est$method, rho = est$rho)
#' @export
sde <- function(x){

  X <- standardise_X(x, centre = TRUE, scale = FALSE)

  S <- safe_cov(X, centre = FALSE)

  covx <- diag(diag(S))
  covx <- make_positive_definite(covx)

  return(list(
    covx = covx,
    rho = NA_real_,
    method = "sde"
  ))
}


#' Convex Sum Covariance Estimator
#'
#' Computes a convex combination of the sample covariance matrix and a spherical
#' target matrix.
#'
#' @param x Numeric predictor matrix.
#' @param rho Convex shrinkage intensity between 0 and 1.
#'
#' @return A list with covariance matrix `covx`, shrinkage intensity `rho`,
#' and method label.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(120), nrow = 20, ncol = 6)
#' est <- cse(X, rho = 0.5)
#' dim(est$covx)
#' @export
cse <- function(x, rho = 0.5){

  if(!is.numeric(rho) || length(rho) != 1){
    stop("rho must be a single numeric value.")
  }

  if(rho < 0 || rho > 1){
    stop("rho must lie between 0 and 1.")
  }

  X <- standardise_X(x, centre = TRUE, scale = FALSE)

  p <- ncol(X)

  S <- safe_cov(X, centre = FALSE)

  target <- mean(diag(S)) * diag(p)

  covx <- rho * S + (1 - rho) * target
  covx <- make_positive_definite(covx)

  return(list(
    covx = covx,
    rho = rho,
    method = "cse"
  ))
}


#' Maximum Entropy Covariance Estimator
#'
#' Computes the Maximum Entropy Covariance estimator for sufficient dimension
#' reduction. The response is first sliced. Slice-specific covariance matrices
#' are computed, the entropy-maximising slice covariance is selected, and an
#' entropy-guided covariance estimate is stabilised by a convex shrinkage step.
#'
#' @param x Numeric predictor matrix.
#' @param y Response vector.
#' @param ytype Character. One of `"auto"`, `"continuous"`, or `"categorical"`.
#' @param nslices Number of slices for continuous responses.
#' @param eps Positive numerical floor for eigenvalues.
#'
#' @return A list with covariance matrix `covx`, shrinkage intensity `rho`,
#' method label, response type, number of slices used, and entropy values.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(180), nrow = 30, ncol = 6)
#' y <- X[, 1] - 0.7 * X[, 2] + rnorm(30, sd = 0.5)
#' est <- mec(X, y, ytype = 'continuous', nslices = 5)
#' c(method = est$method, rho = est$rho)
#' @export
mec <- function(
    x,
    y,
    ytype = c("auto", "continuous", "categorical"),
    nslices = 5,
    eps = 1e-8
){

  ytype <- match.arg(ytype)

  X <- standardise_X(x, centre = TRUE, scale = FALSE)

  n <- nrow(X)
  p <- ncol(X)

  slices <- make_slices(
    y = y,
    ytype = ytype,
    nslices = nslices
  )

  h <- slices$nslices
  slice_indicator <- slices$slice_indicator
  labels <- slices$labels

  if(n - h <= 0){
    stop(
      "n - h must be positive for pooled slice covariance. ",
      "Reduce nslices or increase sample size."
    )
  }

  S_list <- vector("list", length = h)
  ng <- numeric(h)

  for(i in seq_len(h)){

    Xi <- X[slice_indicator == labels[i], , drop = FALSE]

    ni <- nrow(Xi)
    ng[i] <- ni

    if(ni <= 1){
      S_list[[i]] <- eps * diag(p)
    } else {
      S_list[[i]] <- safe_cov(Xi, centre = TRUE)
      S_list[[i]] <- make_positive_definite(S_list[[i]], eps = eps)
    }
  }

  Sp <- Reduce(
    "+",
    lapply(seq_len(h), function(j){
      (ng[j] - 1) * S_list[[j]]
    })
  ) / (n - h)

  Sp <- make_positive_definite(Sp, eps = eps)

  entropy_values <- numeric(h)

  for(i in seq_len(h)){

    eig_values <- eigen(
      S_list[[i]],
      symmetric = TRUE,
      only.values = TRUE
    )$values

    eig_values <- pmax(eig_values, eps)

    entropy_values[i] <- p / 2 +
      0.5 * sum(log(eig_values)) +
      p / 2 * log(2 * pi)
  }

  entropy_slice <- which.max(entropy_values)

  mix <- S_list[[entropy_slice]] + Sp
  mix <- make_positive_definite(mix, eps = eps)

  eig_mix <- eigen(mix, symmetric = TRUE)

  me <- eig_mix$vectors

  Z <- diag(t(me) %*% S_list[[entropy_slice]] %*% me)
  Zp <- diag(t(me) %*% Sp %*% me)

  Zme <- diag(
    pmax(mean(Z), Zp),
    nrow = p
  )

  covxx <- me %*% Zme %*% t(me)
  covxx <- make_positive_definite(covxx, eps = eps)

  target <- mean(diag(covxx)) * diag(p)

  beta <- (sum(diag(covxx)))^2 / sum(diag(covxx %*% covxx))

  if(abs(p - beta) < .Machine$double.eps){
    alpha <- 0
  } else {
    alpha <- 2 * (p * (1 + beta) - 2) / (p - beta)
  }

  mm <- 0.5 * alpha

  rho <- n / (n + mm)
  rho <- min(max(rho, 0), 1)

  covx <- rho * covxx + (1 - rho) * target
  covx <- make_positive_definite(covx, eps = eps)

  return(list(
    covx = covx,
    rho = rho,
    method = "mec",
    response_type = slices$ytype,
    nslices_used = h,
    entropy_values = entropy_values,
    entropy_slice = entropy_slice
  ))
}


#' Estimate a Regularised Covariance Matrix
#'
#' Unified interface for covariance estimators used in shrinkage principal
#' fitted components regression.
#'
#' @param X Numeric predictor matrix.
#' @param y Optional response vector. Required when `method = "mec"`.
#' @param ytype Character. One of `"auto"`, `"continuous"`, or `"categorical"`.
#' @param method Covariance estimator. One of `"mec"`, `"oas"`, `"sre"`,
#' `"sde"`, or `"cse"`.
#' @param nslices Number of slices for MEC with continuous responses.
#' @param gamma Ridge constant for SRE.
#' @param rho Convex shrinkage intensity for CSE.
#' @param eps Positive numerical floor for eigenvalues.
#'
#' @return A list containing the covariance estimate, shrinkage intensity,
#' and method diagnostics.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(120), nrow = 20, ncol = 6)
#' est <- estimate_covariance(X, method = 'oas')
#' c(method = est$method, rho = est$rho)
#' @export
estimate_covariance <- function(
    X,
    y = NULL,
    ytype = c("auto", "continuous", "categorical"),
    method = c("mec", "oas", "sre", "sde", "cse"),
    nslices = 5,
    gamma = 0.1,
    rho = 0.5,
    eps = 1e-8
){

  ytype <- match.arg(ytype)
  method <- match.arg(method)

  if(method == "mec"){

    if(is.null(y)){
      stop("y must be supplied when method = 'mec'.")
    }

    return(mec(
      x = X,
      y = y,
      ytype = ytype,
      nslices = nslices,
      eps = eps
    ))
  }

  if(method == "oas"){
    return(oas(x = X))
  }

  if(method == "sre"){
    return(sre(x = X, gamma = gamma))
  }

  if(method == "sde"){
    return(sde(x = X))
  }

  if(method == "cse"){
    return(cse(x = X, rho = rho))
  }

  stop("Unknown covariance estimation method.")
}

#' Fast Diagonal Covariance Estimator for High-dimensional SPFC
#'
#' Computes a diagonal covariance representation for ultra-high-dimensional
#' settings where forming a full p by p covariance matrix is computationally
#' expensive.
#'
#' @param X Numeric predictor matrix.
#' @param y Optional response vector.
#' @param ytype Response type.
#' @param method Covariance method.
#' @param nslices Number of slices for MEC-style response-guided shrinkage.
#' @param gamma Ridge constant.
#' @param rho Convex shrinkage intensity.
#' @param eps Numerical floor.
#'
#' @return A list containing diagonal covariance entries and diagnostics.
#'
#' @noRd
estimate_covariance_diag_fast <- function(
    X,
    y = NULL,
    ytype = c("auto", "continuous", "categorical"),
    method = c("mec", "oas", "sre", "sde", "cse"),
    nslices = 5,
    gamma = 0.1,
    rho = 0.5,
    eps = 1e-8
){

  ytype <- match.arg(ytype)
  method <- match.arg(method)

  X <- as.matrix(X)

  n <- nrow(X)
  p <- ncol(X)

  diag_sample <- apply(
    X,
    2,
    stats::var
  )

  diag_sample[is.na(diag_sample) | diag_sample <= eps] <- eps

  target <- mean(diag_sample)

  diag_target <- rep(
    target,
    p
  )

  rho_used <- NA_real_
  nslices_used <- NA_integer_

  if(method == "sde"){

    diag_cov <- diag_sample

  } else if(method == "sre"){

    diag_cov <- diag_sample + gamma
    rho_used <- gamma

  } else if(method == "cse"){

    rho_used <- rho

    diag_cov <- (1 - rho_used) * diag_sample +
      rho_used * diag_target

  } else if(method == "oas"){

    # Lightweight OAS-style shrinkage intensity.
    # This avoids forming the full covariance matrix.
    trace_s <- sum(diag_sample)
    trace_s2_diag <- sum(diag_sample^2)

    numerator <- (1 - 2 / p) * trace_s2_diag + trace_s^2
    denominator <- (n + 1 - 2 / p) *
      (trace_s2_diag - trace_s^2 / p)

    if(is.finite(denominator) && denominator > eps){
      rho_used <- min(
        1,
        max(
          0,
          numerator / denominator
        )
      )
    } else {
      rho_used <- 1
    }

    diag_cov <- (1 - rho_used) * diag_sample +
      rho_used * diag_target

  } else if(method == "mec"){

    if(is.null(y)){
      stop("y must be supplied for method = 'mec'.")
    }

    detected_ytype <- detect_y_type(
      y = y,
      ytype = ytype
    )

    slices <- make_slices(
      y = y,
      ytype = detected_ytype,
      nslices = nslices
    )

    slice_id <- slices$slice_indicator
    nslices_used <- slices$nslices

    diag_within <- rep(
      0,
      p
    )

    for(h in sort(unique(slice_id))){

      idx <- which(slice_id == h)

      if(length(idx) > 1){

        wh <- length(idx) / n

        diag_h <- apply(
          X[idx, , drop = FALSE],
          2,
          stats::var
        )

        diag_h[is.na(diag_h) | diag_h <= eps] <- eps

        diag_within <- diag_within + wh * diag_h
      }
    }

    diag_within[diag_within <= eps] <- eps

    # Response-guided shrinkage intensity.
    diff_norm <- sum(
      (diag_sample - diag_within)^2
    )

    base_norm <- sum(
      (diag_sample - diag_target)^2
    )

    if(is.finite(base_norm) && base_norm > eps){
      rho_used <- min(
        1,
        max(
          0,
          diff_norm / (diff_norm + base_norm)
        )
      )
    } else {
      rho_used <- 0.5
    }

    diag_cov <- (1 - rho_used) * diag_sample +
      rho_used * diag_within
  }

  diag_cov[is.na(diag_cov) | diag_cov <= eps] <- eps

  out <- list(
    method = method,
    representation = "diagonal",
    diag_cov = diag_cov,
    rho = rho_used,
    nslices_used = nslices_used,
    n = n,
    p = p
  )

  return(out)
}
