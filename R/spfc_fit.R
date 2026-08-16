#' Solve a Symmetric Generalised Eigenvalue Problem
#'
#' Solves the generalised eigenvalue problem
#' \deqn{Bv = \lambda Av}
#' for symmetric matrices, where `A` is regularised to be positive definite.
#'
#' @param A Numeric positive definite matrix, usually a covariance estimate.
#' @param B Numeric symmetric matrix, usually a fitted covariance matrix.
#' @param d Number of leading eigenvectors to return.
#' @param eps Positive numerical floor for eigenvalues.
#'
#' @return A list containing eigenvalues and eigenvectors.
#'
#' @noRd
generalised_eigen <- function(
    A,
    B,
    d,
    eps = 1e-8
){

  if(!is.matrix(A)){
    A <- as.matrix(A)
  }

  if(!is.matrix(B)){
    B <- as.matrix(B)
  }

  if(nrow(A) != ncol(A)){
    stop("A must be a square matrix.")
  }

  if(nrow(B) != ncol(B)){
    stop("B must be a square matrix.")
  }

  if(nrow(A) != nrow(B)){
    stop("A and B must have the same dimensions.")
  }

  p <- nrow(A)

  if(d < 1 || d > p){
    stop("d must lie between 1 and ncol(A).")
  }

  A <- make_positive_definite(
    A,
    eps = eps
  )

  B <- (B + t(B)) / 2

  chol_A <- chol(A)

  A_inv_half <- backsolve(
    chol_A,
    diag(p)
  )

  M <- t(A_inv_half) %*% B %*% A_inv_half
  M <- (M + t(M)) / 2

  eig <- eigen(
    M,
    symmetric = TRUE
  )

  ord <- order(
    eig$values,
    decreasing = TRUE
  )

  eigenvalues <- eig$values[ord]

  eigenvectors_std <- eig$vectors[
    ,
    ord,
    drop = FALSE
  ]

  eigenvectors <- A_inv_half %*% eigenvectors_std

  eigenvectors <- eigenvectors[
    ,
    seq_len(d),
    drop = FALSE
  ]

  eigenvectors <- qr.Q(
    qr(eigenvectors)
  )

  return(list(
    values = eigenvalues,
    vectors = eigenvectors
  ))
}

#' Fast Low-rank Generalised Eigen Solver for SPFC
#'
#' Solves the SPFC generalised eigenproblem using the low-rank
#' structure of the fitted covariance matrix.
#'
#' @param A Covariance matrix.
#' @param X Centred/scaled predictor matrix.
#' @param F Response basis matrix.
#' @param d Number of directions.
#' @param eps Numerical floor.
#'
#' @return A list with eigenvalues and eigenvectors.
#'
#' @noRd
generalised_eigen_pfc_fast <- function(
    A,
    X,
    F,
    d,
    eps = 1e-8
){

  A <- as.matrix(A)
  X <- as.matrix(X)
  F <- as.matrix(F)

  p <- ncol(X)
  n <- nrow(X)

  if(d < 1 || d > p){
    stop("d must lie between 1 and ncol(X).")
  }

  A <- make_positive_definite(
    A,
    eps = eps
  )

  qrf <- qr(F)

  rank_f <- qrf$rank

  if(rank_f < 1){
    stop("The response basis matrix F has rank zero.")
  }

  Q <- qr.Q(
    qrf,
    complete = FALSE
  )

  Q <- Q[
    ,
    seq_len(rank_f),
    drop = FALSE
  ]

  L <- crossprod(
    X,
    Q
  ) / sqrt(n)

  R <- chol(A)

  C <- forwardsolve(
    t(R),
    L
  )

  sv <- svd(
    C,
    nu = min(d, ncol(C), nrow(C)),
    nv = 0
  )

  available_d <- ncol(sv$u)

  use_d <- min(
    d,
    available_d
  )

  U <- sv$u[
    ,
    seq_len(use_d),
    drop = FALSE
  ]

  eigenvalues <- rep(
    0,
    p
  )

  eigenvalues[
    seq_len(length(sv$d))
  ] <- sv$d^2

  eigenvectors_use <- backsolve(
    R,
    U
  )

  if(use_d < d){

    padding <- matrix(
      0,
      nrow = p,
      ncol = d - use_d
    )

    eigenvectors_use <- cbind(
      eigenvectors_use,
      padding
    )
  }

  eigenvectors <- eigenvectors_use[
    ,
    seq_len(d),
    drop = FALSE
  ]

  eigenvectors <- qr.Q(
    qr(eigenvectors)
  )

  colnames(eigenvectors) <- paste0(
    "dir",
    seq_len(ncol(eigenvectors))
  )

  return(
    list(
      values = eigenvalues,
      vectors = eigenvectors
    )
  )
}


#' Fast Diagonal-dual Generalised Eigen Solver for SPFC
#'
#' Uses a diagonal covariance representation and the low-rank fitted
#' covariance structure to avoid forming or factorising a p by p matrix.
#'
#' @param diag_A Positive covariance diagonal.
#' @param X Centred/scaled predictor matrix.
#' @param F Response basis matrix.
#' @param d Number of directions.
#' @param eps Numerical floor.
#'
#' @return A list with eigenvalues and eigenvectors.
#'
#' @noRd
generalised_eigen_pfc_dual_diag <- function(
    diag_A,
    X,
    F,
    d,
    eps = 1e-8
){

  X <- as.matrix(X)
  F <- as.matrix(F)

  p <- ncol(X)
  n <- nrow(X)

  if(length(diag_A) != p){
    stop("length(diag_A) must equal ncol(X).")
  }

  if(d < 1 || d > p){
    stop("d must lie between 1 and ncol(X).")
  }

  diag_A[is.na(diag_A) | diag_A <= eps] <- eps

  qrf <- qr(F)
  rank_f <- qrf$rank

  if(rank_f < 1){
    stop("The response basis matrix F has rank zero.")
  }

  Q <- qr.Q(
    qrf,
    complete = FALSE
  )

  Q <- Q[
    ,
    seq_len(rank_f),
    drop = FALSE
  ]

  L <- crossprod(
    X,
    Q
  ) / sqrt(n)

  C <- sweep(
    L,
    1,
    sqrt(diag_A),
    FUN = "/"
  )

  sv <- svd(
    C,
    nu = min(d, ncol(C), nrow(C)),
    nv = 0
  )

  available_d <- ncol(sv$u)

  use_d <- min(
    d,
    available_d
  )

  U <- sv$u[
    ,
    seq_len(use_d),
    drop = FALSE
  ]

  eigenvectors_use <- sweep(
    U,
    1,
    sqrt(diag_A),
    FUN = "/"
  )

  if(use_d < d){

    padding <- matrix(
      0,
      nrow = p,
      ncol = d - use_d
    )

    eigenvectors_use <- cbind(
      eigenvectors_use,
      padding
    )
  }

  eigenvectors <- eigenvectors_use[
    ,
    seq_len(d),
    drop = FALSE
  ]

  eigenvectors <- qr.Q(
    qr(eigenvectors)
  )

  eigenvalues <- rep(
    0,
    p
  )

  eigenvalues[
    seq_len(length(sv$d))
  ] <- sv$d^2

  colnames(eigenvectors) <- paste0(
    "dir",
    seq_len(ncol(eigenvectors))
  )

  return(
    list(
      values = eigenvalues,
      vectors = eigenvectors
    )
  )
}


#' Fit Shrinkage Principal Fitted Components
#'
#' Fits a shrinkage principal fitted components model by replacing the unstable
#' covariance object in the PFC reduction step with a selected shrinkage
#' covariance estimator.
#'
#' @param X Numeric predictor matrix with observations in rows and variables in columns.
#' @param y Response vector.
#' @param d Number of reduction directions.
#' @param fy Optional user-supplied response basis matrix.
#' @param ytype Character. One of `"auto"`, `"continuous"`, or `"categorical"`.
#' @param cov_method Covariance estimator. One of `"mec"`, `"oas"`, `"sre"`,
#' `"sde"`, or `"cse"`.
#' @param nslices Number of response slices for continuous responses.
#' @param poly_degree Polynomial degree for continuous response basis construction.
#' @param gamma Ridge constant for SRE.
#' @param rho Convex shrinkage intensity for CSE.
#' @param centre_x Logical. If `TRUE`, centre predictor columns.
#' @param scale_x Logical. If `TRUE`, scale predictor columns.
#' @param eps Positive numerical floor for eigenvalues.
#'
#' @return A `spfc_fit` object containing reduced scores, directions, fitted
#' covariance objects, response basis, preprocessing metadata, and diagnostics.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(180), nrow = 30, ncol = 6)
#' y <- X[, 1] - 0.7 * X[, 2] + rnorm(30, sd = 0.5)
#' fit <- spfc_fit(
#'   X, y, d = 1, ytype = 'continuous',
#'   cov_method = 'oas', poly_degree = 2
#' )
#' coef(fit)
#' @export
spfc_fit <- function(
    X,
    y,
    d,
    fy = NULL,
    ytype = c("auto", "continuous", "categorical"),
    cov_method = c("mec", "oas", "sre", "sde", "cse"),
    nslices = 5,
    poly_degree = 1,
    gamma = 0.1,
    rho = 0.5,
    centre_x = TRUE,
    scale_x = FALSE,
    eps = 1e-8
){

  ytype <- match.arg(ytype)
  cov_method <- match.arg(cov_method)

  X_original <- as.matrix(X)

  x_means <- colMeans(
    X_original,
    na.rm = TRUE
  )

  x_sds <- apply(
    X_original,
    2,
    stats::sd,
    na.rm = TRUE
  )

  x_sds[is.na(x_sds) | x_sds == 0] <- 1

  X <- standardise_X(
    X_original,
    centre = centre_x,
    scale = scale_x
  )

  n <- nrow(X)
  p <- ncol(X)

  if(length(y) != n){
    stop("length(y) must equal nrow(X).")
  }

  if(d < 1 || d > p){
    stop("d must lie between 1 and ncol(X).")
  }

  detected_ytype <- detect_y_type(
    y = y,
    ytype = ytype
  )

  F <- build_fy(
    y = y,
    ytype = detected_ytype,
    fy = fy,
    nslices = nslices,
    poly_degree = poly_degree,
    centre = TRUE
  )

  if(nrow(F) != n){
    stop("nrow(fy) must equal nrow(X).")
  }

  FtF <- crossprod(F)

  if(qr(FtF)$rank < ncol(FtF)){
    FtF <- FtF + eps * diag(ncol(FtF))
  }

  PF <- F %*% solve(FtF) %*% t(F)

  S_fit <- crossprod(
    X,
    PF %*% X
  ) / n

  S_fit <- (S_fit + t(S_fit)) / 2

  use_fast_diag <- p > n

  use_fast_diag <- p > n

  if(isTRUE(use_fast_diag)){

    cov_est <- estimate_covariance_diag_fast(
      X = X,
      y = y,
      ytype = detected_ytype,
      method = cov_method,
      nslices = nslices,
      gamma = gamma,
      rho = rho,
      eps = eps
    )

    Sigma_hat <- NULL

    eig <- generalised_eigen_pfc_dual_diag(
      diag_A = cov_est$diag_cov,
      X = X,
      F = F,
      d = d,
      eps = eps
    )

  } else {

    cov_est <- estimate_covariance(
      X = X,
      y = y,
      ytype = detected_ytype,
      method = cov_method,
      nslices = nslices,
      gamma = gamma,
      rho = rho,
      eps = eps
    )

    Sigma_hat <- cov_est$covx

    eig <- generalised_eigen_pfc_fast(
      A = Sigma_hat,
      X = X,
      F = F,
      d = d,
      eps = eps
    )
  }

  V <- eig$vectors
  Z <- X %*% V

  colnames(V) <- paste0(
    "dir",
    seq_len(ncol(V))
  )

  colnames(Z) <- paste0(
    "SPFC",
    seq_len(ncol(Z))
  )

  out <- list(
    call = match.call(),
    X = X,
    X_original = X_original,
    y = y,
    ytype = detected_ytype,
    F = F,
    PF = PF,
    S_fit = S_fit,
    Sigma = Sigma_hat,
    cov_method = cov_est$method,
    covariance = cov_est,
    solver = if(isTRUE(use_fast_diag)){
      "diagonal_dual"
    } else {
      "low_rank_primal"
    },
    V = V,
    Z = Z,
    eigvals = eig$values,
    d = d,
    rho = cov_est$rho,
    nslices_used = if(!is.null(cov_est$nslices_used)){
      cov_est$nslices_used
    } else {
      NA_integer_
    },
    preprocessing = list(
      centred = centre_x,
      scaled = scale_x,
      means = x_means,
      sds = x_sds
    )
  )

  class(out) <- c("spfc_fit", "list")

  return(out)
}
