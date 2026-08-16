#' Check Whether a Matrix is Positive Definite
#'
#' Checks whether a square numeric matrix is positive definite by inspecting
#' the eigenvalues of its symmetrised form.
#'
#' @param A A numeric square matrix.
#' @param tol Numeric tolerance. Eigenvalues greater than `tol` are treated as positive.
#'
#' @return A logical value. `TRUE` if the matrix is positive definite, otherwise `FALSE`.
#'
#' @noRd
is_positive_definite <- function(A, tol = 1e-8){

  if(!is.matrix(A)){
    A <- as.matrix(A)
  }

  if(nrow(A) != ncol(A)){
    return(FALSE)
  }

  A <- (A + t(A)) / 2

  eig_values <- eigen(
    A,
    symmetric = TRUE,
    only.values = TRUE
  )$values

  return(all(eig_values > tol))
}


#' Force a Symmetric Matrix to be Positive Definite
#'
#' Converts a symmetric matrix to a numerically positive definite matrix by
#' flooring its eigenvalues at a small positive value.
#'
#' @param A A numeric square matrix.
#' @param eps Positive numerical floor for eigenvalues.
#'
#' @return A symmetric positive definite matrix.
#'
#' @noRd
make_positive_definite <- function(A, eps = 1e-8){

  if(!is.matrix(A)){
    A <- as.matrix(A)
  }

  if(nrow(A) != ncol(A)){
    stop("A must be a square matrix.")
  }

  A <- (A + t(A)) / 2

  eig <- eigen(
    A,
    symmetric = TRUE
  )

  eig_values <- pmax(eig$values, eps)

  A_pd <- eig$vectors %*%
    diag(eig_values, nrow = length(eig_values)) %*%
    t(eig$vectors)

  A_pd <- (A_pd + t(A_pd)) / 2

  return(A_pd)
}


#' Compute a Safe Maximum Likelihood-Type Covariance Matrix
#'
#' Computes the maximum likelihood-type covariance matrix using divisor `n`
#' rather than `n - 1`.
#'
#' @param X A numeric matrix with observations in rows and variables in columns.
#' @param centre Logical. If `TRUE`, variables are centred before covariance computation.
#'
#' @return A numeric covariance matrix.
#'
#' @noRd
safe_cov <- function(X, centre = TRUE){

  if(!is.matrix(X)){
    X <- as.matrix(X)
  }

  if(!is.numeric(X)){
    stop("X must be numeric.")
  }

  n <- nrow(X)

  if(n <= 1){
    stop("At least two observations are required to compute covariance.")
  }

  if(centre){
    X <- scale(X, center = TRUE, scale = FALSE)
  }

  S <- crossprod(X) / n

  S <- (S + t(S)) / 2

  return(S)
}


#' Compute the Condition Number of a Matrix
#'
#' Computes the spectral condition number of a symmetric matrix using its
#' positive eigenvalues.
#'
#' @param A A numeric square matrix.
#' @param eps Small positive threshold below which eigenvalues are ignored.
#'
#' @return A numeric condition number.
#'
#' @noRd
condition_number <- function(A, eps = 1e-10){

  if(!is.matrix(A)){
    A <- as.matrix(A)
  }

  if(nrow(A) != ncol(A)){
    stop("A must be a square matrix.")
  }

  A <- (A + t(A)) / 2

  eig_values <- eigen(
    A,
    symmetric = TRUE,
    only.values = TRUE
  )$values

  eig_values <- eig_values[eig_values > eps]

  if(length(eig_values) == 0){
    return(Inf)
  }

  return(max(eig_values) / min(eig_values))
}


#' Standardise a Numeric Predictor Matrix
#'
#' Centres and optionally scales a numeric predictor matrix.
#'
#' @param X A numeric matrix.
#' @param centre Logical. If `TRUE`, centre columns.
#' @param scale Logical. If `TRUE`, scale columns to unit variance.
#'
#' @return A numeric matrix.
#'
#' @noRd
standardise_X <- function(X, centre = TRUE, scale = FALSE){

  if(!is.matrix(X)){
    X <- as.matrix(X)
  }

  if(!is.numeric(X)){
    stop("X must be numeric.")
  }

  Xs <- scale(
    X,
    center = centre,
    scale = scale
  )

  return(as.matrix(Xs))
}
