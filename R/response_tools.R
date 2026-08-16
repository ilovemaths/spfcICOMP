#' Detect Response Type
#'
#' Detects whether a response variable should be treated as continuous or
#' categorical. The user may override automatic detection by setting `ytype`
#' explicitly.
#'
#' @param y Response vector.
#' @param ytype Character. One of `"auto"`, `"continuous"`, or `"categorical"`.
#' @param unique_cutoff Integer cutoff for treating numeric responses with few
#' unique values as categorical.
#' @param prop_unique_cutoff Proportion cutoff for treating numeric responses
#' with few unique values as categorical.
#'
#' @return A character string: `"continuous"` or `"categorical"`.
#'
#' @noRd
detect_y_type <- function(
    y,
    ytype = c("auto", "continuous", "categorical"),
    unique_cutoff = 10,
    prop_unique_cutoff = 0.1
){

  ytype <- match.arg(ytype)

  if(ytype != "auto"){
    return(ytype)
  }

  if(is.factor(y) || is.character(y) || is.logical(y)){
    return("categorical")
  }

  if(is.numeric(y) || is.integer(y)){

    n <- length(y)
    nunique <- length(unique(y))
    prop_unique <- nunique / n

    if(nunique <= unique_cutoff && prop_unique <= prop_unique_cutoff){
      return("categorical")
    }

    return("continuous")
  }

  stop("Unable to determine response type. Please set ytype explicitly.")
}


#' Create Response Slices
#'
#' Constructs response slices for continuous or categorical responses. For
#' categorical responses, the observed classes are used as slices. For continuous
#' responses, empirical quantile slices are used.
#'
#' @param y Response vector.
#' @param ytype Character. One of `"auto"`, `"continuous"`, or `"categorical"`.
#' @param nslices Number of slices for continuous responses.
#'
#' @return A list containing response type, number of slices, slice indicators,
#' and slice labels.
#'
#' @noRd
make_slices <- function(
    y,
    ytype = c("auto", "continuous", "categorical"),
    nslices = 5
){

  detected_ytype <- detect_y_type(y, ytype = ytype)

  if(detected_ytype == "categorical"){

    yf <- as.factor(y)

    slice_indicator <- as.integer(yf)

    return(list(
      ytype = "categorical",
      nslices = length(levels(yf)),
      slice_indicator = slice_indicator,
      labels = seq_along(levels(yf)),
      levels = levels(yf)
    ))
  }

  y_numeric <- as.numeric(y)

  if(anyNA(y_numeric)){
    stop("Continuous response contains missing or non-numeric values.")
  }

  if(length(unique(y_numeric)) < 2){
    stop("Continuous response must contain at least two unique values.")
  }

  probs <- seq(0, 1, length.out = nslices + 1)

  breaks <- unique(
    as.numeric(
      stats::quantile(
        y_numeric,
        probs = probs,
        na.rm = TRUE,
        type = 7
      )
    )
  )

  if(length(breaks) <= 2){
    stop(
      "Continuous response cannot be sliced meaningfully. ",
      "Reduce nslices or check the variation in y."
    )
  }

  slice_indicator <- cut(
    y_numeric,
    breaks = breaks,
    include.lowest = TRUE,
    labels = FALSE
  )

  slice_indicator <- as.integer(slice_indicator)

  labels <- sort(unique(slice_indicator))

  return(list(
    ytype = "continuous",
    nslices = length(labels),
    slice_indicator = slice_indicator,
    labels = labels,
    levels = NULL
  ))
}


#' Build PFC Basis Functions for the Response
#'
#' Builds the response basis matrix used in principal fitted components. For
#' categorical responses, indicator basis functions are used. For continuous
#' responses, polynomial basis functions are used unless `fy` is supplied.
#'
#' @param y Response vector.
#' @param ytype Character. One of `"auto"`, `"continuous"`, or `"categorical"`.
#' @param fy Optional user-supplied basis matrix.
#' @param nslices Number of slices. Currently retained for interface consistency.
#' @param poly_degree Polynomial degree for continuous responses.
#' @param centre Logical. If `TRUE`, centre the basis columns.
#'
#' @return A numeric matrix of response basis functions.
#'
#' @noRd
build_fy <- function(
    y,
    ytype = c("auto", "continuous", "categorical"),
    fy = NULL,
    nslices = 5,
    poly_degree = 1,
    centre = TRUE
){

  if(!is.null(fy)){

    fy <- as.matrix(fy)

    if(!is.numeric(fy)){
      stop("fy must be numeric.")
    }

    return(fy)
  }

  detected_ytype <- detect_y_type(y, ytype = ytype)

  if(detected_ytype == "categorical"){

    yf <- as.factor(y)

    F <- stats::model.matrix(~ yf - 1)

    if(ncol(F) >= 2){
      F <- F[, -1, drop = FALSE]
    }

    if(centre){
      F <- scale(F, center = TRUE, scale = FALSE)
    }

    return(as.matrix(F))
  }

  y_numeric <- as.numeric(y)

  if(poly_degree < 1){
    stop("poly_degree must be at least 1.")
  }

  basis_list <- vector("list", length = poly_degree)

  for(k in seq_len(poly_degree)){
    basis_list[[k]] <- y_numeric^k
  }

  F <- do.call(cbind, basis_list)

  colnames(F) <- paste0("y_power_", seq_len(poly_degree))

  if(centre){
    F <- scale(F, center = TRUE, scale = FALSE)
  }

  return(as.matrix(F))
}
