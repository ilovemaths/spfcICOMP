#' Compute Adaptive Weighted L1 Weights
#'
#' Computes adaptive weighted L1-type sparsity weights from the estimated
#' SPFC loading matrix. This is the thesis-safe adaptive sparsity route.
#'
#' @param V Numeric loading or direction matrix, usually `spfc_fit()$V`.
#' @param gamma Positive adaptivity exponent.
#' @param eps Small positive constant to avoid division by zero.
#'
#' @return A numeric vector of normalised adaptive weights.
#'
#' @noRd
adaptive_weighted_l1 <- function(
    V,
    gamma = 1,
    eps = 1e-6
){

  V <- as.matrix(V)

  if(!is.numeric(V)){
    stop("V must be numeric.")
  }

  if(gamma <= 0){
    stop("gamma must be positive.")
  }

  loading_strength <- rowSums(V^2)

  weights <- 1 / ((abs(loading_strength)^gamma) + eps)

  weights <- weights / sum(weights)

  return(weights)
}


#' Compute Scale-Invariant C1F Complexity
#'
#' Computes Bozdogan's scale-invariant C1F covariance complexity measure from
#' positive eigenvalues.
#'
#' @param evals Numeric vector of eigenvalues.
#' @param eps Small positive threshold for retaining eigenvalues.
#'
#' @return A numeric C1F value.
#'
#' @noRd
c1f_scale_invariant <- function(
    evals,
    eps = 1e-10
){

  evals <- as.numeric(evals)

  evals <- evals[
    is.finite(evals) &
      evals > eps
  ]

  if(length(evals) == 0){
    stop("No positive eigenvalues supplied.")
  }

  lambda_bar <- mean(evals)

  c1f <- sum((evals - lambda_bar)^2) /
    (4 * lambda_bar^2)

  return(as.numeric(c1f))
}


#' Compute C1F-Informed Extension Weights
#'
#' Computes an experimental C1F-informed adaptive sparsity weighting vector from
#' SPFC directions and associated eigenvalues. This is the research-extension
#' route and should be interpreted as a complexity-informed weighting scheme,
#' not as the original scalar C1F criterion itself.
#'
#' @param V Numeric SPFC direction matrix, usually `spfc_fit()$V`.
#' @param eigvals Numeric vector of eigenvalues, usually `spfc_fit()$eigvals`.
#' @param eps Small positive constant to avoid zero weights.
#'
#' @return A numeric vector of normalised C1F-informed weights.
#'
#' @noRd
c1f_extension <- function(
    V,
    eigvals,
    eps = 1e-10
){

  V <- as.matrix(V)

  if(!is.numeric(V)){
    stop("V must be numeric.")
  }

  eigvals <- as.numeric(eigvals)

  d <- ncol(V)

  if(length(eigvals) < d){
    stop("eigvals must contain at least ncol(V) values.")
  }

  eig_use <- eigvals[seq_len(d)]

  eig_use <- pmax(eig_use, eps)

  contribution <- rowSums(
    sweep(
      V^2,
      MARGIN = 2,
      STATS = eig_use,
      FUN = "*"
    )
  )

  contribution <- pmax(contribution, eps)

  weights <- contribution / sum(contribution)

  return(weights)
}


#' Select Variables from SPFC Directions
#'
#' Selects variables using either the thesis-safe adaptive weighted L1 route or
#' the experimental C1F-informed extension route.
#'
#' @param fit Object returned by `spfc_fit()`.
#' @param method Character. One of `"adaptive_weighted_l1"` or `"c1f_extension"`.
#' @param selection_rule Character. One of `"c1f"`, `"quantile"`, or `"fixed"`.
#' @param reduced_model Fitted reduced model used to obtain the model-parameter covariance when `selection_rule = "c1f"`.
#' @param threshold Numeric threshold applied to variable importance scores.
#' If `NULL`, a quantile-based threshold is used.
#' @param quantile_cut Numeric quantile used when `threshold = NULL`.
#' @param gamma Positive adaptivity exponent for `adaptive_weighted_l1`.
#' @param c1f_calibration Character. `"icomp_hd_floor"` is the frozen thesis rule that combines a robust empirical loading scale, the high-dimensional multiplicity floor `sqrt(2 log(p) / n)`, and the bounded C1F multiplier `1 + C1F/(1 + C1F)`. `"robust_universal"` retains the version 0.0.0.9004 calibration for reproducibility, and `"raw"` retains the historical uncalibrated C1F route.
#' @param eps Small positive numerical constant.
#'
#' @return A data frame containing variable indices, importance scores, weights,
#' and selected status.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(180), nrow = 30, ncol = 6)
#' y <- X[, 1] - 0.7 * X[, 2] + rnorm(30, sd = 0.5)
#' fit <- spfc_fit(X, y, d = 1, ytype = 'continuous',
#'                 cov_method = 'oas', poly_degree = 2)
#' reduced_model <- fit_reduced_model(
#'   fit$Z, y, ytype = 'continuous'
#' )
#' sel <- spfc_select_variables(
#'   fit, method = 'c1f_extension', selection_rule = 'c1f',
#'   reduced_model = reduced_model
#' )
#' sel[, c('variable', 'shrunk_importance', 'selected')]
#' @export
spfc_select_variables <- function(
    fit,
    method = c("adaptive_weighted_l1", "c1f_extension"),
    selection_rule = c("c1f", "quantile", "fixed"),
    reduced_model = NULL,
    threshold = NULL,
    quantile_cut = 0.75,
    gamma = 1,
    c1f_calibration = c("icomp_hd_floor", "robust_universal", "raw"),
    eps = 1e-10
){

  method <- match.arg(method)
  selection_rule <- match.arg(selection_rule)
  c1f_calibration <- match.arg(c1f_calibration)

  if(is.null(fit$V)){
    stop("fit must contain V, the SPFC direction matrix.")
  }

  V <- as.matrix(fit$V)
  p <- nrow(V)
  importance <- sqrt(rowSums(V^2))

  if(method == "adaptive_weighted_l1"){
    weights <- adaptive_weighted_l1(
      V = V,
      gamma = gamma,
      eps = eps
    )
  } else {
    if(is.null(fit$eigvals)){
      stop("fit must contain eigvals for method = 'c1f_extension'.")
    }

    contribution_weights <- c1f_extension(
      V = V,
      eigvals = fit$eigvals,
      eps = eps
    )

    # Convert contribution weights to penalty weights: variables with stronger
    # fitted contributions receive smaller adaptive L1 penalties.
    weights <- 1 / pmax(contribution_weights, eps)
    weights <- weights / sum(weights)
  }

  c1f_value <- NA_real_
  penalty <- rep(NA_real_, p)
  shrunk_importance <- importance

  if(selection_rule == "c1f"){
    if(is.null(reduced_model)){
      stop(
        "selection_rule = 'c1f' requires reduced_model so that C1F can be ",
        "computed from the fitted model covariance."
      )
    }

    Sigma_model <- tryCatch(
      extract_ifim(reduced_model, eps = eps),
      error = function(e) NULL
    )

    if(is.null(Sigma_model)){
      stop(
        "Could not obtain the fitted-model covariance required for C1F-based ",
        "variable selection."
      )
    }

    c1f_value <- compute_c1f(Sigma_model, eps = eps)

    if(c1f_calibration == "raw"){
      # Legacy research route: use the raw scalar C1F with the historical
      # sum-to-one adaptive weights. Retained for reproducibility only.
      penalty <- c1f_value * weights
      complexity_fraction <- NA_real_
      loading_scale <- NA_real_
      hd_factor <- NA_real_
      complexity_multiplier <- NA_real_
      global_penalty <- NA_real_
    } else {
      # Both calibrated routes use relative adaptive weights with median one.
      # This removes the 1/p attenuation created by sum-to-one normalisation.
      positive_w <- weights[is.finite(weights) & weights > eps]
      weight_median <- if(length(positive_w)) stats::median(positive_w) else 1
      relative_weights <- weights / max(weight_median, eps)

      # Robust empirical loading scale.
      positive_i <- importance[is.finite(importance) & importance > eps]
      loading_scale <- if(length(positive_i)) stats::median(positive_i) else 0

      complexity_fraction <- c1f_value / (1 + c1f_value)

      if(c1f_calibration == "robust_universal"){
        # Version 0.0.0.9004 calibration retained unchanged for reproducibility.
        hd_factor <- 2 * log(max(p, 2L))
        complexity_multiplier <- complexity_fraction
        global_penalty <- loading_scale * hd_factor * complexity_multiplier
      } else {
        # Frozen thesis calibration (ICOMP-HD floor):
        #
        # lambda_IC-HD = s_Gamma * sqrt(2 log(p) / n) *
        #                [1 + C1F/(1 + C1F)].
        #
        # The sqrt(2 log(p) / n) term supplies a deterministic
        # high-dimensional multiplicity/noise floor. The bounded C1F term
        # strengthens that floor according to fitted-model covariance
        # complexity. No cross-validation or user-chosen LASSO lambda is used.
        n_obs <- if(!is.null(fit$X)) nrow(as.matrix(fit$X)) else NA_integer_
        if(!is.finite(n_obs) || n_obs < 2L){
          stop("The ICOMP-HD C1F calibration requires at least two observations in fit$X.")
        }
        hd_factor <- sqrt(2 * log(max(p, 2L)) / n_obs)
        complexity_multiplier <- 1 + complexity_fraction
        global_penalty <- loading_scale * hd_factor * complexity_multiplier
      }

      penalty <- global_penalty * relative_weights
    }

    shrunk_importance <- pmax(importance - penalty, 0)
    selected <- shrunk_importance > eps
    threshold_used <- NA_real_

  } else if(selection_rule == "quantile"){
    if(is.null(threshold)){
      threshold <- as.numeric(
        stats::quantile(
          importance,
          probs = quantile_cut,
          na.rm = TRUE
        )
      )
    }
    selected <- importance >= threshold
    threshold_used <- threshold

  } else {
    if(is.null(threshold) || length(threshold) != 1L || !is.finite(threshold)){
      stop("selection_rule = 'fixed' requires one finite numeric threshold.")
    }
    selected <- importance >= threshold
    threshold_used <- threshold
  }

  out <- data.frame(
    variable = seq_len(p),
    importance = importance,
    weight = weights,
    penalty = penalty,
    shrunk_importance = shrunk_importance,
    selected = selected
  )

  out <- out[order(out$shrunk_importance, out$importance, decreasing = TRUE), ]
  rownames(out) <- NULL

  attr(out, "method") <- method
  attr(out, "selection_rule") <- selection_rule
  attr(out, "threshold") <- threshold_used
  attr(out, "quantile_cut") <- if(selection_rule == "quantile") quantile_cut else NA_real_
  attr(out, "C1F") <- c1f_value
  attr(out, "c1f_calibration") <- if(selection_rule == "c1f") c1f_calibration else NA_character_
  if(selection_rule == "c1f" && c1f_calibration != "raw"){
    attr(out, "c1f_complexity_fraction") <- complexity_fraction
    attr(out, "loading_scale") <- loading_scale
    attr(out, "high_dimensional_factor") <- hd_factor
    attr(out, "c1f_complexity_multiplier") <- complexity_multiplier
    attr(out, "global_penalty") <- global_penalty
  } else {
    attr(out, "c1f_complexity_fraction") <- NA_real_
    attr(out, "loading_scale") <- NA_real_
    attr(out, "high_dimensional_factor") <- NA_real_
    attr(out, "c1f_complexity_multiplier") <- NA_real_
    attr(out, "global_penalty") <- NA_real_
  }

  return(out)
}
