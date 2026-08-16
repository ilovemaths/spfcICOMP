#' Generate an AR(1) Covariance Matrix
#'
#' @param p Number of variables.
#' @param rho Correlation parameter.
#'
#' @return A p by p covariance matrix.
#'
#' @noRd
make_ar1_cov <- function(p, rho = 0.5){

  if(p < 1){
    stop("p must be positive.")
  }

  if(abs(rho) >= 1){
    stop("rho must lie between -1 and 1.")
  }

  idx <- seq_len(p)

  Sigma <- outer(
    idx,
    idx,
    function(i, j) rho^abs(i - j)
  )

  return(Sigma)
}


#' Generate Sparse True Directions
#'
#' @param p Number of predictors.
#' @param d Structural dimension.
#' @param s Number of active variables.
#' @param seed Optional random seed.
#'
#' @return A list containing true direction matrix and active set.
#'
#' @noRd
generate_sparse_directions <- function(
    p,
    d = 1,
    s = 5,
    seed = NULL
){

  if(!is.null(seed)){
    set.seed(seed)
  }

  if(s > p){
    stop("s must not exceed p.")
  }

  if(d > s){
    stop("d must not exceed s for sparse direction generation.")
  }

  active_set <- seq_len(s)

  B <- matrix(
    0,
    nrow = p,
    ncol = d
  )

  B[active_set, ] <- matrix(
    stats::rnorm(s * d),
    nrow = s,
    ncol = d
  )

  B <- qr.Q(qr(B))

  return(list(
    B = B[, seq_len(d), drop = FALSE],
    active_set = active_set
  ))
}


#' Simulate Continuous High-Dimensional PFC Data
#'
#' Generates data directly from a principal fitted components inverse model,
#' so the nominal structural dimension is identifiable by construction. For
#' `d = 1`, the inverse mean uses `f(Y) = Y`. For `d = 2`, it uses the two
#' population-orthonormal basis functions `Y` and `(Y^2 - 1)/sqrt(2)` with
#' `Y ~ N(0,1)`. The inverse-mean loading matrix is constructed so that the true
#' sufficient reduction subspace is the sparse matrix returned as `B_true`,
#' even when the residual predictor covariance is correlated.
#'
#' An independent test sample is generated from the same inverse model for
#' honest out-of-sample prediction assessment.
#'
#' @param n Training sample size.
#' @param p Number of predictors.
#' @param d Structural dimension. The validated thesis DGP supports 1 or 2.
#' @param s Number of active variables.
#' @param rho_x AR(1) correlation parameter of the inverse-model residual
#' covariance.
#' @param snr Inverse-model signal-to-noise ratio in the sufficient coordinates.
#' @param n_test Independent test-sample size. Defaults to `n`.
#' @param seed Optional random seed.
#'
#' @return A list containing training and test data, true PFC quantities, sparse
#' reduction directions, active set, and residual covariance.
#'
#' @examples
#' sim <- simulate_spfc_continuous(
#'   n = 30, p = 8, d = 1, s = 3, rho_x = 0.3,
#'   snr = 2, n_test = 10, seed = 123
#' )
#' c(train_n = length(sim$y), test_n = length(sim$y_test))
#' @export
simulate_spfc_continuous <- function(
    n = 100,
    p = 50,
    d = 1,
    s = 5,
    rho_x = 0.5,
    snr = 1,
    n_test = n,
    seed = NULL
){

  if(!is.null(seed)){
    set.seed(seed)
  }

  if(n < 2 || n_test < 1){
    stop("n must be at least 2 and n_test must be positive.")
  }

  if(!d %in% c(1L, 2L)){
    stop("The validated continuous thesis DGP currently supports d = 1 or d = 2.")
  }

  if(!is.finite(snr) || snr <= 0){
    stop("snr must be positive and finite.")
  }

  Sigma <- make_ar1_cov(
    p = p,
    rho = rho_x
  )

  dir <- generate_sparse_directions(
    p = p,
    d = d,
    s = s
  )

  # B is the true sparse sufficient-reduction basis. In the PFC inverse model
  # X|Y = mu + Gamma A f(Y) + epsilon, the reduction is span(Delta^{-1}Gamma).
  # Setting Gamma below proportional to Sigma B therefore makes that reduction
  # subspace exactly span(B), while retaining correlated residual noise.
  B <- dir$B
  G <- crossprod(B, Sigma %*% B)
  G <- (G + t(G)) / 2

  eg <- eigen(G, symmetric = TRUE)
  vals <- pmax(eg$values, .Machine$double.eps)
  G_inv_sqrt <- eg$vectors %*% diag(vals^(-0.5), nrow = d) %*% t(eg$vectors)

  Gamma_true <- Sigma %*% B %*% G_inv_sqrt

  # With these orthonormal response basis functions, Cov{f(Y)} = I_d in the
  # population. Multiplying the inverse mean by sqrt(snr) gives an exact SNR of
  # snr in the sufficient coordinates B'X:
  # Var[B'E(X|Y)] = snr * G and Var(B'epsilon) = G.
  make_f <- function(y){
    if(d == 1L){
      matrix(as.numeric(y), ncol = 1L)
    } else {
      cbind(
        y_power_1 = as.numeric(y),
        y_hermite_2 = (as.numeric(y)^2 - 1) / sqrt(2)
      )
    }
  }

  y <- stats::rnorm(n)
  y_test <- stats::rnorm(n_test)
  F_true <- make_f(y)
  F_test_true <- make_f(y_test)

  mean_X <- sqrt(snr) * F_true %*% t(Gamma_true)
  mean_X_test <- sqrt(snr) * F_test_true %*% t(Gamma_true)

  eps_all <- MASS::mvrnorm(
    n = n + n_test,
    mu = rep(0, p),
    Sigma = Sigma
  )

  epsilon <- eps_all[seq_len(n), , drop = FALSE]
  epsilon_test <- eps_all[n + seq_len(n_test), , drop = FALSE]

  X <- mean_X + epsilon
  X_test <- mean_X_test + epsilon_test

  # Exact population SNR in the sufficient coordinates, retained for audit.
  noise_cov_reduced <- G
  signal_cov_reduced <- snr * G
  snr_check <- sum(diag(signal_cov_reduced)) / sum(diag(noise_cov_reduced))

  return(list(
    X = X,
    y = as.numeric(y),
    X_test = X_test,
    y_test = as.numeric(y_test),
    F_true = F_true,
    F_test_true = F_test_true,
    inverse_mean = mean_X,
    inverse_mean_test = mean_X_test,
    Gamma_true = Gamma_true,
    B_true = B,
    active_set = dir$active_set,
    Sigma = Sigma,
    reduced_noise_covariance = noise_cov_reduced,
    reduced_signal_covariance = signal_cov_reduced,
    snr_check = as.numeric(snr_check),
    ytype = "continuous",
    n = n,
    n_test = n_test,
    p = p,
    d = d,
    s = s,
    rho_x = rho_x,
    snr = snr,
    dgp = if(d == 1L) "pfc_inverse_rank1" else "pfc_inverse_rank2"
  ))
}


#' Simulate Binary High-Dimensional Classification Data
#'
#' @param n Sample size.
#' @param p Number of predictors.
#' @param d Structural dimension.
#' @param s Number of active variables.
#' @param rho_x Predictor correlation.
#' @param signal_strength Signal strength in the logistic model.
#' @param seed Optional random seed.
#'
#' @return A list containing X, y, true directions, active set, and covariance.
#'
#' @examples
#' sim <- simulate_spfc_binary(
#'   n = 30, p = 8, d = 1, s = 3, rho_x = 0.3,
#'   signal_strength = 1, seed = 123
#' )
#' table(sim$y)
#' @export
simulate_spfc_binary <- function(
    n = 100,
    p = 50,
    d = 1,
    s = 5,
    rho_x = 0.5,
    signal_strength = 1,
    seed = NULL
){

  if(!is.null(seed)){
    set.seed(seed)
  }

  if(d != 1L){
    stop(
      "The current binary PFC simulation supports d = 1 only. ",
      "A binary response supplies a rank-one indicator response basis; ",
      "use a validated multiclass simulation before studying d > 1 classification."
    )
  }

  Sigma <- make_ar1_cov(
    p = p,
    rho = rho_x
  )

  X <- MASS::mvrnorm(
    n = n,
    mu = rep(0, p),
    Sigma = Sigma
  )

  dir <- generate_sparse_directions(
    p = p,
    d = d,
    s = s
  )

  B <- dir$B

  eta <- as.numeric(
    rowSums(X %*% B)
  )

  linpred <- signal_strength * eta

  prob <- 1 / (1 + exp(-linpred))

  y_num <- stats::rbinom(
    n = n,
    size = 1,
    prob = prob
  )

  y <- factor(
    y_num,
    levels = c(0, 1)
  )

  return(list(
    X = X,
    y = y,
    probability = prob,
    B_true = B,
    active_set = dir$active_set,
    Sigma = Sigma,
    ytype = "categorical",
    n = n,
    p = p,
    d = d,
    s = s,
    rho_x = rho_x,
    signal_strength = signal_strength
  ))
}


#' Compute Subspace Distance
#'
#' Computes a sine-angle distance between two subspaces.
#'
#' @param B_true True basis matrix.
#' @param B_hat Estimated basis matrix.
#'
#' @return Numeric subspace distance.
#'
#' @examples
#' B1 <- matrix(c(1, 0, 0), ncol = 1)
#' B2 <- matrix(c(0.9, 0.1, 0), ncol = 1)
#' subspace_distance(B1, B2)
#' @export
subspace_distance <- function(B_true, B_hat){

  B_true <- as.matrix(B_true)
  B_hat <- as.matrix(B_hat)

  Qt <- qr.Q(qr(B_true))
  Qh <- qr.Q(qr(B_hat))

  svals <- svd(
    t(Qt) %*% Qh,
    nu = 0,
    nv = 0
  )$d

  svals <- pmin(
    pmax(svals, 0),
    1
  )

  dist <- sqrt(
    sum(1 - svals^2)
  )

  return(as.numeric(dist))
}


#' Compute Variable-Selection Performance
#'
#' @param selected Integer vector of selected variable indices.
#' @param active_set Integer vector of truly active variable indices.
#' @param p Total number of variables.
#'
#' @return A list of precision, recall, F1, false positives, and false negatives.
#'
#' @examples
#' variable_selection_metrics(
#'   selected = c(1, 2, 5), active_set = c(1, 2, 3), p = 8
#' )
#' @export
variable_selection_metrics <- function(
    selected,
    active_set,
    p
){

  selected <- unique(as.integer(selected))
  active_set <- unique(as.integer(active_set))

  tp <- length(intersect(selected, active_set))
  fp <- length(setdiff(selected, active_set))
  fn <- length(setdiff(active_set, selected))
  tn <- p - tp - fp - fn

  precision <- ifelse((tp + fp) == 0, NA_real_, tp / (tp + fp))
  recall <- ifelse((tp + fn) == 0, NA_real_, tp / (tp + fn))

  f1 <- ifelse(
    is.na(precision) || is.na(recall) || (precision + recall) == 0,
    NA_real_,
    2 * precision * recall / (precision + recall)
  )

  return(list(
    tp = tp,
    fp = fp,
    fn = fn,
    tn = tn,
    precision = as.numeric(precision),
    recall = as.numeric(recall),
    f1 = as.numeric(f1),
    n_selected = length(selected)
  ))
}


#' Run One SPFC Simulation Replicate
#'
#' @param sim_data Simulated data object from `simulate_spfc_continuous()` or
#' `simulate_spfc_binary()`.
#' @param cov_method Covariance estimator.
#' @param d_fit Optional fitted structural dimension. When `NULL`, the dimension is selected over `d_grid` separately for each requested criterion.
#' @param d_grid Candidate structural dimensions used when `d_fit = NULL`.
#' @param criteria Character vector of information criteria used when `d_fit = NULL`.
#' @param variable_method Variable-selection method.
#' @param selection_rule Character screening rule: `"c1f"`, `"quantile"`, or `"fixed"`.
#' @param threshold Optional finite numeric screening threshold. It is required for `selection_rule = "fixed"`; for `"quantile"`, `NULL` uses `quantile_cut`; it is ignored by the C1F-calibrated rule.
#' @param quantile_cut Quantile cut for variable selection.
#' @param nslices Number of slices.
#' @param poly_degree Polynomial degree for continuous response.
#' @param gamma Ridge constant for SRE.
#' @param rho Convex shrinkage intensity for CSE.
#' @param verbose Logical.
#'
#' @return A one-row data frame with performance metrics.
#'
#' @noRd
run_spfc_simulation_replicate <- function(
    sim_data,
    cov_method = c("mec", "oas", "sre", "sde", "cse"),
    d_fit = NULL,
    d_grid = NULL,
    criteria = c("AIC", "BIC", "CAIC", "ICOMP_IFIM", "CICOMP"),
    variable_method = c("adaptive_weighted_l1", "c1f_extension"),
    selection_rule = c("c1f", "quantile", "fixed"),
    threshold = NULL,
    quantile_cut = 0.75,
    nslices = 5,
    poly_degree = 2,
    gamma = 0.1,
    rho = 0.5,
    verbose = FALSE
){

  cov_method <- match.arg(cov_method)
  variable_method <- match.arg(variable_method)
  selection_rule <- match.arg(selection_rule)

  X <- sim_data$X
  y <- sim_data$y

  allowed_criteria <- c("AIC", "BIC", "CAIC", "ICOMP_IFIM", "ICOMP_MISSPEC", "CICOMP")
  if(any(!criteria %in% allowed_criteria)){
    stop("Unknown criterion supplied: ", paste(setdiff(criteria, allowed_criteria), collapse = ", "))
  }

  # Legacy/known-d mode remains available when d_fit is supplied. The thesis
  # simulation uses d_fit = NULL and estimates d independently by each criterion.
  if(!is.null(d_fit)){
    d_grid <- d_fit
    criteria <- "KNOWN_D"
  } else {
    if(is.null(d_grid)){
      d_grid <- seq_len(min(5L, ncol(X), max(1L, nrow(X) - 2L)))
    }
    d_grid <- sort(unique(as.integer(d_grid[d_grid >= 1])))
  }

  t1 <- Sys.time()

  if(identical(criteria, "KNOWN_D")){
    fit <- spfc_fit(
      X = X, y = y, d = d_fit, ytype = sim_data$ytype,
      cov_method = cov_method, nslices = nslices,
      poly_degree = poly_degree, gamma = gamma, rho = rho,
      centre_x = TRUE, scale_x = FALSE
    )
    reduced_fit <- fit_reduced_model(
      Z = fit$Z, y = y, ytype = sim_data$ytype, classifier = "auto"
    )
    candidate <- list(KNOWN_D = list(
      criterion = "KNOWN_D", selected_d = d_fit,
      spfc_fit = fit, reduced_model = reduced_fit
    ))
  } else {
    dsel <- spfc_select_dimension(
      X = X,
      y = y,
      d_grid = d_grid,
      cov_method = cov_method,
      ytype = sim_data$ytype,
      nslices = nslices,
      poly_degree = poly_degree,
      gamma = gamma,
      rho = rho,
      centre_x = TRUE,
      scale_x = FALSE,
      verbose = verbose
    )

    candidate <- lapply(criteria, function(cr){
      sr <- dsel$selected[dsel$selected$criterion == cr, , drop = FALSE]
      if(nrow(sr) != 1L) stop("Could not select dimension for criterion: ", cr)
      d_hat <- sr$selected_d[1]
      fi <- dsel$fits[[paste0("d", d_hat)]]
      list(
        criterion = cr,
        selected_d = d_hat,
        spfc_fit = fi$spfc_fit,
        reduced_model = fi$reduced_model
      )
    })
    names(candidate) <- criteria
  }

  t2 <- Sys.time()
  runtime_total <- as.numeric(difftime(t2, t1, units = "secs"))

  rows <- lapply(candidate, function(can){
    fit <- can$spfc_fit
    reduced_fit <- can$reduced_model

    # Predictive performance must be evaluated out of sample whenever an
    # independent test sample is supplied by the simulation generator. The
    # SPFC basis, dimension choice, variable selection, and downstream model are
    # all estimated from training data only.
    if(!is.null(sim_data$X_test) && !is.null(sim_data$y_test)){
      Z_eval <- stats::predict(fit, newdata = sim_data$X_test)
      y_eval <- sim_data$y_test
      prediction_sample <- "independent_test"
      n_eval <- nrow(sim_data$X_test)
    } else {
      Z_eval <- fit$Z
      y_eval <- y
      prediction_sample <- "training"
      n_eval <- nrow(fit$Z)
    }

    eval <- evaluate_reduced_model(
      object = reduced_fit,
      Z = Z_eval,
      y = y_eval
    )

    vsel <- spfc_select_variables(
      fit = fit,
      method = variable_method,
      selection_rule = selection_rule,
      reduced_model = reduced_fit,
      threshold = threshold,
      quantile_cut = quantile_cut
    )

    selected_vars <- vsel$variable[vsel$selected]
    vmetrics <- variable_selection_metrics(
      selected = selected_vars,
      active_set = sim_data$active_set,
      p = sim_data$p
    )

    # Subspace distance is reported only when the structural dimension is
    # correctly recovered; dimension-recovery probability is reported
    # separately. This avoids disguising under/over-selection by truncating
    # unequal-dimensional subspaces.
    sdist <- if(can$selected_d == sim_data$d){
      subspace_distance(
        B_true = sim_data$B_true,
        B_hat = fit$V[, seq_len(sim_data$d), drop = FALSE]
      )
    } else {
      NA_real_
    }

    common <- list(
      ytype = sim_data$ytype,
      cov_method = cov_method,
      criterion = can$criterion,
      variable_method = variable_method,
      selection_rule = selection_rule,
      n = sim_data$n,
      p = sim_data$p,
      true_d = sim_data$d,
      fitted_d = can$selected_d,
      dimension_correct = as.integer(can$selected_d == sim_data$d),
      s = sim_data$s,
      rho_x = sim_data$rho_x,
      shrinkage_rho = fit$rho,
      c1f_l1 = attr(vsel, "C1F"),
      c1f_complexity_fraction = attr(vsel, "c1f_complexity_fraction"),
      c1f_loading_scale = attr(vsel, "loading_scale"),
      c1f_hd_factor = attr(vsel, "high_dimensional_factor"),
      c1f_complexity_multiplier = attr(vsel, "c1f_complexity_multiplier"),
      c1f_global_penalty = attr(vsel, "global_penalty"),
      runtime_sec = runtime_total / length(candidate),
      prediction_sample = prediction_sample,
      n_prediction = n_eval,
      subspace_distance = sdist,
      precision = vmetrics$precision,
      recall = vmetrics$recall,
      f1_variable = vmetrics$f1,
      tp = vmetrics$tp,
      fp = vmetrics$fp,
      fn = vmetrics$fn,
      tn = vmetrics$tn,
      n_selected = vmetrics$n_selected
    )

    if(sim_data$ytype == "continuous"){
      data.frame(
        common,
        snr = sim_data$snr,
        signal_strength = NA_real_,
        rmse = eval$rmse,
        mae = eval$mae,
        accuracy = NA_real_,
        sensitivity = NA_real_,
        specificity = NA_real_,
        f1_classification = NA_real_,
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        common,
        snr = NA_real_,
        signal_strength = sim_data$signal_strength,
        rmse = NA_real_,
        mae = NA_real_,
        accuracy = eval$accuracy,
        sensitivity = eval$sensitivity,
        specificity = eval$specificity,
        f1_classification = eval$f1,
        stringsAsFactors = FALSE
      )
    }
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out$error <- NA_character_
  return(out)
}


#' Run a Small SPFC Simulation Study
#'
#' @param response_type Character. `"continuous"` or `"categorical"`.
#' @param nrep Number of Monte Carlo replicates.
#' @param n Sample size.
#' @param p Number of predictors.
#' @param d Structural dimension.
#' @param s Number of active variables.
#' @param rho_x Predictor correlation.
#' @param snr Signal-to-noise ratio for continuous response.
#' @param signal_strength Signal strength for categorical response.
#' @param cov_methods Covariance estimators to compare.
#' @param variable_methods Variable-selection methods to compare.
#' @param seed Random seed.
#'
#' @return A data frame of simulation results.
#'
#' @examples
#' \donttest{
#' ans <- run_spfc_simulation(
#'   response_type = 'continuous', nrep = 1, n = 30, p = 8,
#'   d = 1, s = 3, rho_x = 0.3, snr = 2,
#'   cov_methods = c('oas', 'sre'),
#'   variable_methods = 'adaptive_weighted_l1', seed = 123
#' )
#' ans[, c('cov_method', 'rmse', 'subspace_distance')]
#' }
#' @export
run_spfc_simulation <- function(
    response_type = c("continuous", "categorical"),
    nrep = 10,
    n = 100,
    p = 50,
    d = 1,
    s = 5,
    rho_x = 0.5,
    snr = 1,
    signal_strength = 1,
    cov_methods = c("mec", "oas", "sre", "sde", "cse"),
    variable_methods = c("adaptive_weighted_l1", "c1f_extension"),
    seed = 123
){

  response_type <- match.arg(response_type)

  results <- list()
  counter <- 1

  for(rep in seq_len(nrep)){

    seed_rep <- seed + rep

    if(response_type == "continuous"){

      sim_data <- simulate_spfc_continuous(
        n = n,
        p = p,
        d = d,
        s = s,
        rho_x = rho_x,
        snr = snr,
        seed = seed_rep
      )

    } else {

      sim_data <- simulate_spfc_binary(
        n = n,
        p = p,
        d = d,
        s = s,
        rho_x = rho_x,
        signal_strength = signal_strength,
        seed = seed_rep
      )
    }

    for(cm in cov_methods){

      for(vm in variable_methods){

        results[[counter]] <- run_spfc_simulation_replicate(
          sim_data = sim_data,
          cov_method = cm,
          d_fit = d,
          variable_method = vm
        )

        results[[counter]]$replicate <- rep

        counter <- counter + 1
      }
    }
  }

  out <- do.call(rbind, results)

  rownames(out) <- NULL

  return(out)
}
