#' Create an SPFC Simulation Design Grid
#'
#' Creates a factorial simulation design grid for SPFC simulation studies.
#'
#' @param response_type Character vector. One or both of `"continuous"` and `"categorical"`.
#' @param n_values Sample sizes.
#' @param p_values Predictor dimensions.
#' @param d_values True structural dimensions.
#' @param s_values Sparsity levels.
#' @param rho_x_values Predictor correlation values.
#' @param snr_values Signal-to-noise ratios for continuous responses.
#' @param signal_strength_values Signal strengths for categorical responses.
#' @param cov_methods Covariance estimators.
#' @param variable_methods Variable-selection methods.
#' @param nrep Number of Monte Carlo replicates per design setting.
#'
#' @return A data frame containing the full simulation design.
#'
#' @examples
#' create_spfc_design_grid(
#'   response_type = 'continuous', n_values = 30, p_values = 8,
#'   d_values = 1, s_values = 3, rho_x_values = 0.3,
#'   snr_values = 2, cov_methods = 'oas',
#'   variable_methods = 'adaptive_weighted_l1', nrep = 1
#' )
#' @export
create_spfc_design_grid <- function(
    response_type = c("continuous", "categorical"),
    n_values = c(50, 100),
    p_values = c(20, 50),
    d_values = c(1, 2),
    s_values = c(5),
    rho_x_values = c(0.2, 0.5, 0.8),
    snr_values = c(1, 2),
    signal_strength_values = c(0.8, 1.2),
    cov_methods = c("mec", "oas", "sre", "sde", "cse"),
    variable_methods = c("adaptive_weighted_l1", "c1f_extension"),
    nrep = 10
){

  response_type <- match.arg(
    response_type,
    choices = c("continuous", "categorical"),
    several.ok = TRUE
  )

  design_list <- list()
  counter <- 1

  if("continuous" %in% response_type){

    design_list[[counter]] <- expand.grid(
      response_type = "continuous",
      n = n_values,
      p = p_values,
      d = d_values,
      s = s_values,
      rho_x = rho_x_values,
      snr = snr_values,
      signal_strength = NA_real_,
      cov_method = cov_methods,
      variable_method = variable_methods,
      replicate = seq_len(nrep),
      stringsAsFactors = FALSE
    )

    counter <- counter + 1
  }

  if("categorical" %in% response_type){

    design_list[[counter]] <- expand.grid(
      response_type = "categorical",
      n = n_values,
      p = p_values,
      d = d_values,
      s = s_values,
      rho_x = rho_x_values,
      snr = NA_real_,
      signal_strength = signal_strength_values,
      cov_method = cov_methods,
      variable_method = variable_methods,
      replicate = seq_len(nrep),
      stringsAsFactors = FALSE
    )
  }

  design <- do.call(rbind, design_list)

  design$design_id <- seq_len(nrow(design))

  rownames(design) <- NULL

  return(design)
}


#' Run One Row of an SPFC Simulation Design
#'
#' Runs one simulation design row produced by `create_spfc_design_grid()`.
#'
#' @param design_row One-row data frame from the simulation design grid.
#' @param base_seed Integer base seed.
#' @param d_grid Candidate structural dimensions evaluated when the dimension is selected.
#' @param criteria Character vector of information criteria used for structural-dimension selection.
#' @param selection_rule Character screening rule passed to `spfc_select_variables()`.
#' @param nslices Number of response slices for continuous responses.
#' @param poly_degree Polynomial degree for continuous response basis.
#' @param gamma Ridge constant for SRE.
#' @param rho Convex shrinkage intensity for CSE.
#'
#' @return A one-row data frame of simulation results.
#'
#' @noRd
run_spfc_design_row <- function(
    design_row,
    base_seed = 123,
    d_grid = 1:5,
    criteria = c("AIC", "BIC", "CAIC", "ICOMP_IFIM", "CICOMP"),
    selection_rule = "c1f",
    nslices = 5,
    poly_degree = 2,
    gamma = 0.1,
    rho = 0.5
){

  if(nrow(design_row) != 1){
    stop("design_row must contain exactly one row.")
  }

  seed_i <- base_seed + design_row$design_id

  if(design_row$response_type == "continuous"){
    sim_data <- simulate_spfc_continuous(
      n = design_row$n, p = design_row$p, d = design_row$d,
      s = design_row$s, rho_x = design_row$rho_x,
      snr = design_row$snr, seed = seed_i
    )
  } else {
    sim_data <- simulate_spfc_binary(
      n = design_row$n, p = design_row$p, d = design_row$d,
      s = design_row$s, rho_x = design_row$rho_x,
      signal_strength = design_row$signal_strength, seed = seed_i
    )
  }

  res <- run_spfc_simulation_replicate(
    sim_data = sim_data,
    cov_method = design_row$cov_method,
    d_fit = NULL,
    d_grid = d_grid,
    criteria = criteria,
    variable_method = design_row$variable_method,
    selection_rule = selection_rule,
    nslices = nslices,
    poly_degree = poly_degree,
    gamma = gamma,
    rho = rho
  )

  res$design_id <- design_row$design_id
  res$replicate <- design_row$replicate
  return(res)
}


#' Run an SPFC Simulation Design
#'
#' Runs an SPFC simulation design grid sequentially.
#'
#' @param design Data frame returned by `create_spfc_design_grid()`.
#' @param base_seed Integer base seed.
#' @param d_grid Candidate structural dimensions evaluated for each simulation design row.
#' @param criteria Character vector of information criteria used for structural-dimension selection.
#' @param selection_rule Character screening rule passed to `run_spfc_design_row()`.
#' @param nslices Number of slices.
#' @param poly_degree Polynomial degree.
#' @param gamma Ridge constant for SRE.
#' @param rho Convex shrinkage intensity for CSE.
#' @param save_every Optional integer. If supplied, intermediate results are saved every `save_every` rows.
#' @param save_path Optional path for intermediate RDS output.
#' @param verbose Logical. If `TRUE`, print progress.
#'
#' @return A data frame of simulation results.
#'
#' @examples
#' \donttest{
#' design <- create_spfc_design_grid(
#'   response_type = 'continuous', n_values = 30, p_values = 8,
#'   d_values = 1, s_values = 3, rho_x_values = 0.3,
#'   snr_values = 2, cov_methods = 'oas',
#'   variable_methods = 'adaptive_weighted_l1', nrep = 1
#' )
#' ans <- run_spfc_design(
#'   design, base_seed = 123, d_grid = 1, criteria = 'AIC',
#'   selection_rule = 'c1f', poly_degree = 2, verbose = FALSE
#' )
#' ans[, c('cov_method', 'criterion', 'fitted_d', 'rmse')]
#' }
#' @export
run_spfc_design <- function(
    design,
    base_seed = 123,
    d_grid = 1:5,
    criteria = c("AIC", "BIC", "CAIC", "ICOMP_IFIM", "CICOMP"),
    selection_rule = "c1f",
    nslices = 5,
    poly_degree = 2,
    gamma = 0.1,
    rho = 0.5,
    save_every = NULL,
    save_path = "spfc_simulation_partial.rds",
    verbose = TRUE
){

  if(!is.data.frame(design)){
    stop("design must be a data frame.")
  }

  results <- vector("list", nrow(design))

  for(i in seq_len(nrow(design))){

    if(verbose){
      message(
        "Running design row ",
        i,
        " of ",
        nrow(design),
        " | method = ",
        design$cov_method[i],
        " | response = ",
        design$response_type[i],
        " | n = ",
        design$n[i],
        " | p = ",
        design$p[i]
      )
    }

    results[[i]] <- tryCatch(
      run_spfc_design_row(
        design_row = design[i, , drop = FALSE],
        base_seed = base_seed,
        d_grid = d_grid,
        criteria = criteria,
        selection_rule = selection_rule,
        nslices = nslices,
        poly_degree = poly_degree,
        gamma = gamma,
        rho = rho
      ),
      error = function(e){

        data.frame(
          ytype = design$response_type[i],
          cov_method = design$cov_method[i],
          variable_method = design$variable_method[i],
          n = design$n[i],
          p = design$p[i],
          true_d = design$d[i],
          criterion = NA_character_,
          selection_rule = selection_rule,
          fitted_d = NA_integer_,
          dimension_correct = NA_integer_,
          s = design$s[i],
          rho_x = design$rho_x[i],
          snr = design$snr[i],
          signal_strength = design$signal_strength[i],
          shrinkage_rho = NA_real_,
          c1f_l1 = NA_real_,
          c1f_complexity_fraction = NA_real_,
          c1f_loading_scale = NA_real_,
          c1f_hd_factor = NA_real_,
          c1f_complexity_multiplier = NA_real_,
          c1f_global_penalty = NA_real_,
          runtime_sec = NA_real_,
          prediction_sample = NA_character_,
          n_prediction = NA_integer_,
          subspace_distance = NA_real_,
          rmse = NA_real_,
          mae = NA_real_,
          accuracy = NA_real_,
          sensitivity = NA_real_,
          specificity = NA_real_,
          f1_classification = NA_real_,
          precision = NA_real_,
          recall = NA_real_,
          f1_variable = NA_real_,
          tp = NA_real_, fp = NA_real_, fn = NA_real_, tn = NA_real_,
          n_selected = NA_real_,
          design_id = design$design_id[i],
          replicate = design$replicate[i],
          error = conditionMessage(e),
          stringsAsFactors = FALSE
        )
      }
    )

    if(!is.null(save_every)){
      if(i %% save_every == 0){
        partial <- do.call(rbind, results[seq_len(i)])
        saveRDS(partial, save_path)
      }
    }
  }

  out <- do.call(rbind, results)

  rownames(out) <- NULL

  return(out)
}


#' Create a Small Pilot Design for SPFC
#'
#' Creates a small simulation design suitable for checking that the simulation
#' machinery works before running a larger study.
#'
#' @return A small design data frame.
#'
#' @noRd
create_spfc_pilot_design <- function(){

  create_spfc_design_grid(
    response_type = c("continuous", "categorical"),
    n_values = c(50),
    p_values = c(20),
    d_values = c(1),
    s_values = c(5),
    rho_x_values = c(0.5),
    snr_values = c(2),
    signal_strength_values = c(1),
    cov_methods = c("mec", "oas"),
    variable_methods = c("adaptive_weighted_l1"),
    nrep = 2
  )
}
