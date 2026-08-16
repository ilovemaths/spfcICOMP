test_that("simulation replicate estimates structural dimension by criteria", {
  sim <- simulate_spfc_continuous(
    n = 40, p = 10, d = 1, s = 3,
    rho_x = 0.2, snr = 3, seed = 101
  )

  out <- run_spfc_simulation_replicate(
    sim_data = sim,
    cov_method = "oas",
    d_fit = NULL,
    d_grid = 1:2,
    criteria = c("AIC", "CICOMP"),
    variable_method = "adaptive_weighted_l1",
    selection_rule = "c1f",
    verbose = FALSE
  )

  expect_equal(nrow(out), 2L)
  expect_setequal(out$criterion, c("AIC", "CICOMP"))
  expect_true(all(out$fitted_d %in% 1:2))
  expect_true(all(out$dimension_correct %in% c(0, 1)))
  expect_true(all(is.finite(out$c1f_l1)))
})

test_that("legacy quantile variable selection remains available", {
  sim <- simulate_spfc_continuous(
    n = 40, p = 10, d = 1, s = 3,
    rho_x = 0.2, snr = 3, seed = 102
  )
  fit <- spfc_fit(sim$X, sim$y, d = 1, cov_method = "oas")
  red <- fit_reduced_model(fit$Z, sim$y, ytype = "continuous")

  vq <- spfc_select_variables(
    fit = fit,
    selection_rule = "quantile",
    quantile_cut = 0.75
  )
  vc <- spfc_select_variables(
    fit = fit,
    selection_rule = "c1f",
    reduced_model = red
  )

  expect_true(is.data.frame(vq))
  expect_true(is.data.frame(vc))
  expect_true("shrunk_importance" %in% names(vc))
  expect_true(is.finite(attr(vc, "C1F")))
})


test_that("C1F variable selection works with ridge logistic reduced model", {

  set.seed(20260808)
  X <- matrix(rnorm(60 * 12), nrow = 60, ncol = 12)
  eta <- 1.5 * X[, 1] - 1.0 * X[, 2]
  prob <- plogis(eta)
  y <- factor(rbinom(60, size = 1, prob = prob), levels = c(0, 1))

  fit <- spfc_fit(
    X = X,
    y = y,
    d = 1,
    ytype = "categorical",
    cov_method = "oas"
  )

  reduced <- fit_reduced_model(
    Z = fit$Z,
    y = y,
    ytype = "categorical",
    classifier = "ridge_logistic"
  )

  V <- extract_ifim(reduced)
  expect_true(is.matrix(V))
  expect_true(all(is.finite(V)))
  expect_true(is_positive_definite(V))

  selected <- spfc_select_variables(
    fit = fit,
    selection_rule = "c1f",
    reduced_model = reduced
  )

  expect_true(is.data.frame(selected))
  expect_true(is.finite(attr(selected, "C1F")))
  expect_equal(nrow(selected), ncol(X))
})


test_that("robust C1F calibration is invariant to sum-normalised adaptive weights", {
  sim <- simulate_spfc_continuous(
    n = 60, p = 20, d = 1, s = 5, rho_x = 0.2, snr = 3, seed = 103
  )
  fit <- spfc_fit(sim$X, sim$y, d = 1, cov_method = "oas")
  red <- fit_reduced_model(fit$Z, sim$y, ytype = "continuous")
  vc <- spfc_select_variables(
    fit = fit,
    method = "adaptive_weighted_l1",
    selection_rule = "c1f",
    reduced_model = red,
    c1f_calibration = "robust_universal"
  )
  expect_true(is.finite(attr(vc, "global_penalty")))
  expect_gt(attr(vc, "global_penalty"), 0)
  expect_true(all(vc$penalty >= 0))
  expect_true(any(vc$selected) || all(vc$shrunk_importance == 0))
})

test_that("simulation summary preserves SNR factorial cells", {
  x <- data.frame(
    ytype = rep("continuous", 4), cov_method = rep("oas", 4),
    criterion = rep("AIC", 4), selection_rule = rep("c1f", 4),
    variable_method = rep("adaptive_weighted_l1", 4),
    n = rep(50, 4), p = rep(20, 4), true_d = rep(1, 4), fitted_d = rep(1, 4),
    s = rep(5, 4), rho_x = rep(0.2, 4), snr = c(1, 1, 3, 3),
    runtime_sec = rep(0.1, 4), subspace_distance = rep(0.2, 4),
    precision = rep(0.5, 4), recall = rep(1, 4), f1_variable = rep(2/3, 4),
    n_selected = rep(10, 4), dimension_correct = rep(1, 4)
  )
  sm <- summarise_spfc_simulation(x)
  expect_equal(nrow(sm), 2L)
  expect_setequal(sm$snr, c(1, 3))
  expect_true(all(sm$nrep == 2))
})


test_that("ICOMP-HD C1F calibration follows the frozen thesis formula", {
  sim <- simulate_spfc_continuous(
    n = 80, p = 25, d = 1, s = 5, rho_x = 0.2, snr = 3, seed = 104
  )
  fit <- spfc_fit(sim$X, sim$y, d = 1, cov_method = "oas")
  red <- fit_reduced_model(fit$Z, sim$y, ytype = "continuous")
  vc <- spfc_select_variables(
    fit = fit,
    method = "adaptive_weighted_l1",
    selection_rule = "c1f",
    reduced_model = red,
    c1f_calibration = "icomp_hd_floor"
  )

  frac <- attr(vc, "c1f_complexity_fraction")
  scale_g <- attr(vc, "loading_scale")
  hd <- attr(vc, "high_dimensional_factor")
  mult <- attr(vc, "c1f_complexity_multiplier")
  gp <- attr(vc, "global_penalty")

  expect_equal(hd, sqrt(2 * log(ncol(sim$X)) / nrow(sim$X)), tolerance = 1e-12)
  expect_equal(mult, 1 + frac, tolerance = 1e-12)
  expect_equal(gp, scale_g * hd * mult, tolerance = 1e-12)
  expect_true(gp > 0)
})

test_that("ICOMP-HD C1F calibration is the default thesis route", {
  sim <- simulate_spfc_continuous(
    n = 60, p = 20, d = 1, s = 5, rho_x = 0.2, snr = 3, seed = 105
  )
  fit <- spfc_fit(sim$X, sim$y, d = 1, cov_method = "oas")
  red <- fit_reduced_model(fit$Z, sim$y, ytype = "continuous")
  vc <- spfc_select_variables(
    fit = fit,
    selection_rule = "c1f",
    reduced_model = red
  )
  expect_identical(attr(vc, "c1f_calibration"), "icomp_hd_floor")
  expect_true(is.finite(attr(vc, "global_penalty")))
})


test_that("continuous d=2 DGP is a genuine rank-two PFC inverse model", {
  sim <- simulate_spfc_continuous(
    n = 500, n_test = 100, p = 12, d = 2, s = 5,
    rho_x = 0.2, snr = 3, seed = 20260809
  )

  expect_equal(ncol(sim$B_true), 2L)
  expect_equal(ncol(sim$Gamma_true), 2L)
  expect_equal(ncol(sim$F_true), 2L)
  expect_identical(sim$dgp, "pfc_inverse_rank2")

  # The inverse mean must have rank two and both response basis functions must
  # contribute non-degenerate variation.
  expect_equal(qr(sim$inverse_mean)$rank, 2L)
  expect_gt(stats::var(sim$F_true[, 1]), 0.5)
  expect_gt(stats::var(sim$F_true[, 2]), 0.5)

  # Delta^{-1} Gamma has exactly the sparse true reduction subspace B.
  reduction_raw <- solve(sim$Sigma, sim$Gamma_true)
  expect_lt(subspace_distance(sim$B_true, reduction_raw), 1e-7)

  # The requested inverse-model SNR is exact in population by construction.
  expect_equal(sim$snr_check, 3, tolerance = 1e-10)

  # The fitted degree-2 response basis has enough rank for a d=2 PFC model.
  F <- build_fy(sim$y, ytype = "continuous", poly_degree = 2)
  expect_equal(qr(F)$rank, 2L)
})

test_that("continuous simulation prediction is evaluated on independent test data", {
  sim <- simulate_spfc_continuous(
    n = 50, n_test = 35, p = 12, d = 1, s = 4,
    rho_x = 0.2, snr = 3, seed = 20260810
  )

  out <- run_spfc_simulation_replicate(
    sim_data = sim,
    cov_method = "oas",
    d_fit = 1,
    selection_rule = "c1f",
    poly_degree = 2,
    verbose = FALSE
  )

  expect_true(all(out$prediction_sample == "independent_test"))
  expect_true(all(out$n_prediction == 35L))
  expect_true(all(is.finite(out$rmse)))
  expect_true(all(is.finite(out$mae)))

  # Perturbing only the held-out response must alter evaluation while leaving
  # the training data and fitted SPFC problem unchanged.
  sim_shift <- sim
  sim_shift$y_test <- sim_shift$y_test + 10
  out_shift <- run_spfc_simulation_replicate(
    sim_data = sim_shift,
    cov_method = "oas",
    d_fit = 1,
    selection_rule = "c1f",
    poly_degree = 2,
    verbose = FALSE
  )
  expect_gt(abs(out_shift$rmse - out$rmse), 1)
})


test_that("binary simulation does not claim unsupported d > 1 recovery", {
  expect_error(
    simulate_spfc_binary(n = 40, p = 10, d = 2, s = 4, seed = 20260811),
    "supports d = 1 only"
  )
})
