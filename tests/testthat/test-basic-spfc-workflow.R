test_that("basic SPFC workflow runs for continuous response", {

  set.seed(123)

  X <- matrix(rnorm(40 * 10), nrow = 40, ncol = 10)
  y <- X[, 1] - 0.5 * X[, 2] + rnorm(40)

  fit <- spfc_fit(
    X = X,
    y = y,
    d = 1,
    ytype = "continuous",
    cov_method = "mec",
    nslices = 5
  )

  expect_equal(nrow(fit$Z), 40)
  expect_equal(ncol(fit$Z), 1)
  expect_equal(fit$ytype, "continuous")
  expect_equal(fit$cov_method, "mec")
})


test_that("basic SPFC workflow runs for categorical response", {

  set.seed(123)

  X <- matrix(rnorm(40 * 10), nrow = 40, ncol = 10)
  y <- factor(rep(c("A", "B"), each = 20))

  fit <- spfc_fit(
    X = X,
    y = y,
    d = 1,
    ytype = "categorical",
    cov_method = "mec"
  )

  mod <- fit_reduced_model(
    Z = fit$Z,
    y = y,
    ytype = "categorical"
  )

  eval <- evaluate_reduced_model(
    object = mod,
    Z = fit$Z,
    y = y
  )

  expect_equal(nrow(fit$Z), 40)
  expect_equal(ncol(fit$Z), 1)
  expect_equal(fit$ytype, "categorical")
  expect_true(is.numeric(eval$accuracy))
})
