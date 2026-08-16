test_that("spfc_fit S3 methods work for continuous response", {

  set.seed(123)

  X_test <- matrix(
    rnorm(40 * 10),
    nrow = 40,
    ncol = 10
  )

  y_cont <- X_test[, 1] - 0.5 * X_test[, 2] + rnorm(40)

  fit_cont <- spfc_fit(
    X = X_test,
    y = y_cont,
    d = 2,
    ytype = "continuous",
    cov_method = "mec",
    nslices = 5,
    poly_degree = 2,
    centre_x = TRUE,
    scale_x = FALSE
  )

  expect_s3_class(fit_cont, "spfc_fit")

  expect_output(
    print(fit_cont),
    "Shrinkage Principal Fitted Components Fit"
  )

  fit_summary <- summary(fit_cont)

  expect_s3_class(fit_summary, "summary.spfc_fit")

  expect_equal(
    fit_summary$d,
    2
  )

  expect_equal(
    fit_summary$cov_method,
    "mec"
  )

  expect_true(
    is.matrix(coef(fit_cont))
  )

  expect_equal(
    dim(coef(fit_cont)),
    c(10, 2)
  )

  expect_true(
    is.matrix(fitted(fit_cont))
  )

  expect_equal(
    dim(fitted(fit_cont)),
    c(40, 2)
  )

  Z_pred <- predict(
    fit_cont,
    newdata = X_test
  )

  expect_true(
    is.matrix(Z_pred)
  )

  expect_equal(
    dim(Z_pred),
    c(40, 2)
  )

  expect_equal(
    as.matrix(Z_pred),
    as.matrix(fit_cont$Z),
    tolerance = 1e-8
  )

  expect_true(
    is.list(fit_cont$preprocessing)
  )

  expect_true(
    fit_cont$preprocessing$centred
  )

  expect_false(
    fit_cont$preprocessing$scaled
  )
})


test_that("spfc_fit predict rejects wrong number of columns", {

  set.seed(123)

  X_test <- matrix(
    rnorm(40 * 10),
    nrow = 40,
    ncol = 10
  )

  y_cont <- X_test[, 1] - 0.5 * X_test[, 2] + rnorm(40)

  fit_cont <- spfc_fit(
    X = X_test,
    y = y_cont,
    d = 2,
    ytype = "continuous",
    cov_method = "mec",
    nslices = 5,
    poly_degree = 2
  )

  X_bad <- matrix(
    rnorm(40 * 9),
    nrow = 40,
    ncol = 9
  )

  expect_error(
    predict(
      fit_cont,
      newdata = X_bad
    ),
    "expects 10 columns"
  )
})


test_that("spfc_fit S3 methods work for categorical response", {

  set.seed(123)

  X_test <- matrix(
    rnorm(40 * 10),
    nrow = 40,
    ncol = 10
  )

  y_cat <- factor(
    rep(c("A", "B"), each = 20)
  )

  fit_cat <- spfc_fit(
    X = X_test,
    y = y_cat,
    d = 1,
    ytype = "categorical",
    cov_method = "mec",
    centre_x = TRUE,
    scale_x = FALSE
  )

  expect_s3_class(fit_cat, "spfc_fit")

  expect_output(
    print(fit_cat),
    "Shrinkage Principal Fitted Components Fit"
  )

  fit_summary <- summary(fit_cat)

  expect_s3_class(fit_summary, "summary.spfc_fit")

  expect_equal(
    fit_summary$d,
    1
  )

  expect_equal(
    fit_summary$ytype,
    "categorical"
  )

  expect_true(
    is.matrix(coef(fit_cat))
  )

  expect_equal(
    dim(coef(fit_cat)),
    c(10, 1)
  )

  expect_true(
    is.matrix(fitted(fit_cat))
  )

  expect_equal(
    dim(fitted(fit_cat)),
    c(40, 1)
  )

  Z_pred <- predict(
    fit_cat,
    newdata = X_test
  )

  expect_equal(
    dim(Z_pred),
    c(40, 1)
  )

  expect_equal(
    as.matrix(Z_pred),
    as.matrix(fit_cat$Z),
    tolerance = 1e-8
  )
})
