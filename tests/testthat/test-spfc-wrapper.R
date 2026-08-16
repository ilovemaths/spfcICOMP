test_that("spfc complete workflow works for continuous response with selected dimension", {

  set.seed(123)

  X_test <- matrix(
    rnorm(40 * 10),
    nrow = 40,
    ncol = 10
  )

  y_cont <- X_test[, 1] - 0.5 * X_test[, 2] + rnorm(40)

  fit_full <- spfc(
    X = X_test,
    y = y_cont,
    d = NULL,
    d_grid = 1:3,
    ytype = "continuous",
    cov_method = "mec",
    selector = "BIC",
    nslices = 5,
    poly_degree = 2,
    verbose = FALSE
  )

  expect_s3_class(
    fit_full,
    "spfc"
  )

  expect_equal(
    fit_full$ytype,
    "continuous"
  )

  expect_equal(
    fit_full$cov_method,
    "mec"
  )

  expect_equal(
    fit_full$selector,
    "BIC"
  )

  expect_true(
    fit_full$d %in% 1:3
  )

  expect_s3_class(
    fit_full$fit,
    "spfc_fit"
  )

  expect_true(
    is.list(fit_full$dimension_selection)
  )

  expect_true(
    is.list(fit_full$reduced_model)
  )

  expect_true(
    is.list(fit_full$evaluation)
  )

  expect_true(
    is.data.frame(fit_full$variable_selection)
  )

  expect_true(
    is.finite(fit_full$evaluation$rmse)
  )

  expect_true(
    is.finite(fit_full$evaluation$mae)
  )

  expect_equal(
    dim(coef(fit_full)),
    c(10, fit_full$d)
  )

  expect_equal(
    dim(fitted(fit_full)),
    c(40, fit_full$d)
  )

  Z_pred <- predict(
    fit_full,
    newdata = X_test
  )

  expect_equal(
    dim(Z_pred),
    c(40, fit_full$d)
  )

  expect_equal(
    as.matrix(Z_pred),
    as.matrix(fit_full$fit$Z),
    tolerance = 1e-8
  )

  fit_summary <- summary(fit_full)

  expect_s3_class(
    fit_summary,
    "summary.spfc"
  )

  expect_equal(
    fit_summary$ytype,
    "continuous"
  )

  expect_true(
    is.data.frame(fit_summary$evaluation)
  )
})


test_that("spfc complete workflow works for categorical response with selected dimension", {

  set.seed(123)

  X_test <- matrix(
    rnorm(40 * 10),
    nrow = 40,
    ncol = 10
  )

  y_cat <- factor(
    rep(c("A", "B"), each = 20)
  )

  fit_full <- spfc(
    X = X_test,
    y = y_cat,
    d = NULL,
    d_grid = 1:3,
    ytype = "categorical",
    cov_method = "mec",
    selector = "BIC",
    verbose = FALSE
  )

  expect_s3_class(
    fit_full,
    "spfc"
  )

  expect_equal(
    fit_full$ytype,
    "categorical"
  )

  expect_equal(
    fit_full$cov_method,
    "mec"
  )

  expect_true(
    fit_full$d %in% 1:3
  )

  expect_s3_class(
    fit_full$fit,
    "spfc_fit"
  )

  expect_true(
    is.list(fit_full$reduced_model)
  )

  expect_true(
    is.list(fit_full$evaluation)
  )

  expect_true(
    is.data.frame(fit_full$variable_selection)
  )

  expect_true(
    is.finite(fit_full$evaluation$accuracy)
  )

  expect_true(
    is.finite(fit_full$evaluation$f1)
  )

  expect_equal(
    dim(coef(fit_full)),
    c(10, fit_full$d)
  )

  expect_equal(
    dim(fitted(fit_full)),
    c(40, fit_full$d)
  )

  Z_pred <- predict(
    fit_full,
    newdata = X_test
  )

  expect_equal(
    dim(Z_pred),
    c(40, fit_full$d)
  )

  fit_summary <- summary(fit_full)

  expect_s3_class(
    fit_summary,
    "summary.spfc"
  )

  expect_equal(
    fit_summary$ytype,
    "categorical"
  )

  expect_true(
    is.data.frame(fit_summary$evaluation)
  )
})


test_that("spfc complete workflow works with fixed structural dimension", {

  set.seed(123)

  X_test <- matrix(
    rnorm(40 * 10),
    nrow = 40,
    ncol = 10
  )

  y_cont <- X_test[, 1] - 0.5 * X_test[, 2] + rnorm(40)

  fit_full <- spfc(
    X = X_test,
    y = y_cont,
    d = 2,
    ytype = "continuous",
    cov_method = "oas",
    verbose = FALSE
  )

  expect_s3_class(
    fit_full,
    "spfc"
  )

  expect_null(
    fit_full$dimension_selection
  )

  expect_equal(
    fit_full$d,
    2
  )

  expect_equal(
    fit_full$cov_method,
    "oas"
  )

  expect_equal(
    dim(fitted(fit_full)),
    c(40, 2)
  )
})


test_that("print and summary methods for spfc produce expected output", {

  set.seed(123)

  X_test <- matrix(
    rnorm(40 * 10),
    nrow = 40,
    ncol = 10
  )

  y_cont <- X_test[, 1] - 0.5 * X_test[, 2] + rnorm(40)

  fit_full <- spfc(
    X = X_test,
    y = y_cont,
    d = 1,
    ytype = "continuous",
    cov_method = "mec",
    verbose = FALSE
  )

  expect_output(
    print(fit_full),
    "Complete SPFC Workflow"
  )

  fit_summary <- summary(fit_full)

  expect_output(
    print(fit_summary),
    "Summary of Complete SPFC Workflow"
  )
})
