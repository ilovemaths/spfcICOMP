test_that("benchmark_spfc works for continuous response", {

  set.seed(123)

  X_test <- matrix(
    rnorm(40 * 10),
    nrow = 40,
    ncol = 10
  )

  y_cont <- X_test[, 1] - 0.5 * X_test[, 2] + rnorm(40)

  bench_cont <- benchmark_spfc(
    X = X_test,
    y = y_cont,
    d = 2,
    ytype = "continuous",
    methods = c("mec", "oas"),
    nslices = 5,
    poly_degree = 2,
    verbose = FALSE
  )

  expect_s3_class(
    bench_cont,
    "spfc_benchmark"
  )

  expect_true(
    is.data.frame(bench_cont$summary)
  )

  expect_equal(
    nrow(bench_cont$summary),
    2
  )

  expect_true(
    all(c("mec", "oas") %in% bench_cont$summary$cov_method)
  )

  expect_true(
    all(c("rmse", "mae", "runtime_sec", "n_selected") %in% names(bench_cont$summary))
  )

  expect_true(
    all(is.finite(bench_cont$summary$runtime_sec))
  )

  expect_true(
    all(is.na(bench_cont$summary$accuracy))
  )

  expect_true(
    all(!is.na(bench_cont$summary$rmse))
  )

  bench_summary <- summary(bench_cont)

  expect_s3_class(
    bench_summary,
    "summary.spfc_benchmark"
  )

  expect_equal(
    bench_summary$metric,
    "rmse"
  )

  expect_true(
    bench_summary$smaller_is_better
  )

  expect_true(
    bench_summary$best_method %in% c("mec", "oas")
  )
})


test_that("benchmark_spfc works for categorical response", {

  set.seed(123)

  X_test <- matrix(
    rnorm(40 * 10),
    nrow = 40,
    ncol = 10
  )

  y_cat <- factor(
    rep(c("A", "B"), each = 20)
  )

  bench_cat <- benchmark_spfc(
    X = X_test,
    y = y_cat,
    d = 1,
    ytype = "categorical",
    methods = c("mec", "oas"),
    verbose = FALSE
  )

  expect_s3_class(
    bench_cat,
    "spfc_benchmark"
  )

  expect_true(
    is.data.frame(bench_cat$summary)
  )

  expect_equal(
    nrow(bench_cat$summary),
    2
  )

  expect_true(
    all(c("mec", "oas") %in% bench_cat$summary$cov_method)
  )

  expect_true(
    all(c("accuracy", "sensitivity", "specificity", "f1") %in% names(bench_cat$summary))
  )

  expect_true(
    all(is.na(bench_cat$summary$rmse))
  )

  expect_true(
    all(!is.na(bench_cat$summary$accuracy))
  )

  bench_summary <- summary(bench_cat)

  expect_s3_class(
    bench_summary,
    "summary.spfc_benchmark"
  )

  expect_equal(
    bench_summary$metric,
    "accuracy"
  )

  expect_false(
    bench_summary$smaller_is_better
  )

  expect_true(
    bench_summary$best_method %in% c("mec", "oas")
  )
})


test_that("summary.spfc_benchmark rejects unknown metric", {

  set.seed(123)

  X_test <- matrix(
    rnorm(40 * 10),
    nrow = 40,
    ncol = 10
  )

  y_cont <- X_test[, 1] - 0.5 * X_test[, 2] + rnorm(40)

  bench_cont <- benchmark_spfc(
    X = X_test,
    y = y_cont,
    d = 2,
    ytype = "continuous",
    methods = c("mec", "oas"),
    verbose = FALSE
  )

  expect_error(
    summary(
      bench_cont,
      metric = "not_a_metric"
    ),
    "metric not found"
  )
})
