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

  expect_equal(
    names(bench_summary$display_table),
    c(
      "cov_method",
      "d",
      "reduced_model",
      "runtime_sec",
      "rmse",
      "mae",
      "n_selected"
    )
  )

  expect_false(
    anyNA(bench_summary$display_table)
  )

  expect_false(
    anyNA(bench_summary$covariance_details)
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

  expect_equal(
    names(bench_summary$display_table),
    c(
      "cov_method",
      "d",
      "reduced_model",
      "runtime_sec",
      "accuracy",
      "sensitivity",
      "specificity",
      "precision",
      "f1",
      "balanced_accuracy",
      "n_selected"
    )
  )

  expect_false(
    anyNA(bench_summary$covariance_details)
  )
})


test_that("benchmark display selects cross-validation metrics", {

  object <- structure(
    list(
      ytype = "continuous",
      validation = "cv",
      summary = data.frame(
        cov_method = c("mec", "sde"),
        d = c(1, 1),
        rho = c(0.4, NA_real_),
        nslices_used = c(5L, NA_integer_),
        reduced_model = c("lm", "lm"),
        runtime_sec = c(0.1, 0.2),
        rmse = c(NA_real_, NA_real_),
        mae = c(NA_real_, NA_real_),
        mean_rmse = c(0.5, 0.6),
        sd_rmse = c(0.05, 0.06),
        mean_mae = c(0.4, 0.5),
        sd_mae = c(0.04, 0.05),
        accuracy = c(NA_real_, NA_real_),
        n_selected = c(4L, 5L)
      )
    ),
    class = c("spfc_benchmark", "list")
  )

  display <- benchmark_display_table(object)
  details <- benchmark_covariance_details(object)

  expect_equal(
    names(display),
    c(
      "cov_method",
      "d",
      "reduced_model",
      "runtime_sec",
      "mean_rmse",
      "sd_rmse",
      "mean_mae",
      "sd_mae",
      "n_selected"
    )
  )

  expect_false(anyNA(display))
  expect_false(anyNA(details))
  expect_equal(nrow(details), 2L)
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
