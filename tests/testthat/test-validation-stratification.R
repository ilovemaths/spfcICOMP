test_that("numeric binary train-test splitting is stratified", {

  y <- c(
    rep(0, 30),
    rep(1, 10)
  )

  split <- create_train_test_split(
    y = y,
    test_fraction = 0.20,
    seed = 123,
    stratify = TRUE
  )

  test_counts <- table(
    factor(
      y[split$test],
      levels = c(0, 1)
    )
  )

  train_counts <- table(
    factor(
      y[split$train],
      levels = c(0, 1)
    )
  )

  expect_equal(
    as.integer(test_counts),
    c(6L, 2L)
  )

  expect_equal(
    as.integer(train_counts),
    c(24L, 8L)
  )
})


test_that("numeric binary cross-validation folding is stratified", {

  y <- c(
    rep(0, 12),
    rep(1, 6)
  )

  fold_id <- create_cv_folds(
    y = y,
    folds = 3,
    seed = 123,
    stratify = TRUE
  )

  counts <- table(
    factor(
      y,
      levels = c(0, 1)
    ),
    factor(
      fold_id,
      levels = 1:3
    )
  )

  expect_equal(
    as.integer(counts[1, ]),
    rep(4L, 3)
  )

  expect_equal(
    as.integer(counts[2, ]),
    rep(2L, 3)
  )
})
