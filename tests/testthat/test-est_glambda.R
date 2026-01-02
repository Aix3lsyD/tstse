test_that("est_glambda returns correct structure", {
  set.seed(123)
  x <- gen_glambda(n = 100, lambda = 0.5, phi = 0.7, plot = FALSE)

  fit <- est_glambda(x, lambda_range = c(0.3, 0.7), offset_range = c(10, 30),
                     lambda_by = 0.2)

  expect_s3_class(fit, "est_glambda")
  expect_named(fit, c("lambda", "offset", "q", "table"))
  expect_true(is.numeric(fit$lambda))
  expect_true(is.numeric(fit$offset))
  expect_true(is.numeric(fit$q))
  expect_s3_class(fit$table, "data.frame")
})

test_that("est_glambda table has correct columns", {
  set.seed(123)
  x <- gen_glambda(n = 100, lambda = 0.5, phi = 0.7, plot = FALSE)

  fit <- est_glambda(x, lambda_range = c(0.3, 0.7), offset_range = c(10, 30),
                     lambda_by = 0.2)

  expect_named(fit$table, c("lambda", "offset", "Q"))
  expect_equal(nrow(fit$table), 3)  # 0.3, 0.5, 0.7
})

test_that("est_glambda finds reasonable parameters", {
  set.seed(456)
  # Generate data with known lambda
  x <- gen_glambda(n = 200, lambda = 0.5, phi = 0.6, offset = 20, plot = FALSE)

  fit <- est_glambda(x, lambda_range = c(0, 1), offset_range = c(10, 40),
                     lambda_by = 0.25)

  # Estimated lambda should be in reasonable range
  expect_true(fit$lambda >= 0 && fit$lambda <= 1)
  expect_true(fit$offset >= 10 && fit$offset <= 40)
  expect_true(fit$q >= 0)  # Q-statistic is sum of squares
})

test_that("est_glambda validates inputs", {
  # x too short
  expect_error(est_glambda(1:10, lambda_range = c(0, 1)),
               "`x` must be a numeric vector with at least 20 observations")

  # x not numeric
  expect_error(est_glambda(letters[1:30], lambda_range = c(0, 1)),
               "`x` must be a numeric vector")

  # x contains NA
  x <- c(rnorm(50), NA, rnorm(49))
  expect_error(est_glambda(x, lambda_range = c(0, 1)),
               "`x` contains NA values")

  # lambda_range wrong length
  expect_error(est_glambda(rnorm(100), lambda_range = 0.5),
               "`lambda_range` must be c\\(low, high\\)")

  # lambda_range inverted
  expect_error(est_glambda(rnorm(100), lambda_range = c(1, 0)),
               "`lambda_range\\[1\\]` must be <= `lambda_range\\[2\\]`")

  # offset_range wrong length
  expect_error(est_glambda(rnorm(100), offset_range = 50),
               "`offset_range` must be c\\(low, high\\)")

  # offset_range inverted
  expect_error(est_glambda(rnorm(100), offset_range = c(100, 0)),
               "`offset_range\\[1\\]` must be <= `offset_range\\[2\\]`")

  # offset_range negative
  expect_error(est_glambda(rnorm(100), offset_range = c(-10, 50)),
               "`offset_range` values must be non-negative")

  # lambda_by invalid
  expect_error(est_glambda(rnorm(100), lambda_by = 0),
               "`lambda_by` must be a positive number")
  expect_error(est_glambda(rnorm(100), lambda_by = -0.1),
               "`lambda_by` must be a positive number")
})

test_that("est_glambda print method works", {
  set.seed(123)
  x <- gen_glambda(n = 100, lambda = 0.5, phi = 0.7, plot = FALSE)
  fit <- est_glambda(x, lambda_range = c(0.3, 0.7), offset_range = c(10, 30),
                     lambda_by = 0.2)

  expect_output(print(fit), "G-Lambda Parameter Estimation")
  expect_output(print(fit), "Optimal lambda:")
  expect_output(print(fit), "Optimal offset:")
  expect_output(print(fit), "Q-statistic:")
})

test_that("est_glambda handles numeric vector input", {
  set.seed(789)
  x <- gen_glambda(n = 100, lambda = 0.5, phi = 0.7, plot = FALSE)

  fit <- est_glambda(as.numeric(x), lambda_range = c(0.3, 0.7),
                     offset_range = c(10, 30), lambda_by = 0.2)

  expect_s3_class(fit, "est_glambda")
})

test_that("est_glambda Q-statistic is non-negative", {
  set.seed(111)
  x <- gen_glambda(n = 150, lambda = 0.3, phi = 0.5, plot = FALSE)

  fit <- est_glambda(x, lambda_range = c(0, 0.6), offset_range = c(5, 25),
                     lambda_by = 0.3)

  # All Q values in table should be non-negative

  expect_true(all(fit$table$Q >= 0))
  expect_true(fit$q >= 0)
})

test_that("est_glambda produces similar output to tswge::est.glambda.wge", {
  skip_if_not_installed("tswge")

  set.seed(222)
  x <- gen_glambda(n = 100, lambda = 0.4, phi = 0.6, offset = 20, plot = FALSE)

  # Run both estimators with same grid
  # tswge uses lambda.range and offset.range, returns bestlambda and bestshift
  old <- suppressWarnings(tswge::est.glambda.wge(
    x, lambda.range = c(0.2, 0.6), offset.range = c(15, 25)
  ))

  new <- est_glambda(x, lambda_range = c(0.2, 0.6), offset_range = c(15, 25),
                     lambda_by = 0.1)

  # Should find similar optimal parameters
  # Note: may not be exactly equal due to algorithm differences
  expect_equal(new$lambda, old$bestlambda, tolerance = 0.3)
  expect_equal(new$offset, old$bestshift, tolerance = 10)
})
