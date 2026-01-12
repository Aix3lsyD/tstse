# Tests for co_tas and co_tas_boot functions

# =============================================================================
# co_tas tests
# =============================================================================

test_that("co_tas returns correct structure", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.5), n = 100)

  result <- co_tas(x, maxp = 5)

  expect_s3_class(result, "co_tas")
  expect_true("pvalue" %in% names(result))
  expect_true("tco" %in% names(result))
  expect_true("phi" %in% names(result))
  expect_true("n_a" %in% names(result))
  expect_true("p" %in% names(result))
})

test_that("co_tas p-value is between 0 and 1", {
  set.seed(456)
  x <- arima.sim(list(ar = 0.7), n = 100)

  result <- co_tas(x, maxp = 5)

  expect_gte(result$pvalue, 0)
  expect_lte(result$pvalue, 1)
})

test_that("co_tas detects trend in series with trend", {
  set.seed(789)
  n <- 200
  t <- 1:n
  z <- arima.sim(list(ar = 0.5), n = n)
  x <- 5 + 0.2 * t + z  # Strong trend

  result <- co_tas(x, maxp = 5)

  # Should detect trend (low p-value)
  expect_lt(result$pvalue, 0.05)
})

test_that("co_tas handles white noise", {
  set.seed(101)
  x <- rnorm(100)

  result <- co_tas(x, maxp = 5)

  expect_s3_class(result, "co_tas")
  expect_true(is.numeric(result$pvalue))
  expect_true(is.numeric(result$tco))
})

test_that("co_tas effective sample size is reasonable", {
  set.seed(202)
  n <- 100
  x <- arima.sim(list(ar = 0.9), n = n)  # High autocorrelation

  result <- co_tas(x, maxp = 5)

  # Effective sample size should be less than actual n for correlated data
  expect_lt(result$n_a, n)
  expect_gt(result$n_a, 0)
})

test_that("co_tas input validation works", {
  # Non-numeric input
  expect_error(co_tas("not numeric"), "must be a non-empty numeric vector")

  # Empty vector
  expect_error(co_tas(numeric(0)), "must be a non-empty numeric vector")

  # Invalid maxp
  expect_error(co_tas(rnorm(100), maxp = -1), "must be a positive integer")
  expect_error(co_tas(rnorm(100), maxp = 0), "must be a positive integer")
})

test_that("co_tas respects type argument", {
  set.seed(303)
  x <- arima.sim(list(ar = 0.6), n = 100)

  result_aic <- co_tas(x, maxp = 5, type = "aic")
  result_bic <- co_tas(x, maxp = 5, type = "bic")

  # Both should return valid results (may have same or different AR orders)
  expect_s3_class(result_aic, "co_tas")
  expect_s3_class(result_bic, "co_tas")
})

test_that("co_tas print method works", {
  set.seed(404)
  x <- arima.sim(list(ar = 0.5), n = 100)
  result <- co_tas(x, maxp = 5)

  expect_output(print(result), "Cochrane-Orcutt Trend Test")
  expect_output(print(result), "t-statistic:")
  expect_output(print(result), "p-value:")
  expect_output(print(result), "Effective sample size:")
})

# =============================================================================
# co_tas_boot tests
# =============================================================================

test_that("co_tas_boot returns correct structure", {
  set.seed(505)
  x <- arima.sim(list(ar = 0.5), n = 50)

  result <- co_tas_boot(x, nb = 19, maxp = 3, seed = 123)

  expect_s3_class(result, "co_tas_boot")
  expect_true("pvalue" %in% names(result))
  expect_true("tco" %in% names(result))
  expect_true("pvalue_asymptotic" %in% names(result))
  expect_true("phi" %in% names(result))
  expect_true("nb" %in% names(result))
})

test_that("co_tas_boot p-value is between 0 and 1", {
  set.seed(606)
  x <- arima.sim(list(ar = 0.5), n = 50)

  result <- co_tas_boot(x, nb = 19, seed = 42)

  expect_gte(result$pvalue, 0)
  expect_lte(result$pvalue, 1)
})

test_that("co_tas_boot is reproducible with seed", {
  set.seed(707)
  x <- arima.sim(list(ar = 0.5), n = 50)

  result1 <- co_tas_boot(x, nb = 19, seed = 999)
  result2 <- co_tas_boot(x, nb = 19, seed = 999)

  expect_equal(result1$pvalue, result2$pvalue)
  expect_equal(result1$tco, result2$tco)
})

test_that("co_tas_boot input validation works", {
  # Non-numeric input
  expect_error(co_tas_boot("not numeric"), "must be a non-empty numeric vector")

  # Invalid nb
  expect_error(co_tas_boot(rnorm(50), nb = -1), "must be a positive integer")
  expect_error(co_tas_boot(rnorm(50), nb = 0), "must be a positive integer")
})

test_that("co_tas_boot print method works", {
  set.seed(808)
  x <- arima.sim(list(ar = 0.5), n = 50)
  result <- co_tas_boot(x, nb = 19, seed = 42)

  expect_output(print(result), "Bootstrap Cochrane-Orcutt Trend Test")
  expect_output(print(result), "Bootstrap p-value:")
  expect_output(print(result), "Asymptotic p-value:")
  expect_output(print(result), "Bootstrap replicates:")
})

test_that("co_tas_boot detects trend with sufficient power", {
  skip_on_cran()  # Skip on CRAN due to computation time

  set.seed(909)
  n <- 100
  t <- 1:n
  z <- arima.sim(list(ar = 0.5), n = n)
  x <- 5 + 0.15 * t + z  # Moderate trend

  result <- co_tas_boot(x, nb = 99, seed = 42)

  # Both bootstrap and asymptotic should detect trend
  expect_lt(result$pvalue, 0.10)
})
