# Tests for co_tas_fast and co_tas_boot_fast functions

# =============================================================================
# co_tas_fast tests
# =============================================================================

test_that("co_tas_fast returns correct structure", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.5), n = 100)

  result <- co_tas_fast(x, maxp = 5)

  expect_s3_class(result, "co_tas")
  expect_true("pvalue" %in% names(result))
  expect_true("tco" %in% names(result))
  expect_true("phi" %in% names(result))
  expect_true("n_a" %in% names(result))
  expect_true("p" %in% names(result))
})

test_that("co_tas_fast p-value is between 0 and 1", {
  set.seed(456)
  x <- arima.sim(list(ar = 0.7), n = 100)

  result <- co_tas_fast(x, maxp = 5)

  expect_gte(result$pvalue, 0)
  expect_lte(result$pvalue, 1)
})

test_that("co_tas_fast detects trend in series with trend", {
  set.seed(789)
  n <- 200
  t <- 1:n
  z <- arima.sim(list(ar = 0.5), n = n)
  x <- 5 + 0.2 * t + z  # Strong trend

  result <- co_tas_fast(x, maxp = 5)

  # Should detect trend (low p-value)
  expect_lt(result$pvalue, 0.05)
})

test_that("co_tas_fast handles white noise", {
  set.seed(101)
  x <- rnorm(100)

  result <- co_tas_fast(x, maxp = 5)

  expect_s3_class(result, "co_tas")
  expect_true(is.numeric(result$pvalue))
  expect_true(is.numeric(result$tco))
})

test_that("co_tas_fast effective sample size is reasonable", {
  set.seed(202)
  n <- 100
  x <- arima.sim(list(ar = 0.9), n = n)  # High autocorrelation

  result <- co_tas_fast(x, maxp = 5)

  # Effective sample size should be less than actual n for correlated data
  expect_lt(result$n_a, n)
  expect_gt(result$n_a, 0)
})

test_that("co_tas_fast input validation works", {
  # Non-numeric input
  expect_error(co_tas_fast("not numeric"), "must be a non-empty numeric vector")

  # Empty vector
  expect_error(co_tas_fast(numeric(0)), "must be a non-empty numeric vector")

  # Invalid maxp
  expect_error(co_tas_fast(rnorm(100), maxp = -1), "must be a positive integer")
  expect_error(co_tas_fast(rnorm(100), maxp = 0), "must be a positive integer")

  # Series too short
  expect_error(co_tas_fast(rnorm(5)), "at least 10 observations")
})

test_that("co_tas_fast respects type argument", {
  set.seed(303)
  x <- arima.sim(list(ar = 0.6), n = 100)

  result_aic <- co_tas_fast(x, maxp = 5, type = "aic")
  result_bic <- co_tas_fast(x, maxp = 5, type = "bic")

  # Both should return valid results
  expect_s3_class(result_aic, "co_tas")
  expect_s3_class(result_bic, "co_tas")
})

test_that("co_tas_fast matches co_tas output", {
  set.seed(404)
  x <- arima.sim(list(ar = 0.6), n = 100)

  result_fast <- co_tas_fast(x, maxp = 5)
  result_orig <- co_tas(x, maxp = 5)

  # Results should be very close (allowing for minor numerical differences)
  expect_equal(result_fast$p, result_orig$p)
  expect_equal(result_fast$phi, result_orig$phi, tolerance = 1e-6)
  expect_equal(result_fast$n_a, result_orig$n_a, tolerance = 0.1)
  expect_equal(result_fast$tco, result_orig$tco, tolerance = 0.01)
  # p-values may differ slightly due to different approximation methods
  expect_equal(result_fast$pvalue, result_orig$pvalue, tolerance = 0.05)
})

# =============================================================================
# co_tas_boot_fast tests
# =============================================================================

test_that("co_tas_boot_fast returns correct structure", {
  set.seed(505)
  x <- arima.sim(list(ar = 0.5), n = 50)

  result <- co_tas_boot_fast(x, nb = 19, maxp = 3, seed = 123)

  expect_s3_class(result, "co_tas_boot")
  expect_true("pvalue" %in% names(result))
  expect_true("tco" %in% names(result))
  expect_true("pvalue_asymptotic" %in% names(result))
  expect_true("phi" %in% names(result))
  expect_true("nb" %in% names(result))
  expect_true("boot_pvals" %in% names(result))
  expect_true("boot_seeds" %in% names(result))
})

test_that("co_tas_boot_fast p-value is between 0 and 1", {
  set.seed(606)
  x <- arima.sim(list(ar = 0.5), n = 50)

  result <- co_tas_boot_fast(x, nb = 19, seed = 42)

  expect_gte(result$pvalue, 0)
  expect_lte(result$pvalue, 1)
})

test_that("co_tas_boot_fast is reproducible with seed", {
  set.seed(707)
  x <- arima.sim(list(ar = 0.5), n = 50)

  result1 <- co_tas_boot_fast(x, nb = 19, seed = 999)
  result2 <- co_tas_boot_fast(x, nb = 19, seed = 999)

  expect_equal(result1$pvalue, result2$pvalue)
  expect_equal(result1$tco, result2$tco)
  expect_equal(result1$boot_pvals, result2$boot_pvals)
})

test_that("co_tas_boot_fast input validation works", {
  # Non-numeric input
  expect_error(co_tas_boot_fast("not numeric"), "must be a non-empty numeric vector")

  # Invalid nb
  expect_error(co_tas_boot_fast(rnorm(50), nb = -1), "must be a positive integer")
  expect_error(co_tas_boot_fast(rnorm(50), nb = 0), "must be a positive integer")

  # maxp too large
  expect_error(co_tas_boot_fast(rnorm(100), maxp = 25), "must be <= 20")
})

test_that("co_tas_boot_fast has correct number of bootstrap samples", {
  set.seed(808)
  x <- arima.sim(list(ar = 0.5), n = 50)
  nb <- 37

  result <- co_tas_boot_fast(x, nb = nb, seed = 42)

  expect_equal(length(result$boot_pvals), nb)
  expect_equal(length(result$boot_seeds), nb)
  expect_equal(result$nb, nb)
})

test_that("co_tas_boot_fast uses print method from co_tas_boot", {
  set.seed(909)
  x <- arima.sim(list(ar = 0.5), n = 50)
  result <- co_tas_boot_fast(x, nb = 19, seed = 42)

  # Should use the same print method as co_tas_boot
  expect_output(print(result), "Bootstrap Cochrane-Orcutt Trend Test")
  expect_output(print(result), "Bootstrap p-value:")
  expect_output(print(result), "Asymptotic p-value:")
  expect_output(print(result), "Bootstrap replicates:")
})

test_that("co_tas_boot_fast detects trend with sufficient power", {
  skip_on_cran()  # Skip on CRAN due to computation time

  set.seed(1001)
  n <- 100
  t <- 1:n
  z <- arima.sim(list(ar = 0.5), n = n)
  x <- 5 + 0.15 * t + z  # Moderate trend

  result <- co_tas_boot_fast(x, nb = 99, seed = 42)

  # Both bootstrap and asymptotic should detect trend
  expect_lt(result$pvalue, 0.10)
})

test_that("co_tas_boot_fast handles AR(0) case", {
  set.seed(1100)
  # White noise - should result in AR(0) or very low order
  x <- rnorm(100)

  result <- co_tas_boot_fast(x, nb = 19, maxp = 5, seed = 42)

  expect_s3_class(result, "co_tas_boot")
  expect_true(is.numeric(result$pvalue))
  expect_true(is.numeric(result$tco))
})
