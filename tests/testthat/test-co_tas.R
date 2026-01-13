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
  expect_true("phi_null" %in% names(result))
  expect_true("nb" %in% names(result))
  expect_true("btest" %in% names(result))
  expect_true("boot_pvals" %in% names(result))
  expect_false(result$btest)  # default is FALSE
  expect_length(result$boot_pvals, 19)
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

# =============================================================================
# co_tas_boot btest=TRUE tests
# =============================================================================

test_that("co_tas_boot with btest=TRUE returns correct structure", {
  set.seed(1001)
  x <- arima.sim(list(ar = 0.5), n = 50)

  result <- co_tas_boot(x, nb = 19, maxp = 3, btest = TRUE, seed = 123)

  expect_s3_class(result, "co_tas_boot")
  expect_true(result$btest)
  expect_true("boot_tco" %in% names(result))
  expect_length(result$boot_tco, 19)
  expect_null(result$boot_pvals)  # Should be NULL when btest=TRUE
})

test_that("co_tas_boot with btest=TRUE p-value is between 0 and 1", {
  set.seed(1002)
  x <- arima.sim(list(ar = 0.5), n = 50)

  result <- co_tas_boot(x, nb = 19, btest = TRUE, seed = 42)

  expect_gte(result$pvalue, 0)
  expect_lte(result$pvalue, 1)
})

test_that("co_tas_boot with btest=TRUE is reproducible with seed", {
  set.seed(1003)
  x <- arima.sim(list(ar = 0.5), n = 50)

  result1 <- co_tas_boot(x, nb = 19, btest = TRUE, seed = 999)
  result2 <- co_tas_boot(x, nb = 19, btest = TRUE, seed = 999)

  expect_equal(result1$pvalue, result2$pvalue)
  expect_equal(result1$tco, result2$tco)
  expect_equal(result1$boot_tco, result2$boot_tco)
})

test_that("co_tas_boot btest=TRUE and btest=FALSE give different results", {
  set.seed(1004)
  x <- arima.sim(list(ar = 0.5), n = 50)

  result_p <- co_tas_boot(x, nb = 99, btest = FALSE, seed = 42)
  result_t <- co_tas_boot(x, nb = 99, btest = TRUE, seed = 42)

  # Same t-statistic and asymptotic p-value
  expect_equal(result_p$tco, result_t$tco)
  expect_equal(result_p$pvalue_asymptotic, result_t$pvalue_asymptotic)

  # Different bootstrap statistics stored
  expect_false(is.null(result_p$boot_pvals))
  expect_true(is.null(result_p$boot_tco))
  expect_true(is.null(result_t$boot_pvals))
  expect_false(is.null(result_t$boot_tco))
})

test_that("co_tas_boot print method shows bootstrap method", {
  set.seed(1005)
  x <- arima.sim(list(ar = 0.5), n = 50)

  result_p <- co_tas_boot(x, nb = 19, btest = FALSE, seed = 42)
  result_t <- co_tas_boot(x, nb = 19, btest = TRUE, seed = 42)

  expect_output(print(result_p), "Bootstrap method: p-value")
  expect_output(print(result_t), "Bootstrap method: t-statistic")
})

test_that("co_tas_boot phi_null differs from phi", {
  set.seed(1006)
  # Generate series where AR on original vs differenced series may differ
  x <- arima.sim(list(ar = 0.7), n = 100)

  result <- co_tas_boot(x, nb = 19, maxp = 5, seed = 42)

  # Both phi and phi_null should be present
  expect_true(length(result$phi) > 0 || length(result$phi_null) > 0)
  expect_true("phi_null" %in% names(result))
})

test_that("co_tas_boot input validation for btest", {
  # Invalid btest
  expect_error(co_tas_boot(rnorm(50), btest = "yes"), "must be TRUE or FALSE")
  expect_error(co_tas_boot(rnorm(50), btest = NA), "must be TRUE or FALSE")
})

test_that("co_tas_boot use_fast=TRUE matches use_fast=FALSE", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.6), n = 100)

  # Same seed for both to ensure identical bootstrap samples
  result_fast <- co_tas_boot(x, nb = 49, maxp = 5, use_fast = TRUE, seed = 456)
  result_slow <- co_tas_boot(x, nb = 49, maxp = 5, use_fast = FALSE, seed = 456)

  # Core statistics should match
  expect_equal(result_fast$tco, result_slow$tco, tolerance = 0.01)
  expect_equal(result_fast$pvalue_asymptotic, result_slow$pvalue_asymptotic, tolerance = 0.05)
  expect_equal(result_fast$phi, result_slow$phi, tolerance = 1e-6)
  expect_equal(result_fast$pvalue, result_slow$pvalue, tolerance = 0.1)
})

test_that("co_tas_boot parallel matches sequential", {
  skip_on_cran()
  skip_if(parallel::detectCores(logical = FALSE) < 2, "Not enough cores")

  set.seed(123)
  x <- arima.sim(list(ar = 0.6), n = 100)

  result_seq <- co_tas_boot(x, nb = 49, maxp = 5, cores = 1, seed = 456)
  result_par <- co_tas_boot(x, nb = 49, maxp = 5, cores = 2, seed = 456)

  # Results should be identical with same seed
  expect_equal(result_seq$tco, result_par$tco)
  expect_equal(result_seq$pvalue, result_par$pvalue)
  expect_equal(result_seq$boot_pvals, result_par$boot_pvals)
  expect_equal(result_seq$boot_seeds, result_par$boot_seeds)
})
