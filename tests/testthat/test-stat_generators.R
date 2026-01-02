# Tests for stat_generators (statistic factory functions)

# ==============================================================================
# Factory function return type tests
# ==============================================================================

test_that("make_stat_co returns a function", {
  stat_fn <- make_stat_co()
  expect_true(is.function(stat_fn))
})

test_that("make_stat_ols_t returns a function", {
  stat_fn <- make_stat_ols_t()
  expect_true(is.function(stat_fn))
})

test_that("make_stat_ols_slope returns a function", {
  stat_fn <- make_stat_ols_slope()
  expect_true(is.function(stat_fn))
})

test_that("make_stat_spearman returns a function", {
  stat_fn <- make_stat_spearman()
  expect_true(is.function(stat_fn))
})

test_that("make_stat_sen returns a function", {
  stat_fn <- make_stat_sen()
  expect_true(is.function(stat_fn))
})

test_that("make_stat_bn returns a function", {
  stat_fn <- make_stat_bn()
  expect_true(is.function(stat_fn))
})

test_that("make_stat_lr returns a function", {
  stat_fn <- make_stat_lr()
  expect_true(is.function(stat_fn))
})


# ==============================================================================
# Statistic output tests - returns single numeric
# ==============================================================================

test_that("make_stat_co returns single numeric", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 100)
  stat_fn <- make_stat_co()
  result <- stat_fn(x)

  expect_type(result, "double")
  expect_length(result, 1)
  expect_false(is.na(result))
})

test_that("make_stat_ols_t returns single numeric", {
  set.seed(123)
  x <- rnorm(100)
  stat_fn <- make_stat_ols_t()
  result <- stat_fn(x)

  expect_type(result, "double")
  expect_length(result, 1)
  expect_false(is.na(result))
})

test_that("make_stat_ols_slope returns single numeric", {
  set.seed(123)
  x <- rnorm(100)
  stat_fn <- make_stat_ols_slope()
  result <- stat_fn(x)

  expect_type(result, "double")
  expect_length(result, 1)
  expect_false(is.na(result))
})

test_that("make_stat_spearman returns single numeric", {
  set.seed(123)
  x <- rnorm(100)
  stat_fn <- make_stat_spearman()
  result <- stat_fn(x)

  expect_type(result, "double")
  expect_length(result, 1)
  expect_false(is.na(result))
})

test_that("make_stat_sen returns single numeric", {
  set.seed(123)
  x <- rnorm(100)
  stat_fn <- make_stat_sen()
  result <- stat_fn(x)

  expect_type(result, "double")
  expect_length(result, 1)
  expect_false(is.na(result))
})

test_that("make_stat_bn returns single numeric", {
  set.seed(123)
  x <- rnorm(100)
  stat_fn <- make_stat_bn()
  result <- stat_fn(x)

  expect_type(result, "double")
  expect_length(result, 1)
  expect_false(is.na(result))
})

test_that("make_stat_lr returns single numeric", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.5), n = 100)
  stat_fn <- make_stat_lr()
  result <- stat_fn(x)

  expect_type(result, "double")
  expect_length(result, 1)
  expect_false(is.na(result))
})


# ==============================================================================
# Trend detection tests - statistics should be larger for trending data
# ==============================================================================

test_that("make_stat_ols_t detects trend", {
  stat_fn <- make_stat_ols_t()

  # No trend
  set.seed(123)
  x_no_trend <- rnorm(200)
  t_no_trend <- stat_fn(x_no_trend)

  # Strong trend
  set.seed(123)
  x_trend <- 0.1 * (1:200) + rnorm(200)
  t_trend <- stat_fn(x_trend)

  # Trending series should have larger |t|
  expect_gt(abs(t_trend), abs(t_no_trend))
  expect_gt(abs(t_trend), 2)  # Should be significant
})

test_that("make_stat_ols_slope detects trend direction", {
  stat_fn <- make_stat_ols_slope()

  # Upward trend
  set.seed(123)
  x_up <- 0.1 * (1:100) + rnorm(100)
  slope_up <- stat_fn(x_up)
  expect_gt(slope_up, 0)

  # Downward trend
  set.seed(123)
  x_down <- -0.1 * (1:100) + rnorm(100)
  slope_down <- stat_fn(x_down)
  expect_lt(slope_down, 0)
})

test_that("make_stat_spearman detects monotonic trend", {
  stat_fn <- make_stat_spearman()

  # Strong upward trend
  set.seed(123)
  x_trend <- 0.2 * (1:100) + rnorm(100, sd = 0.5)
  rho <- stat_fn(x_trend)

  expect_gt(rho, 0.5)  # Should be strongly positive
})

test_that("make_stat_sen estimates slope correctly", {
  stat_fn <- make_stat_sen()

  # Known slope
  set.seed(123)
  true_slope <- 0.1
  x <- true_slope * (1:200) + rnorm(200, sd = 0.5)
  estimated_slope <- stat_fn(x)

  expect_equal(estimated_slope, true_slope, tolerance = 0.05)
})


# ==============================================================================
# Parameter passing tests
# ==============================================================================

test_that("make_stat_co respects ar_method parameter", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 100)

  stat_burg <- make_stat_co(ar_method = "burg")
  stat_mle <- make_stat_co(ar_method = "mle")

  # Both should return valid results (may differ slightly)
  result_burg <- stat_burg(x)
  result_mle <- stat_mle(x)

  expect_type(result_burg, "double")
  expect_type(result_mle, "double")
})

test_that("make_stat_co respects criterion parameter", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 100)

  stat_aic <- make_stat_co(criterion = "aic")
  stat_bic <- make_stat_co(criterion = "bic")

  # Both should return valid results
  result_aic <- stat_aic(x)
  result_bic <- stat_bic(x)

  expect_type(result_aic, "double")
  expect_type(result_bic, "double")
})

test_that("make_stat_bn respects order.max parameter", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 200)

  stat_fn1 <- make_stat_bn(order.max = 5)
  stat_fn2 <- make_stat_bn(order.max = 10)

  # Both should work
  result1 <- stat_fn1(x)
  result2 <- stat_fn2(x)

  expect_type(result1, "double")
  expect_type(result2, "double")
})

test_that("make_stat_lr respects order parameter", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 200)

  stat_ar1 <- make_stat_lr(order = 1)
  stat_ar2 <- make_stat_lr(order = 2)

  result1 <- stat_ar1(x)
  result2 <- stat_ar2(x)

  expect_type(result1, "double")
  expect_type(result2, "double")
})


# ==============================================================================
# Optional package tests
# ==============================================================================

test_that("make_stat_mk requires Kendall package", {
  skip_if_not_installed("Kendall")

  stat_fn <- make_stat_mk()
  expect_true(is.function(stat_fn))

  set.seed(123)
  x <- 0.1 * (1:100) + rnorm(100)
  result <- stat_fn(x)

  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("make_stat_hac requires sandwich and lmtest packages", {
  skip_if_not_installed("sandwich")
  skip_if_not_installed("lmtest")

  stat_fn <- make_stat_hac()
  expect_true(is.function(stat_fn))

  set.seed(123)
  x <- 0.1 * (1:100) + rnorm(100)
  result <- stat_fn(x)

  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("make_stat_gls requires nlme package", {
  skip_if_not_installed("nlme")

  stat_fn <- make_stat_gls()
  expect_true(is.function(stat_fn))

  set.seed(123)
  x <- 0.1 * (1:100) + arima.sim(list(ar = 0.5), n = 100)
  result <- stat_fn(x)

  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("make_stat_mk errors without Kendall package", {
  skip_if(requireNamespace("Kendall", quietly = TRUE))
  expect_error(make_stat_mk(), "Kendall")
})

test_that("make_stat_hac errors without sandwich package", {
  skip_if(requireNamespace("sandwich", quietly = TRUE))
  expect_error(make_stat_hac(), "sandwich")
})

test_that("make_stat_gls errors without nlme package", {
  skip_if(requireNamespace("nlme", quietly = TRUE))
  expect_error(make_stat_gls(), "nlme")
})


# ==============================================================================
# Sen's slope performance test (verify no O(n^2) memory issue)
# ==============================================================================

test_that("make_stat_sen handles moderately large n", {
  stat_fn <- make_stat_sen()

  # This should complete quickly without memory issues
  set.seed(123)
  x <- 0.01 * (1:500) + rnorm(500)

  # Should complete in reasonable time
  result <- stat_fn(x)
  expect_type(result, "double")
  expect_false(is.na(result))
})


# ==============================================================================
# Closure capture tests
# ==============================================================================

test_that("make_stat_co captures parameters correctly", {
  # Create two different stat functions
  stat5 <- make_stat_co(maxp = 5)
  stat3 <- make_stat_co(maxp = 3)

  # They should be different functions that use different maxp
  set.seed(123)
  x <- arima.sim(list(ar = c(0.5, 0.2, 0.1, 0.05)), n = 200)

  # Both should work
  r5 <- stat5(x)
  r3 <- stat3(x)

  expect_type(r5, "double")
  expect_type(r3, "double")
})


# ==============================================================================
# Edge case tests
# ==============================================================================

test_that("stat functions handle short series", {
  set.seed(123)
  x <- rnorm(20)

  # OLS-based should work
  expect_no_error(make_stat_ols_t()(x))
  expect_no_error(make_stat_ols_slope()(x))
  expect_no_error(make_stat_spearman()(x))
  expect_no_error(make_stat_sen()(x))
})

test_that("stat functions handle constant series", {
  x <- rep(5, 50)

  # OLS slope should be 0
  expect_equal(make_stat_ols_slope()(x), 0, ignore_attr = TRUE)

  # Spearman should be NA or 0 for constant
  # (depends on implementation handling of ties)
  result <- make_stat_spearman()(x)
  expect_true(is.na(result) || result == 0)
})
