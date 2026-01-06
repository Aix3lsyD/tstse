# ==============================================================================
# Tests for jt() - Jacob Turner Trend Test
# ==============================================================================

test_that("jt returns correct structure", {
  set.seed(42)
  x <- arima.sim(list(ar = 0.7), n = 200)
  result <- jt(x)

  # Check all expected components are present
  expected_names <- c("x", "z_x", "b0hat", "b1hat", "z_order", "z_phi",
                      "x_order", "x_phi", "vif", "rel_n", "est_df",
                      "pvalue", "tco")
  expect_true(all(expected_names %in% names(result)))

  # Check dimensions
  expect_length(result$x, 200)
  expect_length(result$z_x, 200)
})

test_that("jt VIF values are reasonable", {
  set.seed(42)
  x <- arima.sim(list(ar = 0.7), n = 200)
  result <- jt(x)

  # VIF should be >= 1 (more autocorrelation = higher VIF)
  expect_gte(result$vif, 1)

  # Effective n should be <= actual n - 2
  expect_lte(result$rel_n, 200 - 2)

  # df should be positive
  expect_gt(result$est_df, 0)

  # p-value should be in valid range
  expect_gte(result$pvalue, 0)
  expect_lte(result$pvalue, 1)
})

test_that("jt handles white noise correctly", {
  set.seed(42)
  x <- rnorm(200)
  result <- jt(x)

  # For white noise, VIF should be close to 1
  # (some sampling variability expected)
  expect_lt(abs(result$vif - 1), 0.5)

  # Effective n should be close to actual n - 2
  expect_gt(result$rel_n, 150)  # Should be near 198
})

test_that("jt detects trend correctly", {
  set.seed(42)
  n <- 200
  # Series with clear trend
  x_trend <- 0.1 * (1:n) + rnorm(n)
  result <- jt(x_trend)

  # Should detect significant trend
  expect_lt(result$pvalue, 0.05)
  expect_gt(abs(result$tco), 2)
})

test_that("jt handles no trend correctly", {
  set.seed(42)
  # Pure AR(1) without trend
  x <- arima.sim(list(ar = 0.5), n = 200)
  result <- jt(x)

  # Should not detect significant trend (most of the time)
  # Using a relaxed threshold since this is probabilistic
  expect_gt(result$pvalue, 0.01)
})

test_that("jt validates input correctly", {
  # Non-numeric input
  expect_error(jt("not numeric"), "`x` must be a non-empty numeric vector")

  # Empty vector
  expect_error(jt(numeric(0)), "`x` must be a non-empty numeric vector")

  # Invalid maxp
  expect_error(jt(rnorm(100), maxp = -1), "`maxp` must be a positive integer")
  expect_error(jt(rnorm(100), maxp = "five"), "`maxp` must be a positive integer")

  # Series too short
  expect_error(jt(rnorm(10), maxp = 5), "Series too short")
})

test_that("jt works with different type options", {
  set.seed(42)
  x <- arima.sim(list(ar = 0.7), n = 100)

  result_aic <- jt(x, type = "aic")
  result_bic <- jt(x, type = "bic")
  result_aicc <- jt(x, type = "aicc")

  # All should return valid results
  expect_true(is.numeric(result_aic$pvalue))
  expect_true(is.numeric(result_bic$pvalue))
  expect_true(is.numeric(result_aicc$pvalue))
})

test_that("jt AR orders are within bounds", {
  set.seed(42)
  x <- arima.sim(list(ar = c(0.5, 0.3)), n = 200)
  result <- jt(x, maxp = 10)

  # AR orders should be within bounds
 expect_gte(result$z_order, 0)
  expect_lte(result$z_order, 10)
  expect_gte(result$x_order, 0)
  expect_lte(result$x_order, 10)
})

test_that("jt phi coefficients match AR order", {
  set.seed(42)
  x <- arima.sim(list(ar = 0.8), n = 200)
  result <- jt(x)

  # Length of phi should match order
  if (result$z_order > 0) {
    expect_length(result$z_phi, result$z_order)
  }
  if (result$x_order > 0) {
    expect_length(result$x_phi, result$x_order)
  }
})

test_that("jt z_x transformation is correct", {
  set.seed(42)
  x <- arima.sim(list(ar = 0.7), n = 100)
  result <- jt(x)

  # z_x should start with x[1]
  expect_equal(result$z_x[1], x[1])

  # z_x should have same length as x
  expect_length(result$z_x, length(x))
})

test_that("jt highly autocorrelated series has high VIF", {
  set.seed(42)
  # Very persistent AR(1)
  x <- arima.sim(list(ar = 0.95), n = 300)
  result <- jt(x)

  # Should have high VIF due to strong autocorrelation
  expect_gt(result$vif, 5)

  # Effective n should be much lower than actual n
  expect_lt(result$rel_n, 100)
})
