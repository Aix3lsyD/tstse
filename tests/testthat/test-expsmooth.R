# Tests for expsmooth()
# Simple exponential smoothing wrapper around HoltWinters

test_that("expsmooth matches tswge::expsmooth.wge alpha estimation", {
  skip_if_not_installed("tswge")

  set.seed(123)
  x <- cumsum(rnorm(100))

  result <- expsmooth(x, plot = FALSE)
  expected <- suppressWarnings(tswge::expsmooth.wge(x, plot = FALSE))

  # Alpha should match
  expect_equal(result$alpha, expected$alpha, tolerance = 1e-6)
})


test_that("expsmooth matches tswge::expsmooth.wge fitted values", {
  skip_if_not_installed("tswge")

  set.seed(456)
  x <- cumsum(rnorm(50))

  result <- expsmooth(x, alpha = 0.3, plot = FALSE)
  expected <- suppressWarnings(tswge::expsmooth.wge(x, alpha = 0.3, plot = FALSE))

  # Fitted values should match (tswge calls it 'u')
  expect_equal(result$fitted, expected$u, tolerance = 1e-10)
})


test_that("expsmooth matches tswge::expsmooth.wge with forecasting", {
  skip_if_not_installed("tswge")

  set.seed(789)
  x <- cumsum(rnorm(80))

  result <- expsmooth(x, alpha = 0.5, n_ahead = 10, plot = FALSE)
  expected <- suppressWarnings(tswge::expsmooth.wge(x, alpha = 0.5, n.ahead = 10, plot = FALSE))

  # Fitted + forecast should match tswge's u vector
  combined <- c(result$fitted, result$forecast)
  expect_equal(combined, expected$u, tolerance = 1e-10)
})


test_that("expsmooth returns correct structure", {
  set.seed(111)
  x <- rnorm(50)

  result <- expsmooth(x, plot = FALSE)

  expect_s3_class(result, "expsmooth")
  expect_named(result, c("alpha", "fitted", "forecast", "hw"))
  expect_true(is.numeric(result$alpha))
  expect_true(is.numeric(result$fitted))
  expect_null(result$forecast)  # No forecasting
  expect_s3_class(result$hw, "HoltWinters")
})


test_that("expsmooth returns correct structure with forecasting", {
  set.seed(222)
  x <- rnorm(50)

  result <- expsmooth(x, n_ahead = 5, plot = FALSE)

  expect_s3_class(result, "expsmooth")
  expect_length(result$fitted, 50)
  expect_length(result$forecast, 5)
})


test_that("expsmooth respects fixed alpha", {
  set.seed(333)
  x <- cumsum(rnorm(100))

  result <- expsmooth(x, alpha = 0.2, plot = FALSE)
  expect_equal(result$alpha, 0.2)

  result2 <- expsmooth(x, alpha = 0.8, plot = FALSE)
  expect_equal(result2$alpha, 0.8)
})


test_that("expsmooth estimates alpha when NULL", {
  set.seed(444)
  x <- cumsum(rnorm(100))

  result <- expsmooth(x, alpha = NULL, plot = FALSE)

  expect_true(result$alpha > 0)
  expect_true(result$alpha < 1)
})


test_that("expsmooth fitted values have correct length", {
  x <- rnorm(100)
  result <- expsmooth(x, plot = FALSE)
  expect_length(result$fitted, 100)

  x2 <- rnorm(50)
  result2 <- expsmooth(x2, plot = FALSE)
  expect_length(result2$fitted, 50)
})


test_that("expsmooth first fitted value equals first observation", {
  set.seed(555)
  x <- rnorm(50)
  result <- expsmooth(x, plot = FALSE)

  # First fitted value should be x[1] (convention)
  expect_equal(result$fitted[1], x[1])
})


test_that("expsmooth forecasts are constant (simple exponential smoothing)", {
  # For simple exponential smoothing, all forecasts equal the last smoothed value
  set.seed(666)
  x <- cumsum(rnorm(100))

  result <- expsmooth(x, n_ahead = 10, plot = FALSE)

  # All forecasts should be equal (flat forecast)
  expect_true(all(result$forecast == result$forecast[1]))

  # Forecast should equal the last level from HoltWinters
  last_level <- as.numeric(result$hw$coefficients["a"])
  expect_equal(result$forecast[1], last_level)
})


test_that("expsmooth handles ts objects", {
  x <- ts(rnorm(60), frequency = 12, start = c(2020, 1))
  result <- expsmooth(x, plot = FALSE)

  expect_length(result$fitted, 60)
  expect_s3_class(result, "expsmooth")
})


test_that("expsmooth validates x parameter", {
  expect_error(expsmooth("abc", plot = FALSE), "`x` must be a numeric")
  expect_error(expsmooth(1, plot = FALSE), "at least 2 observations")
  expect_error(expsmooth(c(1, NA, 3), plot = FALSE), "contains NA")
  expect_error(expsmooth(numeric(0), plot = FALSE), "at least 2 observations")
})


test_that("expsmooth validates alpha parameter", {
  x <- rnorm(50)
  expect_error(expsmooth(x, alpha = "high", plot = FALSE), "`alpha` must be")
  expect_error(expsmooth(x, alpha = c(0.3, 0.5), plot = FALSE), "`alpha` must be")
  expect_error(expsmooth(x, alpha = NA, plot = FALSE), "`alpha` must be")
  expect_error(expsmooth(x, alpha = 0, plot = FALSE), "between 0 and 1")
  expect_error(expsmooth(x, alpha = 1, plot = FALSE), "between 0 and 1")
  expect_error(expsmooth(x, alpha = -0.5, plot = FALSE), "between 0 and 1")
  expect_error(expsmooth(x, alpha = 1.5, plot = FALSE), "between 0 and 1")
})


test_that("expsmooth validates n_ahead parameter", {
  x <- rnorm(50)
  expect_error(expsmooth(x, n_ahead = "five", plot = FALSE), "`n_ahead` must be")
  expect_error(expsmooth(x, n_ahead = NA, plot = FALSE), "`n_ahead` must be")
  expect_error(expsmooth(x, n_ahead = c(5, 10), plot = FALSE), "`n_ahead` must be")
  expect_error(expsmooth(x, n_ahead = -5, plot = FALSE), "non-negative")
})


test_that("expsmooth validates plot parameter", {
  x <- rnorm(50)
  expect_error(expsmooth(x, plot = "yes"), "`plot` must be TRUE or FALSE")
  expect_error(expsmooth(x, plot = NA), "`plot` must be TRUE or FALSE")
  expect_error(expsmooth(x, plot = c(TRUE, FALSE)), "`plot` must be TRUE or FALSE")
})


test_that("expsmooth plot = FALSE produces no plot", {
  x <- rnorm(50)
  # This should not open any graphics device or error
  expect_silent(result <- expsmooth(x, plot = FALSE))
})


test_that("expsmooth plot = TRUE produces plot without error", {
  skip_if_not(capabilities("png"))

  x <- rnorm(50)

  # Use a null device to avoid opening graphics window
  png(tempfile())
  on.exit(dev.off())

  expect_no_error(expsmooth(x, plot = TRUE))
})


test_that("expsmooth plot with forecasting works", {
  skip_if_not(capabilities("png"))

  x <- rnorm(50)

  png(tempfile())
  on.exit(dev.off())

  expect_no_error(expsmooth(x, n_ahead = 10, plot = TRUE))
})


test_that("expsmooth print method works", {
  x <- rnorm(50)
  result <- expsmooth(x, plot = FALSE)

  expect_output(print(result), "Simple Exponential Smoothing")
  expect_output(print(result), "Alpha:")
  expect_output(print(result), "Observations: 50")
})


test_that("expsmooth print method shows forecast info", {
  x <- rnorm(50)
  result <- expsmooth(x, n_ahead = 10, plot = FALSE)

  expect_output(print(result), "Forecasts: 10 periods ahead")
})


test_that("expsmooth returns invisible result", {
  x <- rnorm(50)

  png(tempfile())
  on.exit(dev.off())

  # invisible() means assignment works but print doesn't happen automatically
  expect_invisible(expsmooth(x, plot = TRUE))
})


test_that("expsmooth handles short series", {
  # Minimum length is 2
  x <- c(1, 2)
  result <- expsmooth(x, plot = FALSE)
  expect_length(result$fitted, 2)

  x <- c(1, 2, 3)
  result <- expsmooth(x, plot = FALSE)
  expect_length(result$fitted, 3)
})


test_that("expsmooth higher alpha gives more responsive smoothing", {
  set.seed(777)
  # Create a series with a jump
  x <- c(rep(10, 20), rep(20, 20))

  result_low <- expsmooth(x, alpha = 0.1, plot = FALSE)
  result_high <- expsmooth(x, alpha = 0.9, plot = FALSE)

  # After the jump (at position 21), high alpha should adapt faster
  # So fitted[25] should be closer to 20 for high alpha
  expect_true(result_high$fitted[25] > result_low$fitted[25])
})


test_that("expsmooth n_ahead = 0 gives no forecasts", {
  x <- rnorm(50)
  result <- expsmooth(x, n_ahead = 0, plot = FALSE)
  expect_null(result$forecast)
})


test_that("expsmooth underlying HoltWinters object is accessible", {
  x <- rnorm(50)
  result <- expsmooth(x, alpha = 0.4, plot = FALSE)

  # Should be able to use the HoltWinters object directly
  expect_s3_class(result$hw, "HoltWinters")
  expect_false(result$hw$beta)
  expect_false(result$hw$gamma)
})
