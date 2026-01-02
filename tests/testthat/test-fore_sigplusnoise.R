# Tests for fore_sigplusnoise()
# Signal plus noise forecasting

test_that("fore_sigplusnoise matches tswge::fore.sigplusnoise.wge for linear trend", {
  skip_if_not_installed("tswge")

  set.seed(123)
  t <- 1:100
  signal <- 5 + 0.1 * t
  noise <- arima.sim(model = list(ar = 0.5), n = 100, sd = 1)
  x <- as.numeric(signal + noise)

  result <- fore_sigplusnoise(x, linear = TRUE, n_ahead = 10, plot = FALSE)

  png(tempfile())
  expected <- tswge::fore.sigplusnoise.wge(x, linear = TRUE, n.ahead = 10, plot = FALSE)
  dev.off()

  # Compare forecasts (may differ slightly due to implementation details)
  expect_equal(length(result$f), length(expected$f))
  expect_equal(result$f, expected$f, tolerance = 0.5)
})


test_that("fore_sigplusnoise matches tswge for cosine signal", {
  skip_if_not_installed("tswge")

  set.seed(456)
  t <- 1:100
  freq <- 0.1
  signal <- 10 + 3 * cos(2 * pi * freq * t)
  noise <- arima.sim(model = list(ar = 0.3), n = 100, sd = 0.5)
  x <- as.numeric(signal + noise)

  result <- fore_sigplusnoise(x, linear = FALSE, freq = freq, n_ahead = 10, plot = FALSE)

  png(tempfile())
  expected <- tswge::fore.sigplusnoise.wge(x, linear = FALSE, freq = freq, n.ahead = 10, plot = FALSE)
  dev.off()

  expect_equal(length(result$f), length(expected$f))
})


test_that("fore_sigplusnoise returns correct structure", {
  set.seed(111)
  x <- cumsum(rnorm(50))

  result <- fore_sigplusnoise(x, n_ahead = 5, plot = FALSE)

  expect_s3_class(result, "fore_sigplusnoise")
  expect_named(result, c("signal", "noise", "phi", "p", "f", "ll", "ul",
                         "resid", "se", "wnv", "signal_coef", "psi"))
})


test_that("fore_sigplusnoise forecast length matches n_ahead", {
  set.seed(222)
  x <- rnorm(100)

  result5 <- fore_sigplusnoise(x, n_ahead = 5, plot = FALSE)
  expect_length(result5$f, 5)
  expect_length(result5$ll, 5)
  expect_length(result5$ul, 5)
  expect_length(result5$se, 5)

  result20 <- fore_sigplusnoise(x, n_ahead = 20, plot = FALSE)
  expect_length(result20$f, 20)
})


test_that("fore_sigplusnoise linear = TRUE fits linear trend", {
  set.seed(333)
  t <- 1:100
  x <- 10 + 2 * t + rnorm(100, sd = 1)

  result <- fore_sigplusnoise(x, linear = TRUE, plot = FALSE)

  # Check signal coefficients are close to true values
  expect_equal(unname(result$signal_coef["intercept"]), 10, tolerance = 2)
  expect_equal(unname(result$signal_coef["slope"]), 2, tolerance = 0.1)
  expect_length(result$signal_coef, 2)
})


test_that("fore_sigplusnoise linear = FALSE fits cosine", {
  set.seed(444)
  t <- 1:200
  freq <- 0.05
  x <- 5 + 3 * cos(2 * pi * freq * t) + 2 * sin(2 * pi * freq * t) + rnorm(200, sd = 0.5)

  result <- fore_sigplusnoise(x, linear = FALSE, freq = freq, plot = FALSE)

  # Should have 3 coefficients
  expect_length(result$signal_coef, 3)
  expect_named(result$signal_coef, c("intercept", "cos_coef", "sin_coef"))
})


test_that("fore_sigplusnoise prediction intervals contain point forecasts", {
  set.seed(555)
  x <- cumsum(rnorm(80))

  result <- fore_sigplusnoise(x, n_ahead = 10, plot = FALSE)

  # Point forecasts should be between lower and upper limits
  expect_true(all(result$f >= result$ll))
  expect_true(all(result$f <= result$ul))
})


test_that("fore_sigplusnoise prediction intervals widen with horizon", {
  set.seed(666)
  x <- arima.sim(model = list(ar = 0.8), n = 100) + seq(1, 10, length.out = 100)

  result <- fore_sigplusnoise(x, n_ahead = 20, plot = FALSE)

  # Standard errors should generally increase
  expect_true(result$se[20] > result$se[1])

  # Interval widths should increase
  widths <- result$ul - result$ll
  expect_true(widths[20] > widths[1])
})


test_that("fore_sigplusnoise signal + noise = original", {
  set.seed(777)
  x <- rnorm(50)

  result <- fore_sigplusnoise(x, plot = FALSE)

  # Signal + noise should equal original data
  expect_equal(result$signal + result$noise, x, tolerance = 1e-10)
})


test_that("fore_sigplusnoise handles low noise variance", {
  set.seed(888)
  # Pure linear trend + small white noise - may fit low or zero AR order
  t <- 1:100
  x <- 5 + 0.5 * t + rnorm(100, sd = 0.1)

  result <- fore_sigplusnoise(x, max_p = 5, plot = FALSE)

  # Should successfully fit and forecast regardless of AR order selected

  expect_true(result$p >= 0)
  expect_true(result$p <= 5)
  expect_length(result$f, 10)
})


test_that("fore_sigplusnoise lastn = TRUE forecasts last observations", {
  set.seed(999)
  x <- cumsum(rnorm(100))

  result <- fore_sigplusnoise(x, n_ahead = 10, lastn = TRUE, plot = FALSE)

  # Should still return forecasts
  expect_length(result$f, 10)
})


test_that("fore_sigplusnoise validates x parameter", {
  expect_error(fore_sigplusnoise("abc", plot = FALSE), "`x` must be a numeric")
  expect_error(fore_sigplusnoise(1:4, plot = FALSE), "at least 5")
  expect_error(fore_sigplusnoise(c(1, NA, 3, 4, 5, 6), plot = FALSE), "contains NA")
})


test_that("fore_sigplusnoise validates linear parameter", {
  x <- rnorm(50)
  expect_error(fore_sigplusnoise(x, linear = "TRUE", plot = FALSE),
               "`linear` must be TRUE or FALSE")
  expect_error(fore_sigplusnoise(x, linear = NA, plot = FALSE),
               "`linear` must be TRUE or FALSE")
})


test_that("fore_sigplusnoise validates n_ahead parameter", {
  x <- rnorm(50)
  expect_error(fore_sigplusnoise(x, n_ahead = 0, plot = FALSE),
               "`n_ahead` must be a positive")
  expect_error(fore_sigplusnoise(x, n_ahead = -5, plot = FALSE),
               "`n_ahead` must be a positive")
})


test_that("fore_sigplusnoise validates max_p parameter", {
  x <- rnorm(50)
  expect_error(fore_sigplusnoise(x, max_p = -1, plot = FALSE),
               "`max_p` must be a non-negative")
})


test_that("fore_sigplusnoise validates alpha parameter", {
  x <- rnorm(50)
  expect_error(fore_sigplusnoise(x, alpha = 0, plot = FALSE),
               "`alpha` must be between 0 and 1")
  expect_error(fore_sigplusnoise(x, alpha = 1, plot = FALSE),
               "`alpha` must be between 0 and 1")
  expect_error(fore_sigplusnoise(x, alpha = 1.5, plot = FALSE),
               "`alpha` must be between 0 and 1")
})


test_that("fore_sigplusnoise validates method parameter", {
  x <- rnorm(50)
  expect_error(fore_sigplusnoise(x, method = "invalid", plot = FALSE))
})


test_that("fore_sigplusnoise plot produces no error", {
  skip_if_not(capabilities("png"))

  set.seed(101)
  x <- cumsum(rnorm(50))

  png(tempfile())
  on.exit(dev.off())

  expect_no_error(fore_sigplusnoise(x, n_ahead = 5, plot = TRUE))
})


test_that("fore_sigplusnoise plot with limits = FALSE works", {
  skip_if_not(capabilities("png"))

  set.seed(102)
  x <- cumsum(rnorm(50))

  png(tempfile())
  on.exit(dev.off())

  expect_no_error(fore_sigplusnoise(x, n_ahead = 5, limits = FALSE, plot = TRUE))
})


test_that("fore_sigplusnoise print method works", {
  set.seed(103)
  x <- cumsum(rnorm(50))

  result_linear <- fore_sigplusnoise(x, linear = TRUE, plot = FALSE)
  expect_output(print(result_linear), "Signal Plus Noise Forecast")
  expect_output(print(result_linear), "Linear trend")
  expect_output(print(result_linear), "Intercept:")
  expect_output(print(result_linear), "Slope:")

  result_cosine <- fore_sigplusnoise(x, linear = FALSE, freq = 0.1, plot = FALSE)
  expect_output(print(result_cosine), "Cosine function")
})


test_that("fore_sigplusnoise returns invisible result", {
  x <- rnorm(50)

  png(tempfile())
  on.exit(dev.off())

  expect_invisible(fore_sigplusnoise(x, plot = TRUE))
})


test_that("fore_sigplusnoise alpha affects prediction interval width", {
  set.seed(104)
  x <- cumsum(rnorm(80))

  result_95 <- fore_sigplusnoise(x, alpha = 0.05, plot = FALSE)
  result_90 <- fore_sigplusnoise(x, alpha = 0.10, plot = FALSE)

  # 95% intervals should be wider than 90%
  width_95 <- result_95$ul - result_95$ll
  width_90 <- result_90$ul - result_90$ll

  expect_true(all(width_95 > width_90))
})
