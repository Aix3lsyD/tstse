library(tswge)

skip_if_not_installed("signal")


test_that("butterworth returns correct structure", {

  set.seed(123)
  x <- rnorm(100)

  result <- butterworth(x, order = 4, type = "low", cutoff = 0.2, plot = FALSE)

  expect_type(result, "list")
  expect_named(result, "x_filt")
  expect_length(result$x_filt, 100)
  expect_true(is.numeric(result$x_filt))
})


test_that("butterworth matches tswge for lowpass filter", {

  set.seed(456)
  x <- cumsum(rnorm(200))

  old <- tswge::butterworth.wge(x, order = 4, type = "low", cutoff = 0.1, plot = FALSE)
  new <- butterworth(x, order = 4, type = "low", cutoff = 0.1, plot = FALSE)

  expect_equal(new$x_filt, old$x.filt, tolerance = 1e-10)
})


test_that("butterworth matches tswge for highpass filter", {

  set.seed(789)
  x <- sin(seq(0, 10 * pi, length.out = 200)) + rnorm(200, sd = 0.1)

  old <- tswge::butterworth.wge(x, order = 3, type = "high", cutoff = 0.3, plot = FALSE)
  new <- butterworth(x, order = 3, type = "high", cutoff = 0.3, plot = FALSE)

  expect_equal(new$x_filt, old$x.filt, tolerance = 1e-10)
})


test_that("butterworth matches tswge for bandpass filter", {

  set.seed(111)
  x <- rnorm(200)

  old <- tswge::butterworth.wge(x, order = 2, type = "pass", cutoff = c(0.1, 0.3), plot = FALSE)
  new <- butterworth(x, order = 2, type = "pass", cutoff = c(0.1, 0.3), plot = FALSE)

  expect_equal(new$x_filt, old$x.filt, tolerance = 1e-10)
})


test_that("butterworth matches tswge for bandstop filter", {

  set.seed(222)
  x <- rnorm(200)

  old <- tswge::butterworth.wge(x, order = 2, type = "stop", cutoff = c(0.2, 0.4), plot = FALSE)
  new <- butterworth(x, order = 2, type = "stop", cutoff = c(0.2, 0.4), plot = FALSE)

  expect_equal(new$x_filt, old$x.filt, tolerance = 1e-10)
})


test_that("butterworth validates input", {

  x <- rnorm(100)

  expect_error(butterworth("not numeric", type = "low", cutoff = 0.2),
               "`x` must be a non-empty numeric vector")

  expect_error(butterworth(x, order = -1, type = "low", cutoff = 0.2),
               "`order` must be a positive integer")

  expect_error(butterworth(x, type = "invalid", cutoff = 0.2),
               "'arg' should be one of")

  expect_error(butterworth(x, type = "low", cutoff = "bad"),
               "`cutoff` must be numeric")
})


test_that("butterworth validates cutoff for lowpass/highpass", {

  x <- rnorm(100)

  expect_error(butterworth(x, type = "low", cutoff = c(0.1, 0.2)),
               "must be a single value")

  expect_error(butterworth(x, type = "high", cutoff = 0),
               "must be between 0 and 0.5")

  expect_error(butterworth(x, type = "low", cutoff = 0.5),
               "must be between 0 and 0.5")

  expect_error(butterworth(x, type = "low", cutoff = -0.1),
               "must be between 0 and 0.5")
})


test_that("butterworth validates cutoff for bandpass/bandstop", {

  x <- rnorm(100)

  expect_error(butterworth(x, type = "pass", cutoff = 0.2),
               "must be a 2-element vector")

  expect_error(butterworth(x, type = "stop", cutoff = c(0.3, 0.2)),
               "must be less than")

  expect_error(butterworth(x, type = "pass", cutoff = c(0, 0.3)),
               "must be between 0 and 0.5")

  expect_error(butterworth(x, type = "stop", cutoff = c(0.2, 0.5)),
               "must be between 0 and 0.5")
})


test_that("butterworth lowpass reduces high frequency content", {

  # Create signal with low and high frequency components
  t <- seq(0, 1, length.out = 500)
  x <- sin(2 * pi * 5 * t) + sin(2 * pi * 100 * t)  # 5 Hz + 100 Hz

  result <- butterworth(x, order = 4, type = "low", cutoff = 0.1, plot = FALSE)

  # Variance of filtered signal should be less (high freq removed)
  expect_lt(var(result$x_filt), var(x))

  # The low frequency component should be preserved approximately
  low_freq <- sin(2 * pi * 5 * t)
  expect_equal(cor(result$x_filt, low_freq), 1, tolerance = 0.1)
})


test_that("butterworth plotting works", {

  set.seed(333)
  x <- rnorm(100)

  expect_no_error(butterworth(x, order = 4, type = "low", cutoff = 0.2, plot = TRUE))
})


test_that("butterworth with different filter orders", {

  set.seed(444)
  x <- rnorm(200)

  # Different orders should produce different results
  result2 <- butterworth(x, order = 2, type = "low", cutoff = 0.2, plot = FALSE)
  result4 <- butterworth(x, order = 4, type = "low", cutoff = 0.2, plot = FALSE)
  result8 <- butterworth(x, order = 8, type = "low", cutoff = 0.2, plot = FALSE)

  # Higher order = sharper cutoff, but all should be correlated
  expect_gt(cor(result2$x_filt, result4$x_filt), 0.9)
  expect_gt(cor(result4$x_filt, result8$x_filt), 0.9)

  # But not identical
  expect_false(all(result2$x_filt == result4$x_filt))
})
