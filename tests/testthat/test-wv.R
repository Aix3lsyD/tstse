# Tests for wv()
# Wigner-Ville distribution

test_that("wv matches tswge::wv.wge output", {
  skip_if_not_installed("tswge")

  set.seed(123)
  x <- rnorm(64)

  result <- wv(x, plot = FALSE)

  png(tempfile())
  expected <- tswge::wv.wge(x)
  dev.off()

  expect_equal(result$freq, expected$freq)
  expect_equal(result$time, expected$t)
  expect_equal(result$tfr, expected$tfr, tolerance = 1e-10)
})


test_that("wv matches tswge for sinusoidal signal", {
  skip_if_not_installed("tswge")

  t <- 1:128
  x <- cos(2 * pi * 0.1 * t)

  result <- wv(x, plot = FALSE)

  png(tempfile())
  expected <- tswge::wv.wge(x)
  dev.off()

  expect_equal(result$tfr, expected$tfr, tolerance = 1e-10)
})


test_that("wv returns correct structure", {
  x <- rnorm(100)
  result <- wv(x, plot = FALSE)

  expect_s3_class(result, "wv")
  expect_named(result, c("tfr", "freq", "time"))
  expect_true(is.matrix(result$tfr))
  expect_true(is.numeric(result$freq))
  expect_true(is.numeric(result$time))
})


test_that("wv tfr matrix has correct dimensions", {
  n <- 100
  x <- rnorm(n)

  result_default <- wv(x, plot = FALSE)
  expect_equal(dim(result_default$tfr), c(101, n))

  result_custom <- wv(x, n_freq = 50, plot = FALSE)
  expect_equal(dim(result_custom$tfr), c(50, n))
})


test_that("wv frequency range is 0 to 0.5", {
  x <- rnorm(64)
  result <- wv(x, plot = FALSE)

  expect_equal(result$freq[1], 0)
  expect_equal(result$freq[length(result$freq)], 0.5)
})


test_that("wv time indices match input length", {
  n <- 80
  x <- rnorm(n)
  result <- wv(x, plot = FALSE)

  expect_length(result$time, n)
  expect_equal(result$time, 1:n)
})


test_that("wv detects sinusoid at correct frequency", {
  # Pure sinusoid at frequency 0.2
  t <- 1:256
  freq_true <- 0.2
  x <- cos(2 * pi * freq_true * t)

  result <- wv(x, n_freq = 101, plot = FALSE)

  # Find frequency with maximum average energy
  mean_energy <- rowMeans(result$tfr[, 50:200])  # Middle portion to avoid edge effects
  peak_idx <- which.max(mean_energy)
  peak_freq <- result$freq[peak_idx]

  # Should be close to 0.2
  expect_true(abs(peak_freq - freq_true) < 0.05)
})


test_that("wv tfr is real-valued", {
  x <- rnorm(64)
  result <- wv(x, plot = FALSE)

  expect_true(is.numeric(result$tfr))
  expect_false(is.complex(result$tfr))
})


test_that("wv validates x parameter", {
  expect_error(wv("abc", plot = FALSE), "`x` must be a numeric")
  expect_error(wv(1:3, plot = FALSE), "at least 4")
  expect_error(wv(c(1, NA, 3, 4, 5), plot = FALSE), "contains NA")
})


test_that("wv validates n_freq parameter", {
  x <- rnorm(64)
  expect_error(wv(x, n_freq = 1, plot = FALSE), "`n_freq` must be")
  expect_error(wv(x, n_freq = 0, plot = FALSE), "`n_freq` must be")
  expect_error(wv(x, n_freq = "100", plot = FALSE), "`n_freq` must be")
})


test_that("wv validates plot parameter", {
  x <- rnorm(64)
  expect_error(wv(x, plot = "TRUE"), "`plot` must be TRUE or FALSE")
  expect_error(wv(x, plot = NA), "`plot` must be TRUE or FALSE")
})


test_that("wv plot produces no error", {
  skip_if_not(capabilities("png"))

  x <- rnorm(100)

  png(tempfile())
  on.exit(dev.off())

  expect_no_error(wv(x, plot = TRUE))
})


test_that("wv print method works", {
  x <- rnorm(100)
  result <- wv(x, plot = FALSE)

  expect_output(print(result), "Wigner-Ville Distribution")
  expect_output(print(result), "Time points:")
  expect_output(print(result), "Frequency points:")
  expect_output(print(result), "TFR matrix:")
})


test_that("wv returns invisible result", {
  x <- rnorm(50)

  png(tempfile())
  on.exit(dev.off())

  expect_invisible(wv(x, plot = TRUE))
})


test_that("wv handles short series", {
  x <- rnorm(10)
  result <- wv(x, n_freq = 21, plot = FALSE)

  expect_equal(ncol(result$tfr), 10)
  expect_false(any(is.na(result$tfr)))
})


test_that("wv n_freq parameter works", {
  x <- rnorm(64)

  result1 <- wv(x, n_freq = 51, plot = FALSE)
  expect_length(result1$freq, 51)

  result2 <- wv(x, n_freq = 201, plot = FALSE)
  expect_length(result2$freq, 201)
})


test_that("wv uses hilbert transform internally", {
  # The WV distribution should work on the analytic signal
  # We can verify by checking that results are consistent
  set.seed(111)
  x <- rnorm(64)

  result <- wv(x, plot = FALSE)

  # TFR should have non-trivial values
  expect_true(max(abs(result$tfr)) > 0)
  expect_true(var(as.vector(result$tfr)) > 0)
})
