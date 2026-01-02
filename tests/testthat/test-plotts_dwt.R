# Tests for plotts_dwt()
# Discrete Wavelet Transform plotting

test_that("plotts_dwt matches tswge::plotts.dwt.wge output", {
  skip_if_not_installed("tswge")
  skip_if_not_installed("waveslim")

  set.seed(123)
  x <- rnorm(128)

  result <- plotts_dwt(x, n_levels = 4, wavelet = "S8", plot = FALSE)

  png(tempfile())
  expected <- tswge::plotts.dwt.wge(x, n.levels = 4, type = "S8")
  dev.off()

  # Compare DWT coefficients
  for (j in 1:5) {  # 4 detail + 1 smooth
    expect_equal(result$dwt[[j]], expected[[j]], tolerance = 1e-10)
  }
})


test_that("plotts_dwt matches tswge with Haar wavelet", {
  skip_if_not_installed("tswge")
  skip_if_not_installed("waveslim")

  set.seed(456)
  x <- rnorm(64)

  result <- plotts_dwt(x, n_levels = 3, wavelet = "haar", plot = FALSE)

  png(tempfile())
  expected <- tswge::plotts.dwt.wge(x, n.levels = 3, type = "haar")
  dev.off()

  for (j in 1:4) {
    expect_equal(result$dwt[[j]], expected[[j]], tolerance = 1e-10)
  }
})


test_that("plotts_dwt returns correct structure", {
  skip_if_not_installed("waveslim")

  x <- rnorm(256)
  result <- plotts_dwt(x, n_levels = 4, plot = FALSE)

  expect_s3_class(result, "plotts_dwt")
  expect_named(result, c("dwt", "n_levels", "wavelet", "n"))
  expect_equal(result$n_levels, 4)
  expect_equal(result$wavelet, "S8")
  expect_equal(result$n, 256)
})


test_that("plotts_dwt dwt object has correct number of components", {
  skip_if_not_installed("waveslim")

  x <- rnorm(128)

  result3 <- plotts_dwt(x, n_levels = 3, plot = FALSE)
  expect_length(result3$dwt, 4)  # d1, d2, d3, s3

  result5 <- plotts_dwt(x, n_levels = 5, plot = FALSE)
  expect_length(result5$dwt, 6)  # d1, d2, d3, d4, d5, s5
})


test_that("plotts_dwt detail coefficients have correct lengths", {
  skip_if_not_installed("waveslim")

  n <- 256
  x <- rnorm(n)
  result <- plotts_dwt(x, n_levels = 4, plot = FALSE)

  # Detail coefficients halve at each level
  expect_length(result$dwt[[1]], n / 2)    # d1: 128
  expect_length(result$dwt[[2]], n / 4)    # d2: 64
  expect_length(result$dwt[[3]], n / 8)    # d3: 32
  expect_length(result$dwt[[4]], n / 16)   # d4: 16
  expect_length(result$dwt[[5]], n / 16)   # s4: 16
})


test_that("plotts_dwt all wavelet types work", {
  skip_if_not_installed("waveslim")

  x <- rnorm(64)

  expect_no_error(plotts_dwt(x, n_levels = 2, wavelet = "S8", plot = FALSE))
  expect_no_error(plotts_dwt(x, n_levels = 2, wavelet = "haar", plot = FALSE))
  expect_no_error(plotts_dwt(x, n_levels = 2, wavelet = "D4", plot = FALSE))
  expect_no_error(plotts_dwt(x, n_levels = 2, wavelet = "D6", plot = FALSE))
  expect_no_error(plotts_dwt(x, n_levels = 2, wavelet = "D8", plot = FALSE))
})


test_that("plotts_dwt wavelet attribute is set correctly", {
  skip_if_not_installed("waveslim")

  x <- rnorm(64)

  result <- plotts_dwt(x, n_levels = 2, wavelet = "D4", plot = FALSE)
  expect_equal(result$wavelet, "D4")
  expect_equal(attr(result$dwt, "wavelet"), "D4")
})


test_that("plotts_dwt validates x parameter", {
  skip_if_not_installed("waveslim")

  expect_error(plotts_dwt("abc", plot = FALSE), "`x` must be a numeric")
  expect_error(plotts_dwt(1:3, plot = FALSE), "at least 4")
  expect_error(plotts_dwt(c(1, NA, 3, 4), plot = FALSE), "contains NA")
})


test_that("plotts_dwt validates n_levels parameter", {
  skip_if_not_installed("waveslim")

  x <- rnorm(64)
  expect_error(plotts_dwt(x, n_levels = 0, plot = FALSE),
               "`n_levels` must be a positive")
  expect_error(plotts_dwt(x, n_levels = -1, plot = FALSE),
               "`n_levels` must be a positive")
})


test_that("plotts_dwt validates sample size divisibility", {
  skip_if_not_installed("waveslim")

  x <- rnorm(100)  # Not divisible by 16 (2^4)
  expect_error(plotts_dwt(x, n_levels = 4, plot = FALSE),
               "must be divisible")
})


test_that("plotts_dwt validates wavelet parameter", {
  skip_if_not_installed("waveslim")

  x <- rnorm(64)
  expect_error(plotts_dwt(x, wavelet = "invalid", plot = FALSE))
})


test_that("plotts_dwt validates plot parameter", {
  skip_if_not_installed("waveslim")

  x <- rnorm(64)
  expect_error(plotts_dwt(x, plot = "TRUE"), "`plot` must be TRUE or FALSE")
  expect_error(plotts_dwt(x, plot = NA), "`plot` must be TRUE or FALSE")
})


test_that("plotts_dwt plot produces no error", {
  skip_if_not_installed("waveslim")
  skip_if_not(capabilities("png"))

  x <- rnorm(128)

  png(tempfile())
  on.exit(dev.off())

  expect_no_error(plotts_dwt(x, n_levels = 4, plot = TRUE))
})


test_that("plotts_dwt print method works", {
  skip_if_not_installed("waveslim")

  x <- rnorm(256)
  result <- plotts_dwt(x, n_levels = 4, plot = FALSE)

  expect_output(print(result), "Discrete Wavelet Transform")
  expect_output(print(result), "Sample size: 256")
  expect_output(print(result), "Levels: 4")
  expect_output(print(result), "Wavelet: S8")
  expect_output(print(result), "d1:")
  expect_output(print(result), "s4:")
})


test_that("plotts_dwt returns invisible result", {
  skip_if_not_installed("waveslim")

  x <- rnorm(64)

  png(tempfile())
  on.exit(dev.off())

  expect_invisible(plotts_dwt(x, n_levels = 2, plot = TRUE))
})

test_that("plotts_dwt handles minimum valid input", {
  skip_if_not_installed("waveslim")

  # Minimum: 4 observations, 1 level
  x <- rnorm(4)
  result <- plotts_dwt(x, n_levels = 1, plot = FALSE)

  expect_equal(result$n, 4)
  expect_equal(result$n_levels, 1)
  expect_length(result$dwt[[1]], 2)  # d1
  expect_length(result$dwt[[2]], 2)  # s1
})


test_that("plotts_dwt DWT is energy preserving", {
  skip_if_not_installed("waveslim")

  set.seed(789)
  x <- rnorm(128)
  result <- plotts_dwt(x, n_levels = 4, plot = FALSE)

  # Energy in original signal
  energy_orig <- sum(x^2)

  # Energy in DWT coefficients (Parseval's theorem)
  energy_dwt <- 0
  for (j in seq_len(result$n_levels + 1)) {
    energy_dwt <- energy_dwt + sum(result$dwt[[j]]^2)
  }

  expect_equal(energy_dwt, energy_orig, tolerance = 1e-10)
})
