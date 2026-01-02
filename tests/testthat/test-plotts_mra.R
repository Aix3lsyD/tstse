# Tests for plotts_mra()
# Multiresolution Analysis plotting

test_that("plotts_mra matches tswge::plotts.mra.wge output", {
  skip_if_not_installed("tswge")
  skip_if_not_installed("waveslim")

  set.seed(123)
  x <- rnorm(128)

  result <- plotts_mra(x, n_levels = 4, wavelet = "S8", plot = FALSE)

  png(tempfile())
  expected <- tswge::plotts.mra.wge(x, n.levels = 4, type = "S8")
  dev.off()

  # Compare MRA components
  for (j in 1:5) {  # 4 detail + 1 smooth
    expect_equal(result$mra[[j]], expected[[j]], tolerance = 1e-10)
  }
})


test_that("plotts_mra matches tswge with Haar wavelet", {
  skip_if_not_installed("tswge")
  skip_if_not_installed("waveslim")

  set.seed(456)
  x <- rnorm(64)

  result <- plotts_mra(x, n_levels = 3, wavelet = "haar", plot = FALSE)

  png(tempfile())
  expected <- tswge::plotts.mra.wge(x, n.levels = 3, type = "haar")
  dev.off()

  for (j in 1:4) {
    expect_equal(result$mra[[j]], expected[[j]], tolerance = 1e-10)
  }
})


test_that("plotts_mra returns correct structure", {
  skip_if_not_installed("waveslim")

  x <- rnorm(256)
  result <- plotts_mra(x, n_levels = 4, plot = FALSE)

  expect_s3_class(result, "plotts_mra")
  expect_named(result, c("mra", "smooth", "n_levels", "wavelet", "n"))
  expect_equal(result$n_levels, 4)
  expect_equal(result$wavelet, "S8")
  expect_equal(result$n, 256)
})


test_that("plotts_mra smooth matrix has correct dimensions", {
  skip_if_not_installed("waveslim")

  n <- 128
  x <- rnorm(n)

  result <- plotts_mra(x, n_levels = 4, plot = FALSE)

  # smooth matrix: (n_levels + 1) rows x n columns
  expect_equal(dim(result$smooth), c(5, n))
})


test_that("plotts_mra smooth approximations are cumulative", {
  skip_if_not_installed("waveslim")

  x <- rnorm(64)
  result <- plotts_mra(x, n_levels = 3, plot = FALSE)

  # S0 = S1 + d1
  expect_equal(result$smooth[1, ],
               result$smooth[2, ] + result$mra[[1]],
               tolerance = 1e-10)

  # S1 = S2 + d2
  expect_equal(result$smooth[2, ],
               result$smooth[3, ] + result$mra[[2]],
               tolerance = 1e-10)

  # S2 = S3 + d3
  expect_equal(result$smooth[3, ],
               result$smooth[4, ] + result$mra[[3]],
               tolerance = 1e-10)

  # S3 = s3 (smooth component)
  expect_equal(result$smooth[4, ], result$mra[[4]], tolerance = 1e-10)
})


test_that("plotts_mra S0 reconstructs original data", {
  skip_if_not_installed("waveslim")

  set.seed(789)
  x <- rnorm(128)
  result <- plotts_mra(x, n_levels = 4, plot = FALSE)

  # S0 should equal original data (sum of all MRA components)
  reconstructed <- Reduce(`+`, result$mra)
  expect_equal(result$smooth[1, ], reconstructed, tolerance = 1e-10)
  expect_equal(result$smooth[1, ], x, tolerance = 1e-10)
})


test_that("plotts_mra all wavelet types work", {
  skip_if_not_installed("waveslim")

  x <- rnorm(64)

  expect_no_error(plotts_mra(x, n_levels = 2, wavelet = "S8", plot = FALSE))
  expect_no_error(plotts_mra(x, n_levels = 2, wavelet = "haar", plot = FALSE))
  expect_no_error(plotts_mra(x, n_levels = 2, wavelet = "D4", plot = FALSE))
  expect_no_error(plotts_mra(x, n_levels = 2, wavelet = "D6", plot = FALSE))
  expect_no_error(plotts_mra(x, n_levels = 2, wavelet = "D8", plot = FALSE))
})


test_that("plotts_mra validates x parameter", {
  skip_if_not_installed("waveslim")

  expect_error(plotts_mra("abc", plot = FALSE), "`x` must be a numeric")
  expect_error(plotts_mra(1:3, plot = FALSE), "at least 4")
  expect_error(plotts_mra(c(1, NA, 3, 4), plot = FALSE), "contains NA")
})


test_that("plotts_mra validates n_levels parameter", {
  skip_if_not_installed("waveslim")

  x <- rnorm(64)
  expect_error(plotts_mra(x, n_levels = 0, plot = FALSE),
               "`n_levels` must be a positive")
  expect_error(plotts_mra(x, n_levels = -1, plot = FALSE),
               "`n_levels` must be a positive")
})


test_that("plotts_mra validates sample size divisibility", {
  skip_if_not_installed("waveslim")

  x <- rnorm(100)  # Not divisible by 16 (2^4)
  expect_error(plotts_mra(x, n_levels = 4, plot = FALSE),
               "must be divisible")
})


test_that("plotts_mra validates wavelet parameter", {
  skip_if_not_installed("waveslim")

  x <- rnorm(64)
  expect_error(plotts_mra(x, wavelet = "invalid", plot = FALSE))
})


test_that("plotts_mra validates plot parameter", {
  skip_if_not_installed("waveslim")

  x <- rnorm(64)
  expect_error(plotts_mra(x, plot = "TRUE"), "`plot` must be TRUE or FALSE")
  expect_error(plotts_mra(x, plot = NA), "`plot` must be TRUE or FALSE")
})


test_that("plotts_mra plot produces no error", {
  skip_if_not_installed("waveslim")
  skip_if_not(capabilities("png"))

  x <- rnorm(128)

  png(tempfile())
  on.exit(dev.off())

  expect_no_error(plotts_mra(x, n_levels = 4, plot = TRUE))
})


test_that("plotts_mra print method works", {
  skip_if_not_installed("waveslim")

  x <- rnorm(256)
  result <- plotts_mra(x, n_levels = 4, plot = FALSE)

  expect_output(print(result), "Multiresolution Analysis")
  expect_output(print(result), "Sample size: 256")
  expect_output(print(result), "Levels: 4")
  expect_output(print(result), "Wavelet: S8")
  expect_output(print(result), "S0:")
})


test_that("plotts_mra returns invisible result", {
  skip_if_not_installed("waveslim")

  x <- rnorm(64)

  png(tempfile())
  on.exit(dev.off())

  expect_invisible(plotts_mra(x, n_levels = 2, plot = TRUE))
})


test_that("plotts_mra smooth approximations get smoother at higher levels", {
  skip_if_not_installed("waveslim")

  # Generate signal with high-frequency noise
  set.seed(111)
  t <- 1:256
  x <- sin(2 * pi * t / 64) + rnorm(256, sd = 0.5)

  result <- plotts_mra(x, n_levels = 4, plot = FALSE)

  # Each successive smooth should have lower variance (be smoother)
  vars <- apply(result$smooth, 1, var)

  # S1 variance < S0 variance (removed finest detail)
  expect_true(vars[2] < vars[1])

  # Generally decreasing trend (may not be strictly monotonic)
  expect_true(vars[5] < vars[1])  # SJ much smoother than S0
})
