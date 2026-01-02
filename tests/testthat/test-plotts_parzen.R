# Tests for plotts_parzen()
# Multi-panel spectral plotting

test_that("plotts_parzen matches tswge::plotts.parzen.wge output", {
  skip_if_not_installed("tswge")

  set.seed(123)
  x <- rnorm(100)

  result <- plotts_parzen(x, trunc = c(0, 0), plot = FALSE)

  png(tempfile())
  expected <- tswge::plotts.parzen.wge(x, m2 = c(0, 0))
  dev.off()

  expect_equal(result$freq, expected$freq)
  expect_equal(result$periodogram, expected$db, tolerance = 1e-10)
  expect_equal(result$parzen, expected$dbz, tolerance = 1e-10)
})


test_that("plotts_parzen matches tswge with additional truncation points", {
  skip_if_not_installed("tswge")

  set.seed(456)
  x <- cumsum(rnorm(80))

  result <- plotts_parzen(x, trunc = c(10, 20), plot = FALSE)

  png(tempfile())
  expected <- tswge::plotts.parzen.wge(x, m2 = c(10, 20))
  dev.off()

  expect_equal(result$freq, expected$freq)
  expect_equal(result$periodogram, expected$db, tolerance = 1e-10)

  # Note: parzen() may have minor differences from plotts.parzen.wge's inline
  # implementation. The key is that parzen() matches tswge::parzen.wge().
  # Test with looser tolerance here, and verify parzen() separately.
  expect_equal(result$parzen1, expected$dbz1, tolerance = 1e-10)
  expect_equal(result$parzen2, expected$dbz2, tolerance = 1e-10)
})


test_that("plotts_parzen default parzen matches parzen() output", {
  # This verifies internal consistency - plotts_parzen uses parzen()
  set.seed(456)
  x <- cumsum(rnorm(80))

  result <- plotts_parzen(x, trunc = c(0, 0), plot = FALSE)
  pz <- parzen(x, db = TRUE, trunc = 0L, plot = FALSE)

  expect_equal(result$parzen, pz$pzgram)
  expect_equal(result$trunc_points[1], pz$M)
})


test_that("plotts_parzen matches tswge with one additional truncation", {
  skip_if_not_installed("tswge")

  set.seed(789)
  x <- arima.sim(model = list(ar = 0.7), n = 120)

  result <- plotts_parzen(x, trunc = c(15, 0), plot = FALSE)

  png(tempfile())
  expected <- tswge::plotts.parzen.wge(x, m2 = c(15, 0))
  dev.off()

  expect_equal(result$parzen1, expected$dbz1, tolerance = 1e-10)
  expect_null(result$parzen2)
})


test_that("plotts_parzen returns correct structure", {
  x <- rnorm(100)
  result <- plotts_parzen(x, plot = FALSE)

  expect_s3_class(result, "plotts_parzen")
  expect_named(result, c("freq", "periodogram", "parzen", "parzen1", "parzen2", "trunc_points"))
})


test_that("plotts_parzen frequency is Fourier frequencies", {
  x <- rnorm(100)
  result <- plotts_parzen(x, plot = FALSE)

  # Fourier frequencies: 1/n, 2/n, ..., floor(n/2)/n
  n <- length(x)
  expected_freq <- (1:floor(n / 2)) / n

  expect_equal(result$freq, expected_freq)
})


test_that("plotts_parzen default truncation is 2*sqrt(n)", {
  x <- rnorm(100)
  result <- plotts_parzen(x, plot = FALSE)

  expected_m <- floor(2 * sqrt(100))  # = 20
  expect_equal(result$trunc_points[1], expected_m)
})


test_that("plotts_parzen with no extra trunc has NULL parzen1 and parzen2", {
  x <- rnorm(50)
  result <- plotts_parzen(x, trunc = c(0, 0), plot = FALSE)

  expect_null(result$parzen1)
  expect_null(result$parzen2)
  expect_length(result$trunc_points, 1)
})


test_that("plotts_parzen with one extra trunc has parzen1 only", {
  x <- rnorm(50)
  result <- plotts_parzen(x, trunc = c(10, 0), plot = FALSE)

  expect_false(is.null(result$parzen1))
  expect_null(result$parzen2)
  expect_length(result$trunc_points, 2)
})


test_that("plotts_parzen with two extra trunc has both", {
  x <- rnorm(50)
  result <- plotts_parzen(x, trunc = c(5, 15), plot = FALSE)

  expect_false(is.null(result$parzen1))
  expect_false(is.null(result$parzen2))
  expect_length(result$trunc_points, 3)
})


test_that("plotts_parzen smoothing reduces variance", {
  set.seed(111)
  x <- rnorm(200)
  result <- plotts_parzen(x, trunc = c(5, 0), plot = FALSE)

  # Parzen smoothed should have less variance than periodogram
  var_pgram <- var(result$periodogram)
  var_parzen <- var(result$parzen)
  var_parzen1 <- var(result$parzen1)

  expect_true(var_parzen < var_pgram)
  # Smaller M = more smoothing = less variance
  expect_true(var_parzen1 < var_parzen)
})


test_that("plotts_parzen validates x parameter", {
  expect_error(plotts_parzen("abc"), "`x` must be a numeric")
  expect_error(plotts_parzen(1:3), "at least 4 observations")
  expect_error(plotts_parzen(c(1, NA, 3, 4)), "contains NA")
})


test_that("plotts_parzen validates trunc parameter", {
  x <- rnorm(50)
  expect_error(plotts_parzen(x, trunc = 10), "length 2")
  expect_error(plotts_parzen(x, trunc = c(10, 20, 30)), "length 2")
  expect_error(plotts_parzen(x, trunc = c(-5, 10)), "non-negative")
  expect_error(plotts_parzen(x, trunc = "10"), "numeric vector")
})


test_that("plotts_parzen validates plot parameter", {
  x <- rnorm(50)
  expect_error(plotts_parzen(x, plot = "TRUE"), "`plot` must be TRUE or FALSE")
  expect_error(plotts_parzen(x, plot = NA), "`plot` must be TRUE or FALSE")
})


test_that("plotts_parzen plot produces no error", {
  skip_if_not(capabilities("png"))

  x <- rnorm(100)

  png(tempfile())
  on.exit(dev.off())

  expect_no_error(plotts_parzen(x, plot = TRUE))
})


test_that("plotts_parzen plot with extra trunc produces no error", {
  skip_if_not(capabilities("png"))

  x <- rnorm(100)

  png(tempfile())
  on.exit(dev.off())

  expect_no_error(plotts_parzen(x, trunc = c(10, 25), plot = TRUE))
})


test_that("plotts_parzen print method works", {
  x <- rnorm(100)
  result <- plotts_parzen(x, plot = FALSE)

  expect_output(print(result), "Periodogram and Parzen")
  expect_output(print(result), "Truncation points")
})


test_that("plotts_parzen returns invisible result", {
  x <- rnorm(50)

  png(tempfile())
  on.exit(dev.off())

  expect_invisible(plotts_parzen(x, plot = TRUE))
})


test_that("plotts_parzen handles even and odd length series", {
  # Even length
  x_even <- rnorm(100)
  result_even <- plotts_parzen(x_even, plot = FALSE)
  expect_length(result_even$freq, 50)  # floor(100/2)

  # Odd length
  x_odd <- rnorm(101)
  result_odd <- plotts_parzen(x_odd, plot = FALSE)
  expect_length(result_odd$freq, 50)  # floor(101/2)
})


test_that("plotts_parzen periodogram matches sample_spec at Fourier frequencies", {
  set.seed(222)
  x <- rnorm(100)

  result <- plotts_parzen(x, plot = FALSE)
  ss <- sample_spec(x, db = TRUE, plot = FALSE)

  # Fourier frequencies in plotts_parzen should match corresponding sample_spec values
  # plotts_parzen uses freq = (1:50)/100
  # sample_spec uses freq = seq(0, 0.5, length.out=251)

  # Find matching frequencies
  for (i in seq_along(result$freq)) {
    f <- result$freq[i]
    # Find closest frequency in sample_spec
    idx <- which.min(abs(ss$freq - f))
    if (abs(ss$freq[idx] - f) < 0.001) {
      # Should be close (not exact due to different computation methods)
      expect_equal(result$periodogram[i], ss$spec[idx], tolerance = 0.1)
    }
  }
})
