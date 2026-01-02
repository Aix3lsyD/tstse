# Tests for is_glambda()
# Theoretical instantaneous spectrum (G-Lambda)

test_that("is_glambda runs without error for basic case", {
  expect_no_error(
    result <- is_glambda(n = 100, phi = 0.9, lambda = 1, offset = 10, plot = FALSE)
  )
  expect_s3_class(result, "is_glambda")
})


test_that("is_glambda returns correct structure", {
  result <- is_glambda(n = 50, phi = 0.5, lambda = 1, offset = 10, plot = FALSE)

  expect_s3_class(result, "is_glambda")
  expect_named(result, c("time", "freq", "gfreq", "spec", "phi", "lambda", "offset"))
  expect_equal(result$lambda, 1)
  expect_equal(result$offset, 10)
  expect_equal(result$phi, 0.5)
})


test_that("is_glambda time ranges from 1 to n", {
  result <- is_glambda(n = 50, phi = 0.5, lambda = 1, offset = 10, plot = FALSE)

  expect_true(min(result$time) == 1)
  expect_true(max(result$time) == 50)
})


test_that("is_glambda handles white noise (phi = 0)", {
  result <- is_glambda(n = 50, phi = 0, lambda = 1, offset = 10, plot = FALSE)

  expect_length(result$phi, 0)

  # White noise should have flat spectrum (all values equal after normalization)
  # Since spec starts at 0 (normalized), all values should be 0
  expect_true(all(result$spec == 0))
})


test_that("is_glambda handles AR(1)", {
  result <- is_glambda(n = 50, phi = 0.9, lambda = 1, offset = 10, plot = FALSE)

  expect_length(result$phi, 1)
  expect_equal(result$phi, 0.9)
})


test_that("is_glambda handles AR(2)", {
  result <- is_glambda(n = 50, phi = c(1.5, -0.75), lambda = 1, offset = 10, plot = FALSE)

  expect_length(result$phi, 2)
  expect_equal(result$phi, c(1.5, -0.75))
})


test_that("is_glambda works with lambda = 0", {
  expect_no_error(
    result <- is_glambda(n = 50, phi = 0.5, lambda = 0, offset = 10, plot = FALSE)
  )
  expect_equal(result$lambda, 0)
})


test_that("is_glambda works with lambda = 1", {
  expect_no_error(
    result <- is_glambda(n = 50, phi = 0.5, lambda = 1, offset = 10, plot = FALSE)
  )
  expect_equal(result$lambda, 1)
})


test_that("is_glambda works with lambda = 2", {
  expect_no_error(
    result <- is_glambda(n = 50, phi = 0.5, lambda = 2, offset = 10, plot = FALSE)
  )
  expect_equal(result$lambda, 2)
})


test_that("is_glambda spec values are non-negative", {
  result <- is_glambda(n = 100, phi = 0.9, lambda = 1, offset = 10, plot = FALSE)

  # Spec should be normalized to start at 0
  expect_true(min(result$spec) >= 0)
})


test_that("is_glambda n_freq parameter works", {
  result1 <- is_glambda(n = 20, phi = 0.5, lambda = 1, offset = 10,
                        n_freq = 100, plot = FALSE)
  result2 <- is_glambda(n = 20, phi = 0.5, lambda = 1, offset = 10,
                        n_freq = 200, plot = FALSE)

  # Different n_freq should give different total lengths
  # Total length = n * n_freq
  expect_equal(length(result1$freq), 20 * 100)
  expect_equal(length(result2$freq), 20 * 200)
})


test_that("is_glambda validates n parameter", {
  expect_error(is_glambda(n = 0, phi = 0.5, lambda = 1, offset = 10),
               "`n` must be a positive")
  expect_error(is_glambda(n = -10, phi = 0.5, lambda = 1, offset = 10),
               "`n` must be a positive")
})


test_that("is_glambda validates phi parameter", {
  expect_error(is_glambda(n = 50, phi = "0.5", lambda = 1, offset = 10),
               "`phi` must be numeric")
})


test_that("is_glambda validates sigma2 parameter", {
  expect_error(is_glambda(n = 50, phi = 0.5, sigma2 = 0, lambda = 1, offset = 10),
               "`sigma2` must be a positive")
  expect_error(is_glambda(n = 50, phi = 0.5, sigma2 = -1, lambda = 1, offset = 10),
               "`sigma2` must be a positive")
})


test_that("is_glambda validates lambda parameter", {
  expect_error(is_glambda(n = 50, phi = 0.5, lambda = "1", offset = 10),
               "`lambda` must be")
  expect_error(is_glambda(n = 50, phi = 0.5, lambda = c(1, 2), offset = 10),
               "`lambda` must be")
})


test_that("is_glambda validates offset parameter", {
  expect_error(is_glambda(n = 50, phi = 0.5, lambda = 1, offset = 0),
               "`offset` must be a positive")
  expect_error(is_glambda(n = 50, phi = 0.5, lambda = 1, offset = -5),
               "`offset` must be a positive")
})


test_that("is_glambda validates plot parameter", {
  expect_error(is_glambda(n = 50, phi = 0.5, lambda = 1, offset = 10, plot = "TRUE"),
               "`plot` must be TRUE or FALSE")
})


test_that("is_glambda plot produces no error", {
  skip_if_not(capabilities("png"))

  png(tempfile())
  on.exit(dev.off())

  expect_no_error(is_glambda(n = 50, phi = 0.9, lambda = 1, offset = 10, plot = TRUE))
})


test_that("is_glambda print method works", {
  result <- is_glambda(n = 50, phi = c(1.5, -0.75), lambda = 1, offset = 10, plot = FALSE)

  expect_output(print(result), "Theoretical Instantaneous Spectrum")
  expect_output(print(result), "AR\\(2\\)")
  expect_output(print(result), "Lambda:")
  expect_output(print(result), "Phi:")
})


test_that("is_glambda print method handles white noise", {
  result <- is_glambda(n = 50, phi = 0, lambda = 1, offset = 10, plot = FALSE)

  expect_output(print(result), "White noise")
})


test_that("is_glambda returns invisible result", {
  png(tempfile())
  on.exit(dev.off())

  expect_invisible(is_glambda(n = 50, phi = 0.5, lambda = 1, offset = 10, plot = TRUE))
})


test_that("is_glambda AR(1) spectrum peaks at low frequency for positive phi", {
  result <- is_glambda(n = 100, phi = 0.9, lambda = 1, offset = 10,
                       n_freq = 100, plot = FALSE)

  # For AR(1) with positive phi, spectrum peaks at low frequencies
  # Get unique frequencies and their spectra
  unique_freq <- unique(result$gfreq)
  unique_spec <- result$spec[1:length(unique_freq)]

  # Max should be at low frequency end
  max_idx <- which.max(unique_spec)
  expect_true(max_idx < length(unique_freq) / 2)
})
