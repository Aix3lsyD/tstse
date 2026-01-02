# Tests for gen_ar_fast function

# ==============================================================================
# Basic functionality tests
# ==============================================================================

test_that("gen_ar_fast returns correct length", {
  x <- gen_ar_fast(100, phi = 0.7, seed = 123)
  expect_length(x, 100)

  x <- gen_ar_fast(500, phi = c(0.5, 0.3), seed = 123)
  expect_length(x, 500)
})

test_that("gen_ar_fast returns numeric vector", {
  x <- gen_ar_fast(100, phi = 0.7, seed = 123)
  expect_type(x, "double")
  expect_true(is.numeric(x))
})

test_that("gen_ar_fast produces no NA values", {
  x <- gen_ar_fast(100, phi = 0.7, seed = 123)
  expect_false(any(is.na(x)))

  x <- gen_ar_fast(100, phi = c(0.5, 0.3), seed = 123)
  expect_false(any(is.na(x)))
})


# ==============================================================================
# Reproducibility tests
# ==============================================================================

test_that("gen_ar_fast is reproducible with seed", {
  x1 <- gen_ar_fast(100, phi = 0.7, seed = 123)
  x2 <- gen_ar_fast(100, phi = 0.7, seed = 123)
  expect_identical(x1, x2)
})

test_that("gen_ar_fast produces different results with different seeds", {
  x1 <- gen_ar_fast(100, phi = 0.7, seed = 123)
  x2 <- gen_ar_fast(100, phi = 0.7, seed = 456)
  expect_false(identical(x1, x2))
})

test_that("gen_ar_fast produces different results without seed", {
  set.seed(111)
  x1 <- gen_ar_fast(100, phi = 0.7)
  set.seed(222)
  x2 <- gen_ar_fast(100, phi = 0.7)
  expect_false(identical(x1, x2))
})


# ==============================================================================
# White noise case tests
# ==============================================================================

test_that("gen_ar_fast handles white noise (empty phi)", {
  x <- gen_ar_fast(100, phi = numeric(0), seed = 123)
  expect_length(x, 100)
  expect_false(any(is.na(x)))
})

test_that("gen_ar_fast handles white noise (phi = 0)", {
  x <- gen_ar_fast(100, phi = 0, seed = 123)
  expect_length(x, 100)
  expect_false(any(is.na(x)))
})

test_that("gen_ar_fast handles white noise (all zeros)", {
  x <- gen_ar_fast(100, phi = c(0, 0, 0), seed = 123)
  expect_length(x, 100)
  expect_false(any(is.na(x)))
})

test_that("gen_ar_fast white noise has approximately correct variance", {
  # With vara = 1, variance should be approximately 1
  x <- gen_ar_fast(10000, phi = numeric(0), vara = 1, seed = 123)
  expect_equal(var(x), 1, tolerance = 0.1)

  # With vara = 4, variance should be approximately 4
  x <- gen_ar_fast(10000, phi = numeric(0), vara = 4, seed = 123)
  expect_equal(var(x), 4, tolerance = 0.2)
})


# ==============================================================================
# AR(1) tests
# ==============================================================================

test_that("gen_ar_fast AR(1) has correct autocorrelation structure", {
  # AR(1) with phi = 0.7 should have lag-1 ACF near 0.7
  x <- gen_ar_fast(5000, phi = 0.7, seed = 123)
  acf_vals <- acf(x, lag.max = 1, plot = FALSE)$acf
  expect_equal(acf_vals[2], 0.7, tolerance = 0.1)
})

test_that("gen_ar_fast AR(1) variance is correct", {
  # AR(1) theoretical variance = vara / (1 - phi^2)
  phi <- 0.7
  vara <- 1
  theoretical_var <- vara / (1 - phi^2)

  x <- gen_ar_fast(10000, phi = phi, vara = vara, seed = 123)
  expect_equal(var(x), theoretical_var, tolerance = 0.2)
})


# ==============================================================================
# Near unit root tests
# ==============================================================================

test_that("gen_ar_fast handles near unit root (phi = 0.99)", {
  # Should not error and should produce valid output
  x <- gen_ar_fast(200, phi = 0.99, seed = 123)
  expect_length(x, 200)
  expect_false(any(is.na(x)))
  expect_false(any(is.infinite(x)))
})

test_that("gen_ar_fast handles phi = 0.999", {
  x <- gen_ar_fast(200, phi = 0.999, seed = 123)
  expect_length(x, 200)
  expect_false(any(is.na(x)))
  expect_false(any(is.infinite(x)))
})


# ==============================================================================
# AR(2) tests
# ==============================================================================

test_that("gen_ar_fast handles AR(2)", {
  x <- gen_ar_fast(500, phi = c(0.5, 0.3), seed = 123)
  expect_length(x, 500)
  expect_false(any(is.na(x)))
})

test_that("gen_ar_fast AR(2) produces stationary output for stationary phi", {
  # Stationary AR(2): phi = (0.5, 0.3)
  x <- gen_ar_fast(5000, phi = c(0.5, 0.3), seed = 123)

  # Should have finite variance (not exploding)
  expect_true(is.finite(var(x)))
  expect_lt(var(x), 100)  # Reasonable bound
})


# ==============================================================================
# Variance parameter tests
# ==============================================================================

test_that("gen_ar_fast respects vara parameter", {
  # Different vara should produce different scales
  x1 <- gen_ar_fast(5000, phi = 0.5, vara = 1, seed = 123)
  x2 <- gen_ar_fast(5000, phi = 0.5, vara = 4, seed = 123)

  # Ratio of variances should be approximately 4
  ratio <- var(x2) / var(x1)
  expect_equal(ratio, 4, tolerance = 0.5)
})


# ==============================================================================
# Input validation tests
# ==============================================================================

test_that("gen_ar_fast validates n", {
  expect_error(gen_ar_fast(0, phi = 0.7), "positive integer")
  expect_error(gen_ar_fast(-1, phi = 0.7), "positive integer")
  expect_error(gen_ar_fast("abc", phi = 0.7), "positive integer")
  expect_error(gen_ar_fast(c(10, 20), phi = 0.7), "positive integer")
})

test_that("gen_ar_fast validates phi", {
  expect_error(gen_ar_fast(100, phi = "abc"), "numeric vector")
})

test_that("gen_ar_fast validates vara", {
  expect_error(gen_ar_fast(100, phi = 0.7, vara = 0), "positive number")
  expect_error(gen_ar_fast(100, phi = 0.7, vara = -1), "positive number")
  expect_error(gen_ar_fast(100, phi = 0.7, vara = "abc"), "positive number")
})


# ==============================================================================
# Comparison with gen_arma tests
# ==============================================================================

test_that("gen_ar_fast produces statistically similar output to gen_arma", {
  # Both should produce AR(1) with phi = 0.7
  # They use different burn-in strategies, so won't be identical
  # But statistical properties should be similar

  set.seed(123)
  x_fast <- gen_ar_fast(5000, phi = 0.7)

  set.seed(123)
  x_arma <- gen_arma(5000, phi = 0.7, theta = 0, plot = FALSE)

  # Means should be similar (both near 0)
  expect_equal(mean(x_fast), 0, tolerance = 0.2)
  expect_equal(mean(x_arma), 0, tolerance = 0.2)

  # Variances should be similar
  # Theoretical variance = 1 / (1 - 0.7^2) = 1.96
  expect_equal(var(x_fast), var(x_arma), tolerance = 0.5)
})


# ==============================================================================
# Edge case tests
# ==============================================================================

test_that("gen_ar_fast handles n = 1", {
  x <- gen_ar_fast(1, phi = 0.7, seed = 123)
  expect_length(x, 1)
  expect_false(is.na(x))
})

test_that("gen_ar_fast handles high order AR", {
  # AR(5)
  phi <- c(0.2, 0.15, 0.1, 0.08, 0.05)
  x <- gen_ar_fast(500, phi = phi, seed = 123)
  expect_length(x, 500)
  expect_false(any(is.na(x)))
})

test_that("gen_ar_fast handles negative phi", {
  # AR(1) with negative coefficient (oscillating)
  x <- gen_ar_fast(500, phi = -0.7, seed = 123)
  expect_length(x, 500)
  expect_false(any(is.na(x)))

  # Should have negative lag-1 ACF
  acf_vals <- acf(x, lag.max = 1, plot = FALSE)$acf
  expect_lt(acf_vals[2], 0)
})
