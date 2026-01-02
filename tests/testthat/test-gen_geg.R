# Tests for gen_geg()
# Generate Gegenbauer process realizations

test_that("gen_geg matches tswge::gen.geg.wge output with seed", {
  skip_if_not_installed("tswge")

  # Use same seed and parameters
  u <- 0.8
  d <- 0.3
  n <- 100
  n_trunc <- 1000
  seed <- 123

  result <- gen_geg(n = n, u = u, d = d, n_trunc = n_trunc, seed = seed)
  expected <- tswge::gen.geg.wge(n = n, u = u, lambda = d, trun = n_trunc, sn = seed)

  # Allow some tolerance due to potential differences in gegenb implementation
  expect_equal(result, expected, tolerance = 0.1)
})


test_that("gen_geg matches tswge with different parameters", {
  skip_if_not_installed("tswge")

  u <- 0.5
  d <- 0.4
  n <- 150
  n_trunc <- 500
  seed <- 456

  result <- gen_geg(n = n, u = u, d = d, n_trunc = n_trunc, seed = seed)
  expected <- tswge::gen.geg.wge(n = n, u = u, lambda = d, trun = n_trunc, sn = seed)

  # Allow some tolerance due to potential differences in gegenb implementation
  expect_equal(result, expected, tolerance = 0.1)
})


test_that("gen_geg returns correct length", {
  result <- gen_geg(n = 100, u = 0.8, d = 0.3, seed = 123)
  expect_length(result, 100)

  result2 <- gen_geg(n = 500, u = 0.5, d = 0.2, seed = 456)
  expect_length(result2, 500)
})


test_that("gen_geg returns numeric vector", {
  result <- gen_geg(n = 50, u = 0.8, d = 0.3, seed = 123)

  expect_true(is.numeric(result))
  expect_true(is.vector(result))
  expect_false(any(is.na(result)))
})


test_that("gen_geg seed produces reproducible results", {
  result1 <- gen_geg(n = 100, u = 0.8, d = 0.3, seed = 123)
  result2 <- gen_geg(n = 100, u = 0.8, d = 0.3, seed = 123)

  expect_equal(result1, result2)
})


test_that("gen_geg NULL seed produces different results", {
  # Without seed, results should differ (with high probability)
  set.seed(NULL)
  result1 <- gen_geg(n = 100, u = 0.8, d = 0.3, seed = NULL)
  result2 <- gen_geg(n = 100, u = 0.8, d = 0.3, seed = NULL)

  expect_false(isTRUE(all.equal(result1, result2)))
})


test_that("gen_geg var_a affects variance", {
  # Same seed, different variance
  result1 <- gen_geg(n = 500, u = 0.8, d = 0.3, var_a = 1, seed = 123)
  result2 <- gen_geg(n = 500, u = 0.8, d = 0.3, var_a = 4, seed = 123)

  # var_a = 4 should produce higher variance output
  expect_true(var(result2) > var(result1))
})


test_that("gen_geg u = 1 gives long memory at low frequency", {
  # u = 1 means spectral peak at f = 0
  set.seed(123)
  x <- gen_geg(n = 500, u = 1, d = 0.3, n_trunc = 500)

  # ACF should decay slowly (long memory characteristic)
  # Check that ACF is positive at early lags
  acf_vals <- acf(x, lag.max = 20, plot = FALSE)$acf
  # Should have positive correlation at lag 1-5
  expect_true(all(acf_vals[2:6] > 0))  # acf_vals[1] is lag 0
})


test_that("gen_geg u = -1 gives long memory at high frequency", {
  # u = -1 means spectral peak at f = 0.5 (Nyquist)
  set.seed(123)
  x <- gen_geg(n = 500, u = -1, d = 0.3, n_trunc = 500)

  # Should have alternating pattern
  # ACF at odd lags should be negative
  acf_vals <- acf(x, lag.max = 10, plot = FALSE)$acf
  expect_true(acf_vals[2] < 0)  # Lag 1 should be negative
})


test_that("gen_geg handles u = 0 (frequency 0.25)", {
  set.seed(123)
  result <- gen_geg(n = 200, u = 0, d = 0.3, n_trunc = 500)

  expect_length(result, 200)
  expect_false(any(is.na(result)))
})


test_that("gen_geg warns for d >= 0.5", {
  expect_warning(
    gen_geg(n = 50, u = 0.8, d = 0.5, seed = 123),
    "non-stationary"
  )
})


test_that("gen_geg warns for d < 0", {
  expect_warning(
    gen_geg(n = 50, u = 0.8, d = -0.1, seed = 123),
    "anti-persistent"
  )
})


test_that("gen_geg validates n parameter", {
  expect_error(gen_geg(n = 0, u = 0.8, d = 0.3), "`n` must be a positive")
  expect_error(gen_geg(n = -10, u = 0.8, d = 0.3), "`n` must be a positive")
  expect_error(gen_geg(n = "100", u = 0.8, d = 0.3), "`n` must be a positive")
})


test_that("gen_geg validates u parameter", {
  expect_error(gen_geg(n = 100, u = "0.8", d = 0.3), "`u` must be numeric")
})


test_that("gen_geg validates d parameter", {
  expect_error(gen_geg(n = 100, u = 0.8, d = "0.3"), "`d` must be numeric")
})


test_that("gen_geg validates u and d length match", {
  expect_error(
    gen_geg(n = 100, u = c(0.8, 0.5), d = 0.3),
    "same length"
  )
})


test_that("gen_geg validates n_trunc parameter", {
  expect_error(gen_geg(n = 100, u = 0.8, d = 0.3, n_trunc = 0),
               "`n_trunc` must be a positive")
  expect_error(gen_geg(n = 100, u = 0.8, d = 0.3, n_trunc = -100),
               "`n_trunc` must be a positive")
})


test_that("gen_geg validates var_a parameter", {
  expect_error(gen_geg(n = 100, u = 0.8, d = 0.3, var_a = 0),
               "`var_a` must be a positive")
  expect_error(gen_geg(n = 100, u = 0.8, d = 0.3, var_a = -1),
               "`var_a` must be a positive")
})


test_that("gen_geg validates seed parameter", {
  expect_error(gen_geg(n = 100, u = 0.8, d = 0.3, seed = "abc"),
               "`seed` must be")
})


test_that("gen_geg larger n_trunc gives better approximation", {
  # With same seed, larger truncation should give similar but slightly
  # different results (converging as n_trunc -> infinity)
  result_small <- gen_geg(n = 100, u = 0.8, d = 0.3, n_trunc = 100, seed = 123)
  result_large <- gen_geg(n = 100, u = 0.8, d = 0.3, n_trunc = 1000, seed = 123)

  # Should be different (due to different truncation)
  expect_false(isTRUE(all.equal(result_small, result_large)))
})
