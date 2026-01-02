# tests/testthat/test-gen_arch.R

test_that("gen_arch produces correct length output", {
  x <- gen_arch(n = 100, alpha0 = 0.1, alpha = 0.3, plot = FALSE, seed = 42)
  expect_length(x, 100)
})

test_that("gen_arch is reproducible with seed", {
  x1 <- gen_arch(n = 50, alpha0 = 0.1, alpha = 0.3, plot = FALSE, seed = 123)
  x2 <- gen_arch(n = 50, alpha0 = 0.1, alpha = 0.3, plot = FALSE, seed = 123)
  expect_equal(x1, x2)
})

test_that("gen_arch handles ARCH(q) with multiple alpha coefficients", {
  x <- gen_arch(n = 100, alpha0 = 0.1, alpha = c(0.2, 0.15, 0.1),
                plot = FALSE, seed = 42)
  expect_length(x, 100)
  expect_true(is.numeric(x))
})

test_that("gen_arch validates inputs", {
  expect_error(gen_arch(n = -1, alpha0 = 0.1, alpha = 0.3),
               "`n` must be a positive integer")
  expect_error(gen_arch(n = 100, alpha0 = -0.1, alpha = 0.3),
               "`alpha0` must be a positive number")
  expect_error(gen_arch(n = 100, alpha0 = 0.1, alpha = numeric(0)),
               "`alpha` must be a non-empty numeric vector")
  expect_error(gen_arch(n = 100, alpha0 = 0.1, alpha = -0.3),
               "all `alpha` coefficients must be non-negative")
})

test_that("gen_arch returns invisibly when plot = TRUE", {
  expect_invisible(gen_arch(n = 50, alpha0 = 0.1, alpha = 0.3, seed = 42))
})

test_that("gen_arch matches tswge::gen.arch.wge output with same seed", {
  skip_if_not_installed("tswge")

  # ARCH(1) test
  tstse_result <- gen_arch(n = 200, alpha0 = 0.1, alpha = 0.3,
                           plot = FALSE, seed = 42)
  tswge_result <- tswge::gen.arch.wge(n = 200, alpha0 = 0.1, alpha = 0.3,
                                      plot = FALSE, sn = 42)
  expect_equal(tstse_result, tswge_result, tolerance = 1e-10)

  # ARCH(2) test
  tstse_result2 <- gen_arch(n = 200, alpha0 = 0.1, alpha = c(0.2, 0.15),
                            plot = FALSE, seed = 123)
  tswge_result2 <- tswge::gen.arch.wge(n = 200, alpha0 = 0.1, alpha = c(0.2, 0.15),
                                       plot = FALSE, sn = 123)
  expect_equal(tstse_result2, tswge_result2, tolerance = 1e-10)
})

test_that("gen_arch produces volatility clustering", {
  set.seed(999)
  x <- gen_arch(n = 1000, alpha0 = 0.1, alpha = 0.7, plot = FALSE, seed = 999)

  # ARCH processes should show autocorrelation in squared values
  acf_sq <- acf(x^2, lag.max = 10, plot = FALSE)$acf[-1]

  # With alpha = 0.7, we expect significant positive autocorrelation in squares

  expect_true(any(acf_sq > 0.1))
})

test_that("gen_arch produces zero-mean series", {
  x <- gen_arch(n = 5000, alpha0 = 0.1, alpha = 0.3, plot = FALSE, seed = 42)

  # Mean should be approximately 0 for large samples
  expect_lt(abs(mean(x)), 0.1)
})

test_that("gen_arch respects spin parameter", {
  # Different spin values should produce different results
  # (even with same seed) since burn-in affects the state
  x1 <- gen_arch(n = 100, alpha0 = 0.1, alpha = 0.3, spin = 500,
                 plot = FALSE, seed = 42)
  x2 <- gen_arch(n = 100, alpha0 = 0.1, alpha = 0.3, spin = 2000,
                 plot = FALSE, seed = 42)

  # With longer spin, more randomness is "burned in"
  # They should NOT be equal

  expect_false(isTRUE(all.equal(x1, x2)))
})
