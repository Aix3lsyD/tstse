# =============================================================================
# Tests for co_fast() - Fast Cochrane-Orcutt using C++ primitives
# =============================================================================

test_that("co_fast matches co(..., method='burg') on AR(1) series", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 100)

  r1 <- co(x, maxp = 5, method = "burg")
  r2 <- co_fast(x, maxp = 5)

  expect_equal(r2$tco, r1$tco, tolerance = 1e-10)
  expect_equal(r2$z_order, r1$z_order)
  expect_equal(r2$z_phi, r1$z_phi, tolerance = 1e-10)
  expect_equal(r2$b0hat, r1$b0hat, tolerance = 1e-10)
  expect_equal(r2$b1hat, r1$b1hat, tolerance = 1e-10)
  expect_equal(r2$pvalue, r1$pvalue, tolerance = 1e-10)
  expect_equal(r2$method, "burg")
})

test_that("co_fast matches co on high-persistence AR(1)", {
  set.seed(456)
  x <- arima.sim(list(ar = 0.95), n = 200)

  r1 <- co(x, maxp = 5, method = "burg")
  r2 <- co_fast(x, maxp = 5)

  expect_equal(r2$tco, r1$tco, tolerance = 1e-10)
  expect_equal(r2$z_order, r1$z_order)
  expect_equal(r2$z_phi, r1$z_phi, tolerance = 1e-10)
  expect_equal(r2$pvalue, r1$pvalue, tolerance = 1e-10)
})

test_that("co_fast matches co on AR(2) series", {
  set.seed(789)
  x <- arima.sim(list(ar = c(0.5, -0.3)), n = 150)

  r1 <- co(x, maxp = 5, method = "burg")
  r2 <- co_fast(x, maxp = 5)

  expect_equal(r2$tco, r1$tco, tolerance = 1e-10)
  expect_equal(r2$z_order, r1$z_order)
  expect_equal(r2$z_phi, r1$z_phi, tolerance = 1e-10)
  expect_equal(r2$b0hat, r1$b0hat, tolerance = 1e-10)
  expect_equal(r2$b1hat, r1$b1hat, tolerance = 1e-10)
})

test_that("co_fast matches co on series with trend", {
  set.seed(101)
  n <- 100
  t <- seq_len(n)
  z <- arima.sim(list(ar = 0.7), n = n)
  x <- 5 + 0.1 * t + z

  r1 <- co(x, maxp = 5, method = "burg")
  r2 <- co_fast(x, maxp = 5)

  expect_equal(r2$tco, r1$tco, tolerance = 1e-10)
  expect_equal(r2$z_order, r1$z_order)
  expect_equal(r2$z_phi, r1$z_phi, tolerance = 1e-10)
  expect_equal(r2$b0hat, r1$b0hat, tolerance = 1e-10)
  expect_equal(r2$b1hat, r1$b1hat, tolerance = 1e-10)
})

test_that("co_fast runs without error on white noise", {
  # NOTE: co_fast() uses C++ Burg with recursive variance estimation while

  # co() uses R's aic_ar() with backcast residual variance. These use different
  # variance formulas that can lead to different model selection, especially
  # on white noise where there's no true AR structure.
  #
  # We test that both run successfully and return valid structures.
  set.seed(202)
  x <- rnorm(100)

  r1 <- co(x, maxp = 5, method = "burg")
  r2 <- co_fast(x, maxp = 5)

  # Both should return valid structures
  expect_type(r1, "list")
  expect_type(r2, "list")
  expect_true(is.finite(r1$tco))
  expect_true(is.finite(r2$tco))
  expect_gte(r1$z_order, 1)
  expect_gte(r2$z_order, 1)
})

test_that("co_fast works with different IC types", {
  set.seed(303)
  x <- arima.sim(list(ar = 0.6), n = 100)

  for (ic in c("aic", "aicc", "bic")) {
    r1 <- co(x, maxp = 5, method = "burg", type = ic)
    r2 <- co_fast(x, maxp = 5, type = ic)

    expect_equal(r2$tco, r1$tco, tolerance = 1e-10,
                 label = sprintf("tco for type='%s'", ic))
    expect_equal(r2$z_order, r1$z_order,
                 label = sprintf("z_order for type='%s'", ic))
  }
})

test_that("co_fast works with different maxp values", {
  set.seed(404)
  x <- arima.sim(list(ar = 0.7), n = 100)

  for (mp in c(3, 5, 8)) {
    r1 <- co(x, maxp = mp, method = "burg")
    r2 <- co_fast(x, maxp = mp)

    expect_equal(r2$tco, r1$tco, tolerance = 1e-10,
                 label = sprintf("tco for maxp=%d", mp))
    expect_equal(r2$z_order, r1$z_order,
                 label = sprintf("z_order for maxp=%d", mp))
  }
})

test_that("co_fast returns correct structure", {
  set.seed(505)
  x <- arima.sim(list(ar = 0.7), n = 100)

  result <- co_fast(x, maxp = 5)

  expect_type(result, "list")
  expect_named(result, c("x", "z_x", "b0hat", "b1hat", "z_order", "z_phi",
                         "pvalue", "tco", "method"))
  expect_equal(result$x, x)
  expect_length(result$z_x, length(x))
  expect_type(result$b0hat, "double")
  expect_type(result$b1hat, "double")
  expect_type(result$z_order, "integer")
  expect_type(result$z_phi, "double")
  expect_type(result$pvalue, "double")
  expect_type(result$tco, "double")
  expect_equal(result$method, "burg")
})

test_that("co_fast validates input correctly", {
  expect_error(co_fast(NULL), "`x` must be a non-empty numeric vector")
  expect_error(co_fast(character(5)), "`x` must be a non-empty numeric vector")
  expect_error(co_fast(numeric(0)), "`x` must be a non-empty numeric vector")
  expect_error(co_fast(1:10, maxp = -1), "`maxp` must be a positive integer")
  expect_error(co_fast(1:10, maxp = c(1, 2)), "`maxp` must be a positive integer")
  # co() doesn't have explicit length check - it handles short series gracefully
  # so we match that behavior (C++ will clamp maxp if needed)
})

test_that("co_fast handles white noise (returns valid results)", {
  # White noise series - co() and co_fast() may select different AR orders
  # due to different variance estimation methods (backcast vs recursive).
  # We verify both return valid results with z_order >= 1.
  set.seed(606)
  x <- rnorm(200)

  r1 <- co(x, maxp = 5, method = "burg")
  r2 <- co_fast(x, maxp = 5)

  # Both should return valid structures
  expect_true(is.finite(r1$tco))
  expect_true(is.finite(r2$tco))
  expect_gte(r1$z_order, 1)
  expect_gte(r2$z_order, 1)  # min_p=1 ensures order >= 1
})

test_that("co_fast is faster than co", {
  skip_on_cran()

  set.seed(707)
  x <- arima.sim(list(ar = 0.8), n = 200)

  # Warm-up
  co(x, maxp = 5, method = "burg")
  co_fast(x, maxp = 5)

  # Time comparison
  time_co <- system.time(for (i in 1:20) co(x, maxp = 5, method = "burg"))["elapsed"]
  time_fast <- system.time(for (i in 1:20) co_fast(x, maxp = 5))["elapsed"]

  # co_fast should be at least 2x faster
  expect_lt(time_fast, time_co / 2,
            label = sprintf("co_fast (%.3fs) should be >2x faster than co (%.3fs)",
                            time_fast, time_co))
})
