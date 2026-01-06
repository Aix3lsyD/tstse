# ==============================================================================
# Tests for jt_fast() - Optimized Jacob Turner Trend Test
# ==============================================================================

# ------------------------------------------------------------------------------
# jt_fast() must produce identical results to jt()
# ------------------------------------------------------------------------------

test_that("jt_fast matches jt on AR(1) low persistence", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.3), n = 100)

  r1 <- jt(x, maxp = 5)
  r2 <- jt_fast(x, maxp = 5)

  expect_equal(r2$tco, r1$tco, tolerance = 1e-10)
  expect_equal(r2$pvalue, r1$pvalue, tolerance = 1e-10)
  expect_equal(r2$vif, r1$vif, tolerance = 1e-10)
  expect_equal(r2$b0hat, r1$b0hat, tolerance = 1e-10)
  expect_equal(r2$b1hat, r1$b1hat, tolerance = 1e-10)
  expect_equal(r2$z_order, r1$z_order)
  expect_equal(r2$x_order, r1$x_order)
  expect_equal(r2$rel_n, r1$rel_n, tolerance = 1e-10)
  expect_equal(r2$est_df, r1$est_df, tolerance = 1e-10)
})

test_that("jt_fast matches jt on AR(1) high persistence", {
 set.seed(123)
  x <- arima.sim(list(ar = 0.9), n = 100)

  r1 <- jt(x, maxp = 5)
  r2 <- jt_fast(x, maxp = 5)

  expect_equal(r2$tco, r1$tco, tolerance = 1e-10)
  expect_equal(r2$pvalue, r1$pvalue, tolerance = 1e-10)
  expect_equal(r2$vif, r1$vif, tolerance = 1e-10)
  expect_equal(r2$z_order, r1$z_order)
  expect_equal(r2$x_order, r1$x_order)
})

test_that("jt_fast matches jt on AR(2)", {
  set.seed(123)
  x <- arima.sim(list(ar = c(0.5, -0.3)), n = 100)

  r1 <- jt(x, maxp = 5)
  r2 <- jt_fast(x, maxp = 5)

  expect_equal(r2$tco, r1$tco, tolerance = 1e-10)
  expect_equal(r2$pvalue, r1$pvalue, tolerance = 1e-10)
  expect_equal(r2$vif, r1$vif, tolerance = 1e-10)
  expect_equal(r2$z_order, r1$z_order)
  expect_equal(r2$x_order, r1$x_order)
})

test_that("jt_fast matches jt on white noise", {
  set.seed(123)
  x <- rnorm(100)

  r1 <- jt(x, maxp = 5)
  r2 <- jt_fast(x, maxp = 5)

  expect_equal(r2$tco, r1$tco, tolerance = 1e-10)
  expect_equal(r2$pvalue, r1$pvalue, tolerance = 1e-10)
  expect_equal(r2$vif, r1$vif, tolerance = 1e-10)
  expect_equal(r2$z_order, r1$z_order)
  expect_equal(r2$x_order, r1$x_order)
})

test_that("jt_fast matches jt on series with trend", {
  set.seed(123)
  n <- 100
  x <- 0.1 * (1:n) + arima.sim(list(ar = 0.7), n = n)

  r1 <- jt(x, maxp = 5)
  r2 <- jt_fast(x, maxp = 5)

  expect_equal(r2$tco, r1$tco, tolerance = 1e-10)
  expect_equal(r2$pvalue, r1$pvalue, tolerance = 1e-10)
  expect_equal(r2$vif, r1$vif, tolerance = 1e-10)
  expect_equal(r2$z_order, r1$z_order)
  expect_equal(r2$x_order, r1$x_order)
})

test_that("jt_fast matches jt with different maxp values", {
  set.seed(42)
  x <- arima.sim(list(ar = 0.8), n = 150)

  for (maxp in c(3, 5, 10)) {
    r1 <- jt(x, maxp = maxp)
    r2 <- jt_fast(x, maxp = maxp)

    expect_equal(r2$tco, r1$tco, tolerance = 1e-10,
                 info = paste("maxp =", maxp))
    expect_equal(r2$pvalue, r1$pvalue, tolerance = 1e-10,
                 info = paste("maxp =", maxp))
    expect_equal(r2$vif, r1$vif, tolerance = 1e-10,
                 info = paste("maxp =", maxp))
  }
})

test_that("jt_fast matches jt with different type options", {
  set.seed(42)
  x <- arima.sim(list(ar = 0.7), n = 100)

  for (type in c("aic", "aicc", "bic")) {
    r1 <- jt(x, maxp = 5, type = type)
    r2 <- jt_fast(x, maxp = 5, type = type)

    expect_equal(r2$tco, r1$tco, tolerance = 1e-10,
                 info = paste("type =", type))
    expect_equal(r2$pvalue, r1$pvalue, tolerance = 1e-10,
                 info = paste("type =", type))
    expect_equal(r2$vif, r1$vif, tolerance = 1e-10,
                 info = paste("type =", type))
  }
})

test_that("jt_fast matches jt on longer series", {
  set.seed(42)
  x <- arima.sim(list(ar = 0.8), n = 500)

  r1 <- jt(x, maxp = 5)
  r2 <- jt_fast(x, maxp = 5)

  expect_equal(r2$tco, r1$tco, tolerance = 1e-10)
  expect_equal(r2$pvalue, r1$pvalue, tolerance = 1e-10)
  expect_equal(r2$vif, r1$vif, tolerance = 1e-10)
})

test_that("jt_fast matches jt on very persistent series", {
  set.seed(42)
  x <- arima.sim(list(ar = 0.95), n = 200)

  r1 <- jt(x, maxp = 5)
  r2 <- jt_fast(x, maxp = 5)

  expect_equal(r2$tco, r1$tco, tolerance = 1e-10)
  expect_equal(r2$pvalue, r1$pvalue, tolerance = 1e-10)
  expect_equal(r2$vif, r1$vif, tolerance = 1e-10)
})

# ------------------------------------------------------------------------------
# jt_fast() should have same structure as jt()
# ------------------------------------------------------------------------------

test_that("jt_fast returns correct structure", {
  set.seed(42)
  x <- arima.sim(list(ar = 0.7), n = 200)
  result <- jt_fast(x)

  # Check all expected components are present
  expected_names <- c("x", "z_x", "b0hat", "b1hat", "z_order", "z_phi",
                      "x_order", "x_phi", "vif", "rel_n", "est_df",
                      "pvalue", "tco")
  expect_true(all(expected_names %in% names(result)))

  # Check dimensions
  expect_length(result$x, 200)
  expect_length(result$z_x, 200)
})

# ------------------------------------------------------------------------------
# jt_fast() input validation should match jt()
# ------------------------------------------------------------------------------

test_that("jt_fast validates input correctly", {
  # Non-numeric input
  expect_error(jt_fast("not numeric"), "`x` must be a non-empty numeric vector")

  # Empty vector
  expect_error(jt_fast(numeric(0)), "`x` must be a non-empty numeric vector")

  # Invalid maxp
  expect_error(jt_fast(rnorm(100), maxp = -1), "`maxp` must be a positive integer")
  expect_error(jt_fast(rnorm(100), maxp = "five"), "`maxp` must be a positive integer")

  # Series too short
  expect_error(jt_fast(rnorm(10), maxp = 5), "Series too short")
})

# ------------------------------------------------------------------------------
# Randomized equivalence tests
# ------------------------------------------------------------------------------

test_that("jt_fast matches jt on random AR(1) series", {
  for (i in 1:5) {
    set.seed(1000 + i)
    phi <- runif(1, 0.1, 0.95)
    n <- sample(50:200, 1)
    x <- arima.sim(list(ar = phi), n = n)

    r1 <- jt(x, maxp = 5)
    r2 <- jt_fast(x, maxp = 5)

    expect_equal(r2$tco, r1$tco, tolerance = 1e-10,
                 info = paste("seed =", 1000 + i, ", phi =", round(phi, 3), ", n =", n))
    expect_equal(r2$pvalue, r1$pvalue, tolerance = 1e-10,
                 info = paste("seed =", 1000 + i))
    expect_equal(r2$vif, r1$vif, tolerance = 1e-10,
                 info = paste("seed =", 1000 + i))
  }
})
