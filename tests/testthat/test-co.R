# Tests for co() - Cochrane-Orcutt estimation

test_that("co returns expected structure", {
  set.seed(123)
  x <- rnorm(100)

  result <- co(x, maxp = 3)

  expect_type(result, "list")
  expect_named(result, c("x", "z_x", "b0hat", "b1hat", "z_order", "z_phi",
                         "pvalue", "tco", "method"))
  expect_equal(result$x, x)
  expect_length(result$z_x, 100)
  expect_true(is.numeric(result$b0hat))
  expect_true(is.numeric(result$b1hat))
  expect_true(result$z_order >= 0)
  expect_true(result$pvalue >= 0 && result$pvalue <= 1)
})

test_that("co detects known trend", {
  set.seed(456)
  n <- 200
  t <- seq_len(n)
  true_slope <- 0.1
  x <- 5 + true_slope * t + rnorm(n)

  result <- co(x, maxp = 5)

  # Slope should be close to true value
  expect_true(abs(result$b1hat - true_slope) < 0.05)
  # Should be significantly different from zero
  expect_true(result$pvalue < 0.05)
})

test_that("co works with AR(1) residuals", {
  set.seed(789)
  n <- 200
  t <- seq_len(n)
  true_phi <- 0.7
  z <- arima.sim(list(ar = true_phi), n = n)
  x <- 10 + 0.05 * t + z

  result <- co(x, maxp = 5)

  # Should detect AR structure
  expect_true(result$z_order >= 1)
  # If AR(1) selected, phi should be reasonable
  if (result$z_order == 1) {
    expect_true(abs(result$z_phi[1] - true_phi) < 0.3)
  }
})

test_that("co works with different methods", {
  set.seed(111)
  x <- cumsum(rnorm(100)) + seq_len(100) * 0.1

  result_mle <- co(x, method = "mle")
  result_burg <- co(x, method = "burg")
  result_yw <- co(x, method = "yw")

  # All should return valid results
  expect_true(is.finite(result_mle$b1hat))
  expect_true(is.finite(result_burg$b1hat))
  expect_true(is.finite(result_yw$b1hat))

  # Methods should be recorded correctly
  expect_equal(result_mle$method, "mle")
  expect_equal(result_burg$method, "burg")
  expect_equal(result_yw$method, "yw")
})

test_that("co works with different criteria", {
  set.seed(222)
  x <- rnorm(100) + seq_len(100) * 0.05

  result_aic <- co(x, type = "aic")
  result_aicc <- co(x, type = "aicc")
  result_bic <- co(x, type = "bic")

  # All should return valid results
  expect_true(is.finite(result_aic$b1hat))
  expect_true(is.finite(result_aicc$b1hat))
  expect_true(is.finite(result_bic$b1hat))
})

test_that("co handles white noise (no trend)", {
  set.seed(333)
  x <- rnorm(100)  # Pure white noise

  result <- co(x, maxp = 5)

  # Slope should be close to zero
  expect_true(abs(result$b1hat) < 0.1)
  # p-value should be large (not significant)
  expect_true(result$pvalue > 0.05)
})

test_that("co validates input", {
  expect_error(co(character(10)), "`x` must be a non-empty numeric vector")
  expect_error(co(numeric(0)), "`x` must be a non-empty numeric vector")
  expect_error(co(1:10, maxp = -1), "`maxp` must be a positive integer")
  expect_error(co(1:10, maxp = c(1, 2)), "`maxp` must be a positive integer")
})

test_that("co works with short series", {
  set.seed(444)
  x <- rnorm(20) + seq_len(20) * 0.1

  # Should work with short series
  result <- co(x, maxp = 3)
  expect_true(is.finite(result$b1hat))
  expect_true(result$z_order <= 3)
})

test_that("co residuals have correct length", {
  set.seed(555)
  n <- 50
  x <- rnorm(n)

  result <- co(x, maxp = 3)

  # OLS residuals should be length n
  expect_length(result$z_x, n)
})

test_that("co matches tswge for basic case", {
  skip_if_not_installed("tswge")

  set.seed(666)
  x <- rnorm(100) + seq_len(100) * 0.05

  result_tstse <- co(x, maxp = 5, method = "burg")

  # Compare with tswge (if available)
  result_tswge <- tswge::co.wge(x, maxp = 5)

  # Results should be similar (not exact due to implementation differences)
  expect_true(abs(result_tstse$b1hat - result_tswge$b1hat) < 0.01)
  expect_equal(result_tstse$z_order, result_tswge$z.order)
})
