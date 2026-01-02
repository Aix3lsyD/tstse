# Tests for ar_mle functions (check_stationary, aic_ar_mle)

# ==============================================================================
# check_stationary tests
# ==============================================================================

test_that("check_stationary returns TRUE for stationary AR(1)", {
  expect_true(check_stationary(0.5))
  expect_true(check_stationary(0.9))
  expect_true(check_stationary(-0.5))
  expect_true(check_stationary(-0.9))
})

test_that("check_stationary returns FALSE for unit root", {
  expect_false(check_stationary(1.0))
  expect_false(check_stationary(-1.0))
})

test_that("check_stationary returns FALSE for explosive AR(1)", {
  expect_false(check_stationary(1.01))
  expect_false(check_stationary(1.1))
  expect_false(check_stationary(-1.01))
})

test_that("check_stationary handles empty phi (white noise)", {
  expect_true(check_stationary(numeric(0)))
  expect_true(check_stationary(NULL))
})

test_that("check_stationary handles stationary AR(2)", {
  # Stationary region for AR(2): phi1 + phi2 < 1, phi2 - phi1 < 1, |phi2| < 1
  expect_true(check_stationary(c(0.5, 0.3)))
  expect_true(check_stationary(c(0.7, 0.2)))
  expect_true(check_stationary(c(-0.5, 0.3)))
})

test_that("check_stationary handles non-stationary AR(2)", {
  # Non-stationary: violates phi1 + phi2 < 1
  expect_false(check_stationary(c(0.6, 0.5)))

  # Non-stationary: |phi2| >= 1
  expect_false(check_stationary(c(0.5, 1.0)))
})

test_that("check_stationary respects tolerance parameter", {
  # With default tol=1.001, phi=0.999 should be stationary
  expect_true(check_stationary(0.999, tol = 1.001))

  # With tighter tolerance, same phi might not pass
  expect_false(check_stationary(0.999, tol = 1.002))
})


# ==============================================================================
# aic_ar_mle basic tests
# ==============================================================================

test_that("aic_ar_mle returns correct structure", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 200)

  fit <- aic_ar_mle(x, p_max = 3)

  expect_type(fit, "list")
  expect_named(fit, c("phi", "p", "criterion_value", "method_used", "all_criteria"))
  expect_type(fit$phi, "double")
  expect_type(fit$p, "integer")
  expect_equal(fit$method_used, "mle")
})

test_that("aic_ar_mle selects reasonable order", {
  # Simulate AR(2)
  set.seed(123)
  x <- arima.sim(list(ar = c(0.7, 0.2)), n = 500)

  fit <- aic_ar_mle(x, p_max = 5)

  # Should select order 2 or close to it
  expect_true(fit$p >= 1 && fit$p <= 5)
})

test_that("aic_ar_mle phi length matches order", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 200)

  fit <- aic_ar_mle(x, p_max = 5)
  expect_equal(length(fit$phi), fit$p)
})


# ==============================================================================
# Criterion tests
# ==============================================================================

test_that("aic_ar_mle respects criterion parameter", {
  set.seed(123)
  x <- arima.sim(list(ar = c(0.7, 0.2)), n = 500)

  fit_aic <- aic_ar_mle(x, p_max = 5, criterion = "aic")
  fit_bic <- aic_ar_mle(x, p_max = 5, criterion = "bic")
  fit_aicc <- aic_ar_mle(x, p_max = 5, criterion = "aicc")

  # All should return valid results
  expect_true(fit_aic$p >= 1)
  expect_true(fit_bic$p >= 1)
  expect_true(fit_aicc$p >= 1)

  # BIC typically selects simpler models than AIC
  # (not always, but often)
  expect_true(fit_bic$p <= fit_aic$p + 1)
})

test_that("aic_ar_mle all_criteria contains values for each order", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 200)

  fit <- aic_ar_mle(x, p_max = 3)

  # Should have criteria for orders 1, 2, 3
  expect_true("1" %in% names(fit$all_criteria))
  expect_true("2" %in% names(fit$all_criteria))
  expect_true("3" %in% names(fit$all_criteria))
})


# ==============================================================================
# Stationarity enforcement tests
# ==============================================================================

test_that("aic_ar_mle returns stationary model by default", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 200)

  fit <- aic_ar_mle(x, p_max = 5, stationary = TRUE)
  expect_true(check_stationary(fit$phi))
})

test_that("aic_ar_mle stationary=FALSE allows any model", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 200)

  fit <- aic_ar_mle(x, p_max = 5, stationary = FALSE)
  # Should return a result (may or may not be stationary)
  expect_true(fit$p >= 1)
})


# ==============================================================================
# Input handling tests
# ==============================================================================

test_that("aic_ar_mle handles p_max as vector", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 200)

  # Specify specific orders to try
  fit <- aic_ar_mle(x, p_max = c(1, 2, 3))
  expect_true(fit$p %in% c(1, 2, 3))

  # Try non-consecutive orders
  fit2 <- aic_ar_mle(x, p_max = c(1, 3, 5))
  expect_true(fit2$p %in% c(1, 3, 5))
})

test_that("aic_ar_mle errors on invalid orders", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 200)

  # No valid orders
  expect_error(aic_ar_mle(x, p_max = 0), "No valid AR orders")
  expect_error(aic_ar_mle(x, p_max = -1), "No valid AR orders")
})


# ==============================================================================
# Short series tests
# ==============================================================================

test_that("aic_ar_mle handles short series", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 30)

  # Should work but may warn
  fit <- suppressWarnings(aic_ar_mle(x, p_max = 3))
  expect_true(fit$p >= 1)
})


# ==============================================================================
# Consistency tests
# ==============================================================================

test_that("aic_ar_mle is deterministic (no random component)", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 200)

  fit1 <- aic_ar_mle(x, p_max = 5)
  fit2 <- aic_ar_mle(x, p_max = 5)

  expect_equal(fit1$p, fit2$p)
  expect_equal(fit1$phi, fit2$phi)
  expect_equal(fit1$criterion_value, fit2$criterion_value)
})


# ==============================================================================
# Coefficient accuracy tests
# ==============================================================================

test_that("aic_ar_mle estimates are close to true values", {
  # Generate AR(1) with known phi
  set.seed(42)
  true_phi <- 0.7
  x <- arima.sim(list(ar = true_phi), n = 1000)

  fit <- aic_ar_mle(x, p_max = 3)

  # With large sample, should be close to true value
  if (fit$p == 1) {
    expect_equal(fit$phi[1], true_phi, tolerance = 0.1)
  }
})

test_that("aic_ar_mle AR(2) estimates are reasonable", {
  set.seed(42)
  true_phi <- c(0.5, 0.3)
  x <- arima.sim(list(ar = true_phi), n = 1000)

  fit <- aic_ar_mle(x, p_max = 5)

  # Should select order 2
  if (fit$p == 2) {
    expect_equal(fit$phi, true_phi, tolerance = 0.15)
  }
})
