# test-cpp-components.R
# Tests for C++ implementations of wbg_boot_fast components
# Validates correctness against R equivalents

test_that("ols_detrend_cpp matches R lm() residuals", {
  set.seed(123)

  # Test various series lengths
  for (n in c(50, 100, 500, 1000)) {
    x <- 5 + 0.1 * (1:n) + rnorm(n)

    resid_cpp <- as.numeric(ols_detrend_cpp(x))
    resid_r <- as.numeric(resid(lm(x ~ seq_along(x))))

    expect_equal(length(resid_cpp), length(resid_r))
    expect_equal(resid_cpp, resid_r, tolerance = 1e-10,
                 info = paste("n =", n))
  }
})

test_that("ols_detrend_cpp handles edge cases", {
  # Exact linear trend (residuals should be ~0)
  x <- 1:100
  resid_cpp <- as.numeric(ols_detrend_cpp(x))
  expect_true(all(abs(resid_cpp) < 1e-10))

  # Constant series
  x <- rep(5, 100)
  resid_cpp <- as.numeric(ols_detrend_cpp(x))
  expect_true(all(abs(resid_cpp) < 1e-10))
})

test_that("burg_fit_cpp matches ar.burg() coefficients", {
  set.seed(456)

  # Test with known AR processes (stationary)
  test_cases <- list(
    list(phi = 0.7, p = 1),
    list(phi = c(0.5, 0.2), p = 2),
    list(phi = c(0.4, 0.2, 0.1), p = 3)
  )

  for (tc in test_cases) {
    # Generate AR(p) data
    x <- arima.sim(list(ar = tc$phi), n = 500)

    phi_cpp <- as.numeric(burg_fit_cpp(x, tc$p))
    fit_r <- ar.burg(x, order.max = tc$p, aic = FALSE)
    phi_r <- as.numeric(fit_r$ar)

    expect_equal(length(phi_cpp), tc$p)
    expect_equal(phi_cpp, phi_r, tolerance = 1e-8,
                 info = paste("AR order p =", tc$p))
  }
})

test_that("burg_fit_cpp handles AR(1) near unit root", {
  set.seed(789)
  x <- arima.sim(list(ar = 0.95), n = 500)

  phi_cpp <- as.numeric(burg_fit_cpp(x, 1))
  fit_r <- ar.burg(x, order.max = 1, aic = FALSE)

  expect_equal(phi_cpp[1], fit_r$ar[1], tolerance = 1e-8)
})

test_that("burg_fit_full_cpp returns complete output", {
  set.seed(111)
  x <- arima.sim(list(ar = c(0.7, -0.2)), n = 300)

  result <- burg_fit_full_cpp(x, 2)

  expect_true("phi" %in% names(result))
  expect_true("vara" %in% names(result))
  expect_equal(length(result$phi), 2)
  expect_true(result$vara > 0)
})

test_that("burg_aic_select_cpp selects same order as aic_burg()", {
  set.seed(222)

  # AR(2) data
  x <- arima.sim(list(ar = c(0.7, -0.2)), n = 500)

  cpp_result <- burg_aic_select_cpp(x, maxp = 5, criterion = "aic")
  r_result <- aic_burg(x, p = 1:5, type = "aic")

  expect_equal(cpp_result$p, r_result$p,
               info = "Selected AR order should match")
  expect_equal(as.numeric(cpp_result$phi), as.numeric(r_result$phi),
               tolerance = 1e-10,
               info = "AR coefficients should match")
  # Variance may differ slightly due to denominator (n vs n-p)
  expect_equal(cpp_result$vara, r_result$vara, tolerance = 0.01,
               info = "Variance should be close")
})

test_that("burg_aic_select_cpp works with different criteria", {
  set.seed(333)
  x <- arima.sim(list(ar = c(0.5, 0.2)), n = 300)

  for (crit in c("aic", "aicc", "bic")) {
    result <- burg_aic_select_cpp(x, maxp = 5, criterion = crit)

    expect_true(result$p >= 1 && result$p <= 5,
                info = paste("Criterion:", crit))
    expect_equal(length(result$phi), result$p,
                 info = paste("Criterion:", crit))
    expect_true(result$vara > 0,
                info = paste("Criterion:", crit))
  }
})

test_that("ar_transform_cpp matches artrans()", {
  set.seed(444)
  x <- rnorm(200)
  phi <- c(0.7, -0.2)

  w_cpp <- as.numeric(ar_transform_cpp(x, phi))
  w_r <- as.numeric(artrans(x, phi, plot = FALSE))

  expect_equal(length(w_cpp), length(w_r))
  expect_equal(w_cpp, w_r, tolerance = 1e-10)
})

test_that("ar_transform_cpp handles different AR orders", {
  set.seed(555)
  x <- rnorm(100)

  for (p in 1:5) {
    phi <- rep(0.1, p)  # Small coefficients to ensure stationarity

    w_cpp <- as.numeric(ar_transform_cpp(x, phi))
    w_r <- as.numeric(artrans(x, phi, plot = FALSE))

    expect_equal(length(w_cpp), length(x) - p,
                 info = paste("p =", p))
    expect_equal(w_cpp, w_r, tolerance = 1e-10,
                 info = paste("p =", p))
  }
})

test_that("ar_transform_cpp with empty phi returns original", {
  x <- rnorm(100)
  w_cpp <- as.numeric(ar_transform_cpp(x, numeric(0)))
  expect_equal(w_cpp, x)
})

test_that("co_time_transform_cpp matches R implementation", {
  n <- 100
  phi <- c(0.7, -0.2)
  p <- length(phi)

  t_co_cpp <- as.numeric(co_time_transform_cpp(n, phi))

  # R equivalent
  t_co_r <- numeric(n - p)
  for (tt in (p + 1):n) {
    t_co_r[tt - p] <- tt - sum(phi * (tt - seq_len(p)))
  }

  expect_equal(length(t_co_cpp), n - p)
  expect_equal(t_co_cpp, t_co_r, tolerance = 1e-10)
})

test_that("co_time_transform_cpp handles single coefficient", {
  n <- 50
  phi <- 0.8

  t_co_cpp <- as.numeric(co_time_transform_cpp(n, phi))

  # R equivalent for AR(1)
  t_co_r <- (2:n) - phi * (1:(n - 1))

  expect_equal(t_co_cpp, t_co_r, tolerance = 1e-10)
})

test_that("co_tstat_cpp matches co() t-statistic", {
  set.seed(666)

  # Test with trend + AR errors
  x <- arima.sim(list(ar = 0.7), n = 300) + 0.05 * (1:300)

  tstat_cpp <- co_tstat_cpp(x, maxp = 5, criterion = "aic")
  tstat_r <- co(x, maxp = 5, method = "burg", type = "aic")$tco

  expect_equal(tstat_cpp, tstat_r, tolerance = 1e-10)
})

test_that("co_full_cpp returns complete results matching co()", {
  set.seed(777)
  x <- arima.sim(list(ar = c(0.5, -0.1)), n = 400) + 0.02 * (1:400)

  cpp_result <- co_full_cpp(x, maxp = 5, criterion = "aic")
  r_result <- co(x, maxp = 5, method = "burg", type = "aic")

  expect_equal(cpp_result$tco, r_result$tco, tolerance = 1e-10,
               info = "t-statistic should match")
  expect_equal(cpp_result$p, r_result$z_order,
               info = "AR order should match")
  expect_equal(as.numeric(cpp_result$phi), as.numeric(r_result$z_phi),
               tolerance = 1e-10,
               info = "AR coefficients should match")
})

test_that("co_tstat_cpp works with no trend", {
  set.seed(888)

  # Pure AR process (no trend)
  x <- arima.sim(list(ar = 0.6), n = 200)

  tstat_cpp <- co_tstat_cpp(x, maxp = 3)
  tstat_r <- co(x, maxp = 3, method = "burg", type = "aic")$tco

  expect_equal(tstat_cpp, tstat_r, tolerance = 1e-10)
})

test_that("gen_ar_seeded_cpp is reproducible", {
  phi <- c(0.7, -0.2)
  n <- 500

  # Same seed produces identical results
  x1 <- as.numeric(gen_ar_seeded_cpp(n, phi, vara = 1.0, rng_seed = 12345))
  x2 <- as.numeric(gen_ar_seeded_cpp(n, phi, vara = 1.0, rng_seed = 12345))

  expect_identical(x1, x2)

  # Different seeds produce different results
  x3 <- as.numeric(gen_ar_seeded_cpp(n, phi, vara = 1.0, rng_seed = 54321))
  expect_false(identical(x1, x3))
})

test_that("gen_ar_seeded_cpp produces correct length", {
  phi <- c(0.5, 0.2)

  for (n in c(50, 100, 500, 1000)) {
    x <- gen_ar_seeded_cpp(n, phi, vara = 1.0, rng_seed = 42)
    expect_equal(length(x), n)
  }
})

test_that("gen_ar_seeded_cpp produces stationary series", {
  phi <- c(0.7, -0.2)
  n <- 2000

  x <- as.numeric(gen_ar_seeded_cpp(n, phi, vara = 1.0, rng_seed = 999))

  # Check variance is reasonable (not exploding)
  expect_true(var(x) < 100)

  # Check sample ACF roughly matches expected pattern
  acf_vals <- acf(x, lag.max = 5, plot = FALSE)$acf[-1]
  expect_true(acf_vals[1] > 0.5)  # AR(1) component positive
})

test_that("gen_ar_seeded_cpp handles AR(0) (white noise)", {
  n <- 500
  x <- as.numeric(gen_ar_seeded_cpp(n, numeric(0), vara = 1.0, rng_seed = 123))

  expect_equal(length(x), n)
  # Should have near-zero autocorrelation
  acf_vals <- acf(x, lag.max = 5, plot = FALSE)$acf[-1]
  expect_true(all(abs(acf_vals) < 0.15))
})

test_that("gen_ar_seeded_cpp respects variance parameter", {
  phi <- 0.5
  n <- 5000

  # Theoretical variance for AR(1): vara / (1 - phi^2)
  vara <- 2.0
  theoretical_var <- vara / (1 - phi^2)

  x <- as.numeric(gen_ar_seeded_cpp(n, phi, vara = vara, rng_seed = 42))
  sample_var <- var(x)

  # Allow 20% tolerance for sampling variability
  expect_equal(sample_var, theoretical_var, tolerance = 0.2 * theoretical_var)
})

test_that("calc_ar_burnin_cpp returns reasonable values", {
  # Near unit root should have large burn-in
  burn_high <- calc_ar_burnin_cpp(c(0.99), 1000)
  expect_true(burn_high >= 500)

  # Low persistence should have minimal burn-in
  burn_low <- calc_ar_burnin_cpp(c(0.1), 1000)
  expect_equal(burn_low, 50)

  # Normal case
  burn_mid <- calc_ar_burnin_cpp(c(0.7, -0.2), 1000)
  expect_true(burn_mid >= 50 && burn_mid <= 2000)
})

test_that("all C++ functions handle minimum viable input", {
  # Minimum series length
  x <- rnorm(10)

  expect_length(ols_detrend_cpp(x), 10)
  expect_length(burg_fit_cpp(x, 1), 1)
  expect_true(is.list(burg_aic_select_cpp(x, maxp = 2)))
  expect_length(ar_transform_cpp(x, 0.5), 9)
  expect_length(co_time_transform_cpp(10, 0.5), 9)
})

test_that("ols_detrend_full_cpp returns slope and intercept", {
  set.seed(999)
  x <- 10 + 2 * (1:100) + rnorm(100, sd = 0.1)

  result <- ols_detrend_full_cpp(x)

  expect_true("residuals" %in% names(result))
  expect_true("intercept" %in% names(result))
  expect_true("slope" %in% names(result))
  expect_equal(length(result$residuals), 100)
  expect_equal(result$intercept, 10, tolerance = 0.5)
  expect_equal(result$slope, 2, tolerance = 0.1)
})

test_that("ols_tstat_cpp computes correct t-statistic", {
  set.seed(1111)

  # Create data with known slope
  n <- 100
  t_idx <- 1:n
  y <- 5 + 0.5 * t_idx + rnorm(n, sd = 2)

  # C++ t-stat
  tstat_cpp <- ols_tstat_cpp(y, t_idx)

  # R equivalent
  fit_r <- lm(y ~ t_idx)
  tstat_r <- summary(fit_r)$coefficients[2, 3]

  expect_equal(tstat_cpp, tstat_r, tolerance = 1e-10)
})
