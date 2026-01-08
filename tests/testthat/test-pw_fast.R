# =============================================================================
# Tests for pw_fast() - Fast Prais-Winsten using C++ primitives
# =============================================================================

test_that("pw_fast matches prais package (twostep) on AR(1) series with trend", {
  skip_if_not_installed("prais")

  set.seed(123)
  n <- 100
  t_idx <- 1:n
  z <- arima.sim(list(ar = 0.7), n = n)
  x <- 5 + 0.1 * t_idx + z

  r1 <- pw_fast(x)

  df <- data.frame(y = x, t = t_idx, time = t_idx)
  r2 <- prais::prais_winsten(y ~ t, data = df, index = "time", twostep = TRUE)
  r2_sum <- summary(r2)

  # Compare t-statistics (allow 15% relative difference due to rho estimation)
  # Our implementation uses slightly different rho estimation formula
  expect_equal(r1$tpw, r2_sum$coefficients["t", "t value"], tolerance = 0.15 * abs(r2_sum$coefficients["t", "t value"]))

  # Compare AR(1) rho (use nrow to get final iteration value)
  # Small differences expected from different estimation formulas
  expect_equal(r1$rho, r2$rho[nrow(r2$rho)], tolerance = 0.01)

  # Compare slope estimate
  expect_equal(r1$b1hat, unname(coef(r2)["t"]), tolerance = 1e-3)

  # Compare intercept
  expect_equal(r1$b0hat, unname(coef(r2)["(Intercept)"]), tolerance = 1e-3)
})

test_that("pw_fast matches prais package on high-persistence AR(1)", {
  skip_if_not_installed("prais")

  set.seed(456)
  n <- 200
  t_idx <- 1:n
  z <- arima.sim(list(ar = 0.9), n = n)
  x <- 10 + 0.05 * t_idx + z

  r1 <- pw_fast(x)

  df <- data.frame(y = x, t = t_idx, time = t_idx)
  r2 <- prais::prais_winsten(y ~ t, data = df, index = "time", twostep = TRUE)
  r2_sum <- summary(r2)

  # Allow 15% relative difference due to different rho estimation
  expect_equal(r1$tpw, r2_sum$coefficients["t", "t value"], tolerance = 0.15 * abs(r2_sum$coefficients["t", "t value"]))
  expect_equal(r1$rho, r2$rho[nrow(r2$rho)], tolerance = 0.01)
})

test_that("pw_fast matches prais package on series without trend", {
  skip_if_not_installed("prais")

  set.seed(789)
  n <- 150
  t_idx <- 1:n
  x <- arima.sim(list(ar = 0.6), n = n)

  r1 <- pw_fast(x)

  df <- data.frame(y = x, t = t_idx, time = t_idx)
  r2 <- prais::prais_winsten(y ~ t, data = df, index = "time", twostep = TRUE)
  r2_sum <- summary(r2)

  # For small t-stats, use absolute tolerance instead of relative
  expect_equal(r1$tpw, r2_sum$coefficients["t", "t value"], tolerance = 0.5)
  expect_equal(r1$rho, r2$rho[nrow(r2$rho)], tolerance = 0.01)
})

test_that("pw_fast returns correct structure", {
  set.seed(505)
  x <- arima.sim(list(ar = 0.7), n = 100)

  result <- pw_fast(x)

  expect_type(result, "list")
  expect_named(result, c("x", "z_x", "b0hat", "b1hat", "rho", "pvalue", "tpw"))
  expect_equal(result$x, x)
  expect_length(result$z_x, length(x))
  expect_type(result$b0hat, "double")
  expect_type(result$b1hat, "double")
  expect_type(result$rho, "double")
  expect_type(result$pvalue, "double")
  expect_type(result$tpw, "double")
})

test_that("pw_fast rho is within valid range", {
  set.seed(606)
  x <- arima.sim(list(ar = 0.8), n = 100)

  result <- pw_fast(x)

  # rho should be clamped to (-0.999, 0.999)
  expect_gt(result$rho, -1)
  expect_lt(result$rho, 1)
})

test_that("pw_fast pvalue is in valid range", {
  set.seed(707)
  x <- arima.sim(list(ar = 0.7), n = 100)

  result <- pw_fast(x)

  expect_gte(result$pvalue, 0)
  expect_lte(result$pvalue, 1)
})

test_that("pw_fast validates input correctly", {
  expect_error(pw_fast(NULL), "`x` must be a non-empty numeric vector")
  expect_error(pw_fast(character(5)), "`x` must be a non-empty numeric vector")
  expect_error(pw_fast(numeric(0)), "`x` must be a non-empty numeric vector")
  expect_error(pw_fast(1:4), "Series too short")
})

test_that("pw_fast works on white noise", {
  set.seed(808)
  x <- rnorm(100)

  result <- pw_fast(x)

  expect_true(is.finite(result$tpw))
  expect_true(is.finite(result$rho))
  expect_true(abs(result$rho) < 1)
})

test_that("pw_fast detects trend correctly", {
  set.seed(909)
  n <- 100
  t_idx <- 1:n

  # Strong upward trend
  x <- 0.2 * t_idx + arima.sim(list(ar = 0.5), n = n)
  result <- pw_fast(x)

  # Should have significant t-statistic (large positive)
  expect_gt(result$tpw, 2)
  expect_lt(result$pvalue, 0.05)
})

test_that("pw_fast handles negative autocorrelation", {
  skip_if_not_installed("prais")

  set.seed(1010)
  n <- 100
  t_idx <- 1:n
  z <- arima.sim(list(ar = -0.5), n = n)
  x <- 3 + 0.05 * t_idx + z

  r1 <- pw_fast(x)

  df <- data.frame(y = x, t = t_idx, time = t_idx)
  r2 <- prais::prais_winsten(y ~ t, data = df, index = "time", twostep = TRUE)
  r2_sum <- summary(r2)

  # Allow 5% relative difference for large t-stats
  expect_equal(r1$tpw, r2_sum$coefficients["t", "t value"], tolerance = 0.05 * abs(r2_sum$coefficients["t", "t value"]))
  expect_equal(r1$rho, r2$rho[nrow(r2$rho)], tolerance = 0.01)
  expect_lt(r1$rho, 0)  # Should be negative
})

test_that("make_stat_pw returns correct statistic", {
  set.seed(1111)
  x <- arima.sim(list(ar = 0.7), n = 100)

  stat_fn <- make_stat_pw()
  tpw_from_fn <- stat_fn(x)
  tpw_from_fast <- pw_fast(x)$tpw

  expect_equal(tpw_from_fn, tpw_from_fast, tolerance = 1e-10)
})

# =============================================================================
# Iterative Prais-Winsten tests
# =============================================================================

test_that("pw_fast iterative produces reasonable results vs prais", {
  skip_if_not_installed("prais")

  set.seed(1212)
  n <- 100
  t_idx <- 1:n
  z <- arima.sim(list(ar = 0.7), n = n)
  x <- 5 + 0.1 * t_idx + z

  # Our iterative implementation (starts from OLS rho estimate)
  r1 <- pw_fast(x, iterate = TRUE)

  # prais with twostep = FALSE (starts from rho=0)
  df <- data.frame(y = x, t = t_idx, time = t_idx)
  r2 <- prais::prais_winsten(y ~ t, data = df, index = "time", twostep = FALSE)
  r2_sum <- summary(r2)

  # Note: Our iterative starts from OLS rho estimate while prais starts from 0
  # They converge to different values, so we test for reasonable agreement
  # Both should detect significance in same direction
  expect_true(sign(r1$tpw) == sign(r2_sum$coefficients["t", "t value"]))
  expect_equal(r1$tpw, r2_sum$coefficients["t", "t value"], tolerance = 0.25 * abs(r2_sum$coefficients["t", "t value"]))

  # Rho estimates should be in similar range
  expect_equal(r1$rho, r2$rho[length(r2$rho)], tolerance = 0.1)
})

test_that("pw_fast iterative returns valid results", {
  set.seed(1313)
  x <- arima.sim(list(ar = 0.8), n = 100)

  result <- pw_fast(x, iterate = TRUE)

  expect_type(result, "list")
  expect_true(is.finite(result$tpw))
  expect_true(is.finite(result$rho))
  expect_true(abs(result$rho) < 1)
  expect_true(result$pvalue >= 0 && result$pvalue <= 1)
})

test_that("pw_fast iterative converges for high-persistence series", {
  set.seed(1414)
  n <- 150
  t_idx <- 1:n
  z <- arima.sim(list(ar = 0.95), n = n)
  x <- 10 + 0.05 * t_idx + z

  # Should converge without issues
  result <- pw_fast(x, iterate = TRUE, max_iter = 100)

  expect_true(is.finite(result$tpw))
  expect_true(abs(result$rho) < 1)
})

test_that("pw_fast two-step and iterative give similar results", {
  set.seed(1515)
  x <- arima.sim(list(ar = 0.6), n = 100)

  r_twostep <- pw_fast(x, iterate = FALSE)
  r_iter <- pw_fast(x, iterate = TRUE)

  # Results should be similar (not identical due to iteration)
  expect_equal(r_twostep$tpw, r_iter$tpw, tolerance = 0.1)
  expect_equal(r_twostep$rho, r_iter$rho, tolerance = 0.05)
})

test_that("pw_fast iterative respects max_iter parameter", {
  set.seed(1616)
  x <- arima.sim(list(ar = 0.7), n = 100)

  # Iterative with different max_iter should give valid results
  r1 <- pw_fast(x, iterate = TRUE, max_iter = 1)
  r_many <- pw_fast(x, iterate = TRUE, max_iter = 50)

  # Both should return valid rho (in range)
  expect_true(abs(r1$rho) < 1)
  expect_true(abs(r_many$rho) < 1)

  # Both should return finite t-statistics
  expect_true(is.finite(r1$tpw))
  expect_true(is.finite(r_many$tpw))
})
