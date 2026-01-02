# Tests for roll_win_rmse() - Rolling Window RMSE

test_that("roll_win_rmse returns expected structure", {
  set.seed(123)
  x <- rnorm(50)

  result <- roll_win_rmse(x, horizon = 2, phi = 0.5, plot = FALSE)

  expect_type(result, "list")
  expect_named(result, c("rmse", "rmse_all", "n_windows", "training_size",
                         "horizon", "phi", "theta", "d", "s"))
  expect_true(result$rmse > 0)
  expect_true(result$n_windows > 0)
  expect_length(result$rmse_all, result$n_windows)
  expect_equal(result$horizon, 2L)
})

test_that("roll_win_rmse works with AR model", {
  set.seed(123)
  # Generate AR(1) data
  x <- gen_arma(n = 100, phi = 0.7, plot = FALSE)

  result <- roll_win_rmse(x, horizon = 3, phi = 0.7, plot = FALSE)

  # RMSE should be reasonable (not too large for well-specified model)
  expect_true(result$rmse < 3)
  expect_true(result$n_windows > 10)
})

test_that("roll_win_rmse works with ARMA model", {
  set.seed(456)
  x <- gen_arma(n = 100, phi = 0.5, theta = 0.3, plot = FALSE)

  result <- roll_win_rmse(x, horizon = 2, phi = 0.5, theta = 0.3, plot = FALSE)

  expect_true(result$rmse > 0)
  expect_true(is.finite(result$rmse))
})

test_that("roll_win_rmse works with differencing", {
  set.seed(789)
  # Random walk
  x <- cumsum(rnorm(80))

  result <- roll_win_rmse(x, horizon = 2, d = 1, plot = FALSE)

  expect_true(result$rmse > 0)
  expect_equal(result$d, 1L)
})

test_that("roll_win_rmse works with seasonality", {
  set.seed(111)
  # Simple seasonal pattern
  n <- 100
  x <- sin(2 * pi * (1:n) / 12) + rnorm(n, sd = 0.5)

  result <- roll_win_rmse(x, horizon = 3, s = 12, plot = FALSE)

  expect_true(result$rmse > 0)
  expect_equal(result$s, 12L)
})

test_that("roll_win_rmse computes correct number of windows", {
  set.seed(222)
  n <- 50
  x <- rnorm(n)
  horizon <- 5
  phi <- c(0.5, 0.3)  # AR(2)

  result <- roll_win_rmse(x, horizon = horizon, phi = phi, plot = FALSE)

  # For ARMA: training_size = max(p, q) + 1 = 2 + 1 = 3
  # n_windows = n - training_size - horizon + 1 = 50 - 3 - 5 + 1 = 43
  expect_equal(result$training_size, 3L)
  expect_equal(result$n_windows, 43L)
})

test_that("roll_win_rmse errors on short series", {
  x <- rnorm(5)

  expect_error(roll_win_rmse(x, horizon = 10, phi = 0.5, plot = FALSE),
               "Series too short")
})

test_that("roll_win_rmse validates input", {
  expect_error(roll_win_rmse(character(10)), "`x` must be a non-empty numeric")
  expect_error(roll_win_rmse(numeric(0)), "`x` must be a non-empty numeric")
  expect_error(roll_win_rmse(1:50, horizon = -1), "`horizon` must be a positive")
  expect_error(roll_win_rmse(1:50, d = -1), "`d` must be a non-negative")
})

test_that("roll_win_rmse handles pure white noise", {
  set.seed(333)
  x <- rnorm(60)

  # Fit with no AR/MA terms
  result <- roll_win_rmse(x, horizon = 2, plot = FALSE)

  expect_true(result$rmse > 0)
  # For white noise, RMSE should be close to 1 (sd of innovations)
  expect_true(result$rmse < 2)
})

test_that("roll_win_rmse parallel gives same results as sequential", {
  skip_on_cran()
  skip_if(parallel::detectCores() < 2, "Not enough cores")

  set.seed(444)
  x <- gen_arma(n = 80, phi = 0.6, plot = FALSE)

  result_seq <- roll_win_rmse(x, horizon = 3, phi = 0.6,
                               cores = 1, plot = FALSE)
  result_par <- roll_win_rmse(x, horizon = 3, phi = 0.6,
                               cores = 2, plot = FALSE)

  expect_equal(result_seq$rmse, result_par$rmse)
  expect_equal(result_seq$rmse_all, result_par$rmse_all)
})

test_that("roll_win_rmse verbose mode works", {
  set.seed(555)
  x <- rnorm(40)

  expect_message(
    roll_win_rmse(x, horizon = 2, phi = 0.5, verbose = TRUE, plot = FALSE),
    "Rolling Window"
  )
})

test_that("roll_win_rmse matches tswge for basic case", {
  skip_if_not_installed("tswge")

  set.seed(666)
  x <- rnorm(60) + 0.5 * c(0, head(rnorm(60), -1))  # AR(1)-like

  result_tstse <- roll_win_rmse(x, horizon = 2, phi = 0.5, plot = FALSE)

  # Compare with tswge
  result_tswge <- tswge::roll.win.rmse.wge(x, horizon = 2, phi = 0.5)

  # Mean RMSE should be similar
  expect_true(abs(result_tstse$rmse - result_tswge$rwRMSE) < 0.5)
  expect_equal(result_tstse$n_windows, result_tswge$numwindows)
})

test_that("roll_win_rmse handles AR(0) case", {
  set.seed(777)
  x <- rnorm(50)

  # phi = 0 means AR(0)
  result <- roll_win_rmse(x, horizon = 2, phi = 0, theta = 0, plot = FALSE)

  expect_true(result$rmse > 0)
  # Training size should be minimum (max(0, 0) + 1 = 1, but min is 3)
  expect_true(result$training_size >= 1)
})

test_that("roll_win_rmse individual RMSEs are all positive", {
  set.seed(888)
  x <- gen_arma(n = 60, phi = 0.5, plot = FALSE)

  result <- roll_win_rmse(x, horizon = 2, phi = 0.5, plot = FALSE)

  expect_true(all(result$rmse_all > 0))
  expect_true(all(is.finite(result$rmse_all)))
})
