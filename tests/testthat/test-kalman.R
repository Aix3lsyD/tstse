# Tests for kalman()
# Kalman filter and smoother wrapper around astsa::Ksmooth

test_that("kalman matches tswge::kalman.wge output", {
  skip_if_not_installed("tswge")
  skip_if_not_installed("astsa")

  set.seed(123)
  n <- 50
  y <- cumsum(rnorm(n)) + rnorm(n, sd = 0.5)

  # Parameters
  x0 <- 0
  P0 <- 1
  Phi <- 1
  Q <- 1
  A <- 1
  R <- 0.25

  result <- kalman(y = y, x0 = x0, P0 = P0, Phi = Phi, Q = Q, A = A, R = R)

  # tswge uses different parameter names
  expected <- tswge::kalman.wge(
    y = y,
    start = x0,
    gam0 = P0,
    F = Phi,
    gamV = Q,
    G = A,
    gamW = R
  )

  # Compare predictions
  expect_equal(result$predicted$state, expected[, "Prediction"], tolerance = 1e-10)
  expect_equal(result$predicted$var, expected[, "Var_Predict"], tolerance = 1e-10)

  # Compare filtered
  expect_equal(result$filtered$state, expected[, "Filter"], tolerance = 1e-10)
  expect_equal(result$filtered$var, expected[, "Var_Filter"], tolerance = 1e-10)

  # Compare smoothed
  expect_equal(result$smoothed$state, expected[, "Smooth"], tolerance = 1e-10)
  expect_equal(result$smoothed$var, expected[, "Var_Smooth"], tolerance = 1e-10)
})


test_that("kalman matches tswge::kalman.wge with different parameters", {
  skip_if_not_installed("tswge")
  skip_if_not_installed("astsa")

  set.seed(456)
  y <- sin(seq(0, 4 * pi, length.out = 100)) + rnorm(100, sd = 0.3)

  result <- kalman(y = y, x0 = 0, P0 = 10, Phi = 0.9, Q = 0.5, A = 1, R = 0.09)

  expected <- tswge::kalman.wge(
    y = y, start = 0, gam0 = 10, F = 0.9, gamV = 0.5, G = 1, gamW = 0.09
  )

  expect_equal(result$filtered$state, expected[, "Filter"], tolerance = 1e-10)
  expect_equal(result$smoothed$state, expected[, "Smooth"], tolerance = 1e-10)
})


test_that("kalman returns correct structure", {
  skip_if_not_installed("astsa")

  y <- rnorm(30)
  result <- kalman(y = y, x0 = 0, P0 = 1, Phi = 1, Q = 1, A = 1, R = 1)

  expect_s3_class(result, "kalman")
  expect_named(result, c("y", "predicted", "filtered", "smoothed", "result"))

  expect_equal(result$y, y)

  expect_named(result$predicted, c("state", "var"))
  expect_named(result$filtered, c("state", "var"))
  expect_named(result$smoothed, c("state", "var"))

  expect_length(result$predicted$state, 30)
  expect_length(result$predicted$var, 30)
  expect_length(result$filtered$state, 30)
  expect_length(result$filtered$var, 30)
  expect_length(result$smoothed$state, 30)
  expect_length(result$smoothed$var, 30)
})


test_that("kalman filtered variance <= predicted variance", {
  skip_if_not_installed("astsa")

  # Filtering should reduce uncertainty compared to prediction
  set.seed(789)
  y <- cumsum(rnorm(50)) + rnorm(50)

  result <- kalman(y = y, x0 = 0, P0 = 10, Phi = 1, Q = 1, A = 1, R = 1)

  # After first observation, filtered variance should be less than predicted
  expect_true(all(result$filtered$var[-1] <= result$predicted$var[-1] + 1e-10))
})


test_that("kalman smoothed variance <= filtered variance", {
  skip_if_not_installed("astsa")

  # Smoothing uses more information, so should have lower variance
  set.seed(101)
  y <- cumsum(rnorm(50)) + rnorm(50)

  result <- kalman(y = y, x0 = 0, P0 = 10, Phi = 1, Q = 1, A = 1, R = 1)

  # Smoothed variance should be <= filtered variance (using all data)
  expect_true(all(result$smoothed$var <= result$filtered$var + 1e-10))
})


test_that("kalman with zero observation noise gives exact observations", {
  skip_if_not_installed("astsa")

  # If R = 0 (no observation noise), filtered state should equal observations
  y <- c(1, 2, 3, 4, 5)

  result <- kalman(
    y = y, x0 = 0, P0 = 100,
    Phi = 1, Q = 1, A = 1, R = 1e-10  # Near-zero R
  )

  # Filtered should be very close to observations
  expect_equal(result$filtered$state, y, tolerance = 1e-5)
})


test_that("kalman with large initial variance converges", {
  skip_if_not_installed("astsa")

  set.seed(111)
  y <- rep(5, 50) + rnorm(50, sd = 0.1)

  # Very uncertain initial state
  result <- kalman(y = y, x0 = 0, P0 = 1000, Phi = 1, Q = 0.01, A = 1, R = 0.01)

  # Should converge to observations
  expect_true(abs(result$filtered$state[50] - 5) < 0.5)
})


test_that("kalman handles two observations", {
  skip_if_not_installed("astsa")

  # astsa::Ksmooth requires at least 2 observations
  y <- c(3, 5)
  result <- kalman(y = y, x0 = 0, P0 = 1, Phi = 1, Q = 1, A = 1, R = 1)

  expect_length(result$filtered$state, 2)
  expect_length(result$smoothed$state, 2)
})


test_that("kalman errors with single observation", {
  skip_if_not_installed("astsa")

  # astsa::Ksmooth doesn't handle n=1, so we should error gracefully
  # or document this limitation
  y <- 5
  expect_error(kalman(y = y, x0 = 0, P0 = 1, Phi = 1, Q = 1, A = 1, R = 1))
})


test_that("kalman preserves raw astsa result", {
  skip_if_not_installed("astsa")

  y <- rnorm(20)
  result <- kalman(y = y, x0 = 0, P0 = 1, Phi = 1, Q = 1, A = 1, R = 1)

  # Should be able to access raw result
  expect_true(!is.null(result$result))
  expect_true(!is.null(result$result$Xf))
})


test_that("kalman validates y parameter", {
  skip_if_not_installed("astsa")

  expect_error(kalman(y = "abc", x0 = 0, P0 = 1, Phi = 1, Q = 1, A = 1, R = 1),
               "`y` must be a numeric")
  expect_error(kalman(y = numeric(0), x0 = 0, P0 = 1, Phi = 1, Q = 1, A = 1, R = 1),
               "at least 2 observations")
  expect_error(kalman(y = 5, x0 = 0, P0 = 1, Phi = 1, Q = 1, A = 1, R = 1),
               "at least 2 observations")
})


test_that("kalman validates numeric parameters", {
  skip_if_not_installed("astsa")

  y <- rnorm(10)
  expect_error(kalman(y = y, x0 = "a", P0 = 1, Phi = 1, Q = 1, A = 1, R = 1),
               "`x0` must be numeric")
  expect_error(kalman(y = y, x0 = 0, P0 = "b", Phi = 1, Q = 1, A = 1, R = 1),
               "`P0` must be numeric")
  expect_error(kalman(y = y, x0 = 0, P0 = 1, Phi = "c", Q = 1, A = 1, R = 1),
               "`Phi` must be numeric")
  expect_error(kalman(y = y, x0 = 0, P0 = 1, Phi = 1, Q = "d", A = 1, R = 1),
               "`Q` must be numeric")
  expect_error(kalman(y = y, x0 = 0, P0 = 1, Phi = 1, Q = 1, A = "e", R = 1),
               "`A` must be numeric")
  expect_error(kalman(y = y, x0 = 0, P0 = 1, Phi = 1, Q = 1, A = 1, R = "f"),
               "`R` must be numeric")
})


test_that("kalman validates variance parameters are non-negative", {
  skip_if_not_installed("astsa")

  y <- rnorm(10)
  expect_error(kalman(y = y, x0 = 0, P0 = 1, Phi = 1, Q = -1, A = 1, R = 1),
               "`Q`.*non-negative")
  expect_error(kalman(y = y, x0 = 0, P0 = 1, Phi = 1, Q = 1, A = 1, R = -1),
               "`R`.*non-negative")
  expect_error(kalman(y = y, x0 = 0, P0 = -1, Phi = 1, Q = 1, A = 1, R = 1),
               "`P0`.*non-negative")
})

test_that("kalman print method works", {
  skip_if_not_installed("astsa")

  y <- rnorm(25)
  result <- kalman(y = y, x0 = 0, P0 = 1, Phi = 1, Q = 1, A = 1, R = 1)

  expect_output(print(result), "Kalman Filter/Smoother")
  expect_output(print(result), "Observations: 25")
  expect_output(print(result), "predicted")
  expect_output(print(result), "filtered")
  expect_output(print(result), "smoothed")
})


test_that("kalman local level model tracks random walk", {
  skip_if_not_installed("astsa")

  # Generate a true random walk
  set.seed(222)
  n <- 100
  x_true <- cumsum(rnorm(n, sd = 1))  # True state
  y <- x_true + rnorm(n, sd = 2)       # Noisy observations

  # Fit with correct model specification
  result <- kalman(
    y = y,
    x0 = 0,
    P0 = 10,
    Phi = 1,     # Random walk
    Q = 1,       # True state noise variance
    A = 1,       # Direct observation
    R = 4        # True observation noise variance
  )

  # Smoothed estimate should be closer to truth than observations
  obs_mse <- mean((y - x_true)^2)
  smooth_mse <- mean((result$smoothed$state - x_true)^2)

  expect_true(smooth_mse < obs_mse)
})


test_that("kalman AR(1) model works", {
  skip_if_not_installed("astsa")

  # AR(1) state: x_t = 0.8 * x_{t-1} + v_t
  set.seed(333)
  n <- 100
  phi <- 0.8
  x_true <- numeric(n)
  x_true[1] <- rnorm(1)
  for (t in 2:n) {
    x_true[t] <- phi * x_true[t - 1] + rnorm(1, sd = 0.5)
  }
  y <- x_true + rnorm(n, sd = 1)

  result <- kalman(
    y = y,
    x0 = 0,
    P0 = 5,
    Phi = phi,
    Q = 0.25,
    A = 1,
    R = 1
  )

  # Should produce reasonable estimates
  expect_false(any(is.na(result$filtered$state)))
  expect_false(any(is.na(result$smoothed$state)))
})
