# Tests for kalman_miss()
# Kalman filter and smoother with missing observations

test_that("kalman_miss matches tswge::kalman.miss.wge output", {
  skip_if_not_installed("tswge")
  skip_if_not_installed("astsa")

  set.seed(123)
  n <- 50

  # Create data with "missing" observations encoded in time-varying A
  y <- cumsum(rnorm(n)) + rnorm(n, sd = 0.5)

  # Time-varying observation matrix (1x1xn)
  # Set some to 0 to indicate "missing"
  A <- array(1, dim = c(1, 1, n))
  missing_idx <- c(10, 11, 25, 26, 27)
  A[, , missing_idx] <- 0

  # Parameters
  x0 <- 0
  P0 <- 1
  Phi <- 1
  Q <- 1
  R <- 0.25

  result <- kalman_miss(
    y = y,
    x0 = x0,
    P0 = P0,
    Phi = Phi,
    Q = Q,
    A = A,
    R = R
  )

  # tswge uses different parameter names
  expected <- tswge::kalman.miss.wge(
    y = y,
    start = x0,
    gam0 = P0,
    F = Phi,
    gamV = Q,
    Gtmiss = A,
    gamW = R
  )

  # Compare predictions
  expect_equal(result$predicted$state, expected[, "Prediction"],
               tolerance = 1e-10)
  expect_equal(result$filtered$state, expected[, "Filter"],
               tolerance = 1e-10)
  expect_equal(result$smoothed$state, expected[, "Smooth"],
               tolerance = 1e-10)
})


test_that("kalman_miss returns correct structure", {
  skip_if_not_installed("astsa")

  n <- 30
  y <- rnorm(n)
  A <- array(1, dim = c(1, 1, n))

  result <- kalman_miss(
    y = y,
    x0 = 0,
    P0 = 1,
    Phi = 1,
    Q = 1,
    A = A,
    R = 1
  )

  expect_s3_class(result, "kalman")
  expect_named(result, c("y", "predicted", "filtered", "smoothed", "result"))
  expect_named(result$predicted, c("state", "var"))
  expect_named(result$filtered, c("state", "var"))
  expect_named(result$smoothed, c("state", "var"))
})


test_that("kalman_miss handles missing observations correctly", {
  skip_if_not_installed("astsa")

  set.seed(456)
  n <- 50

  # True state
  x_true <- cumsum(rnorm(n, sd = 1))
  y_full <- x_true + rnorm(n, sd = 0.5)

  # Create version with "missing" data
  A_full <- array(1, dim = c(1, 1, n))
  A_miss <- A_full
  missing_idx <- 20:30
  A_miss[, , missing_idx] <- 0

  # Fit both models
  result_full <- kalman_miss(
    y = y_full,
    x0 = 0, P0 = 10, Phi = 1,
    Q = 1, A = A_full, R = 0.25
  )

  result_miss <- kalman_miss(
    y = y_full,
    x0 = 0, P0 = 10, Phi = 1,
    Q = 1, A = A_miss, R = 0.25
  )

  # At missing times, filtered variance should be higher
  var_full <- result_full$filtered$var
  var_miss <- result_miss$filtered$var

  # Variance should be higher when observations are "missing"
  expect_true(all(var_miss[missing_idx] >= var_full[missing_idx] - 1e-10))
})


test_that("kalman_miss with no missing data matches kalman", {
  skip_if_not_installed("astsa")

  set.seed(789)
  n <- 40
  y <- rnorm(n)

  # All observations present
  A <- array(1, dim = c(1, 1, n))

  result_miss <- kalman_miss(
    y = y,
    x0 = 0,
    P0 = 1,
    Phi = 1,
    Q = 1,
    A = A,
    R = 1
  )

  result_reg <- kalman(
    y = y,
    x0 = 0,
    P0 = 1,
    Phi = 1,
    Q = 1,
    A = 1,
    R = 1
  )

  # Results should match
  expect_equal(result_miss$filtered$state,
               result_reg$filtered$state, tolerance = 1e-10)
  expect_equal(result_miss$smoothed$state,
               result_reg$smoothed$state, tolerance = 1e-10)
})


test_that("kalman_miss validates y parameter", {
  skip_if_not_installed("astsa")

  A <- array(1, dim = c(1, 1, 10))
  expect_error(
    kalman_miss(y = "abc", x0 = 0, P0 = 1, Phi = 1,
                Q = 1, A = A, R = 1),
    "`y` must be numeric"
  )
})


test_that("kalman_miss validates minimum observations", {
  skip_if_not_installed("astsa")

  A <- array(1, dim = c(1, 1, 1))
  expect_error(
    kalman_miss(y = 5, x0 = 0, P0 = 1,
                Phi = 1, Q = 1, A = A, R = 1),
    "at least 2"
  )
})


test_that("kalman_miss validates A is 3D array", {
  skip_if_not_installed("astsa")

  n <- 10
  y <- rnorm(n)

  # 2D matrix instead of 3D array
  expect_error(
    kalman_miss(y = y, x0 = 0, P0 = 1, Phi = 1,
                Q = 1, A = matrix(1, 1, n), R = 1),
    "3-dimensional array"
  )
})


test_that("kalman_miss validates A dimension matches y", {
  skip_if_not_installed("astsa")

  n <- 10
  y <- rnorm(n)
  A <- array(1, dim = c(1, 1, 5))  # Wrong length

  expect_error(
    kalman_miss(y = y, x0 = 0, P0 = 1, Phi = 1,
                Q = 1, A = A, R = 1),
    "Third dimension of `A`"
  )
})

test_that("kalman_miss smoothed variance <= filtered variance", {
  skip_if_not_installed("astsa")

  set.seed(101)
  n <- 50
  y <- cumsum(rnorm(n)) + rnorm(n)
  A <- array(1, dim = c(1, 1, n))

  result <- kalman_miss(
    y = y,
    x0 = 0,
    P0 = 10,
    Phi = 1,
    Q = 1,
    A = A,
    R = 1
  )

  var_filt <- result$filtered$var
  var_smooth <- result$smoothed$var

  expect_true(all(var_smooth <= var_filt + 1e-10))
})


test_that("kalman_miss print method works via kalman class", {
  skip_if_not_installed("astsa")

  n <- 20
  y <- rnorm(n)
  A <- array(1, dim = c(1, 1, n))

  result <- kalman_miss(
    y = y, x0 = 0, P0 = 1, Phi = 1,
    Q = 1, A = A, R = 1
  )

  # Uses the kalman class print method
  expect_output(print(result), "Kalman Filter/Smoother")
})
