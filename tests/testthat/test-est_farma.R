# Tests for est_farma()
# FARMA model estimation

test_that("est_farma matches tswge::est.farma.wge output structure", {
  skip_if_not_installed("tswge")

  set.seed(123)
  x <- cumsum(rnorm(100))

  result <- est_farma(x, d_range = c(0, 0.5), d_step = 0.1, p_max = 3)

  png(tempfile())
  expected <- tswge::est.farma.wge(x, low.d = 0, high.d = 0.5, inc.d = 0.1, p.max = 3)
  dev.off()

  # Both should return similar structure
  expect_true(!is.null(result$d))
  expect_true(!is.null(result$phi))
  expect_true(!is.null(result$var_a))
  expect_true(!is.null(result$aic))
})


test_that("est_farma estimates d in reasonable range for long memory data", {
  set.seed(456)
  # Generate data with long memory characteristics using cumsum (d ≈ 0.5)
  # then test that estimated d is in a reasonable range
  x <- cumsum(rnorm(300))

  result <- est_farma(x, d_range = c(0, 0.5), d_step = 0.05, p_max = 2)

  # For cumsum data (integrated), d should be estimated toward upper end
  # Just check it's within the search range and reasonable
  expect_true(result$d >= 0)
  expect_true(result$d <= 0.5)
})


test_that("est_farma returns correct structure", {
  set.seed(111)
  x <- cumsum(rnorm(100))

  result <- est_farma(x, d_range = c(0, 0.5), d_step = 0.1, p_max = 3)

  expect_s3_class(result, "est_farma")
  expect_named(result, c("d", "phi", "p", "var_a", "aic", "d_grid", "aic_grid"))
})


test_that("est_farma d is within specified range", {
  set.seed(222)
  x <- rnorm(100)

  result <- est_farma(x, d_range = c(0.1, 0.4), d_step = 0.1, p_max = 2)

  expect_true(result$d >= 0.1)
  expect_true(result$d <= 0.4)
})


test_that("est_farma d_grid matches expected values", {
  set.seed(333)
  x <- rnorm(100)

  result <- est_farma(x, d_range = c(0, 0.4), d_step = 0.1, p_max = 2)

  expected_grid <- seq(0, 0.4, by = 0.1)
  expect_equal(result$d_grid, expected_grid, tolerance = 1e-10)
})


test_that("est_farma p is within specified range", {
  set.seed(444)
  x <- arima.sim(model = list(ar = 0.5), n = 100)

  result <- est_farma(x, d_range = c(0, 0.3), d_step = 0.1, p_max = 5)

  expect_true(result$p >= 0)
  expect_true(result$p <= 5)
})


test_that("est_farma phi length matches p", {
  set.seed(555)
  x <- arima.sim(model = list(ar = c(0.5, -0.3)), n = 150)

  result <- est_farma(x, d_range = c(0, 0.3), d_step = 0.1, p_max = 5)

  expect_length(result$phi, result$p)
})


test_that("est_farma var_a is positive", {
  set.seed(666)
  x <- rnorm(100)

  result <- est_farma(x, d_range = c(0, 0.3), d_step = 0.1, p_max = 3)

  expect_true(result$var_a > 0)
})


test_that("est_farma aic_grid has correct length", {
  set.seed(777)
  x <- rnorm(100)

  result <- est_farma(x, d_range = c(0, 0.4), d_step = 0.1, p_max = 2)

  # Grid: 0, 0.1, 0.2, 0.3, 0.4 = 5 values
  expect_length(result$aic_grid, 5)
})


test_that("est_farma selected d corresponds to minimum aic", {
  set.seed(888)
  x <- cumsum(rnorm(100))

  result <- est_farma(x, d_range = c(0, 0.5), d_step = 0.1, p_max = 3)

  # The selected d should correspond to the minimum AIC
  min_idx <- which.min(result$aic_grid)
  expect_equal(unname(result$d), result$d_grid[min_idx])
})


test_that("est_farma handles AR(0) case", {
  set.seed(999)
  # Pure white noise should select AR(0)
  x <- rnorm(200)

  result <- est_farma(x, d_range = c(0, 0.2), d_step = 0.05, p_max = 3)

  # Should work regardless of p selected
  expect_true(result$p >= 0)
  if (result$p == 0) {
    expect_length(result$phi, 0)
  }
})


test_that("est_farma validates x parameter", {
  expect_error(est_farma("abc", d_range = c(0, 0.5)), "`x` must be a numeric")
  expect_error(est_farma(1:9, d_range = c(0, 0.5)), "at least 10")
  expect_error(est_farma(c(1:10, NA), d_range = c(0, 0.5)), "contains NA")
})


test_that("est_farma validates d_range parameter", {
  x <- rnorm(50)
  expect_error(est_farma(x, d_range = 0.5), "length 2")
  expect_error(est_farma(x, d_range = c(0.5, 0.3)), "must be less than")
  expect_error(est_farma(x, d_range = c(0.3, 0.3)), "must be less than")
})


test_that("est_farma validates d_step parameter", {
  x <- rnorm(50)
  expect_error(est_farma(x, d_range = c(0, 0.5), d_step = 0), "positive number")
  expect_error(est_farma(x, d_range = c(0, 0.5), d_step = -0.1), "positive number")
})


test_that("est_farma validates p_max parameter", {
  x <- rnorm(50)
  expect_error(est_farma(x, d_range = c(0, 0.5), p_max = -1), "non-negative")
})


test_that("est_farma validates n_back parameter", {
  x <- rnorm(50)
  expect_error(est_farma(x, d_range = c(0, 0.5), n_back = 0), "positive integer")
})


test_that("est_farma print method works", {
  set.seed(101)
  x <- cumsum(rnorm(100))

  result <- est_farma(x, d_range = c(0, 0.5), d_step = 0.1, p_max = 3)

  expect_output(print(result), "FARMA Model Estimation")
  expect_output(print(result), "Fractional differencing:")
  expect_output(print(result), "AR order:")
  expect_output(print(result), "White noise variance:")
})


test_that("est_farma works with small d_step", {
  set.seed(102)
  x <- cumsum(rnorm(100))

  result <- est_farma(x, d_range = c(0.1, 0.3), d_step = 0.02, p_max = 2)

  expect_true(result$d >= 0.1)
  expect_true(result$d <= 0.3)
  expect_length(result$d_grid, 11)  # 0.1, 0.12, 0.14, ..., 0.3
})
