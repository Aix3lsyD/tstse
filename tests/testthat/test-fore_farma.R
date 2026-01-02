# tests/testthat/test-fore_farma.R

test_that("fore_farma returns correct structure", {
  set.seed(123)
  x <- cumsum(rnorm(100))

  result <- fore_farma(x, d = 0.3, n_ahead = 5, plot = FALSE)

  expect_type(result, "list")
  expect_named(result, c("farma_fore", "ar_fore", "ar_fit_order"))
  expect_length(result$farma_fore, 5)
  expect_length(result$ar_fore, 5)
  expect_true(is.numeric(result$ar_fit_order))
})

test_that("fore_farma handles different n_ahead values", {
  set.seed(123)
  x <- cumsum(rnorm(100))

  result5 <- fore_farma(x, d = 0.3, n_ahead = 5, plot = FALSE)
  result10 <- fore_farma(x, d = 0.3, n_ahead = 10, plot = FALSE)
  result20 <- fore_farma(x, d = 0.3, n_ahead = 20, plot = FALSE)

  expect_length(result5$farma_fore, 5)
  expect_length(result10$farma_fore, 10)
  expect_length(result20$farma_fore, 20)
})

test_that("fore_farma handles lastn parameter", {
  set.seed(123)
  x <- cumsum(rnorm(100))

  # Future forecast mode (lastn = FALSE)
  result_future <- fore_farma(x, d = 0.3, n_ahead = 5, lastn = FALSE, plot = FALSE)
  expect_length(result_future$farma_fore, 5)

  # Holdout mode (lastn = TRUE)
  result_holdout <- fore_farma(x, d = 0.3, n_ahead = 5, lastn = TRUE, plot = FALSE)
  expect_length(result_holdout$farma_fore, 5)
})

test_that("fore_farma handles AR coefficients", {
  set.seed(123)
  x <- cumsum(rnorm(100))

  result <- fore_farma(x, d = 0.3, phi = 0.5, n_ahead = 5, plot = FALSE)

  expect_length(result$farma_fore, 5)
  expect_true(all(is.finite(result$farma_fore)))
})

test_that("fore_farma handles MA coefficients", {
  set.seed(123)
  x <- cumsum(rnorm(100))

  result <- fore_farma(x, d = 0.3, theta = 0.3, n_ahead = 5, plot = FALSE)

  expect_length(result$farma_fore, 5)
  expect_true(all(is.finite(result$farma_fore)))
})

test_that("fore_farma handles ARMA coefficients", {
  set.seed(123)
  x <- cumsum(rnorm(100))

  result <- fore_farma(x, d = 0.3, phi = 0.5, theta = 0.3, n_ahead = 5, plot = FALSE)

  expect_length(result$farma_fore, 5)
  expect_true(all(is.finite(result$farma_fore)))
})

test_that("fore_farma validates inputs", {
  set.seed(123)
  x <- cumsum(rnorm(100))

  expect_error(fore_farma(1:5, d = 0.3),
               "`x` must be a numeric vector with at least 10 observations")
  expect_error(fore_farma(x, d = "a"),
               "`d` must be a single numeric value")
  expect_error(fore_farma(x, d = c(0.3, 0.4)),
               "`d` must be a single numeric value")
  expect_error(fore_farma(x, d = 0.3, phi = "a"),
               "`phi` must be numeric")
  expect_error(fore_farma(x, d = 0.3, n_ahead = 0),
               "`n_ahead` must be a positive integer")
})

test_that("fore_farma returns invisibly", {
  set.seed(123)
  x <- cumsum(rnorm(100))

  expect_invisible(fore_farma(x, d = 0.3, n_ahead = 5))
})

test_that("fore_farma matches tswge::fore.farma.wge output", {
  skip_if_not_installed("tswge")

  set.seed(123)
  x <- cumsum(rnorm(100))

  # Test future forecast mode
  set.seed(456)
  tstse_result <- fore_farma(x, d = 0.3, phi = 0.5, n_ahead = 10,
                             lastn = FALSE, plot = FALSE)
  set.seed(456)
  tswge_result <- tswge::fore.farma.wge(x, d = 0.3, phi = 0.5, n.ahead = 10,
                                        lastn = FALSE, plot = FALSE)

  # Compare FARMA forecasts (tolerance accounts for stochastic arima optimization)
  expect_equal(tstse_result$farma_fore, tswge_result$farma.fore, tolerance = 0.01)

  # Compare AR forecasts
  expect_equal(tstse_result$ar_fore, tswge_result$ar.fore, tolerance = 1e-8)

  # Compare AR order
  expect_equal(tstse_result$ar_fit_order, tswge_result$ar.fit.order)
})

test_that("fore_farma lastn mode matches tswge", {
  skip_if_not_installed("tswge")

  set.seed(123)
  x <- cumsum(rnorm(100))

  # Test holdout mode
  set.seed(789)
  tstse_result <- fore_farma(x, d = 0.4, n_ahead = 10, lastn = TRUE, plot = FALSE)
  set.seed(789)
  tswge_result <- tswge::fore.farma.wge(x, d = 0.4, n.ahead = 10, lastn = TRUE, plot = FALSE)

  expect_equal(tstse_result$farma_fore, tswge_result$farma.fore, tolerance = 1e-8)
  expect_equal(tstse_result$ar_fore, tswge_result$ar.fore, tolerance = 1e-8)
})

test_that("fore_farma produces reasonable forecasts for long-memory data", {
  set.seed(42)
  # Generate a simple long-memory-like series
  n <- 200
  x <- numeric(n)
  x[1] <- rnorm(1)
  for (i in 2:n) {
    x[i] <- 0.9 * x[i-1] + rnorm(1, sd = 0.5)
  }

  result <- fore_farma(x, d = 0.3, n_ahead = 10, lastn = FALSE, plot = FALSE)

  # Forecasts should be finite
  expect_true(all(is.finite(result$farma_fore)))
  expect_true(all(is.finite(result$ar_fore)))

  # Forecasts should not be too extreme
  x_range <- diff(range(x))
  expect_true(all(abs(result$farma_fore - mean(x)) < 3 * x_range))
})

test_that("fore_farma handles edge case d = 0", {
  set.seed(123)
  x <- rnorm(100)

  # With d = 0, should behave like regular ARMA forecast
  result <- fore_farma(x, d = 0, phi = 0.5, n_ahead = 5, plot = FALSE)

  expect_length(result$farma_fore, 5)
  expect_true(all(is.finite(result$farma_fore)))
})

test_that("fore_farma n_back parameter works", {
  set.seed(123)
  x <- cumsum(rnorm(100))

  # Different n_back values should work
  result1 <- fore_farma(x, d = 0.3, n_ahead = 5, n_back = 100, plot = FALSE)
  result2 <- fore_farma(x, d = 0.3, n_ahead = 5, n_back = 500, plot = FALSE)

  expect_length(result1$farma_fore, 5)
  expect_length(result2$farma_fore, 5)

  # Results may differ slightly due to backcasting differences
  expect_true(all(is.finite(result1$farma_fore)))
  expect_true(all(is.finite(result2$farma_fore)))
})
