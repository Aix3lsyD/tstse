library(tswge)

skip_if_not_installed("nnfor")


test_that("roll_win_rmse_nn returns correct structure", {

  skip_on_cran()

  set.seed(123)
  x <- rnorm(100)

  # Fit a simple MLP model
  suppressWarnings({
    fit <- nnfor::mlp(ts(x), hd = 2, reps = 1)
  })

  result <- roll_win_rmse_nn(x, horizon = 1, model = fit, plot = FALSE, verbose = FALSE)

  expect_type(result, "list")
  expect_named(result, c("rmse", "mse", "rmse_all", "mse_all",
                         "n_windows", "horizon", "frequency", "diff_order"))
  expect_true(is.numeric(result$rmse))
  expect_true(is.numeric(result$mse))
  expect_equal(length(result$rmse_all), result$n_windows)
  expect_equal(result$horizon, 1L)
})


test_that("roll_win_rmse_nn validates x input", {

  # Create a mock mlp object for validation tests
  mock_model <- structure(list(sdummy = FALSE, lags = 1:3, difforder = 0),
                          class = "mlp")

  expect_error(roll_win_rmse_nn("not numeric", model = mock_model),
               "`x` must be a non-empty numeric vector")

  expect_error(roll_win_rmse_nn(numeric(0), model = mock_model),
               "`x` must be a non-empty numeric vector")
})


test_that("roll_win_rmse_nn validates model input", {

  x <- rnorm(100)

  expect_error(roll_win_rmse_nn(x, model = "not a model"),
               "`model` must be an mlp object")

  expect_error(roll_win_rmse_nn(x, model = list(a = 1)),
               "`model` must be an mlp object")
})


test_that("roll_win_rmse_nn validates horizon", {

  mock_model <- structure(list(sdummy = FALSE, lags = 1:3, difforder = 0),
                          class = "mlp")

  x <- rnorm(100)

  expect_error(roll_win_rmse_nn(x, horizon = 0, model = mock_model),
               "`horizon` must be a positive integer")

  expect_error(roll_win_rmse_nn(x, horizon = -1, model = mock_model),
               "`horizon` must be a positive integer")

  expect_error(roll_win_rmse_nn(x, horizon = "one", model = mock_model),
               "`horizon` must be a positive integer")
})


test_that("roll_win_rmse_nn works with different horizons", {

  skip_on_cran()

  set.seed(456)
  x <- cumsum(rnorm(80))

  suppressWarnings({
    fit <- nnfor::mlp(ts(x), hd = 2, reps = 1)
  })

  result1 <- roll_win_rmse_nn(x, horizon = 1, model = fit, plot = FALSE)
  result5 <- roll_win_rmse_nn(x, horizon = 5, model = fit, plot = FALSE)

  expect_equal(result1$horizon, 1L)
  expect_equal(result5$horizon, 5L)

  # More windows with shorter horizon
  expect_gte(result1$n_windows, result5$n_windows)
})


test_that("roll_win_rmse_nn verbose mode prints messages", {

  skip_on_cran()

  set.seed(789)
  x <- rnorm(60)

  suppressWarnings({
    fit <- nnfor::mlp(ts(x), hd = 2, reps = 1)
  })

  expect_message(
    roll_win_rmse_nn(x, horizon = 1, model = fit, plot = FALSE, verbose = TRUE),
    "Training size"
  )
})


test_that("roll_win_rmse_nn handles short series", {

  set.seed(111)
  x <- rnorm(15)

  suppressWarnings({
    fit <- nnfor::mlp(ts(x), hd = 2, reps = 1, lags = 1:3)
  })

  # With horizon = 10, series might be too short
  expect_error(
    roll_win_rmse_nn(x, horizon = 10, model = fit, plot = FALSE),
    "too short"
  )
})


test_that("roll_win_rmse_nn RMSE is sqrt of MSE", {

  skip_on_cran()

  set.seed(222)
  x <- rnorm(80)

  suppressWarnings({
    fit <- nnfor::mlp(ts(x), hd = 2, reps = 1)
  })

  result <- roll_win_rmse_nn(x, horizon = 1, model = fit, plot = FALSE)

  # Check relationship between RMSE and MSE
  expect_equal(result$rmse_all, sqrt(result$mse_all), tolerance = 1e-10)
  expect_equal(result$rmse, mean(result$rmse_all), tolerance = 1e-10)
  expect_equal(result$mse, mean(result$mse_all), tolerance = 1e-10)
})



test_that("roll_win_rmse_nn plotting works", {

  skip_on_cran()

  set.seed(444)
  x <- rnorm(80)

  suppressWarnings({
    fit <- nnfor::mlp(ts(x), hd = 2, reps = 1)
  })

  # Should not error with plot = TRUE and multiple windows
  expect_no_error(
    roll_win_rmse_nn(x, horizon = 1, model = fit, plot = TRUE, verbose = FALSE)
  )
})


test_that("roll_win_rmse_nn handles ts objects", {

  skip_on_cran()

  set.seed(555)
  x <- ts(rnorm(80), frequency = 4)

  suppressWarnings({
    fit <- nnfor::mlp(x, hd = 2, reps = 1)
  })

  result <- roll_win_rmse_nn(x, horizon = 1, model = fit, plot = FALSE)

  expect_equal(result$frequency, 4)
})
