library(tswge)

test_that("fore_arma returns correct structure", {

  set.seed(123)
  x <- gen_arma(n = 100, phi = 0.8, plot = FALSE)

  result <- fore_arma(x, phi = 0.8, n_ahead = 5, plot = FALSE)

  expect_s3_class(result, "fore_arma")
  expect_s3_class(result, "fore_arima")
  expect_true(all(c("f", "ll", "ul", "resid", "wnv", "xbar", "se", "psi") %in% names(result)))
})


test_that("fore_arma is equivalent to fore_arima with d=0, s=0", {

  set.seed(456)
  x <- gen_arma(n = 100, phi = 0.8, plot = FALSE)

  arma_result <- fore_arma(x, phi = 0.8, n_ahead = 5, plot = FALSE)
  arima_result <- fore_arima(x, phi = 0.8, d = 0, s = 0, n_ahead = 5, plot = FALSE)

  expect_equal(arma_result$f, arima_result$f)
  expect_equal(arma_result$ll, arima_result$ll)
  expect_equal(arma_result$ul, arima_result$ul)
  expect_equal(arma_result$se, arima_result$se)
})


test_that("fore_arma forecasts match tswge for AR(1)", {

  set.seed(789)
  x <- gen_arma(n = 100, phi = 0.8, plot = FALSE)

  old <- tswge::fore.arma.wge(x, phi = 0.8, n.ahead = 5, plot = FALSE)
  new <- fore_arma(x, phi = 0.8, n_ahead = 5, plot = FALSE)

  expect_equal(new$f, old$f, tolerance = 1e-6)
  expect_equal(new$ll, old$ll, tolerance = 1e-6)
  expect_equal(new$ul, old$ul, tolerance = 1e-6)
})


test_that("fore_arma forecasts match tswge for MA(1)", {

  set.seed(111)
  x <- gen_arma(n = 100, theta = 0.6, plot = FALSE)

  old <- tswge::fore.arma.wge(x, theta = 0.6, n.ahead = 5, plot = FALSE)
  new <- fore_arma(x, theta = 0.6, n_ahead = 5, plot = FALSE)

  expect_equal(new$f, old$f, tolerance = 1e-6)
})


test_that("fore_arma forecasts match tswge for ARMA(1,1)", {

  set.seed(222)
  x <- gen_arma(n = 100, phi = 0.7, theta = 0.4, plot = FALSE)

  old <- tswge::fore.arma.wge(x, phi = 0.7, theta = 0.4, n.ahead = 5, plot = FALSE)
  new <- fore_arma(x, phi = 0.7, theta = 0.4, n_ahead = 5, plot = FALSE)

  expect_equal(new$f, old$f, tolerance = 1e-6)
})


test_that("fore_arma computes RMSE and MAD in lastn mode", {

  set.seed(333)
  x <- gen_arma(n = 100, phi = 0.8, plot = FALSE)

  result <- fore_arma(x, phi = 0.8, n_ahead = 10, lastn = TRUE, plot = FALSE)

  expect_true("rmse" %in% names(result))
  expect_true("mad" %in% names(result))
  expect_true(is.numeric(result$rmse))
  expect_true(is.numeric(result$mad))
  expect_true(result$rmse >= 0)
  expect_true(result$mad >= 0)
})


test_that("fore_arma RMSE/MAD match tswge", {

  set.seed(444)
  x <- gen_arma(n = 100, phi = 0.8, plot = FALSE)

  old <- tswge::fore.arma.wge(x, phi = 0.8, n.ahead = 10, lastn = TRUE, plot = FALSE)
  new <- fore_arma(x, phi = 0.8, n_ahead = 10, lastn = TRUE, plot = FALSE)

  expect_equal(new$rmse, old$rmse, tolerance = 1e-6)
  expect_equal(new$mad, old$mad, tolerance = 1e-6)
})


test_that("fore_arma does not include RMSE/MAD when lastn=FALSE", {

  set.seed(555)
  x <- gen_arma(n = 100, phi = 0.8, plot = FALSE)

  result <- fore_arma(x, phi = 0.8, n_ahead = 5, lastn = FALSE, plot = FALSE)

  expect_null(result$rmse)
  expect_null(result$mad)
})


test_that("fore_arma RMSE/MAD are computed correctly", {

  set.seed(666)
  x <- gen_arma(n = 100, phi = 0.8, plot = FALSE)

  n_ahead <- 10
  result <- fore_arma(x, phi = 0.8, n_ahead = n_ahead, lastn = TRUE, plot = FALSE)

  # Manually compute expected values
  n <- length(x)
  actuals <- x[(n - n_ahead + 1):n]
  expected_rmse <- sqrt(mean((result$f - actuals)^2))
  expected_mad <- mean(abs(result$f - actuals))

  expect_equal(result$rmse, expected_rmse)
  expect_equal(result$mad, expected_mad)
})


test_that("fore_arma plotting works", {

  set.seed(777)
  x <- gen_arma(n = 100, phi = 0.8, plot = FALSE)

  expect_no_error(fore_arma(x, phi = 0.8, n_ahead = 5, plot = TRUE))
  expect_no_error(fore_arma(x, phi = 0.8, n_ahead = 5, plot = TRUE, limits = FALSE))
})


test_that("fore_arma print method works", {

  set.seed(888)
  x <- gen_arma(n = 100, phi = 0.8, plot = FALSE)

  result <- fore_arma(x, phi = 0.8, n_ahead = 5, plot = FALSE)
  expect_output(print(result), "ARMA Forecast")

  result_lastn <- fore_arma(x, phi = 0.8, n_ahead = 10, lastn = TRUE, plot = FALSE)
  expect_output(print(result_lastn), "Holdout validation")
  expect_output(print(result_lastn), "RMSE")
  expect_output(print(result_lastn), "MAD")
})


test_that("fore_arma works with AR(2)", {

  set.seed(999)
  x <- gen_arma(n = 100, phi = c(1.5, -0.75), plot = FALSE)

  old <- tswge::fore.arma.wge(x, phi = c(1.5, -0.75), n.ahead = 5, plot = FALSE)
  new <- fore_arma(x, phi = c(1.5, -0.75), n_ahead = 5, plot = FALSE)

  expect_equal(new$f, old$f, tolerance = 1e-6)
})
