library(tswge)

test_that("fore_arima returns correct structure", {

  set.seed(123)
  x <- gen_arma(n = 100, phi = 0.8, plot = FALSE)

  result <- fore_arima(x, phi = 0.8, n_ahead = 5, plot = FALSE)

  expect_s3_class(result, "fore_arima")
  expect_named(result, c("f", "ll", "ul", "resid", "wnv", "xbar", "se", "psi"))
})


test_that("fore_arima returns correct dimensions", {

  set.seed(123)
  x <- gen_arma(n = 100, phi = 0.8, plot = FALSE)

  result <- fore_arima(x, phi = 0.8, n_ahead = 10, plot = FALSE)

  expect_equal(length(result$f), 10)
  expect_equal(length(result$ll), 10)
  expect_equal(length(result$ul), 10)
  expect_equal(length(result$se), 10)
  expect_equal(length(result$resid), 100)
})


test_that("fore_arima forecasts match tswge for AR(1)", {

  set.seed(456)
  x <- gen_arma(n = 100, phi = 0.8, plot = FALSE)

  old <- tswge::fore.arima.wge(x, phi = 0.8, n.ahead = 5, plot = FALSE)
  new <- fore_arima(x, phi = 0.8, n_ahead = 5, plot = FALSE)

  expect_equal(new$f, old$f, tolerance = 1e-6)
  expect_equal(new$ll, old$ll, tolerance = 1e-6)
  expect_equal(new$ul, old$ul, tolerance = 1e-6)
})


test_that("fore_arima forecasts match tswge for AR(2)", {

  set.seed(789)
  x <- gen_arma(n = 100, phi = c(1.5, -0.75), plot = FALSE)

  old <- tswge::fore.arima.wge(x, phi = c(1.5, -0.75), n.ahead = 5, plot = FALSE)
  new <- fore_arima(x, phi = c(1.5, -0.75), n_ahead = 5, plot = FALSE)

  expect_equal(new$f, old$f, tolerance = 1e-6)
})


test_that("fore_arima forecasts match tswge for MA(1)", {

  set.seed(111)
  x <- gen_arma(n = 100, theta = 0.6, plot = FALSE)

  old <- tswge::fore.arima.wge(x, theta = 0.6, n.ahead = 5, plot = FALSE)
  new <- fore_arima(x, theta = 0.6, n_ahead = 5, plot = FALSE)

  expect_equal(new$f, old$f, tolerance = 1e-6)
})


test_that("fore_arima forecasts match tswge for ARMA(1,1)", {

  set.seed(222)
  x <- gen_arma(n = 100, phi = 0.7, theta = 0.4, plot = FALSE)

  old <- tswge::fore.arima.wge(x, phi = 0.7, theta = 0.4, n.ahead = 5, plot = FALSE)
  new <- fore_arima(x, phi = 0.7, theta = 0.4, n_ahead = 5, plot = FALSE)

  expect_equal(new$f, old$f, tolerance = 1e-6)
})


test_that("fore_arima handles differencing d=1", {

  set.seed(333)
  x <- cumsum(rnorm(100))  # Random walk

  old <- tswge::fore.arima.wge(x, phi = 0.5, d = 1, n.ahead = 5, plot = FALSE)
  new <- fore_arima(x, phi = 0.5, d = 1, n_ahead = 5, plot = FALSE)

  expect_equal(new$f, old$f, tolerance = 1e-5)
})


test_that("fore_arima handles differencing d=2", {

  set.seed(444)
  x <- cumsum(cumsum(rnorm(100)))  # Integrated twice

  old <- tswge::fore.arima.wge(x, d = 2, n.ahead = 5, plot = FALSE)
  new <- fore_arima(x, d = 2, n_ahead = 5, plot = FALSE)

  expect_equal(new$f, old$f, tolerance = 1e-4)
})


test_that("fore_arima handles seasonal s=12", {

  # Use a seasonal series
  x <- as.numeric(AirPassengers)

  old <- tswge::fore.arima.wge(x, phi = 0.5, s = 12, n.ahead = 6, plot = FALSE)
  new <- fore_arima(x, phi = 0.5, s = 12, n_ahead = 6, plot = FALSE)

  expect_equal(new$f, old$f, tolerance = 1e-4)
})


test_that("fore_arima lastn mode works", {

  set.seed(555)
  x <- gen_arma(n = 100, phi = 0.8, plot = FALSE)

  result <- fore_arima(x, phi = 0.8, n_ahead = 10, lastn = TRUE, plot = FALSE)

  # Forecasts should start from position 90, not 100
  expect_equal(length(result$f), 10)

  # Compare with tswge
  old <- tswge::fore.arima.wge(x, phi = 0.8, n.ahead = 10, lastn = TRUE, plot = FALSE)
  expect_equal(result$f, old$f, tolerance = 1e-6)
})


test_that("fore_arima prediction intervals are ordered correctly", {

  set.seed(666)
  x <- gen_arma(n = 100, phi = 0.8, plot = FALSE)

  result <- fore_arima(x, phi = 0.8, n_ahead = 5, plot = FALSE)

  # ll < f < ul for all horizons

  expect_true(all(result$ll < result$f))
  expect_true(all(result$f < result$ul))
})


test_that("fore_arima prediction intervals widen with horizon", {

  set.seed(777)
  x <- gen_arma(n = 100, phi = 0.8, plot = FALSE)

  result <- fore_arima(x, phi = 0.8, n_ahead = 10, plot = FALSE)

  # Standard errors should increase with horizon
  expect_true(all(diff(result$se) >= 0))
})


test_that("fore_arima validates inputs", {

  set.seed(888)
  x <- gen_arma(n = 100, phi = 0.8, plot = FALSE)

  expect_error(fore_arima(x, d = 4), "`d` must be")
  expect_error(fore_arima(x, d = -1), "`d` must be")
  expect_error(fore_arima(x, s = -1), "`s` must be")
  expect_error(fore_arima(x, n_ahead = 0), "`n_ahead` must be")
  expect_error(fore_arima(x, alpha = 0), "`alpha` must be")
  expect_error(fore_arima(x, alpha = 1), "`alpha` must be")
  expect_error(fore_arima("not numeric"), "`x` must be")
})


test_that("fore_arima plotting works", {

  set.seed(999)
  x <- gen_arma(n = 100, phi = 0.8, plot = FALSE)

  expect_no_error(fore_arima(x, phi = 0.8, n_ahead = 5, plot = TRUE))
  expect_no_error(fore_arima(x, phi = 0.8, n_ahead = 5, plot = TRUE, limits = FALSE))
})


test_that("fore_arima print method works", {

  set.seed(101)
  x <- gen_arma(n = 100, phi = 0.8, plot = FALSE)

  result <- fore_arima(x, phi = 0.8, n_ahead = 5, plot = FALSE)

  expect_output(print(result), "ARIMA Forecast")
  expect_output(print(result), "Horizon: 5")
})


test_that("fore_arima white noise variance matches tswge", {

  set.seed(202)
  x <- gen_arma(n = 100, phi = 0.8, plot = FALSE)

  old <- tswge::fore.arima.wge(x, phi = 0.8, n.ahead = 5, plot = FALSE)
  new <- fore_arima(x, phi = 0.8, n_ahead = 5, plot = FALSE)

  expect_equal(new$wnv, old$wnv, tolerance = 1e-6)
})


test_that("fore_arima standard errors match tswge", {

  set.seed(303)
  x <- gen_arma(n = 100, phi = 0.8, plot = FALSE)

  old <- tswge::fore.arima.wge(x, phi = 0.8, n.ahead = 5, plot = FALSE)
  new <- fore_arima(x, phi = 0.8, n_ahead = 5, plot = FALSE)

  expect_equal(new$se, old$se, tolerance = 1e-6)
})
