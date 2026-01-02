library(tswge)

test_that("fore_aruma returns correct structure", {

  set.seed(123)
  x <- gen_arma(n = 100, phi = 0.8, plot = FALSE)

  result <- fore_aruma(x, phi = 0.8, n_ahead = 5, plot = FALSE)

  expect_s3_class(result, "fore_aruma")
  expect_named(result, c("f", "ll", "ul", "resid", "wnv", "xbar", "se", "psi"))
})


test_that("fore_aruma with lambda=0 matches fore_arima", {

  set.seed(456)
  x <- gen_arma(n = 100, phi = 0.8, plot = FALSE)

  aruma_result <- fore_aruma(x, phi = 0.8, d = 1, lambda = 0, n_ahead = 5, plot = FALSE)
  arima_result <- fore_arima(x, phi = 0.8, d = 1, n_ahead = 5, plot = FALSE)

  expect_equal(aruma_result$f, arima_result$f)
  expect_equal(aruma_result$ll, arima_result$ll)
  expect_equal(aruma_result$ul, arima_result$ul)
  expect_equal(aruma_result$se, arima_result$se)
})


test_that("fore_aruma matches tswge fore.aruma.wge", {

  set.seed(789)
  x <- cumsum(rnorm(100))  # Random walk

  # Test with lambda factor
  old <- tswge::fore.aruma.wge(x, phi = 0.5, d = 1, lambda = c(0.3), n.ahead = 5, plot = FALSE)
  new <- fore_aruma(x, phi = 0.5, d = 1, lambda = c(0.3), n_ahead = 5, plot = FALSE)

  expect_equal(new$f, old$f, tolerance = 1e-5)
  # Prediction limits may differ slightly due to variance estimation
  expect_equal(new$ll, old$ll, tolerance = 1e-2)
  expect_equal(new$ul, old$ul, tolerance = 1e-2)
})


test_that("fore_aruma handles seasonal with lambda", {

  x <- as.numeric(AirPassengers)

  old <- tswge::fore.aruma.wge(x, phi = 0.5, s = 12, lambda = c(0.2), n.ahead = 6, plot = FALSE)
  new <- fore_aruma(x, phi = 0.5, s = 12, lambda = c(0.2), n_ahead = 6, plot = FALSE)

  expect_equal(new$f, old$f, tolerance = 1e-4)
})


test_that("fore_aruma validates lambda input", {

  set.seed(111)
  x <- gen_arma(n = 100, phi = 0.8, plot = FALSE)

  expect_error(fore_aruma(x, lambda = "bad"), "`lambda` must be numeric")
})


test_that("fore_aruma plotting works", {

  set.seed(222)
  x <- gen_arma(n = 100, phi = 0.8, plot = FALSE)

  expect_no_error(fore_aruma(x, phi = 0.8, lambda = c(0.3), n_ahead = 5, plot = TRUE))
})


test_that("fore_aruma print method works", {

  set.seed(333)
  x <- gen_arma(n = 100, phi = 0.8, plot = FALSE)

  result <- fore_aruma(x, phi = 0.8, lambda = c(0.2), n_ahead = 5, plot = FALSE)

  expect_output(print(result), "ARUMA Forecast")
  expect_output(print(result), "Horizon: 5")
})


test_that("fore_aruma with d and s and lambda", {

  set.seed(444)
  x <- cumsum(cumsum(rnorm(100)))

  # All factors combined
  result <- fore_aruma(x, phi = 0.5, d = 1, s = 4, lambda = c(0.2, -0.1),
                       n_ahead = 5, plot = FALSE)

  expect_equal(length(result$f), 5)
  expect_equal(length(result$se), 5)
  expect_true(all(result$ll < result$f))
  expect_true(all(result$f < result$ul))
})


test_that("fore_aruma lastn mode works", {

  set.seed(555)
  x <- gen_arma(n = 100, phi = 0.8, plot = FALSE)

  result <- fore_aruma(x, phi = 0.8, lambda = c(0.3), n_ahead = 10,
                       lastn = TRUE, plot = FALSE)

  expect_equal(length(result$f), 10)
})


test_that("fore_aruma matches tswge with complex lambda", {

  set.seed(666)
  x <- cumsum(rnorm(100))

  # Multi-coefficient lambda
  old <- tswge::fore.aruma.wge(x, phi = c(0.5, -0.3), d = 1, lambda = c(0.4, -0.2),
                               n.ahead = 5, plot = FALSE)
  new <- fore_aruma(x, phi = c(0.5, -0.3), d = 1, lambda = c(0.4, -0.2),
                    n_ahead = 5, plot = FALSE)

  expect_equal(new$f, old$f, tolerance = 1e-5)
})
