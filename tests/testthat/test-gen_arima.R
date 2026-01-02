# tests/testthat/test-gen_arima.R

library(tswge)

test_that("gen_arima matches original gen.arima.wge with seed", {

  # Test case 1: AR(2)
  old1 <- tswge::gen.arima.wge(n = 100, phi = c(1.6, -0.9),
                               plot = FALSE, sn = 123)
  new1 <- gen_arima(n = 100, phi = c(1.6, -0.9),
                    plot = FALSE, seed = 123)

  expect_equal(new1, as.double(old1), tolerance = 1e-10)

  # Test case 2: MA(1)
  old2 <- tswge::gen.arima.wge(n = 100, theta = 0.6,
                               plot = FALSE, sn = 456)
  new2 <- gen_arima(n = 100, theta = 0.6,
                    plot = FALSE, seed = 456)

  expect_equal(new2, as.double(old2), tolerance = 1e-10)

  # Test case 3: ARMA(1,1)
  old3 <- tswge::gen.arima.wge(n = 100, phi = 0.7, theta = 0.3,
                               plot = FALSE, sn = 789)
  new3 <- gen_arima(n = 100, phi = 0.7, theta = 0.3,
                    plot = FALSE, seed = 789)

  expect_equal(new3, as.double(old3), tolerance = 1e-10)

  # Test case 4: ARIMA(1,1,0)
  old4 <- tswge::gen.arima.wge(n = 100, phi = 0.5, d = 1,
                               plot = FALSE, sn = 111)
  new4 <- gen_arima(n = 100, phi = 0.5, d = 1,
                    plot = FALSE, seed = 111)

  expect_equal(new4, as.double(old4), tolerance = 1e-10)

  # Test case 5: With mean
  old5 <- tswge::gen.arima.wge(n = 100, phi = 0.5, mu = 10,
                               plot = FALSE, sn = 222)
  new5 <- gen_arima(n = 100, phi = 0.5, mu = 10,
                    plot = FALSE, seed = 222)

  expect_equal(new5, as.double(old5), tolerance = 1e-10)
})


test_that("gen_arima handles seasonal", {

  old_s <- tswge::gen.arima.wge(n = 100, phi = 0.5, s = 12,
                                plot = FALSE, sn = 333)
  new_s <- gen_arima(n = 100, phi = 0.5, s = 12,
                     plot = FALSE, seed = 333)

  expect_equal(new_s, as.double(old_s), tolerance = 1e-10)
})


test_that("gen_arima returns correct length", {

  x <- gen_arima(n = 150, phi = 0.5, plot = FALSE, seed = 1)
  expect_equal(length(x), 150)
})


test_that("gen_arima validates inputs", {

  expect_error(gen_arima(n = -1, phi = 0.5), "`n` must be")
  expect_error(gen_arima(n = 100, vara = -1), "`vara` must be")
  expect_error(gen_arima(n = 100, d = -1), "`d` must be")
})


test_that("gen_arima is reproducible with seed", {

  x1 <- gen_arima(n = 50, phi = 0.7, theta = 0.3, plot = FALSE, seed = 999)
  x2 <- gen_arima(n = 50, phi = 0.7, theta = 0.3, plot = FALSE, seed = 999)

  expect_identical(x1, x2)
})
