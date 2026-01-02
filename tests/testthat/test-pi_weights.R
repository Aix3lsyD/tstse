# tests/testthat/test-pi_weights.R

library(tswge)

test_that("pi_weights matches original pi.weights.wge", {

  skip_if_not_installed("astsa")

  # Test case 1: MA(1)
  old1 <- tswge::pi.weights.wge(theta = 0.8, lag.max = 10)
  new1 <- pi_weights(theta = 0.8, lag_max = 10)

  expect_equal(new1, old1, tolerance = 1e-10)

  # Test case 2: AR(1)
  old2 <- tswge::pi.weights.wge(phi = 0.6, lag.max = 10)
  new2 <- pi_weights(phi = 0.6, lag_max = 10)

  expect_equal(new2, old2, tolerance = 1e-10)

  # Test case 3: ARMA(1,1)
  old3 <- tswge::pi.weights.wge(phi = 0.7, theta = 0.4, lag.max = 15)
  new3 <- pi_weights(phi = 0.7, theta = 0.4, lag_max = 15)

  expect_equal(new3, old3, tolerance = 1e-10)

  # Test case 4: MA(2)
  old4 <- tswge::pi.weights.wge(theta = c(0.5, -0.3), lag.max = 20)
  new4 <- pi_weights(theta = c(0.5, -0.3), lag_max = 20)

  expect_equal(new4, old4, tolerance = 1e-10)

  # Test case 5: ARMA(2,1)
  old5 <- tswge::pi.weights.wge(phi = c(1.2, -0.6), theta = 0.3, lag.max = 10)
  new5 <- pi_weights(phi = c(1.2, -0.6), theta = 0.3, lag_max = 10)

  expect_equal(new5, old5, tolerance = 1e-10)
})


test_that("pi_weights works without astsa", {

  # Test internal implementation directly
  result <- tstse:::pi_weights_internal(phi = 0.7, theta = 0.4, lag_max = 10)

  expect_equal(length(result), 10)
  expect_true(is.numeric(result))
})


test_that("pi_weights validates inputs", {

  expect_error(pi_weights(phi = "bad"), "`phi` must be")
  expect_error(pi_weights(theta = "bad"), "`theta` must be")
  expect_error(pi_weights(lag_max = 0), "`lag_max` must be")
  expect_error(pi_weights(lag_max = -5), "`lag_max` must be")
})


test_that("pi_weights returns correct length", {

  result <- pi_weights(theta = 0.5, lag_max = 12)
  expect_equal(length(result), 12)
})


test_that("pi_weights internal matches astsa when available", {

  skip_if_not_installed("astsa")

  phi <- c(0.5, -0.3)
  theta <- 0.6
  lag_max <- 15

  with_astsa <- pi_weights(phi = phi, theta = theta, lag_max = lag_max)
  without_astsa <- tstse:::pi_weights_internal(phi = phi, theta = theta, lag_max = lag_max)

  expect_equal(without_astsa, with_astsa, tolerance = 1e-10)
})
