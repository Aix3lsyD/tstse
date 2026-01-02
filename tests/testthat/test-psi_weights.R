# tests/testthat/test-psi_weights.R

library(tswge)

test_that("psi_weights matches original psi.weights.wge", {

  # Test case 1: AR(1)
  old1 <- tswge::psi.weights.wge(phi = 0.8, lag.max = 10)
  new1 <- psi_weights(phi = 0.8, lag_max = 10)

  expect_equal(new1, old1, tolerance = 1e-10)

  # Test case 2: MA(1)
  old2 <- tswge::psi.weights.wge(theta = 0.6, lag.max = 10)
  new2 <- psi_weights(theta = 0.6, lag_max = 10)

  expect_equal(new2, old2, tolerance = 1e-10)

  # Test case 3: ARMA(1,1)
  old3 <- tswge::psi.weights.wge(phi = 0.7, theta = 0.4, lag.max = 15)
  new3 <- psi_weights(phi = 0.7, theta = 0.4, lag_max = 15)

  expect_equal(new3, old3, tolerance = 1e-10)

  # Test case 4: AR(2)
  old4 <- tswge::psi.weights.wge(phi = c(1.5, -0.75), lag.max = 20)
  new4 <- psi_weights(phi = c(1.5, -0.75), lag_max = 20)

  expect_equal(new4, old4, tolerance = 1e-10)

  # Test case 5: ARMA(2,1)
  old5 <- tswge::psi.weights.wge(phi = c(1.2, -0.6), theta = 0.3, lag.max = 10)
  new5 <- psi_weights(phi = c(1.2, -0.6), theta = 0.3, lag_max = 10)

  expect_equal(new5, old5, tolerance = 1e-10)
})


test_that("psi_weights validates inputs", {

  expect_error(psi_weights(phi = "bad"), "`phi` must be")
  expect_error(psi_weights(theta = "bad"), "`theta` must be")
  expect_error(psi_weights(lag_max = 0), "`lag_max` must be")
  expect_error(psi_weights(lag_max = -5), "`lag_max` must be")
})


test_that("psi_weights returns correct length", {

  result <- psi_weights(phi = 0.5, lag_max = 12)
  expect_equal(length(result), 12)
})
