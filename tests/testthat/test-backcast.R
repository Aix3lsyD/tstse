# tests/testthat/test-backcast.R

library(tswge)

test_that("backcast matches original backcast.wge", {

  # Test case 1: AR only
  set.seed(123)
  x1 <- cumsum(rnorm(100))
  phi1 <- 0.7
  theta1 <- 0

  old1 <- tswge::backcast.wge(x1, phi = phi1, theta = theta1, n.back = 50)
  new1 <- backcast(x1, phi = phi1, theta = theta1, n_back = 50)

  expect_equal(new1, old1, tolerance = 1e-10)

  # Test case 2: MA only
  x2 <- rnorm(100)
  phi2 <- 0
  theta2 <- c(0.5, -0.3)

  old2 <- tswge::backcast.wge(x2, phi = phi2, theta = theta2, n.back = 50)
  new2 <- backcast(x2, phi = phi2, theta = theta2, n_back = 50)

  expect_equal(new2, old2, tolerance = 1e-10)

  # Test case 3: ARMA
  set.seed(456)
  x3 <- arima.sim(model = list(ar = 0.6, ma = 0.4), n = 150)
  phi3 <- 0.6
  theta3 <- 0.4

  old3 <- tswge::backcast.wge(x3, phi = phi3, theta = theta3, n.back = 30)
  new3 <- backcast(x3, phi = phi3, theta = theta3, n_back = 30)

  expect_equal(new3, as.double(old3), tolerance = 1e-10)
})


test_that("backcast validates inputs", {

  expect_error(backcast("not numeric", phi = 0.5), "`x` must be")
  expect_error(backcast(numeric(0), phi = 0.5), "`x` must be")
  expect_error(backcast(1:10, phi = "bad"), "`phi` must be")
  expect_error(backcast(1:10, n_back = -1), "`n_back` must be")
})
