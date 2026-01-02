library(tswge)

test_that("true_acf matches original true.arma.aut.wge", {

  # Test case 1: AR(1)
  old1 <- tswge::true.arma.aut.wge(phi = 0.8, lag.max = 10, plot = FALSE)
  new1 <- true_acf(phi = 0.8, lag_max = 10, plot = FALSE)

  expect_equal(new1$acf, old1$acf, tolerance = 1e-10)
  expect_equal(new1$acv, old1$acv, tolerance = 1e-6)

  # Test case 2: AR(2)
  old2 <- tswge::true.arma.aut.wge(phi = c(1.5, -0.75), lag.max = 20, plot = FALSE)
  new2 <- true_acf(phi = c(1.5, -0.75), lag_max = 20, plot = FALSE)

  expect_equal(new2$acf, old2$acf, tolerance = 1e-10)
  expect_equal(new2$acv, old2$acv, tolerance = 1e-6)

  # Test case 3: MA(1)
  old3 <- tswge::true.arma.aut.wge(theta = 0.6, lag.max = 15, plot = FALSE)
  new3 <- true_acf(theta = 0.6, lag_max = 15, plot = FALSE)

  expect_equal(new3$acf, old3$acf, tolerance = 1e-10)
  expect_equal(new3$acv, old3$acv, tolerance = 1e-6)

  # Test case 4: ARMA(1,1)
  old4 <- tswge::true.arma.aut.wge(phi = 0.7, theta = 0.4, lag.max = 15, plot = FALSE)
  new4 <- true_acf(phi = 0.7, theta = 0.4, lag_max = 15, plot = FALSE)

  expect_equal(new4$acf, old4$acf, tolerance = 1e-10)
  expect_equal(new4$acv, old4$acv, tolerance = 1e-6)

  # Test case 5: ARMA(2,1)
  old5 <- tswge::true.arma.aut.wge(phi = c(1.2, -0.6), theta = 0.3, lag.max = 20, plot = FALSE)
  new5 <- true_acf(phi = c(1.2, -0.6), theta = 0.3, lag_max = 20, plot = FALSE)

  expect_equal(new5$acf, old5$acf, tolerance = 1e-10)
  expect_equal(new5$acv, old5$acv, tolerance = 1e-6)
})


test_that("true_acf computes correct AR(1) values", {

  # AR(1): rho(k) = phi^k
  phi <- 0.8
  result <- true_acf(phi = phi, lag_max = 5, plot = FALSE)

  expected_acf <- phi^(0:5)
  expect_equal(result$acf, expected_acf, tolerance = 1e-10)
})


test_that("true_acf handles white noise", {

  result <- true_acf(lag_max = 10, vara = 2, plot = FALSE)

  expect_equal(result$acf[1], 1)
  expect_equal(result$acf[-1], rep(0, 10))
  expect_equal(result$acv[1], 2)
})


test_that("true_acf validates inputs", {

  expect_error(true_acf(phi = "bad"), "`phi` must be")
  expect_error(true_acf(theta = "bad"), "`theta` must be")
  expect_error(true_acf(lag_max = -1), "`lag_max` must be")
  expect_error(true_acf(vara = -1), "`vara` must be")
})


test_that("true_acf returns correct length", {

  result <- true_acf(phi = 0.5, lag_max = 15, plot = FALSE)

  expect_equal(length(result$acf), 16)  # lags 0 to 15
  expect_equal(length(result$acv), 16)
})


test_that("true_acf plot works", {

  expect_no_error(true_acf(phi = 0.8, lag_max = 10, plot = TRUE))
})
