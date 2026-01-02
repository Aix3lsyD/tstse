# tests/testthat/test-gen_arma.R

library(tswge)

test_that("gen_arma matches original gen.arma.wge with seed", {

  # Test case 1: AR(1)
  old1 <- tswge::gen.arma.wge(n = 100, phi = 0.7, plot = FALSE, sn = 123)
  new1 <- gen_arma(n = 100, phi = 0.7, plot = FALSE, seed = 123)

  expect_equal(new1, as.double(old1), tolerance = 1e-10)

  # Test case 2: MA(1)
  old2 <- tswge::gen.arma.wge(n = 100, theta = 0.6, plot = FALSE, sn = 456)
  new2 <- gen_arma(n = 100, theta = 0.6, plot = FALSE, seed = 456)

  expect_equal(new2, as.double(old2), tolerance = 1e-10)

  # Test case 3: ARMA(2,1)
  old3 <- tswge::gen.arma.wge(n = 100, phi = c(1.6, -0.9), theta = 0.3,
                              plot = FALSE, sn = 789)
  new3 <- gen_arma(n = 100, phi = c(1.6, -0.9), theta = 0.3,
                   plot = FALSE, seed = 789)

  expect_equal(new3, as.double(old3), tolerance = 1e-10)

  # Test case 4: With mean and variance
  old4 <- tswge::gen.arma.wge(n = 100, phi = 0.5, mu = 10, vara = 2,
                              plot = FALSE, sn = 111)
  new4 <- gen_arma(n = 100, phi = 0.5, mu = 10, vara = 2,
                   plot = FALSE, seed = 111)

  expect_equal(new4, as.double(old4), tolerance = 1e-10)

  # Test case 5: Pure white noise
  old5 <- tswge::gen.arma.wge(n = 100, plot = FALSE, sn = 222)
  new5 <- gen_arma(n = 100, plot = FALSE, seed = 222)

  expect_equal(new5, as.double(old5), tolerance = 1e-10)
})


test_that("gen_arma returns correct length", {

  x <- gen_arma(n = 150, phi = 0.5, plot = FALSE, seed = 1)
  expect_equal(length(x), 150)
})


test_that("gen_arma is reproducible with seed", {

  x1 <- gen_arma(n = 50, phi = 0.7, theta = 0.3, plot = FALSE, seed = 999)
  x2 <- gen_arma(n = 50, phi = 0.7, theta = 0.3, plot = FALSE, seed = 999)

  expect_identical(x1, x2)
})
