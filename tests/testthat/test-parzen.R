# tests/testthat/test-parzen.R

library(tswge)

test_that("parzen matches original parzen.wge", {

  # Test case 1: Default parameters
  set.seed(123)
  x1 <- rnorm(200)

  old1 <- tswge::parzen.wge(x1, plot = FALSE)
  new1 <- parzen(x1, plot = FALSE)

  expect_equal(new1$freq, old1$freq, tolerance = 1e-10)
  expect_equal(new1$pzgram, old1$pzgram, tolerance = 1e-10)

  # Test case 2: Custom truncation
  set.seed(456)
  x2 <- cumsum(rnorm(150))

  old2 <- tswge::parzen.wge(x2, trunc = 20, plot = FALSE)
  new2 <- parzen(x2, trunc = 20, plot = FALSE)

  expect_equal(new2$freq, old2$freq, tolerance = 1e-10)
  expect_equal(new2$pzgram, old2$pzgram, tolerance = 1e-10)

  # Test case 3: Non-dB scale
  set.seed(789)
  x3 <- arima.sim(model = list(ar = 0.7), n = 100)

  old3 <- tswge::parzen.wge(x3, dbcalc = "FALSE", plot = FALSE)
  new3 <- parzen(x3, db = FALSE, plot = FALSE)

  expect_equal(new3$freq, old3$freq, tolerance = 1e-10)
  expect_equal(new3$pzgram, old3$pzgram, tolerance = 1e-10)
})


test_that("parzen validates inputs", {

  expect_error(parzen("not numeric"), "`x` must be")
  expect_error(parzen(numeric(0)), "`x` must be")
  expect_error(parzen(rnorm(100), trunc = -5), "`trunc` must be")
})


test_that("parzen returns M in output", {

  x <- rnorm(100)

  # Default M
  result1 <- parzen(x, plot = FALSE)
  expect_equal(result1$M, floor(2 * sqrt(100)))

  # Custom M
  result2 <- parzen(x, trunc = 15, plot = FALSE)
  expect_equal(result2$M, 15)
})
