# tests/testthat/test-period.R

library(tswge)

test_that("period matches original period.wge", {

  # Test case 1: Default (dB)
  set.seed(123)
  x1 <- rnorm(200)

  old1 <- tswge::period.wge(x1, dbcalc = "TRUE", plot = FALSE)
  new1 <- period(x1, db = TRUE, plot = FALSE)

  expect_equal(new1$freq, old1$freq, tolerance = 1e-10)
  expect_equal(new1$pgram, old1$pgram, tolerance = 1e-10)

  # Test case 2: Raw scale
  old2 <- tswge::period.wge(x1, dbcalc = "FALSE", plot = FALSE)
  new2 <- period(x1, db = FALSE, plot = FALSE)

  expect_equal(new2$freq, old2$freq, tolerance = 1e-10)
  expect_equal(new2$pgram, old2$pgram, tolerance = 1e-10)

  # Test case 3: AR process
  set.seed(456)
  x2 <- tswge::gen.arma.wge(n = 150, phi = c(1.5, -0.75), plot = FALSE, sn = 456)

  old3 <- tswge::period.wge(x2, plot = FALSE)
  new3 <- period(x2, plot = FALSE)

  expect_equal(new3$freq, old3$freq, tolerance = 1e-10)
  expect_equal(new3$pgram, old3$pgram, tolerance = 1e-10)
})


test_that("period validates inputs", {

  expect_error(period("not numeric"), "`x` must be")
  expect_error(period(numeric(0)), "`x` must be")
})


test_that("period returns correct structure", {

  result <- period(rnorm(100), plot = FALSE)

  expect_true(is.list(result))
  expect_named(result, c("freq", "pgram"))
  expect_equal(length(result$freq), 50)  # floor(100/2)
  expect_equal(length(result$pgram), 50)
})


test_that("period freq range is correct", {

  result <- period(rnorm(100), plot = FALSE)

  expect_true(all(result$freq > 0))
  expect_true(all(result$freq <= 0.5))
})


test_that("period plot works", {

  expect_no_error(period(rnorm(100), plot = TRUE))
  expect_no_error(period(rnorm(100), db = FALSE, plot = TRUE))
})
