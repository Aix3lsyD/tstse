# tests/testthat/test-mult.R

library(tswge)

test_that("mult matches original mult.wge", {

  # Test case 1: Two simple AR(1) factors
  old1 <- tswge::mult.wge(fac1 = c(0.8), fac2 = c(0.5))
  new1 <- mult(c(0.8), c(0.5))

  expect_equal(new1$model_coef, as.double(old1$model.coef), tolerance = 1e-10)

  # Test case 2: AR(2) factor
  old2 <- tswge::mult.wge(fac1 = c(1.6, -0.9))
  new2 <- mult(c(1.6, -0.9))

  expect_equal(new2$model_coef, as.double(old2$model.coef), tolerance = 1e-10)

  # Test case 3: Three factors
  old3 <- tswge::mult.wge(fac1 = c(0.7), fac2 = c(0.5), fac3 = c(1.2, -0.6))
  new3 <- mult(c(0.7), c(0.5), c(1.2, -0.6))

  expect_equal(new3$model_coef, as.double(old3$model.coef), tolerance = 1e-10)

  # Test case 4: Differencing factor (1 - B)
  old4 <- tswge::mult.wge(fac1 = c(1))
  new4 <- mult(c(1))

  expect_equal(new4$model_coef, as.double(old4$model.coef), tolerance = 1e-10)

  # Test case 5: All zeros returns 0
  old5 <- tswge::mult.wge()
  new5 <- mult()

  expect_equal(new5$model_coef, 0)
})


test_that("mult handles more than 6 factors", {

  result <- mult(c(0.5), c(0.5), c(0.5), c(0.5), c(0.5), c(0.5), c(0.5))
  expect_equal(length(result$model_coef), 7)
})


test_that("mult validates inputs", {

  expect_error(mult("not numeric"), "must be numeric")
})


test_that("mult handles seasonal factors", {

  seas4 <- c(0, 0, 0, 1)
  result <- mult(seas4)

  expect_equal(result$model_coef, c(0, 0, 0, 1))
})


test_that("mult returns char_poly when PolynomF available", {

  skip_if_not_installed("PolynomF")

  result <- mult(c(0.8), c(0.5))
  expect_true(!is.null(result$char_poly))
  expect_s3_class(result$char_poly, "polynom")
})


test_that("mult works without PolynomF", {

  # Force convolve path by calling internal function directly
  polys <- list(c(1, -0.8), c(1, -0.5))
  result <- tstse:::.mult_convolve(polys)

  expect_equal(result$model_coef, c(1.3, -0.4), tolerance = 1e-10)
  expect_null(result$char_poly)
})
