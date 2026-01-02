# tests/testthat/test-factor.R

library(tswge)

test_that("factor matches original factor.wge output structure", {

  # Test case 1: AR(2) with complex roots
  phi1 <- c(1.6, -0.9)

  # Capture printed output from both
  old_out <- capture.output(tswge::factor.wge(phi = phi1))
  new_out <- capture.output(result <- factor(phi = phi1))

  # New version should return data
  expect_true(is.data.frame(result$ar))
  expect_equal(nrow(result$ar), 1)  # one conjugate pair = one row
  expect_null(result$ma)

  # Test case 2: AR(1) with real root
  phi2 <- 0.7

  result2 <- factor(phi = phi2, print = FALSE)
  expect_true(is.data.frame(result2$ar))
  expect_equal(result2$ar$is_quadratic, FALSE)

  # Test case 3: ARMA
  phi3 <- c(0.5, -0.3)
  theta3 <- 0.4

  result3 <- factor(phi = phi3, theta = theta3, print = FALSE)
  expect_true(!is.null(result3$ar))
  expect_true(!is.null(result3$ma))
})


test_that("factor computes correct values", {

  # AR(1) with phi = 0.8
  # Root should be 1/0.8 = 1.25, abs_recip = 0.8
  result <- factor(phi = 0.8, print = FALSE)

  expect_equal(result$ar$abs_recip, 0.8, tolerance = 1e-10)
  expect_equal(result$ar$frequency, 0, tolerance = 1e-10)  # real root = freq 0
})


test_that("factor validates inputs", {

  expect_error(factor(phi = "bad"), "`phi` must be")
  expect_error(factor(theta = "bad"), "`theta` must be")
})


test_that("factor returns invisibly when printing", {

  expect_invisible(factor(phi = 0.7))
})
