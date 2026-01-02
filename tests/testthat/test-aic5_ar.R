# tests/testthat/test-aic5_ar.R

library(tswge)

test_that("aic5_ar produces same rankings as original", {

  set.seed(123)
  x <- tswge::gen.arma.wge(n = 200, phi = c(1.5, -0.75), plot = FALSE, sn = 123)

  # Capture output from both
  old_out <- capture.output(tswge::aic5.ar.wge(x, p = 0:5, type = "aic", method = "mle"))
  new_out <- capture.output(result <- aic5_ar(x, p = 0:5, type = "aic", method = "mle"))

  # Check that result is a sorted table
  expect_true(is.data.frame(result))
  expect_true(all(diff(result$aic) >= 0))  # Should be sorted ascending

  # Test BIC
  old_bic <- capture.output(tswge::aic5.ar.wge(x, p = 0:5, type = "bic", method = "mle"))
  new_bic <- capture.output(result_bic <- aic5_ar(x, p = 0:5, type = "bic", method = "mle"))

  expect_true(all(diff(result_bic$bic) >= 0))

  # Test Burg method
  old_burg <- capture.output(tswge::aic5.ar.wge(x, p = 0:5, type = "aic", method = "burg"))
  new_burg <- capture.output(result_burg <- aic5_ar(x, p = 0:5, type = "aic", method = "burg"))

  expect_true(all(diff(result_burg$aic) >= 0))
})


test_that("aic5_ar displays correct output", {

  x <- rnorm(100)

  output <- capture.output(aic5_ar(x, p = 1:5, type = "aic"))

  expect_true(any(grepl("WORKING", output)))
  expect_true(any(grepl("Five Smallest", output)))
  expect_true(any(grepl("Method=", output)))
})


test_that("aic5_ar returns invisible sorted table", {

  x <- rnorm(100)

  capture.output({
    result <- aic5_ar(x, p = 1:5)
  })

  expect_true(is.data.frame(result))
  expect_true("p" %in% names(result))
  expect_true("aic" %in% names(result))
})
