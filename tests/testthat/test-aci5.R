# tests/testthat/test-aic5.R

library(tswge)

test_that("aic5 produces same rankings as original", {

  set.seed(123)
  x <- tswge::gen.arma.wge(n = 200, phi = 0.7, theta = 0.3, plot = FALSE, sn = 123)

  # Capture output from both
  old_out <- capture.output(tswge::aic5.wge(x, p = 0:3, q = 0:2, type = "aic"))
  new_out <- capture.output(result <- aic5(x, p = 0:3, q = 0:2, type = "aic"))

  # Check that result is a sorted table
  expect_true(is.data.frame(result))
  expect_true(all(diff(result$aic) >= 0))  # Should be sorted ascending

  # Test BIC
  old_bic <- capture.output(tswge::aic5.wge(x, p = 0:3, q = 0:2, type = "bic"))
  new_bic <- capture.output(result_bic <- aic5(x, p = 0:3, q = 0:2, type = "bic"))

  expect_true(all(diff(result_bic$bic) >= 0))
})


test_that("aic5 displays correct output", {

  x <- rnorm(100)

  output <- capture.output(aic5(x, p = 0:2, q = 0:1, type = "aic"))

  expect_true(any(grepl("WORKING", output)))
  expect_true(any(grepl("Five Smallest", output)))
})


test_that("aic5 returns invisible sorted table", {

  x <- rnorm(100)

  capture.output({
    result <- aic5(x, p = 0:2, q = 0:1)
  })

  expect_true(is.data.frame(result))
  expect_true(all(c("p", "q", "aic") %in% names(result)))
})


test_that("aic5 works with parallel", {

  skip_on_cran()
  skip_if(parallel::detectCores(logical = FALSE) < 2, "Not enough cores")

  x <- rnorm(100)

  capture.output({
    result_seq <- aic5(x, p = 0:2, q = 0:1, cores = 1)
    result_par <- aic5(x, p = 0:2, q = 0:1, cores = 2)
  })

  expect_equal(result_seq, result_par)
})
