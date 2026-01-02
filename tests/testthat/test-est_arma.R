# tests/testthat/test-est_arma.R

library(tswge)

test_that("est_arma matches original est.arma.wge", {

  # Suppress factor table printing from both
  capture.output({

    # Test case 1: AR(2)
    set.seed(123)
    x1 <- tswge::gen.arma.wge(n = 200, phi = c(1.5, -0.75), plot = FALSE, sn = 123)

    old1 <- tswge::est.arma.wge(x1, p = 2, q = 0, factor = FALSE)
    new1 <- est_arma(x1, p = 2, q = 0, factor = FALSE)

    expect_equal(new1$phi, old1$phi, tolerance = 1e-4)
    expect_equal(new1$theta, old1$theta, tolerance = 1e-4)
    expect_equal(new1$avar, old1$avar, tolerance = 1e-5)
    expect_equal(new1$xbar, old1$xbar, tolerance = 1e-10)
    expect_equal(new1$aic, old1$aic, tolerance = 1e-5)
    expect_equal(new1$bic, old1$bic, tolerance = 1e-5)

    # Test case 2: MA(1)
    set.seed(456)
    x2 <- tswge::gen.arma.wge(n = 200, theta = 0.6, plot = FALSE, sn = 456)

    old2 <- tswge::est.arma.wge(x2, p = 0, q = 1, factor = FALSE)
    new2 <- est_arma(x2, p = 0, q = 1, factor = FALSE)

    expect_equal(new2$phi, old2$phi, tolerance = 1e-4)
    expect_equal(new2$theta, old2$theta, tolerance = 1e-4)
    expect_equal(new2$avar, old2$avar, tolerance = 1e-5)

    # Test case 3: ARMA(1,1)
    set.seed(789)
    x3 <- tswge::gen.arma.wge(n = 200, phi = 0.7, theta = 0.3, plot = FALSE, sn = 789)

    old3 <- tswge::est.arma.wge(x3, p = 1, q = 1, factor = FALSE)
    new3 <- est_arma(x3, p = 1, q = 1, factor = FALSE)

    expect_equal(new3$phi, old3$phi, tolerance = 1e-4)
    expect_equal(new3$theta, old3$theta, tolerance = 1e-2)
    expect_equal(new3$avar, old3$avar, tolerance = 1e-5)
  })
})


test_that("est_arma validates inputs", {

  x <- rnorm(100)

  expect_error(est_arma("not numeric", p = 1), "`x` must be")
  expect_error(est_arma(numeric(0), p = 1), "`x` must be")
  expect_error(est_arma(x, p = -1), "`p` must be")
  expect_error(est_arma(x, q = -1), "`q` must be")
  expect_error(est_arma(x, p = 0, q = 0), "At least one")
})


test_that("est_arma has correct structure", {

  capture.output({
    result <- est_arma(rnorm(100), p = 1, factor = FALSE)
  })

  expect_s3_class(result, "est_arma")
  expect_named(result, c("phi", "theta", "res", "avar", "xbar",
                         "aic", "aicc", "bic", "se_phi", "se_theta"))
})


test_that("est_arma print method works", {

  capture.output({
    result <- est_arma(rnorm(100), p = 1, q = 1, factor = FALSE)
  })

  expect_output(print(result), "ARMA Model Estimation")
})
