# tests/testthat/test-est_ar.R

library(tswge)

test_that("est_ar matches original est.ar.wge", {

  capture.output({

    # Test case 1: MLE method
    set.seed(123)
    x1 <- tswge::gen.arma.wge(n = 200, phi = c(1.5, -0.75), plot = FALSE, sn = 123)

    old1 <- tswge::est.ar.wge(x1, p = 2, factor = FALSE, method = "mle")
    new1 <- est_ar(x1, p = 2, factor = FALSE, method = "mle")

    expect_equal(new1$phi, old1$phi, tolerance = 1e-3)
    expect_equal(new1$avar, old1$avar, tolerance = 1e-4)
    expect_equal(new1$xbar, old1$xbar, tolerance = 1e-10)
    expect_equal(new1$method, old1$method)

    # Test case 2: Burg method
    old2 <- tswge::est.ar.wge(x1, p = 2, factor = FALSE, method = "burg")
    new2 <- est_ar(x1, p = 2, factor = FALSE, method = "burg")

    expect_equal(new2$phi, old2$phi, tolerance = 1e-10)
    expect_equal(new2$avar, old2$avar, tolerance = 1e-4)
    expect_equal(new2$method, old2$method)

    # Test case 3: Yule-Walker method
    old3 <- tswge::est.ar.wge(x1, p = 2, factor = FALSE, method = "yw")
    new3 <- est_ar(x1, p = 2, factor = FALSE, method = "yw")

    expect_equal(new3$phi, old3$phi, tolerance = 1e-10)
    expect_equal(new3$avar, old3$avar, tolerance = 1e-4)
    expect_equal(new3$method, old3$method)
  })
})


test_that("est_ar validates inputs", {

  x <- rnorm(100)

  expect_error(est_ar("not numeric", p = 2), "`x` must be")
  expect_error(est_ar(numeric(0), p = 2), "`x` must be")
  expect_error(est_ar(x, p = 0), "`p` must be")
  expect_error(est_ar(x, p = -1), "`p` must be")
})


test_that("est_ar has correct structure", {

  capture.output({
    result <- est_ar(rnorm(100), p = 2, factor = FALSE)
  })

  expect_s3_class(result, "est_ar")
  expect_named(result, c("method", "phi", "res", "avar", "xbar",
                         "aic", "aicc", "bic"))
})


test_that("est_ar print method works", {

  capture.output({
    result <- est_ar(rnorm(100), p = 2, factor = FALSE)
  })

  expect_output(print(result), "AR Model Estimation")
})


test_that("est_ar method argument works", {

  x <- rnorm(100)

  capture.output({
    r1 <- est_ar(x, p = 2, factor = FALSE, method = "mle")
    r2 <- est_ar(x, p = 2, factor = FALSE, method = "burg")
    r3 <- est_ar(x, p = 2, factor = FALSE, method = "yw")
  })

  expect_equal(r1$method, "mle")
  expect_equal(r2$method, "burg")
  expect_equal(r3$method, "yw")

  # Different methods should give different (but similar) results
  expect_false(identical(r1$phi, r2$phi))
  expect_false(identical(r2$phi, r3$phi))
})
