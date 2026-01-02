# tests/testthat/test-aic_ar.R

library(tswge)

test_that("aic_ar matches original aic.ar.wge", {

  capture.output({

    # Test case 1: AIC selection
    set.seed(123)
    x1 <- tswge::gen.arma.wge(n = 200, phi = c(1.5, -0.75), plot = FALSE, sn = 123)

    old1 <- tswge::aic.ar.wge(x1, p = 1:5, type = "aic", method = "mle")
    new1 <- aic_ar(x1, p = 1:5, type = "aic", method = "mle")

    expect_equal(new1$p, old1$p)
    expect_equal(new1$phi, old1$phi, tolerance = 1e-3)
    expect_equal(new1$value, old1$value, tolerance = 1e-4)
    expect_equal(new1$vara, old1$vara, tolerance = 1e-4)

    # Test case 2: BIC selection
    old2 <- tswge::aic.ar.wge(x1, p = 1:5, type = "bic", method = "mle")
    new2 <- aic_ar(x1, p = 1:5, type = "bic", method = "mle")

    expect_equal(new2$p, old2$p)
    expect_equal(new2$phi, old2$phi, tolerance = 1e-3)
    expect_equal(new2$value, old2$value, tolerance = 1e-4)

    # Test case 3: AICC selection
    old3 <- tswge::aic.ar.wge(x1, p = 1:5, type = "aicc", method = "mle")
    new3 <- aic_ar(x1, p = 1:5, type = "aicc", method = "mle")

    expect_equal(new3$p, old3$p)
    expect_equal(new3$phi, old3$phi, tolerance = 1e-3)

    # Test case 4: Burg method
    old4 <- tswge::aic.ar.wge(x1, p = 1:5, type = "aic", method = "burg")
    new4 <- aic_ar(x1, p = 1:5, type = "aic", method = "burg")

    expect_equal(new4$p, old4$p)
    expect_equal(new4$phi, old4$phi, tolerance = 1e-10)
  })
})


test_that("aic_ar validates inputs", {

  x <- rnorm(100)

  expect_error(aic_ar("not numeric", p = 1:5), "`x` must be")
  expect_error(aic_ar(numeric(0), p = 1:5), "`x` must be")
  expect_error(aic_ar(x, p = -1:5), "`p` must be")
})


test_that("aic_ar returns comparison table", {

  capture.output({
    result <- aic_ar(rnorm(100), p = 1:5)
  })

  expect_true("table" %in% names(result))
  expect_equal(nrow(result$table), 5)
  expect_true(all(c("p", "aic", "aicc", "bic") %in% names(result$table)))
})


test_that("aic_ar has correct structure", {

  capture.output({
    result <- aic_ar(rnorm(100), p = 1:5)
  })

  expect_s3_class(result, "aic_ar")
  expect_named(result, c("type", "method", "value", "p", "phi",
                         "xbar", "vara", "table"))
})


test_that("aic_ar print method works", {

  capture.output({
    result <- aic_ar(rnorm(100), p = 1:5)
  })

  expect_output(print(result), "AR Order Selection")
})
