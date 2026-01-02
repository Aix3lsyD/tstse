# tests/testthat/test-aic_burg.R

library(tswge)

test_that("aic_burg matches original aic.burg.wge", {

  capture.output({

    # Test case 1: AIC selection
    set.seed(123)
    x1 <- tswge::gen.arma.wge(n = 200, phi = c(1.5, -0.75), plot = FALSE, sn = 123)

    old1 <- tswge::aic.burg.wge(x1, p = 1:5, type = "aic")
    new1 <- aic_burg(x1, p = 1:5, type = "aic")

    expect_equal(new1$p, old1$p)
    expect_equal(new1$phi, old1$phi, tolerance = 1e-10)
    expect_equal(new1$value, old1$value, tolerance = 1e-4)
    expect_equal(new1$vara, old1$vara, tolerance = 1e-4)

    # Test case 2: BIC selection
    old2 <- tswge::aic.burg.wge(x1, p = 1:5, type = "bic")
    new2 <- aic_burg(x1, p = 1:5, type = "bic")

    expect_equal(new2$p, old2$p)
    expect_equal(new2$phi, old2$phi, tolerance = 1e-10)

    # Test case 3: AICC selection
    old3 <- tswge::aic.burg.wge(x1, p = 1:5, type = "aicc")
    new3 <- aic_burg(x1, p = 1:5, type = "aicc")

    expect_equal(new3$p, old3$p)
    expect_equal(new3$phi, old3$phi, tolerance = 1e-10)
  })
})


test_that("aic_burg has correct structure", {

  capture.output({
    result <- aic_burg(rnorm(100), p = 1:5)
  })

  expect_s3_class(result, "aic_burg")
  expect_named(result, c("type", "value", "p", "phi", "vara"))
})


test_that("aic_burg print method works", {

  capture.output({
    result <- aic_burg(rnorm(100), p = 1:5)
  })

  expect_output(print(result), "Burg")
})
