# tests/testthat/test-aic.R

library(tswge)

test_that("aic_ts matches original aic.wge", {

  # Test case 1: AIC selection on ARMA data
  set.seed(123)
  x1 <- tswge::gen.arma.wge(n = 200, phi = 0.7, theta = 0.3, plot = FALSE, sn = 123)

  old1 <- tswge::aic.wge(x1, p = 0:3, q = 0:2, type = "aic")
  new1 <- aic_ts(x1, p = 0:3, q = 0:2, type = "aic")

  expect_equal(new1$p, old1$p)
  expect_equal(new1$q, old1$q)
  expect_equal(new1$phi, old1$phi, tolerance = 1e-3)
  expect_equal(new1$theta, old1$theta, tolerance = 1e-3)
  expect_equal(new1$value, old1$value, tolerance = 1e-4)

  # Test case 2: BIC selection
  old2 <- tswge::aic.wge(x1, p = 0:3, q = 0:2, type = "bic")
  new2 <- aic_ts(x1, p = 0:3, q = 0:2, type = "bic")

  expect_equal(new2$p, old2$p)
  expect_equal(new2$q, old2$q)

  # Test case 3: AICC selection
  old3 <- tswge::aic.wge(x1, p = 0:3, q = 0:2, type = "aicc")
  new3 <- aic_ts(x1, p = 0:3, q = 0:2, type = "aicc")

  expect_equal(new3$p, old3$p)
  expect_equal(new3$q, old3$q)

  # Test case 4: AR-only data
  set.seed(456)
  x2 <- tswge::gen.arma.wge(n = 200, phi = c(1.5, -0.75), plot = FALSE, sn = 456)

  old4 <- tswge::aic.wge(x2, p = 0:4, q = 0:2, type = "aic")
  new4 <- aic_ts(x2, p = 0:4, q = 0:2, type = "aic")

  expect_equal(new4$p, old4$p)
  expect_equal(new4$q, old4$q)
})


test_that("aic_ts validates inputs", {

  x <- rnorm(100)

  expect_error(aic_ts("not numeric", p = 0:3, q = 0:2), "`x` must be")
  expect_error(aic_ts(numeric(0), p = 0:3, q = 0:2), "`x` must be")
  expect_error(aic_ts(x, p = -1:3, q = 0:2), "`p` must be")
  expect_error(aic_ts(x, p = 0:3, q = -1:2), "`q` must be")
})


test_that("aic_ts returns comparison table", {

  result <- aic_ts(rnorm(100), p = 0:2, q = 0:1)

  expect_true("table" %in% names(result))
  expect_equal(nrow(result$table), 6)  # 3 * 2 combinations
  expect_true(all(c("p", "q", "aic", "aicc", "bic") %in% names(result$table)))
})


test_that("aic_ts has correct structure", {

  result <- aic_ts(rnorm(100), p = 0:2, q = 0:1)

  expect_s3_class(result, "aic_arma")
  expect_named(result, c("type", "value", "p", "q", "phi", "theta",
                         "xbar", "vara", "table"))
})


test_that("aic_ts print method works", {

  result <- aic_ts(rnorm(100), p = 0:2, q = 0:1)
  expect_output(print(result), "ARMA Order Selection")
})


test_that("aic_ts handles failed models gracefully", {

  # Very short series where some models may fail
  x <- rnorm(30)
  expect_no_error(aic_ts(x, p = 0:2, q = 0:2))
})
