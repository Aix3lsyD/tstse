# tests/testthat/test-pacf_ts.R

library(tswge)

test_that("pacf_ts matches original pacfts.wge", {

  capture.output({

    # Test case 1: Yule-Walker method
    set.seed(123)
    x1 <- tswge::gen.arma.wge(n = 200, phi = c(1.5, -0.75), plot = FALSE, sn = 123)

    old1 <- tswge::pacfts.wge(x1, lag.max = 10, method = "yw", plot = FALSE)
    new1 <- pacf_ts(x1, lag_max = 10, method = "yw", plot = FALSE)

    expect_equal(new1$pacf, old1$pacf_yw, tolerance = 1e-10)

    # Test case 2: Burg method
    old2 <- tswge::pacfts.wge(x1, lag.max = 10, method = "burg", plot = FALSE)
    new2 <- pacf_ts(x1, lag_max = 10, method = "burg", plot = FALSE)

    expect_equal(new2$pacf, old2$pacf_burg, tolerance = 1e-10)

    # Test case 3: MLE method
    old3 <- tswge::pacfts.wge(x1, lag.max = 8, method = "mle", plot = FALSE)
    new3 <- pacf_ts(x1, lag_max = 8, method = "mle", plot = FALSE)

    expect_equal(new3$pacf, old3$pacf_mle, tolerance = 1e-3)  # MLE has numerical variation
  })
})


test_that("pacf_ts validates inputs", {

  expect_error(pacf_ts("not numeric"), "`x` must be")
  expect_error(pacf_ts(numeric(0)), "`x` must be")
  expect_error(pacf_ts(rnorm(100), lag_max = 0), "`lag_max` must be")
})


test_that("pacf_ts returns correct length", {

  capture.output({
    result <- pacf_ts(rnorm(100), lag_max = 15, plot = FALSE)
  })

  expect_equal(length(result$pacf), 15)
})


test_that("pacf_ts adjusts lag_max if too large", {

  capture.output({
    result <- pacf_ts(rnorm(20), lag_max = 100, plot = FALSE)
  })

  expect_equal(length(result$pacf), 19)  # n - 1
})


test_that("pacf_ts works with parallel", {

  skip_on_cran()
  skip_if(parallel::detectCores(logical = FALSE) < 2, "Not enough cores")

  x <- rnorm(100)

  capture.output({
    result_seq <- pacf_ts(x, lag_max = 10, cores = 1, plot = FALSE)
    result_par <- pacf_ts(x, lag_max = 10, cores = 2, plot = FALSE)
  })

  expect_equal(result_seq$pacf, result_par$pacf)
})


test_that("pacf_ts plot works", {

  x <- rnorm(100)

  # Should not error
  expect_no_error(capture.output(pacf_ts(x, lag_max = 10, plot = TRUE)))
  expect_no_error(capture.output(pacf_ts(x, lag_max = 10, plot = TRUE, limits = TRUE)))
})
