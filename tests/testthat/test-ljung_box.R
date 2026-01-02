# tests/testthat/test-ljung_box.R

library(tswge)

test_that("ljung_box matches original ljung.wge", {

  # Suppress the cat() output from original
  capture.output({

    # Test case 1: Default parameters
    set.seed(123)
    x1 <- rnorm(200)

    old1 <- tswge::ljung.wge(x1, K = 24, p = 0, q = 0)
    new1 <- ljung_box(x1, K = 24, p = 0, q = 0)

    expect_equal(new1$chi_square, old1$chi.square, tolerance = 1e-10)
    expect_equal(new1$df, old1$df)
    expect_true(abs(new1$pval - old1$pval) < 1e-10)

    # Test case 2: With p and q adjustment
    set.seed(456)
    x2 <- cumsum(rnorm(150))

    old2 <- tswge::ljung.wge(x2, K = 20, p = 2, q = 1)
    new2 <- ljung_box(x2, K = 20, p = 2, q = 1)

    expect_equal(new2$chi_square, old2$chi.square, tolerance = 1e-10)
    expect_equal(new2$df, old2$df)
    expect_true(abs(new2$pval - old2$pval) < 1e-10)

    # Test case 3: Different K
    set.seed(789)
    x3 <- arima.sim(model = list(ar = 0.7), n = 100)

    old3 <- tswge::ljung.wge(x3, K = 12, p = 0, q = 0)
    new3 <- ljung_box(x3, K = 12, p = 0, q = 0)

    expect_equal(new3$chi_square, old3$chi.square, tolerance = 1e-10)
    expect_true(abs(new3$pval - old3$pval) < 1e-10)
  })
})

test_that("ljung_box validates inputs", {

  x <- rnorm(100)

  expect_error(ljung_box("not numeric"), "`x` must be")
  expect_error(ljung_box(numeric(0)), "`x` must be")
  expect_error(ljung_box(x, K = 0), "`K` must be")
  expect_error(ljung_box(x, K = 150), "`K` must be less than")
  expect_error(ljung_box(x, p = -1), "`p` must be")
  expect_error(ljung_box(x, K = 10, p = 5, q = 6), "Degrees of freedom")
})


test_that("ljung_box has correct structure", {

  result <- ljung_box(rnorm(100))

  expect_s3_class(result, "ljung_box_test")
  expect_named(result, c("test", "K", "chi_square", "df", "pval"))
})


test_that("ljung_box print method works", {

  result <- ljung_box(rnorm(100))
  expect_output(print(result), "Ljung-Box test")
})
