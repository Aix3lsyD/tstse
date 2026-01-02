# tests/testthat/test-plotts.R

test_that("plotts works with base R style", {

  x <- rnorm(100)

  expect_no_error(plotts(x, style = 0))
  expect_no_error(plotts(x, style = 0, main = "Test", col = "blue"))
  expect_no_error(plotts(x, style = 0, xlim = c(0, 50), ylim = c(-2, 2)))
})


test_that("plotts handles short and long series differently", {

  # Short series (n <= 200) - should use points
  x_short <- rnorm(100)
  expect_no_error(plotts(x_short))

  # Long series (n > 200) - should use lines only
  x_long <- rnorm(500)
  expect_no_error(plotts(x_long))
})


test_that("plotts works with ts objects", {

  x_ts <- ts(rnorm(100), start = c(2020, 1), frequency = 12)

  expect_no_error(plotts(x_ts, style = 0))
})


test_that("plotts ggplot2 style works", {

  skip_if_not_installed("ggplot2")

  x <- rnorm(100)

  p <- plotts(x, style = 1)
  expect_s3_class(p, "ggplot")
})


test_that("plotts ggplot2 style works with ts objects", {

  skip_if_not_installed("ggplot2")
  skip_if_not_installed("zoo")

  x_ts <- ts(rnorm(100), start = c(2020, 1), frequency = 12)

  p <- plotts(x_ts, style = 1)
  expect_s3_class(p, "ggplot")
})


test_that("plotts validates inputs", {

  expect_error(plotts("not numeric"), "`x` must be")
  expect_error(plotts(numeric(0)), "`x` must be")
  expect_error(plotts(rnorm(100), style = 2), "`style` must be")
})
