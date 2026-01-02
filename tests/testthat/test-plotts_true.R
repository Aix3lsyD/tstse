library(tswge)

test_that("plotts_true returns correct structure", {

  result <- plotts_true(n = 100, phi = 0.8, plot = FALSE, seed = 123)

  expect_true(is.list(result))
  expect_named(result, c("data", "acf", "acv", "freq", "spec"))
})


test_that("plotts_true returns correct dimensions", {

  result <- plotts_true(n = 100, phi = 0.8, lag_max = 25, plot = FALSE, seed = 123)

  expect_equal(length(result$data), 100)
  expect_equal(length(result$acf), 26)   # lag_max + 1

  expect_equal(length(result$acv), 26)
  expect_equal(length(result$freq), 251)
  expect_equal(length(result$spec), 251)
})


test_that("plotts_true seed produces reproducible results", {

  result1 <- plotts_true(n = 100, phi = 0.8, plot = FALSE, seed = 456)
  result2 <- plotts_true(n = 100, phi = 0.8, plot = FALSE, seed = 456)

  expect_equal(result1$data, result2$data)
  expect_equal(result1$acf, result2$acf)
  expect_equal(result1$spec, result2$spec)
})


test_that("plotts_true ACF matches tswge", {

  # Note: tswge returns aut1, we return acf
  old <- tswge::plotts.true.wge(n = 100, phi = 0.8, sn = 123, plot.data = FALSE)
  new <- plotts_true(n = 100, phi = 0.8, lag_max = 25, plot = FALSE, seed = 123)

  expect_equal(new$acf, old$aut1, tolerance = 1e-10)
})


test_that("plotts_true ACV matches tswge", {

  old <- tswge::plotts.true.wge(n = 100, phi = c(1.5, -0.75), sn = 123, plot.data = FALSE)
  new <- plotts_true(n = 100, phi = c(1.5, -0.75), lag_max = 25, plot = FALSE, seed = 123)

  expect_equal(new$acv, old$acv, tolerance = 1e-6)
})


test_that("plotts_true spec matches tswge", {

  old <- tswge::plotts.true.wge(n = 100, phi = 0.7, theta = 0.4, sn = 123, plot.data = FALSE)
  new <- plotts_true(n = 100, phi = 0.7, theta = 0.4, lag_max = 25, plot = FALSE, seed = 123)

  expect_equal(new$spec, old$spec, tolerance = 1e-6)
})


test_that("plotts_true works with various model types", {

  # AR(1)
  expect_no_error(plotts_true(n = 50, phi = 0.8, plot = FALSE, seed = 1))

  # AR(2)
  expect_no_error(plotts_true(n = 50, phi = c(1.5, -0.75), plot = FALSE, seed = 2))

  # MA(1)
  expect_no_error(plotts_true(n = 50, theta = 0.6, plot = FALSE, seed = 3))

  # ARMA(1,1)
  expect_no_error(plotts_true(n = 50, phi = 0.7, theta = 0.4, plot = FALSE, seed = 4))

  # White noise
  expect_no_error(plotts_true(n = 50, plot = FALSE, seed = 5))
})


test_that("plotts_true plotting works", {

  expect_no_error(plotts_true(n = 50, phi = 0.8, plot = TRUE, seed = 123))
})


test_that("plotts_true validates inputs", {

  expect_error(plotts_true(n = -1), "`n` must be")
  expect_error(plotts_true(n = 100, phi = "bad"), "`phi` must be")
  expect_error(plotts_true(n = 100, theta = "bad"), "`theta` must be")
  expect_error(plotts_true(n = 100, lag_max = -1), "`lag_max` must be")
  expect_error(plotts_true(n = 100, mu = "bad"), "`mu` must be")
  expect_error(plotts_true(n = 100, vara = -1), "`vara` must be")
})


test_that("plotts_true caps lag_max at n-1", {

  result <- plotts_true(n = 20, lag_max = 50, plot = FALSE, seed = 123)

  # lag_max should be capped at 19 (n-1), so acf length = 20
  expect_equal(length(result$acf), 20)
})


test_that("plotts_true freq range is correct", {

  result <- plotts_true(n = 100, phi = 0.5, plot = FALSE, seed = 123)

  expect_equal(result$freq[1], 0)
  expect_equal(result$freq[251], 0.5)
})


test_that("plotts_true always includes data even when plot=FALSE", {

  # This is an intentional deviation from tswge which doesn't generate data

  # when plot.data=FALSE
  result <- plotts_true(n = 100, phi = 0.8, plot = FALSE, seed = 123)

  expect_true("data" %in% names(result))
  expect_equal(length(result$data), 100)
})
