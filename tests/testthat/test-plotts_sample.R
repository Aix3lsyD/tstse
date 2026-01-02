# tests/testthat/test-plotts_sample.R

library(tswge)

test_that("plotts_sample matches original plotts.sample.wge output", {

  set.seed(123)
  x <- tswge::gen.arma.wge(n = 200, phi = c(1.5, -0.75), plot = FALSE, sn = 123)

  # Test without periodogram
  old1 <- tswge::plotts.sample.wge(x, lag.max = 25, trunc = 0, periodogram = FALSE)
  new1 <- plotts_sample(x, lag_max = 25, trunc = 0, periodogram = FALSE)

  expect_equal(new1$xbar, old1$xbar, tolerance = 1e-10)
  expect_equal(new1$acf, old1$autplt, tolerance = 1e-10)
  expect_equal(new1$freq, old1$freq, tolerance = 1e-10)
  expect_equal(new1$parzen_db, old1$dbz, tolerance = 1e-10)

  # Test with periodogram
  old2 <- tswge::plotts.sample.wge(x, lag.max = 25, trunc = 0, periodogram = TRUE)
  new2 <- plotts_sample(x, lag_max = 25, trunc = 0, periodogram = TRUE)

  expect_equal(new2$xbar, old2$xbar, tolerance = 1e-10)
  expect_equal(new2$parzen_db, old2$dbz, tolerance = 1e-10)
  expect_equal(new2$periodogram_db, old2$db, tolerance = 1e-10)
})


test_that("plotts_sample validates inputs", {

  expect_error(plotts_sample("not numeric"), "`x` must be")
  expect_error(plotts_sample(numeric(0)), "`x` must be")
  expect_error(plotts_sample(rnorm(100), lag_max = 0), "`lag_max` must be")
})


test_that("plotts_sample adjusts lag_max if too large", {

  result <- plotts_sample(rnorm(20), lag_max = 100)
  expect_equal(length(result$acf), 20)  # n, since includes lag 0
})


test_that("plotts_sample returns correct structure", {

  # Without periodogram
  result1 <- plotts_sample(rnorm(100), periodogram = FALSE)
  expect_true(all(c("xbar", "acf", "freq", "parzen_db") %in% names(result1)))
  expect_false("periodogram_db" %in% names(result1))

  # With periodogram
  result2 <- plotts_sample(rnorm(100), periodogram = TRUE)
  expect_true("periodogram_db" %in% names(result2))
})


test_that("plotts_sample respects custom trunc", {

  set.seed(123)
  x <- rnorm(100)

  result1 <- plotts_sample(x, trunc = 10)

  set.seed(123)
  x <- rnorm(100)

  result2 <- plotts_sample(x, trunc = 20)

  # Different truncation should give different spectral estimates
  expect_false(identical(result1$parzen_db, result2$parzen_db))
})


test_that("plotts_sample works with arlimits", {

  expect_no_error(plotts_sample(rnorm(100), arlimits = TRUE))
})


test_that("plotts_sample works with custom speclimits", {

  expect_no_error(plotts_sample(rnorm(100), speclimits = c(-20, 20)))
})
