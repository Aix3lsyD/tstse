skip_if_not_installed("MASS")

test_that("true_garma_acf returns correct structure", {
  result <- true_garma_acf(u = 0.5, lambda = 0.3, lag_max = 10, plot = FALSE)

  expect_type(result, "list")
  expect_named(result, c("acf", "acv"))
  expect_length(result$acf, 11)  # lag 0 to lag_max
  expect_length(result$acv, 11)
})


test_that("true_garma_acf first autocorrelation is 1", {
  result <- true_garma_acf(u = 0.3, lambda = 0.4, lag_max = 5, plot = FALSE)

  expect_equal(result$acf[1], 1, tolerance = 1e-10)
})


test_that("true_garma_acf matches tswge for pure Gegenbauer", {
  skip_if_not_installed("tswge")

  old <- suppressWarnings(
    tswge::true.garma.aut.wge(u = 0.5, lambda = 0.3, lag.max = 15, plot = FALSE)
  )
  new <- true_garma_acf(u = 0.5, lambda = 0.3, lag_max = 15, plot = FALSE)

  expect_equal(new$acf, old$acf, tolerance = 1e-6)
  expect_equal(new$acv, old$acv, tolerance = 1e-6)
})


test_that("true_garma_acf matches tswge with AR/MA components", {
  skip_if_not_installed("tswge")

  # With AR and MA
  old <- suppressWarnings(
    tswge::true.garma.aut.wge(u = 0.2, lambda = 0.25, phi = 0.4, theta = 0.3,
                              lag.max = 10, plot = FALSE)
  )
  new <- true_garma_acf(u = 0.2, lambda = 0.25, phi = 0.4, theta = 0.3,
                        lag_max = 10, plot = FALSE)

  expect_equal(new$acf, old$acf, tolerance = 1e-6)
})


test_that("true_garma_acf shows slow decay (long memory)", {
  # Use u close to 1 for standard long memory behavior
  result <- true_garma_acf(u = 0.9, lambda = 0.4, lag_max = 30, plot = FALSE)

  # Long memory should show slow decay - ACF values should remain notable at high lags
  expect_gt(result$acf[16], 0.1)  # Lag 15 should still have notable correlation
})


test_that("true_garma_acf validates u parameter", {
  expect_error(true_garma_acf(u = 1, lambda = 0.3),
               "`u` must be a single numeric value in \\(-1, 1\\)")

  expect_error(true_garma_acf(u = -1, lambda = 0.3),
               "`u` must be a single numeric value in \\(-1, 1\\)")

  expect_error(true_garma_acf(u = 1.5, lambda = 0.3),
               "`u` must be a single numeric value in \\(-1, 1\\)")

  expect_error(true_garma_acf(u = "bad", lambda = 0.3),
               "`u` must be a single numeric value")
})


test_that("true_garma_acf validates lambda parameter", {
  expect_error(true_garma_acf(u = 0.5, lambda = 0),
               "`lambda` must be a positive number")

  expect_error(true_garma_acf(u = 0.5, lambda = -0.5),
               "`lambda` must be a positive number")

  expect_error(true_garma_acf(u = 0.5, lambda = "bad"),
               "`lambda` must be a positive number")
})


test_that("true_garma_acf validates other parameters", {
  expect_error(true_garma_acf(u = 0.5, lambda = 0.3, phi = "bad"),
               "`phi` must be numeric")

  expect_error(true_garma_acf(u = 0.5, lambda = 0.3, theta = "bad"),
               "`theta` must be numeric")

  expect_error(true_garma_acf(u = 0.5, lambda = 0.3, lag_max = -5),
               "`lag_max` must be a non-negative integer")

  expect_error(true_garma_acf(u = 0.5, lambda = 0.3, vara = 0),
               "`vara` must be positive")
})


test_that("true_garma_acf with different vara scales autocovariance", {
  result1 <- true_garma_acf(u = 0.5, lambda = 0.3, vara = 1, lag_max = 5, plot = FALSE)
  result2 <- true_garma_acf(u = 0.5, lambda = 0.3, vara = 4, lag_max = 5, plot = FALSE)

  # ACF should be the same
  expect_equal(result1$acf, result2$acf, tolerance = 1e-5)

  # ACV should scale with vara
  expect_equal(result2$acv / result1$acv, rep(4, length(result1$acv)), tolerance = 1e-4)
})


test_that("true_garma_acf plotting works", {
  expect_no_error({
    pdf(NULL)
    true_garma_acf(u = 0.5, lambda = 0.3, lag_max = 10, plot = TRUE)
    dev.off()
  })
})
