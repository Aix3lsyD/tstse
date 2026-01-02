test_that("true_farma_acf returns correct structure", {
  result <- true_farma_acf(d = 0.3, plot = FALSE)

  expect_type(result, "list")
  expect_named(result, c("acf", "acv"))
  expect_length(result$acf, 51)  # lags 0 to 50
  expect_length(result$acv, 51)
})

test_that("true_farma_acf ACF starts at 1", {
  result <- true_farma_acf(d = 0.3, plot = FALSE)
  expect_equal(result$acf[1], 1)

  result <- true_farma_acf(d = 0.2, phi = 0.5, plot = FALSE)
  expect_equal(result$acf[1], 1)

  result <- true_farma_acf(d = 0.4, theta = 0.3, plot = FALSE)
  expect_equal(result$acf[1], 1)
})

test_that("true_farma_acf reduces to ARMA when d=0", {
  # When d=0, FARMA should equal standard ARMA
  farma_result <- true_farma_acf(d = 0, phi = 0.8, plot = FALSE, lag_max = 20)
  arma_result <- true_acf(phi = 0.8, lag_max = 20, plot = FALSE)

  expect_equal(farma_result$acf, arma_result$acf, tolerance = 1e-10)

  # With MA term
  farma_result <- true_farma_acf(d = 0, phi = 0.5, theta = 0.3, plot = FALSE, lag_max = 20)
  arma_result <- true_acf(phi = 0.5, theta = 0.3, lag_max = 20, plot = FALSE)

  expect_equal(farma_result$acf, arma_result$acf, tolerance = 1e-10)
})

test_that("true_farma_acf shows long memory decay for positive d", {
  # Long memory: ACF decays slowly (hyperbolically)
  result <- true_farma_acf(d = 0.4, plot = FALSE, lag_max = 100)

  # ACF should be positive and slowly decaying
  expect_true(all(result$acf > 0))

  # Compare decay rate - long memory decays slower than exponential
  # At lag 100, should still be relatively high
  expect_gt(result$acf[101], 0.1)  # lag 100 value (index 101)
})

test_that("true_farma_acf shows antipersistence for negative d", {
  # Negative d: antipersistence
  result <- true_farma_acf(d = -0.3, plot = FALSE, lag_max = 50)

  # ACF should have negative values (anti-correlation)
  expect_true(any(result$acf[-1] < 0))
})

test_that("true_farma_acf respects lag_max", {
  result10 <- true_farma_acf(d = 0.3, lag_max = 10, plot = FALSE)
  expect_length(result10$acf, 11)  # lags 0 to 10

  result30 <- true_farma_acf(d = 0.3, lag_max = 30, plot = FALSE)
  expect_length(result30$acf, 31)  # lags 0 to 30

  # First 11 values should be the same
  expect_equal(result10$acf, result30$acf[1:11], tolerance = 1e-10)
})

test_that("true_farma_acf validates inputs", {
  expect_error(true_farma_acf(d = "a", plot = FALSE),
               "`d` must be a single numeric value")
  expect_error(true_farma_acf(d = 0.3, phi = "a", plot = FALSE),
               "`phi` must be a numeric vector")
  expect_error(true_farma_acf(d = 0.3, theta = "a", plot = FALSE),
               "`theta` must be a numeric vector")
  expect_error(true_farma_acf(d = 0.3, lag_max = -1, plot = FALSE),
               "`lag_max` must be a positive integer")
  expect_error(true_farma_acf(d = 0.3, trunc = 0, plot = FALSE),
               "`trunc` must be a positive integer")
  expect_error(true_farma_acf(d = 0.3, vara = 0, plot = FALSE),
               "`vara` must be a positive number")
  expect_error(true_farma_acf(d = 0.3, plot = "yes"),
               "`plot` must be TRUE or FALSE")
})

test_that("true_farma_acf warns for d outside typical range", {
  expect_warning(true_farma_acf(d = 0.6, plot = FALSE),
                 "`d` outside typical range")
  expect_warning(true_farma_acf(d = -0.6, plot = FALSE),
                 "`d` outside typical range")
})

test_that("true_farma_acf plots without error", {
  expect_silent({
    pdf(NULL)
    result <- true_farma_acf(d = 0.3, lag_max = 20)
    dev.off()
  })
})

test_that("true_farma_acf matches tswge::true.farma.aut.wge", {
  skip_if_not_installed("tswge")

  # Pure fractional difference
  old <- suppressWarnings(tswge::true.farma.aut.wge(d = 0.3, phi = 0, theta = 0,
                                                     lag.max = 30, trunc = 1000,
                                                     vara = 1, plot = FALSE))
  new <- true_farma_acf(d = 0.3, phi = 0, theta = 0, lag_max = 30, trunc = 1000,
                        vara = 1, plot = FALSE)
  expect_equal(new$acf, old$acf, tolerance = 1e-10)
  expect_equal(new$acv, old$acv, tolerance = 1e-10)

  # With AR term
  old <- suppressWarnings(tswge::true.farma.aut.wge(d = 0.2, phi = 0.5, theta = 0,
                                                     lag.max = 25, plot = FALSE))
  new <- true_farma_acf(d = 0.2, phi = 0.5, theta = 0, lag_max = 25, plot = FALSE)
  expect_equal(new$acf, old$acf, tolerance = 1e-10)

  # With MA term
  old <- suppressWarnings(tswge::true.farma.aut.wge(d = 0.25, phi = 0, theta = 0.4,
                                                     lag.max = 20, plot = FALSE))
  new <- true_farma_acf(d = 0.25, phi = 0, theta = 0.4, lag_max = 20, plot = FALSE)
  expect_equal(new$acf, old$acf, tolerance = 1e-10)
})
