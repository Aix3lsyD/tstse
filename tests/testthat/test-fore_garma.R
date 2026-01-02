test_that("fore_garma returns correct structure", {
  x <- gen_garma(n = 200, u = 0.8, lambda = 0.3, phi = 0.5,
                 plot = FALSE, seed = 123)

  fit <- fore_garma(x, u = 0.8, lambda = 0.3, phi = 0.5,
                    n_ahead = 10, plot = FALSE)

  expect_s3_class(fit, "fore_garma")
  expect_named(fit, c("f", "ar_f", "ar_order", "original", "n_ahead", "lastn"))
  expect_length(fit$f, 10)
  expect_length(fit$ar_f, 10)
  expect_true(is.numeric(fit$f))
  expect_true(is.numeric(fit$ar_f))
})

test_that("fore_garma holdout mode (lastn = TRUE) works", {
  x <- gen_garma(n = 200, u = 0.8, lambda = 0.3, phi = 0.5,
                 plot = FALSE, seed = 123)

  fit <- fore_garma(x, u = 0.8, lambda = 0.3, phi = 0.5,
                    n_ahead = 10, lastn = TRUE, plot = FALSE)

  expect_length(fit$f, 10)
  expect_true(fit$lastn)
  expect_equal(fit$n_ahead, 10)
})

test_that("fore_garma future mode (lastn = FALSE) works", {
  x <- gen_garma(n = 200, u = 0.8, lambda = 0.3, phi = 0.5,
                 plot = FALSE, seed = 123)

  fit <- fore_garma(x, u = 0.8, lambda = 0.3, phi = 0.5,
                    n_ahead = 20, lastn = FALSE, plot = FALSE)

  expect_length(fit$f, 20)
  expect_false(fit$lastn)
  expect_equal(fit$n_ahead, 20)
})

test_that("fore_garma works with k=1 (single factor)", {
  x <- gen_garma(n = 150, u = 0.5, lambda = 0.3, plot = FALSE, seed = 42)

  # Holdout
  fit1 <- fore_garma(x, u = 0.5, lambda = 0.3, n_ahead = 10,
                     lastn = TRUE, plot = FALSE)
  expect_true(is.numeric(fit1$f))
  expect_false(anyNA(fit1$f))

  # Future
  fit2 <- fore_garma(x, u = 0.5, lambda = 0.3, n_ahead = 15,
                     lastn = FALSE, plot = FALSE)
  expect_true(is.numeric(fit2$f))
  expect_false(anyNA(fit2$f))
})

test_that("fore_garma works with k=2 (two factors)", {
  x <- gen_garma(n = 150, u = c(0.8, 0.3), lambda = c(0.3, 0.2),
                 plot = FALSE, seed = 42)

  # Holdout
  fit1 <- fore_garma(x, u = c(0.8, 0.3), lambda = c(0.3, 0.2),
                     n_ahead = 10, lastn = TRUE, plot = FALSE)
  expect_true(is.numeric(fit1$f))
  expect_length(fit1$f, 10)

  # Future
  fit2 <- fore_garma(x, u = c(0.8, 0.3), lambda = c(0.3, 0.2),
                     n_ahead = 15, lastn = FALSE, plot = FALSE)
  expect_true(is.numeric(fit2$f))
  expect_length(fit2$f, 15)
})

test_that("fore_garma works with AR/MA components", {
  x <- gen_garma(n = 150, u = 0.5, lambda = 0.3, phi = 0.6, theta = 0.3,
                 plot = FALSE, seed = 42)

  fit <- fore_garma(x, u = 0.5, lambda = 0.3, phi = 0.6, theta = 0.3,
                    n_ahead = 10, plot = FALSE)

  expect_true(is.numeric(fit$f))
  expect_length(fit$f, 10)
})

test_that("fore_garma validates inputs", {
  x <- rnorm(100)

  expect_error(fore_garma(1:10, u = 0.5, lambda = 0.3, plot = FALSE),
               "at least 20 observations")
  expect_error(fore_garma(x, u = "a", lambda = 0.3, plot = FALSE),
               "`u` must be numeric")
  expect_error(fore_garma(x, u = 0.5, lambda = "a", plot = FALSE),
               "`lambda` must be numeric")
  expect_error(fore_garma(x, u = c(0.5, 0.3), lambda = 0.3, plot = FALSE),
               "same length")
  expect_error(fore_garma(x, u = c(0.5, 0.3, 0.1), lambda = c(0.3, 0.2, 0.1),
                          plot = FALSE),
               "Only 1 or 2 Gegenbauer factors")
  expect_error(fore_garma(x, u = 0.5, lambda = 0.3, n_ahead = 0, plot = FALSE),
               "`n_ahead` must be a positive integer")
  expect_error(fore_garma(x, u = 0.5, lambda = 0.3, lastn = "yes", plot = FALSE),
               "`lastn` must be TRUE or FALSE")
  expect_error(fore_garma(x, u = 0.5, lambda = 0.3, plot = "yes"),
               "`plot` must be TRUE or FALSE")
})

test_that("fore_garma print method works", {
  x <- gen_garma(n = 100, u = 0.5, lambda = 0.2, plot = FALSE, seed = 123)
  fit <- fore_garma(x, u = 0.5, lambda = 0.2, n_ahead = 5, plot = FALSE)

  output <- capture.output(print(fit))
  expect_true(any(grepl("GARMA Forecast", output)))
  expect_true(any(grepl("Mode", output)))
  expect_true(any(grepl("Steps ahead", output)))
})

test_that("fore_garma plots without error", {
  x <- gen_garma(n = 100, u = 0.5, lambda = 0.3, plot = FALSE, seed = 123)

  # Holdout mode
  expect_silent({
    pdf(NULL)
    fit <- fore_garma(x, u = 0.5, lambda = 0.3, n_ahead = 10,
                      lastn = TRUE, plot = TRUE)
    dev.off()
  })

  # Future mode
  expect_silent({
    pdf(NULL)
    fit <- fore_garma(x, u = 0.5, lambda = 0.3, n_ahead = 10,
                      lastn = FALSE, plot = TRUE)
    dev.off()
  })
})

test_that("fore_garma produces reasonable forecasts", {
  set.seed(123)
  x <- gen_garma(n = 200, u = 0.8, lambda = 0.3, phi = 0.5, plot = FALSE)

  fit <- fore_garma(x, u = 0.8, lambda = 0.3, phi = 0.5,
                    n_ahead = 10, lastn = TRUE, plot = FALSE)

  # Forecasts should be within reasonable range of data
  data_range <- range(x)
  expect_true(all(fit$f > data_range[1] - 3 * sd(x)))
  expect_true(all(fit$f < data_range[2] + 3 * sd(x)))
})

test_that("fore_garma matches tswge::fore.garma.wge", {
  skip_if_not_installed("tswge")

  # Generate data with tswge for exact reproducibility
  x <- tswge::gen.garma.wge(n = 150, u = 0.8, lambda = 0.3, phi = 0.5,
                            plot = FALSE, sn = 456)

  # Compare holdout mode
  old <- suppressWarnings(
    tswge::fore.garma.wge(x, u = 0.8, lambda = 0.3, phi = 0.5,
                          n.ahead = 10, lastn = TRUE, plot = FALSE)
  )
  new <- fore_garma(x, u = 0.8, lambda = 0.3, phi = 0.5,
                    n_ahead = 10, lastn = TRUE, plot = FALSE)

  # Compare GARMA forecasts

  expect_equal(new$f, old$garma.fore, tolerance = 1e-6)

  # Compare AR forecasts
  expect_equal(new$ar_f, old$ar.fore, tolerance = 1e-6)

  # Compare AR order
  expect_equal(new$ar_order, old$ar.fit.order)
})
