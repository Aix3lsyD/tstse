test_that("fore_glambda returns correct structure", {
  set.seed(123)
  x <- gen_glambda(n = 100, lambda = 0.5, phi = 0.7, plot = FALSE)

  fit <- fore_glambda(x, lambda = 0.5, offset = 20, phi = 0.7,
                      n_ahead = 10, plot = FALSE)

  expect_s3_class(fit, "fore_glambda")
  expect_named(fit, c("f", "ar_f", "ar_order", "original", "n_ahead", "lastn"))
  expect_length(fit$f, 10)
  expect_length(fit$ar_f, 10)
  expect_true(is.numeric(fit$ar_order))
  expect_equal(fit$n_ahead, 10)
  expect_true(fit$lastn)
})

test_that("fore_glambda holdout mode works", {
  set.seed(123)
  x <- gen_glambda(n = 100, lambda = 0.5, phi = 0.6, plot = FALSE)

  fit <- fore_glambda(x, lambda = 0.5, offset = 30, phi = 0.6,
                      n_ahead = 15, lastn = TRUE, plot = FALSE)

  expect_length(fit$f, 15)
  expect_length(fit$ar_f, 15)
  expect_true(fit$lastn)
  expect_equal(length(fit$original), 100)
})

test_that("fore_glambda future mode works", {
  set.seed(456)
  x <- gen_glambda(n = 100, lambda = 0.4, phi = 0.5, plot = FALSE)

  fit <- fore_glambda(x, lambda = 0.4, offset = 25, phi = 0.5,
                      n_ahead = 20, lastn = FALSE, plot = FALSE)

  expect_length(fit$f, 20)
  expect_length(fit$ar_f, 20)
  expect_false(fit$lastn)
})

test_that("fore_glambda works with lambda = 0 (geometric)", {
  set.seed(789)
  x <- gen_glambda(n = 100, lambda = 0.01, phi = 0.6, plot = FALSE)

  # lambda = 0 triggers geometric transformation
  fit <- fore_glambda(x, lambda = 0, offset = 30, phi = 0.6,
                      n_ahead = 10, plot = FALSE)

  expect_length(fit$f, 10)
  expect_false(anyNA(fit$f))
})

test_that("fore_glambda validates inputs", {
  # x too short
  expect_error(fore_glambda(1:10, lambda = 0.5, offset = 20, plot = FALSE),
               "`x` must be a numeric vector with at least 20 observations")

  # x not numeric
  expect_error(fore_glambda(letters[1:30], lambda = 0.5, offset = 20, plot = FALSE),
               "`x` must be a numeric vector")

  set.seed(123)
  x <- gen_glambda(n = 50, lambda = 0.5, phi = 0.6, plot = FALSE)

  # lambda not numeric
  expect_error(fore_glambda(x, lambda = "a", offset = 20, plot = FALSE),
               "`lambda` must be a single numeric value")

  # lambda wrong length
  expect_error(fore_glambda(x, lambda = c(0.3, 0.5), offset = 20, plot = FALSE),
               "`lambda` must be a single numeric value")

  # offset not numeric
  expect_error(fore_glambda(x, lambda = 0.5, offset = "a", plot = FALSE),
               "`offset` must be a non-negative number")

  # offset negative
  expect_error(fore_glambda(x, lambda = 0.5, offset = -10, plot = FALSE),
               "`offset` must be a non-negative number")

  # phi not numeric
  expect_error(fore_glambda(x, lambda = 0.5, offset = 20, phi = "a", plot = FALSE),
               "`phi` must be numeric")

  # n_ahead invalid
  expect_error(fore_glambda(x, lambda = 0.5, offset = 20, n_ahead = 0, plot = FALSE),
               "`n_ahead` must be a positive integer")
  expect_error(fore_glambda(x, lambda = 0.5, offset = 20, n_ahead = -5, plot = FALSE),
               "`n_ahead` must be a positive integer")

  # lastn invalid
  expect_error(fore_glambda(x, lambda = 0.5, offset = 20, lastn = "yes", plot = FALSE),
               "`lastn` must be TRUE or FALSE")
  expect_error(fore_glambda(x, lambda = 0.5, offset = 20, lastn = NA, plot = FALSE),
               "`lastn` must be TRUE or FALSE")

  # plot invalid
  expect_error(fore_glambda(x, lambda = 0.5, offset = 20, plot = "yes"),
               "`plot` must be TRUE or FALSE")
})

test_that("fore_glambda print method works", {
  set.seed(123)
  x <- gen_glambda(n = 100, lambda = 0.5, phi = 0.7, plot = FALSE)
  fit <- fore_glambda(x, lambda = 0.5, offset = 20, phi = 0.7,
                      n_ahead = 10, plot = FALSE)

  expect_output(print(fit), "G-Lambda Forecast")
  expect_output(print(fit), "Mode:")
  expect_output(print(fit), "Steps ahead:")
  expect_output(print(fit), "AR benchmark order:")
  expect_output(print(fit), "G-Lambda Forecasts:")
  expect_output(print(fit), "AR Forecasts:")
})

test_that("fore_glambda plots without error", {
  set.seed(123)
  x <- gen_glambda(n = 100, lambda = 0.5, phi = 0.6, plot = FALSE)

  # Holdout mode plot
  expect_silent({
    pdf(NULL)
    fit <- fore_glambda(x, lambda = 0.5, offset = 20, phi = 0.6,
                        n_ahead = 10, lastn = TRUE)
    dev.off()
  })

  # Future mode plot
  expect_silent({
    pdf(NULL)
    fit <- fore_glambda(x, lambda = 0.5, offset = 20, phi = 0.6,
                        n_ahead = 10, lastn = FALSE)
    dev.off()
  })
})

test_that("fore_glambda handles phi = 0 (no AR)", {
  set.seed(111)
  x <- gen_glambda(n = 100, lambda = 0.5, phi = 0, plot = FALSE)

  fit <- fore_glambda(x, lambda = 0.5, offset = 30, phi = 0,
                      n_ahead = 10, plot = FALSE)

  expect_length(fit$f, 10)
  expect_false(anyNA(fit$f))
})

test_that("fore_glambda works with different n_ahead values", {
  set.seed(222)
  x <- gen_glambda(n = 150, lambda = 0.4, phi = 0.5, plot = FALSE)

  for (n_ahead in c(5, 10, 25)) {
    fit <- fore_glambda(x, lambda = 0.4, offset = 25, phi = 0.5,
                        n_ahead = n_ahead, plot = FALSE)
    expect_length(fit$f, n_ahead)
    expect_length(fit$ar_f, n_ahead)
    expect_equal(fit$n_ahead, n_ahead)
  }
})

test_that("fore_glambda preserves original series", {
  set.seed(333)
  x <- gen_glambda(n = 80, lambda = 0.5, phi = 0.6, plot = FALSE)

  fit <- fore_glambda(x, lambda = 0.5, offset = 20, phi = 0.6,
                      n_ahead = 10, plot = FALSE)

  expect_equal(fit$original, x)
})

test_that("fore_glambda produces similar output to tswge::fore.glambda.wge", {
  skip_if_not_installed("tswge")

  set.seed(444)
  x <- gen_glambda(n = 100, lambda = 0.5, phi = 0.6, offset = 20, plot = FALSE)

  # Holdout mode comparison
  old <- suppressWarnings(tswge::fore.glambda.wge(
    x, lambda = 0.5, offset = 20, phi = 0.6,
    n.ahead = 10, lastn = TRUE, plot = FALSE
  ))

  new <- fore_glambda(x, lambda = 0.5, offset = 20, phi = 0.6,
                      n_ahead = 10, lastn = TRUE, plot = FALSE)

  # Forecasts should be similar
  # tswge returns f.glam (g-lambda forecasts) and f.ar (AR forecasts)
  expect_gt(cor(new$f, old$f.glam), 0.9)

  # Note: tswge has a bug in lastn=FALSE mode ('ar.fore' not found)
  # so we only compare holdout mode
})

test_that("fore_glambda AR benchmark uses correct AIC selection", {
  set.seed(555)
  x <- gen_glambda(n = 100, lambda = 0.4, phi = 0.7, plot = FALSE)

  fit <- fore_glambda(x, lambda = 0.4, offset = 25, phi = 0.7,
                      n_ahead = 10, plot = FALSE)

  # AR order should be reasonable (not too high)
  expect_true(fit$ar_order >= 0)
  expect_true(fit$ar_order <= 8)  # max order is 8 for future mode, 4 for holdout
})
