test_that("gen_garma returns correct length", {
  x <- gen_garma(n = 100, u = 0.5, lambda = 0.3, plot = FALSE, seed = 123)
  expect_length(x, 100)

  x <- gen_garma(n = 200, u = 0.8, lambda = 0.4, phi = 0.5, plot = FALSE, seed = 123)
  expect_length(x, 200)
})

test_that("gen_garma is reproducible with seed", {
  x1 <- gen_garma(n = 100, u = 0.5, lambda = 0.3, plot = FALSE, seed = 42)
  x2 <- gen_garma(n = 100, u = 0.5, lambda = 0.3, plot = FALSE, seed = 42)
  expect_equal(x1, x2)

  x3 <- gen_garma(n = 100, u = 0.5, lambda = 0.3, plot = FALSE, seed = 43)
  expect_false(identical(x1, x3))
})

test_that("gen_garma works with one Gegenbauer factor", {
  # Pure Gegenbauer (no AR/MA)
  x <- gen_garma(n = 100, u = 0.8, lambda = 0.4, plot = FALSE, seed = 123)
  expect_true(is.numeric(x))
  expect_false(anyNA(x))

  # With AR
  x <- gen_garma(n = 100, u = 0.5, lambda = 0.3, phi = 0.7, plot = FALSE, seed = 123)
  expect_true(is.numeric(x))
  expect_length(x, 100)

  # With MA
  x <- gen_garma(n = 100, u = 0.5, lambda = 0.3, theta = 0.4, plot = FALSE, seed = 123)
  expect_true(is.numeric(x))

  # With both AR and MA
  x <- gen_garma(n = 100, u = 0.5, lambda = 0.3, phi = 0.6, theta = 0.3,
                 plot = FALSE, seed = 123)
  expect_true(is.numeric(x))
})

test_that("gen_garma works with two Gegenbauer factors", {
  x <- gen_garma(n = 100, u = c(0.8, 0.3), lambda = c(0.3, 0.2),
                 plot = FALSE, seed = 123)
  expect_true(is.numeric(x))
  expect_length(x, 100)
  expect_false(anyNA(x))

  # With AR
  x <- gen_garma(n = 100, u = c(0.8, 0.3), lambda = c(0.3, 0.2),
                 phi = 0.5, plot = FALSE, seed = 123)
  expect_true(is.numeric(x))
})

test_that("gen_garma validates inputs", {
  expect_error(gen_garma(n = -1, u = 0.5, lambda = 0.3, plot = FALSE),
               "`n` must be a positive integer")
  expect_error(gen_garma(n = 100, u = "a", lambda = 0.3, plot = FALSE),
               "`u` must be numeric")
  expect_error(gen_garma(n = 100, u = 0.5, lambda = "a", plot = FALSE),
               "`lambda` must be numeric")
  expect_error(gen_garma(n = 100, u = c(0.5, 0.3), lambda = 0.3, plot = FALSE),
               "`u` and `lambda` must have the same length")
  expect_error(gen_garma(n = 100, u = c(0.5, 0.3, 0.1), lambda = c(0.3, 0.2, 0.1),
                         plot = FALSE),
               "Only 1 or 2 Gegenbauer factors are supported")
  expect_error(gen_garma(n = 100, u = 1.5, lambda = 0.3, plot = FALSE),
               "`u` values must be in")
  expect_error(gen_garma(n = 100, u = 0.5, lambda = 0.3, n_trunc = 0, plot = FALSE),
               "`n_trunc` must be a positive integer")
  expect_error(gen_garma(n = 100, u = 0.5, lambda = 0.3, burn_in = -1, plot = FALSE),
               "`burn_in` must be a non-negative integer")
  expect_error(gen_garma(n = 100, u = 0.5, lambda = 0.3, var_a = 0, plot = FALSE),
               "`var_a` must be a positive number")
  expect_error(gen_garma(n = 100, u = 0.5, lambda = 0.3, plot = "yes"),
               "`plot` must be TRUE or FALSE")
})

test_that("gen_garma plots without error", {
  expect_silent({
    pdf(NULL)
    x <- gen_garma(n = 50, u = 0.5, lambda = 0.3, seed = 123)
    dev.off()
  })
})

test_that("gen_garma produces similar output to tswge::gen.garma.wge", {
  skip_if_not_installed("tswge")

  # Due to subtle differences in MA coefficient computation between

  # tstse's macoef_geg (uses convolve_truncated_cpp) and tswge's version,
  # there are small numerical differences. We test correlation instead
  # of exact equality.

  # Pure Gegenbauer
  old <- suppressWarnings(tswge::gen.garma.wge(n = 200, u = 0.5, lambda = 0.3,
                                                phi = 0, theta = 0,
                                                trun = 1000, burn_in = 600,
                                                vara = 1, plot = FALSE, sn = 123))
  new <- gen_garma(n = 200, u = 0.5, lambda = 0.3, n_trunc = 1000,
                   burn_in = 600, var_a = 1, plot = FALSE, seed = 123)

  # High correlation expected (should be nearly identical)
  expect_gt(cor(new, old), 0.99)

  # Similar variance
  expect_equal(var(new), var(old), tolerance = 0.1)

  # With AR
  old <- suppressWarnings(tswge::gen.garma.wge(n = 200, u = 0.8, lambda = 0.4,
                                                phi = 0.7, theta = 0,
                                                trun = 1000, burn_in = 600,
                                                vara = 1, plot = FALSE, sn = 456))
  new <- gen_garma(n = 200, u = 0.8, lambda = 0.4, phi = 0.7,
                   n_trunc = 1000, burn_in = 600, var_a = 1,
                   plot = FALSE, seed = 456)

  expect_gt(cor(new, old), 0.99)
})
