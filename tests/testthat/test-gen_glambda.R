test_that("gen_glambda returns correct length", {
  x <- gen_glambda(n = 100, lambda = 0, plot = FALSE, seed = 123)
  expect_length(x, 100)

  x <- gen_glambda(n = 200, lambda = 1, plot = FALSE, seed = 123)
  expect_length(x, 200)

  x <- gen_glambda(n = 50, lambda = 2, phi = 0.5, plot = FALSE, seed = 123)
  expect_length(x, 50)
})

test_that("gen_glambda is reproducible with seed", {
  x1 <- gen_glambda(n = 100, lambda = 0, phi = 0.7, plot = FALSE, seed = 42)
  x2 <- gen_glambda(n = 100, lambda = 0, phi = 0.7, plot = FALSE, seed = 42)
  expect_equal(x1, x2)

  x3 <- gen_glambda(n = 100, lambda = 0, phi = 0.7, plot = FALSE, seed = 43)
  expect_false(identical(x1, x3))
})

test_that("gen_glambda works with different lambda values", {
  # lambda = 0 (logarithmic)
  x0 <- gen_glambda(n = 100, lambda = 0, plot = FALSE, seed = 123)
  expect_true(is.numeric(x0))
  expect_false(anyNA(x0))

  # lambda = 0.5
  x05 <- gen_glambda(n = 100, lambda = 0.5, plot = FALSE, seed = 123)
  expect_true(is.numeric(x05))
  expect_false(anyNA(x05))

  # lambda = 1 (linear)
  x1 <- gen_glambda(n = 100, lambda = 1, plot = FALSE, seed = 123)
  expect_true(is.numeric(x1))
  expect_false(anyNA(x1))

  # lambda = 2 (quadratic)
  x2 <- gen_glambda(n = 100, lambda = 2, plot = FALSE, seed = 123)
  expect_true(is.numeric(x2))
  expect_false(anyNA(x2))
})

test_that("gen_glambda works with AR coefficients", {
  # AR(1)
  x <- gen_glambda(n = 100, lambda = 1, phi = 0.9, plot = FALSE, seed = 123)
  expect_true(is.numeric(x))
  expect_length(x, 100)

  # AR(2)
  x <- gen_glambda(n = 100, lambda = 0, phi = c(1.5, -0.75), plot = FALSE, seed = 123)
  expect_true(is.numeric(x))
  expect_length(x, 100)
})

test_that("gen_glambda respects sigma2 parameter", {
  set.seed(999)
  x1 <- gen_glambda(n = 500, lambda = 0, sigma2 = 1, plot = FALSE, seed = 123)
  x2 <- gen_glambda(n = 500, lambda = 0, sigma2 = 4, plot = FALSE, seed = 123)

  # Higher variance should give larger spread
  expect_gt(var(x2), var(x1))
})

test_that("gen_glambda validates inputs", {
  expect_error(gen_glambda(n = -1, lambda = 0, plot = FALSE),
               "`n` must be a positive integer")
  expect_error(gen_glambda(n = 100, lambda = "a", plot = FALSE),
               "`lambda` must be a single numeric value")
  expect_error(gen_glambda(n = 100, lambda = 0, offset = -5, plot = FALSE),
               "`offset` must be a positive number")
  expect_error(gen_glambda(n = 100, lambda = 0, sigma2 = 0, plot = FALSE),
               "`sigma2` must be a positive number")
  expect_error(gen_glambda(n = 100, lambda = 0, plot = "yes"),
               "`plot` must be TRUE or FALSE")
})

test_that("gen_glambda plots without error", {
  expect_silent({
    pdf(NULL)
    x <- gen_glambda(n = 50, lambda = 0, seed = 123)
    dev.off()
  })
})

test_that("gen_glambda matches tswge::gen.glambda.wge", {
  skip_if_not_installed("tswge")

  # Test with lambda = 0
  set.seed(123)
  old <- suppressWarnings(tswge::gen.glambda.wge(n = 100, lambda = 0, phi = 0.7,
                                                  offset = 20, vara = 1,
                                                  plot = FALSE, sn = 123))
  new <- gen_glambda(n = 100, lambda = 0, phi = 0.7, offset = 20, sigma2 = 1,
                     plot = FALSE, seed = 123)
  expect_equal(new, old, tolerance = 1e-10)

  # Test with lambda = 1
  set.seed(456)
  old <- suppressWarnings(tswge::gen.glambda.wge(n = 100, lambda = 1, phi = 0.5,
                                                  offset = 20, vara = 1,
                                                  plot = FALSE, sn = 456))
  new <- gen_glambda(n = 100, lambda = 1, phi = 0.5, offset = 20, sigma2 = 1,
                     plot = FALSE, seed = 456)
  expect_equal(new, old, tolerance = 1e-10)
})
