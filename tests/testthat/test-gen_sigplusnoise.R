# tests/testthat/test-gen_sigplusnoise.R

library(tswge)

test_that("gen_sigplusnoise matches original gen.sigplusnoise.wge", {

  # Test case 1: Linear trend plus white noise
  old1 <- tswge::gen.sigplusnoise.wge(n = 100, b0 = 10, b1 = 0.05,
                                      plot = FALSE, sn = 123)
  new1 <- gen_sigplusnoise(n = 100, b0 = 10, b1 = 0.05,
                           plot = FALSE, seed = 123)

  expect_equal(new1, as.double(old1), tolerance = 1e-10)

  # Test case 2: Single cosine plus AR(1) noise
  old2 <- tswge::gen.sigplusnoise.wge(n = 100, coef = c(5, 0), freq = c(0.1, 0),
                                      phi = 0.7, plot = FALSE, sn = 456)
  new2 <- gen_sigplusnoise(n = 100, coef = c(5, 0), freq = c(0.1, 0),
                           phi = 0.7, plot = FALSE, seed = 456)

  expect_equal(new2, as.double(old2), tolerance = 1e-10)

  # Test case 3: Two cosines with phase shifts
  old3 <- tswge::gen.sigplusnoise.wge(n = 100, coef = c(3, 2), freq = c(0.1, 0.25),
                                      psi = c(0, 1), plot = FALSE, sn = 789)
  new3 <- gen_sigplusnoise(n = 100, coef = c(3, 2), freq = c(0.1, 0.25),
                           psi = c(0, 1), plot = FALSE, seed = 789)

  expect_equal(new3, as.double(old3), tolerance = 1e-10)

  # Test case 4: Full model
  old4 <- tswge::gen.sigplusnoise.wge(n = 150, b0 = 5, b1 = 0.02,
                                      coef = c(2, 1), freq = c(0.05, 0.2),
                                      psi = c(0.5, 1.5), phi = c(0.5, -0.3),
                                      vara = 2, plot = FALSE, sn = 111)
  new4 <- gen_sigplusnoise(n = 150, b0 = 5, b1 = 0.02,
                           coef = c(2, 1), freq = c(0.05, 0.2),
                           psi = c(0.5, 1.5), phi = c(0.5, -0.3),
                           vara = 2, plot = FALSE, seed = 111)

  expect_equal(new4, as.double(old4), tolerance = 1e-10)
})


test_that("gen_sigplusnoise validates inputs", {

  expect_error(gen_sigplusnoise(n = -1), "`n` must be")
  expect_error(gen_sigplusnoise(n = 100, b0 = "bad"), "`b0` must be")
  expect_error(gen_sigplusnoise(n = 100, vara = -1), "`vara` must be")
  expect_error(gen_sigplusnoise(n = 100, coef = c(1, 2), freq = c(0.1)),
               "same length")
})


test_that("gen_sigplusnoise returns correct length", {

  result <- gen_sigplusnoise(n = 150, plot = FALSE, seed = 1)
  expect_equal(length(result), 150)
})


test_that("gen_sigplusnoise is reproducible with seed", {

  x1 <- gen_sigplusnoise(n = 50, b1 = 0.1, phi = 0.5, plot = FALSE, seed = 999)
  x2 <- gen_sigplusnoise(n = 50, b1 = 0.1, phi = 0.5, plot = FALSE, seed = 999)

  expect_identical(x1, x2)
})


test_that("gen_sigplusnoise handles arbitrary number of cosine terms", {

  # More than 2 cosine terms (improvement over original)
  result <- gen_sigplusnoise(n = 100,
                             coef = c(1, 2, 3),
                             freq = c(0.1, 0.2, 0.3),
                             psi = c(0, 0, 0),
                             plot = FALSE, seed = 123)

  expect_equal(length(result), 100)
})
