# tests/testthat/test-gen_garch.R

library(tswge)

test_that("gen_garch matches original gen.garch.wge with seed", {

  # Test case 1: GARCH(1,1) - equal orders
  old1 <- tswge::gen.garch.wge(n = 200, alpha0 = 0.1, alpha = 0.3, beta = 0.5,
                               plot = FALSE, sn = 123)
  new1 <- gen_garch(n = 200, alpha0 = 0.1, alpha = 0.3, beta = 0.5,
                    plot = FALSE, seed = 123)

  expect_equal(new1, as.double(old1), tolerance = 1e-10)

  # Test case 2: GARCH(2,2) - equal orders
  old2 <- tswge::gen.garch.wge(n = 200, alpha0 = 0.1,
                               alpha = c(0.2, 0.1), beta = c(0.3, 0.2),
                               plot = FALSE, sn = 456)
  new2 <- gen_garch(n = 200, alpha0 = 0.1,
                    alpha = c(0.2, 0.1), beta = c(0.3, 0.2),
                    plot = FALSE, seed = 456)

  expect_equal(new2, as.double(old2), tolerance = 1e-10)

  # Test case 3: GARCH(1,2) - q0 > p0 (works in original)
  old3 <- tswge::gen.garch.wge(n = 200, alpha0 = 0.1,
                               alpha = c(0.2, 0.15), beta = 0.5,
                               plot = FALSE, sn = 789)
  new3 <- gen_garch(n = 200, alpha0 = 0.1,
                    alpha = c(0.2, 0.15), beta = 0.5,
                    plot = FALSE, seed = 789)

  expect_equal(new3, as.double(old3), tolerance = 1e-10)
})


test_that("gen_garch handles p0 > q0 (fixed bug from original)", {

  # This case crashes in original tswge::gen.garch.wge
  # Our version handles it correctly
  expect_no_error(
    gen_garch(n = 100, alpha0 = 0.1, alpha = 0.3, beta = c(0.3, 0.2),
              plot = FALSE, seed = 111)
  )
})


test_that("gen_garch returns correct length", {

  result <- gen_garch(n = 100, alpha0 = 0.1, alpha = 0.3, beta = 0.5,
                      plot = FALSE, seed = 1)
  expect_equal(length(result), 100)
})


test_that("gen_garch validates inputs", {

  expect_error(gen_garch(n = -1, alpha0 = 0.1, alpha = 0.3, beta = 0.5),
               "`n` must be")
  expect_error(gen_garch(n = 100, alpha0 = -0.1, alpha = 0.3, beta = 0.5),
               "`alpha0` must be")
  expect_error(gen_garch(n = 100, alpha0 = 0.1, alpha = numeric(0), beta = 0.5),
               "`alpha` must be")
  expect_error(gen_garch(n = 100, alpha0 = 0.1, alpha = 0.3, beta = numeric(0)),
               "`beta` must be")
})


test_that("gen_garch is reproducible with seed", {

  x1 <- gen_garch(n = 50, alpha0 = 0.1, alpha = 0.3, beta = 0.5,
                  plot = FALSE, seed = 999)
  x2 <- gen_garch(n = 50, alpha0 = 0.1, alpha = 0.3, beta = 0.5,
                  plot = FALSE, seed = 999)

  expect_identical(x1, x2)
})
