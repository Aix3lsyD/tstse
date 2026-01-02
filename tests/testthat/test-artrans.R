# tests/testthat/test-artrans.R

library(tswge)

test_that("artrans matches original artrans.wge output", {

  # Test case 1: simple first difference
  x1 <- sunspot.year
  phi1 <- 1

  old1 <- tswge::artrans.wge(x1, phi1, plottr = FALSE)
  new1 <- tstse::artrans(x1, phi1, plot = FALSE)

  expect_equal(new1, as.double(old1), tolerance = 1e-10)

  # Test case 2: second order filter
  x2 <- lh
  phi2 <- c(1.6, -1)

  old2 <- tswge::artrans.wge(x2, phi2, plottr = FALSE)
  new2 <- tstse::artrans(x2, phi2, plot = FALSE)

  expect_equal(new2, as.double(old2), tolerance = 1e-10)

  # Test case 3: higher order with random data
  set.seed(42)
  x3 <- cumsum(rnorm(200))
  phi3 <- c(0.5, 0.3, -0.2)

  old3 <- tswge::artrans.wge(x3, phi3, plottr = FALSE)
  new3 <- tstse::artrans(x3, phi3, plot = FALSE)

  expect_equal(new3, as.double(old3), tolerance = 1e-10)
})


test_that("artrans handles edge cases", {

  set.seed(123)
  x <- rnorm(50)

  # Single coefficient
  expect_no_error(artrans(x, phi = 1, plot = FALSE))

  # Output length is correct
  result <- artrans(x, phi = c(1, 1), plot = FALSE)
  expect_equal(length(result), length(x) - 2)

  # lag_max auto-adjusts when too large
  result2 <- artrans(x, phi = c(1, 1), lag_max = 1000, plot = FALSE)
  expect_equal(length(result2), length(x) - 2)
})


test_that("artrans validates inputs", {

  expect_error(artrans("not numeric", phi = 1), "`x` must be")
  expect_error(artrans(numeric(0), phi = 1), "`x` must be")
  expect_error(artrans(1:10, phi = numeric(0)), "`phi` must be")
  expect_error(artrans(1:10, phi = 1, lag_max = -5), "`lag_max` must be")
})
