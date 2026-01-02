library(tswge)

test_that("true_spec matches original true.arma.spec.wge", {
  # Note: We only compare spec values, not freq.

  # Original has a bug where f is a scalar (0.5) when plot=FALSE.
  # Our implementation correctly returns the full frequency vector.

  # Test case 1: AR(1)
  old1 <- tswge::true.arma.spec.wge(phi = 0.8, plot = FALSE)
  new1 <- true_spec(phi = 0.8, plot = FALSE)
  expect_equal(new1$spec, old1$spec, tolerance = 1e-6)

  # Test case 2: AR(2)
  old2 <- tswge::true.arma.spec.wge(phi = c(1.5, -0.75), plot = FALSE)
  new2 <- true_spec(phi = c(1.5, -0.75), plot = FALSE)
  expect_equal(new2$spec, old2$spec, tolerance = 1e-6)

  # Test case 3: MA(1)
  old3 <- tswge::true.arma.spec.wge(theta = 0.6, plot = FALSE)
  new3 <- true_spec(theta = 0.6, plot = FALSE)
  expect_equal(new3$spec, old3$spec, tolerance = 1e-6)

  # Test case 4: ARMA(1,1)
  old4 <- tswge::true.arma.spec.wge(phi = 0.7, theta = 0.4, plot = FALSE)
  new4 <- true_spec(phi = 0.7, theta = 0.4, plot = FALSE)
  expect_equal(new4$spec, old4$spec, tolerance = 1e-6)

  # Test case 5: ARMA(2,1) with different variance
  old5 <- tswge::true.arma.spec.wge(phi = c(1.2, -0.6), theta = 0.3, vara = 2, plot = FALSE)
  new5 <- true_spec(phi = c(1.2, -0.6), theta = 0.3, vara = 2, plot = FALSE)
  expect_equal(new5$spec, old5$spec, tolerance = 1e-6)
})


test_that("true_spec handles white noise", {

  result <- true_spec(vara = 1, plot = FALSE)

  # White noise: flat spectrum at 0 dB
  expect_true(all(abs(result$spec) < 1e-10))
})


test_that("true_spec validates inputs", {

  expect_error(true_spec(phi = "bad"), "`phi` must be")
  expect_error(true_spec(theta = "bad"), "`theta` must be")
  expect_error(true_spec(vara = -1), "`vara` must be")
  expect_error(true_spec(n_freq = 1), "`n_freq` must be")
})


test_that("true_spec returns correct structure", {

  result <- true_spec(phi = 0.5, plot = FALSE)

  expect_true(is.list(result))
  expect_named(result, c("freq", "spec"))
  expect_equal(length(result$freq), 251)
  expect_equal(length(result$spec), 251)
})


test_that("true_spec freq range is correct", {

  result <- true_spec(phi = 0.5, plot = FALSE)

  expect_equal(result$freq[1], 0)
  expect_equal(result$freq[251], 0.5)
})


test_that("true_spec custom n_freq works", {

  result <- true_spec(phi = 0.5, n_freq = 101, plot = FALSE)

  expect_equal(length(result$freq), 101)
  expect_equal(length(result$spec), 101)
  expect_equal(result$freq[1], 0)
  expect_equal(result$freq[101], 0.2)  # (101-1)/500 = 0.2
})


test_that("true_spec plot works", {

  expect_no_error(true_spec(phi = 0.8, plot = TRUE))
})
