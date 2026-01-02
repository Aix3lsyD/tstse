library(tswge)


test_that("ma_smooth returns correct structure", {

  set.seed(123)
  x <- rnorm(100)

  result <- ma_smooth(x, order = 5, plot = FALSE)

  expect_type(result, "list")
  expect_named(result, c("x", "smooth", "order"))
  expect_equal(result$x, x)
  expect_length(result$smooth, 100)
  expect_equal(result$order, 5)
})


test_that("ma_smooth matches tswge for odd order", {

  set.seed(456)
  x <- cumsum(rnorm(80))

  # Suppress cat() output from tswge
  for (ord in c(3, 5, 7, 11)) {
    old <- capture.output(
      old_result <- tswge::ma.smooth.wge(x, order = ord, plot = FALSE)
    )
    new <- ma_smooth(x, order = ord, plot = FALSE)

    expect_equal(new$smooth, old_result$smooth, tolerance = 1e-10,
                 label = paste("order =", ord))
  }
})


test_that("ma_smooth matches tswge for even order", {

  set.seed(789)
  x <- cumsum(rnorm(80))

  # Even orders use weighted endpoints
  for (ord in c(2, 4, 6, 12)) {
    old <- capture.output(
      old_result <- tswge::ma.smooth.wge(x, order = ord, plot = FALSE)
    )
    new <- ma_smooth(x, order = ord, plot = FALSE)

    expect_equal(new$smooth, old_result$smooth, tolerance = 1e-10,
                 label = paste("order =", ord))
  }
})


test_that("ma_smooth has correct NA pattern for odd order", {

  set.seed(111)
  x <- rnorm(50)

  # For odd order k, first k2 and last k2 values are NA
  # where k2 = floor(k/2)
  result <- ma_smooth(x, order = 7, plot = FALSE)

  k2 <- 3  # floor(7/2)

  # First k2 values should be NA
 expect_true(all(is.na(result$smooth[1:k2])))

  # Last k2 values should be NA
  expect_true(all(is.na(result$smooth[(50 - k2 + 1):50])))

  # Middle values should not be NA
  expect_false(any(is.na(result$smooth[(k2 + 1):(50 - k2)])))
})


test_that("ma_smooth has correct NA pattern for even order", {

  set.seed(222)
  x <- rnorm(50)

  # For even order k, first k2-1 and last k2 values are NA
  result <- ma_smooth(x, order = 6, plot = FALSE)

  k2 <- 3  # 6/2

  # First k2-1 values should be NA (indices 1, 2)
  expect_true(all(is.na(result$smooth[1:(k2 - 1)])))

  # Last k2 values should be NA
  expect_true(all(is.na(result$smooth[(50 - k2 + 1):50])))

  # Check that t_start value is not NA
  t_start <- 6 - 3 + 1  # = 4
  expect_false(is.na(result$smooth[t_start]))
})


test_that("ma_smooth center value is correct for odd order", {

  # Simple test case
  x <- c(1, 2, 3, 4, 5, 6, 7)

  result <- ma_smooth(x, order = 3, plot = FALSE)

  # For order 3, first smoothed value at position 2
  # should be mean(1, 2, 3) = 2
  expect_equal(result$smooth[2], 2, tolerance = 1e-10)

  # Position 3 should be mean(2, 3, 4) = 3
  expect_equal(result$smooth[3], 3, tolerance = 1e-10)

  # Position 4 should be mean(3, 4, 5) = 4
  expect_equal(result$smooth[4], 4, tolerance = 1e-10)
})


test_that("ma_smooth validates input", {

  x <- rnorm(100)

  expect_error(ma_smooth("not numeric", order = 3),
               "`x` must be a non-empty numeric vector")

  expect_error(ma_smooth(numeric(0), order = 3),
               "`x` must be a non-empty numeric vector")

  expect_error(ma_smooth(x, order = -1),
               "`order` must be a positive integer")

  expect_error(ma_smooth(x, order = 0),
               "`order` must be a positive integer")
})


test_that("ma_smooth validates order vs length", {

  x <- rnorm(10)

  expect_error(ma_smooth(x, order = 11),
               "`order` cannot exceed the length of `x`")

  expect_error(ma_smooth(x, order = 20),
               "`order` cannot exceed the length of `x`")
})


test_that("ma_smooth reduces variance", {

  set.seed(333)
  x <- rnorm(200, sd = 2)

  result <- ma_smooth(x, order = 15, plot = FALSE)

  # Remove NAs for comparison
  smooth_valid <- result$smooth[!is.na(result$smooth)]

  # Smoothed series should have lower variance
  expect_lt(var(smooth_valid), var(x))
})


test_that("ma_smooth plotting works", {

  set.seed(444)
  x <- rnorm(50)

  expect_no_error(ma_smooth(x, order = 5, plot = TRUE))
})


test_that("ma_smooth with order = 1 returns original", {

  set.seed(555)
  x <- rnorm(30)

  result <- ma_smooth(x, order = 1, plot = FALSE)

  # With order 1, smoothed should equal original (no NAs)
  expect_equal(result$smooth, x, tolerance = 1e-10)
})


test_that("ma_smooth works with AirPassengers", {

  x <- as.numeric(AirPassengers)

  # 2x12 MA commonly used for monthly data
  old <- capture.output(
    old_result <- tswge::ma.smooth.wge(x, order = 12, plot = FALSE)
  )
  new <- ma_smooth(x, order = 12, plot = FALSE)

  expect_equal(new$smooth, old_result$smooth, tolerance = 1e-10)
})


test_that("ma_smooth preserves trend", {

  # Linear trend plus noise
  set.seed(666)
  t <- 1:100
  x <- 2 * t + rnorm(100, sd = 5)

  result <- ma_smooth(x, order = 11, plot = FALSE)

  # Smoothed values should be close to the true trend
  valid_idx <- which(!is.na(result$smooth))
  correlation <- cor(result$smooth[valid_idx], 2 * t[valid_idx])

  expect_gt(correlation, 0.99)
})
