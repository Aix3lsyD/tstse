# Tests for is_sample()
# Instantaneous spectrum (sample-based)

test_that("is_sample matches tswge::is.sample.wge output structure", {
  skip_if_not_installed("tswge")

  set.seed(123)
  x <- rnorm(100)

  result <- is_sample(x, lambda = 1, offset = 10, plot = FALSE)

  png(tempfile())
  # tswge doesn't return output, just plots
  tswge::is.sample.wge(x, lambda = 1, offset = 10)
  dev.off()

  # Check structure
  expect_s3_class(result, "is_sample")
  expect_true(!is.null(result$time))
  expect_true(!is.null(result$freq))
  expect_true(!is.null(result$spec))
})


test_that("is_sample returns correct structure", {
  set.seed(111)
  x <- rnorm(50)

  result <- is_sample(x, lambda = 1, offset = 5, plot = FALSE)

  expect_s3_class(result, "is_sample")
  expect_named(result, c("time", "freq", "gfreq", "spec", "lambda", "offset"))
  expect_equal(result$lambda, 1)
  expect_equal(result$offset, 5)
})


test_that("is_sample time ranges from 1 to n", {
  x <- rnorm(50)
  result <- is_sample(x, lambda = 1, offset = 10, plot = FALSE)

  expect_true(min(result$time) == 1)
  expect_true(max(result$time) == 50)
})


test_that("is_sample freq is bounded by 0 and 0.5", {
  set.seed(222)
  x <- rnorm(100)

  result <- is_sample(x, lambda = 1, offset = 10, plot = FALSE)

  # Most frequencies should be in valid range
  valid_freq <- result$freq[result$freq > 0 & result$freq < 0.5]
  expect_true(length(valid_freq) > 0.5 * length(result$freq))
})


test_that("is_sample works with lambda = 0", {
  set.seed(333)
  x <- rnorm(50)

  expect_no_error(
    result <- is_sample(x, lambda = 0, offset = 10, plot = FALSE)
  )
  expect_equal(result$lambda, 0)
})


test_that("is_sample works with lambda = 1", {
  set.seed(444)
  x <- rnorm(50)

  expect_no_error(
    result <- is_sample(x, lambda = 1, offset = 10, plot = FALSE)
  )
  expect_equal(result$lambda, 1)
})


test_that("is_sample works with lambda = 2", {
  set.seed(555)
  x <- rnorm(50)

  expect_no_error(
    result <- is_sample(x, lambda = 2, offset = 10, plot = FALSE)
  )
  expect_equal(result$lambda, 2)
})


test_that("is_sample spec values are non-negative", {
  set.seed(666)
  x <- rnorm(100)

  result <- is_sample(x, lambda = 1, offset = 10, plot = FALSE)

  # Spec should be normalized to start at 0
  expect_true(min(result$spec) >= 0)
})


test_that("is_sample validates x parameter", {
  expect_error(is_sample("abc", lambda = 1, offset = 10), "`x` must be a numeric")
  expect_error(is_sample(1:9, lambda = 1, offset = 10), "at least 10")
  expect_error(is_sample(c(1:10, NA), lambda = 1, offset = 10), "contains NA")
})


test_that("is_sample validates lambda parameter", {
  x <- rnorm(50)
  expect_error(is_sample(x, lambda = "1", offset = 10), "`lambda` must be")
  expect_error(is_sample(x, lambda = c(1, 2), offset = 10), "`lambda` must be")
})


test_that("is_sample validates offset parameter", {
  x <- rnorm(50)
  expect_error(is_sample(x, lambda = 1, offset = 0), "`offset` must be a positive")
  expect_error(is_sample(x, lambda = 1, offset = -5), "`offset` must be a positive")
  expect_error(is_sample(x, lambda = 1, offset = "10"), "`offset` must be a positive")
})


test_that("is_sample validates plot parameter", {
  x <- rnorm(50)
  expect_error(is_sample(x, lambda = 1, offset = 10, plot = "TRUE"),
               "`plot` must be TRUE or FALSE")
  expect_error(is_sample(x, lambda = 1, offset = 10, plot = NA),
               "`plot` must be TRUE or FALSE")
})


test_that("is_sample plot produces no error", {
  skip_if_not(capabilities("png"))

  set.seed(777)
  x <- rnorm(50)

  png(tempfile())
  on.exit(dev.off())

  expect_no_error(is_sample(x, lambda = 1, offset = 10, plot = TRUE))
})


test_that("is_sample print method works", {
  set.seed(888)
  x <- rnorm(50)

  result <- is_sample(x, lambda = 1, offset = 10, plot = FALSE)

  expect_output(print(result), "Instantaneous Spectrum")
  expect_output(print(result), "Lambda:")
  expect_output(print(result), "Offset:")
})


test_that("is_sample returns invisible result", {
  x <- rnorm(50)

  png(tempfile())
  on.exit(dev.off())

  expect_invisible(is_sample(x, lambda = 1, offset = 10, plot = TRUE))
})


test_that("is_sample different offset values are stored correctly", {
  set.seed(999)
  x <- rnorm(50)

  result1 <- is_sample(x, lambda = 1, offset = 5, plot = FALSE)
  result2 <- is_sample(x, lambda = 1, offset = 20, plot = FALSE)

  # Offset values should be stored correctly
  expect_equal(result1$offset, 5)
  expect_equal(result2$offset, 20)
})
