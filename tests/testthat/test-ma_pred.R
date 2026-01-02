library(tswge)


test_that("ma_pred returns correct structure", {

  set.seed(123)
  x <- rnorm(100)

  result <- ma_pred(x, order = 5, n_ahead = 3, plot = FALSE)

  expect_type(result, "list")
  expect_named(result, c("x_original", "pred", "order"))
  expect_equal(result$x_original, x)
  expect_length(result$pred, 103)  # n + n_ahead
  expect_equal(result$order, 5)
})


test_that("ma_pred matches tswge for basic case", {

  set.seed(456)
  x <- cumsum(rnorm(50))

  # Suppress cat() output from tswge
  old <- capture.output(
    old_result <- tswge::ma.pred.wge(x, order = 3, n.ahead = 1, plot = FALSE)
  )
  new <- ma_pred(x, order = 3, n_ahead = 1, plot = FALSE)

  expect_equal(new$pred, old_result$pred, tolerance = 1e-10)
  expect_equal(new$order, old_result$order)
})


test_that("ma_pred matches tswge for multi-step forecast", {

  set.seed(789)
  x <- cumsum(rnorm(80))

  old <- capture.output(
    old_result <- tswge::ma.pred.wge(x, order = 5, n.ahead = 10, plot = FALSE)
  )
  new <- ma_pred(x, order = 5, n_ahead = 10, plot = FALSE)

  expect_equal(new$pred, old_result$pred, tolerance = 1e-10)
})


test_that("ma_pred matches tswge with different orders", {

  set.seed(111)
  x <- rnorm(60)

  for (ord in c(2, 4, 7, 10)) {
    old <- capture.output(
      old_result <- tswge::ma.pred.wge(x, order = ord, n.ahead = 5, plot = FALSE)
    )
    new <- ma_pred(x, order = ord, n_ahead = 5, plot = FALSE)

    expect_equal(new$pred, old_result$pred, tolerance = 1e-10,
                 label = paste("order =", ord))
  }
})


test_that("ma_pred first k values are NA", {

  set.seed(222)
  x <- rnorm(50)

  result <- ma_pred(x, order = 7, n_ahead = 1, plot = FALSE)

  # First 'order' values should be NA
  expect_true(all(is.na(result$pred[1:7])))

  # Values from order+1 onward should not be NA
  expect_false(any(is.na(result$pred[8:51])))
})


test_that("ma_pred forecast equals mean of last k values", {

  set.seed(333)
  x <- rnorm(30)

  result <- ma_pred(x, order = 5, n_ahead = 1, plot = FALSE)

  # First forecast should be mean of last 5 observations
  expected_forecast <- mean(x[26:30])
  expect_equal(result$pred[31], expected_forecast, tolerance = 1e-10)
})


test_that("ma_pred recursive forecast works correctly", {

  set.seed(444)
  x <- 1:20  # Simple sequence

  result <- ma_pred(x, order = 3, n_ahead = 3, plot = FALSE)

  # Manual calculation:
  # pred[21] = mean(x[18:20]) = mean(18,19,20) = 19
  # pred[22] = mean(x[19:20], pred[21]) = mean(19,20,19) = 19.333...
  # pred[23] = mean(x[20], pred[21:22]) = mean(20,19,19.333) = 19.444...

  expect_equal(result$pred[21], 19, tolerance = 1e-10)
  expect_equal(result$pred[22], mean(c(19, 20, 19)), tolerance = 1e-10)
})


test_that("ma_pred validates input", {

  x <- rnorm(100)

  expect_error(ma_pred("not numeric", order = 3, n_ahead = 1),
               "`x` must be a non-empty numeric vector")

  expect_error(ma_pred(numeric(0), order = 3, n_ahead = 1),
               "`x` must be a non-empty numeric vector")

  expect_error(ma_pred(x, order = -1, n_ahead = 1),
               "`order` must be a positive integer")

  expect_error(ma_pred(x, order = 0, n_ahead = 1),
               "`order` must be a positive integer")

  expect_error(ma_pred(x, order = 3, n_ahead = 0),
               "`n_ahead` must be a positive integer")

  expect_error(ma_pred(x, order = 3, n_ahead = -1),
               "`n_ahead` must be a positive integer")
})


test_that("ma_pred validates order vs length", {

  x <- rnorm(10)

  expect_error(ma_pred(x, order = 11, n_ahead = 1),
               "`order` cannot exceed the length of `x`")

  expect_error(ma_pred(x, order = 15, n_ahead = 1),
               "`order` cannot exceed the length of `x`")
})


test_that("ma_pred smoothing reduces variance", {

  set.seed(555)
  x <- rnorm(100, sd = 2)

  result <- ma_pred(x, order = 10, n_ahead = 1, plot = FALSE)

  # Remove NAs for comparison
  pred_valid <- result$pred[!is.na(result$pred)]

  # Smoothed series should have lower variance
  expect_lt(var(pred_valid), var(x))
})


test_that("ma_pred plotting works", {

  set.seed(666)
  x <- rnorm(50)

  expect_no_error(ma_pred(x, order = 5, n_ahead = 5, plot = TRUE))
})


test_that("ma_pred with order = 1 returns original values", {

  set.seed(777)
  x <- rnorm(30)

  result <- ma_pred(x, order = 1, n_ahead = 1, plot = FALSE)

  # With order 1, smoothed values should equal original (shifted by 1)
  expect_equal(result$pred[2:31], x, tolerance = 1e-10)
})


test_that("ma_pred works with AirPassengers data", {

  x <- as.numeric(AirPassengers)

  old <- capture.output(
    old_result <- tswge::ma.pred.wge(x, order = 12, n.ahead = 12, plot = FALSE)
  )
  new <- ma_pred(x, order = 12, n_ahead = 12, plot = FALSE)

  expect_equal(new$pred, old_result$pred, tolerance = 1e-10)
})
