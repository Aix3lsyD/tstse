library(tswge)

test_that("slr returns correct structure", {

  set.seed(123)
  x <- rnorm(100)

  result <- slr(x)

  expect_s3_class(result, "slr")
  expect_named(result, c("res", "b0hat", "b1hat", "pvalue", "tstatistic"))
  expect_length(result$res, 100)
  expect_true(is.numeric(result$b0hat))
  expect_true(is.numeric(result$b1hat))
  expect_true(is.numeric(result$pvalue))
  expect_true(is.numeric(result$tstatistic))
})


test_that("slr matches tswge slr.wge", {

  set.seed(456)
  x <- cumsum(rnorm(100))  # Random walk with trend

  old <- tswge::slr.wge(x)
  new <- slr(x)

  expect_equal(new$res, old$res)
  expect_equal(new$b0hat, unname(old$b0hat))
  expect_equal(new$b1hat, unname(old$b1hat))
  expect_equal(new$pvalue, old$pvalue)
  expect_equal(new$tstatistic, old$tstatistic)
})


test_that("slr detects linear trend", {

  # Create series with strong linear trend
  set.seed(789)
  x <- 1:100 * 2 + rnorm(100, sd = 0.5)

  result <- slr(x)

  # Slope should be close to 2
expect_equal(result$b1hat, 2, tolerance = 0.1)

  # P-value should be very small (significant trend)
  expect_lt(result$pvalue, 0.001)
})


test_that("slr residuals have zero mean trend", {

  set.seed(111)
  x <- 1:50 * 0.5 + rnorm(50)

  result <- slr(x)

  # Residuals should have approximately zero trend
  resid_trend <- lm(result$res ~ seq_along(result$res))
  expect_equal(unname(coef(resid_trend)[2]), 0, tolerance = 1e-10)
})


test_that("slr handles no trend case", {

  set.seed(222)
  x <- rnorm(100, mean = 5)  # No trend, just noise around mean

  result <- slr(x)

  # Slope should be near zero
  expect_equal(result$b1hat, 0, tolerance = 0.1)

  # P-value should be large (no significant trend)
  expect_gt(result$pvalue, 0.05)
})


test_that("slr validates input", {

  expect_error(slr("not numeric"), "`x` must be numeric")
  expect_error(slr(1:2), "`x` must have at least 3 observations")
})


test_that("slr print method works", {

  set.seed(333)
  x <- rnorm(100)

  result <- slr(x)

  # Capture output
  out <- capture.output(print(result))
  out_text <- paste(out, collapse = "\n")

  expect_match(out_text, "Simple Linear Regression")
  expect_match(out_text, "Intercept")
  expect_match(out_text, "Slope")
  expect_match(out_text, "t-statistic")
  expect_match(out_text, "p-value")
})


test_that("slr works with AirPassengers data", {

  x <- as.numeric(AirPassengers)

  old <- tswge::slr.wge(x)
  new <- slr(x)

  expect_equal(new$res, old$res)
  expect_equal(new$b0hat, unname(old$b0hat))
  expect_equal(new$b1hat, unname(old$b1hat))
})
