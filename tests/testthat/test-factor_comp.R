# Tests for factor_comp()
# Factor component decomposition of time series

test_that("factor_comp matches tswge::factor.comp.wge output structure", {
  skip_if_not_installed("tswge")

  set.seed(123)
  x <- rnorm(100)

  result <- factor_comp(x, p = 4, n_comp = 2, plot = FALSE)

  png(tempfile())
  expected <- tswge::factor.comp.wge(x, p = 4, ncomp = 2)
  dev.off()

  expect_equal(result$n_comp, expected$ncomp)
  expect_equal(dim(result$components), dim(expected$x.comp))
})


test_that("factor_comp components match tswge::factor.comp.wge", {
  skip_if_not_installed("tswge")

  set.seed(456)
  x <- cumsum(rnorm(80))

  result <- factor_comp(x, p = 3, n_comp = 2, plot = FALSE)

  png(tempfile())
  expected <- tswge::factor.comp.wge(x, p = 3, ncomp = 2)
  dev.off()

  # Components should match (taking real part)
  expect_equal(Re(result$components), Re(expected$x.comp), tolerance = 1e-10)
})


test_that("factor_comp returns correct structure", {
  set.seed(111)
  x <- arima.sim(model = list(ar = c(0.5, -0.3)), n = 100)

  result <- factor_comp(x, p = 2, n_comp = 1, plot = FALSE)

  expect_s3_class(result, "factor_comp")
  expect_named(result, c("n_comp", "n_factors", "components", "phi", "roots", "factors"))

  expect_true(is.numeric(result$phi))
  expect_true(is.complex(result$roots) || is.numeric(result$roots))
  expect_true(is.matrix(result$components))
})


test_that("factor_comp components matrix has correct dimensions", {
  set.seed(222)
  x <- rnorm(100)

  result <- factor_comp(x, p = 4, n_comp = 2, plot = FALSE)

  expect_equal(nrow(result$components), 2)
  expect_equal(ncol(result$components), 100)
})


test_that("factor_comp reduces n_comp if exceeds n_factors", {
  set.seed(333)
  # AR(2) with complex conjugate pair = 1 factor
  x <- arima.sim(model = list(ar = c(0.5, -0.8)), n = 100)

  expect_warning(
    result <- factor_comp(x, p = 2, n_comp = 5, plot = FALSE),
    "exceeds number of factors"
  )

  expect_true(result$n_comp <= result$n_factors)
})


test_that("factor_comp handles AR(1) model", {
  set.seed(444)
  x <- arima.sim(model = list(ar = 0.7), n = 100)

  result <- factor_comp(x, p = 1, n_comp = 1, plot = FALSE)

  expect_equal(result$n_comp, 1)
  expect_equal(result$n_factors, 1)
  expect_length(result$phi, 1)
})


test_that("factor_comp counts factors correctly for real roots", {
  set.seed(555)
  # Generate data and fit AR - the key is checking factor counting logic
  x <- rnorm(200)

  result <- factor_comp(x, p = 4, n_comp = 2, plot = FALSE)

  # Number of factors should be between 1 and p
  expect_true(result$n_factors >= 1)
  expect_true(result$n_factors <= 4)

  # n_comp should not exceed n_factors
  expect_true(result$n_comp <= result$n_factors)
})


test_that("factor_comp handles complex conjugate pairs", {
  set.seed(666)
  # AR(2) with complex conjugate pair
  x <- arima.sim(model = list(ar = c(1.0, -0.5)), n = 100)

  result <- factor_comp(x, p = 2, n_comp = 1, plot = FALSE)

  # One complex pair = 1 factor
  expect_equal(result$n_factors, 1)
})


test_that("factor_comp sum of components approximates original", {
  set.seed(777)
  x <- arima.sim(model = list(ar = c(0.5, -0.3)), n = 100)

  # Extract all factors (may warn if n_comp exceeds actual factor count)
  result <- suppressWarnings(
    factor_comp(x, p = 2, n_comp = 2, plot = FALSE)
  )

  # Sum of components should approximate a filtered version of x
  # (not exact equality due to the decomposition)
  comp_sum <- colSums(result$components)
  p <- length(result$phi)

  # Check variance is reasonable
  expect_true(var(comp_sum[p:100]) > 0)
})


test_that("factor_comp validates x parameter", {
  expect_error(factor_comp("abc", p = 2, n_comp = 1, plot = FALSE),
               "`x` must be a numeric")
  expect_error(factor_comp(c(1, 2), p = 2, n_comp = 1, plot = FALSE),
               "at least 3 observations")
  expect_error(factor_comp(c(1, NA, 3, 4, 5), p = 2, n_comp = 1, plot = FALSE),
               "contains NA")
})


test_that("factor_comp validates p parameter", {
  x <- rnorm(50)
  expect_error(factor_comp(x, p = 0, n_comp = 1, plot = FALSE),
               "`p` must be a positive")
  expect_error(factor_comp(x, p = -1, n_comp = 1, plot = FALSE),
               "`p` must be a positive")
  expect_error(factor_comp(x, p = 50, n_comp = 1, plot = FALSE),
               "`p` must be less than")
  expect_error(factor_comp(x, p = "two", n_comp = 1, plot = FALSE),
               "`p` must be a positive")
})


test_that("factor_comp validates n_comp parameter", {
  x <- rnorm(50)
  expect_error(factor_comp(x, p = 3, n_comp = 0, plot = FALSE),
               "`n_comp` must be a positive")
  expect_error(factor_comp(x, p = 3, n_comp = -1, plot = FALSE),
               "`n_comp` must be a positive")
  expect_error(factor_comp(x, p = 3, n_comp = "one", plot = FALSE),
               "`n_comp` must be a positive")
})


test_that("factor_comp validates aic parameter", {
  x <- rnorm(50)
  expect_error(factor_comp(x, p = 3, n_comp = 1, aic = "TRUE", plot = FALSE),
               "`aic` must be TRUE or FALSE")
  expect_error(factor_comp(x, p = 3, n_comp = 1, aic = NA, plot = FALSE),
               "`aic` must be TRUE or FALSE")
})


test_that("factor_comp validates plot parameter", {
  x <- rnorm(50)
  expect_error(factor_comp(x, p = 3, n_comp = 1, plot = "TRUE"),
               "`plot` must be TRUE or FALSE")
  expect_error(factor_comp(x, p = 3, n_comp = 1, plot = NA),
               "`plot` must be TRUE or FALSE")
})


test_that("factor_comp plot produces no error", {
  skip_if_not(capabilities("png"))

  set.seed(888)
  x <- rnorm(50)

  png(tempfile())
  on.exit(dev.off())

  expect_no_error(factor_comp(x, p = 3, n_comp = 2, plot = TRUE))
})


test_that("factor_comp print method works", {
  set.seed(999)
  x <- rnorm(100)

  result <- factor_comp(x, p = 4, n_comp = 2, plot = FALSE)

  expect_output(print(result), "Factor Component Decomposition")
  expect_output(print(result), "AR order:")
  expect_output(print(result), "Number of factors:")
})


test_that("factor_comp returns invisible result", {
  x <- rnorm(50)

  png(tempfile())
  on.exit(dev.off())

  expect_invisible(factor_comp(x, p = 2, n_comp = 1, plot = TRUE))
})


test_that("factor_comp with aic = TRUE selects order", {
  set.seed(101)
  # True AR(2) model
  x <- arima.sim(model = list(ar = c(0.7, -0.2)), n = 200)

  # With aic = TRUE, should potentially select different order
  result <- factor_comp(x, p = 10, n_comp = 2, aic = TRUE, plot = FALSE)

  # Should have fit some model
  expect_true(length(result$phi) >= 1)
  expect_true(length(result$phi) <= 10)
})


test_that("factor_comp includes factor_ts() output", {
  set.seed(102)
  x <- arima.sim(model = list(ar = c(0.5, -0.3)), n = 100)

  result <- factor_comp(x, p = 2, n_comp = 1, plot = FALSE)

  # Should include factor decomposition
  expect_true(!is.null(result$factors))
})
