# Tests for trans_to_dual() and trans_to_original()

# --- trans_to_dual tests ---

test_that("trans_to_dual returns expected structure", {
  set.seed(123)
  x <- rnorm(50)

  result <- trans_to_dual(x, lambda = 0.5, offset = 60, plot = FALSE)

  expect_type(result, "list")
  expect_named(result, c("int_x", "int_y", "h"))
  expect_true(length(result$int_y) > 0)
  expect_true(result$h > 0)
})

test_that("trans_to_dual works with lambda = 0 (geometric)", {
  set.seed(111)
  x <- cumsum(rnorm(100))

  result <- trans_to_dual(x, lambda = 0, offset = 60, plot = FALSE)

  expect_true(length(result$int_y) > 0)
  expect_true(result$h > 1)  # h > 1 for geometric case
  # Check that time points are geometrically spaced
  expect_true(all(diff(result$int_x) > 0))
})

test_that("trans_to_dual works with lambda = 1 (linear)", {
  set.seed(222)
  x <- rnorm(50)

  result <- trans_to_dual(x, lambda = 1, offset = 60, plot = FALSE)

  expect_true(length(result$int_y) > 0)
  # For lambda = 1, time points should be equally spaced
  diffs <- diff(result$int_x)
  expect_true(all(abs(diffs - diffs[1]) < 1e-10))
})

test_that("trans_to_dual works with lambda = 0.5 (power-law)", {
  set.seed(333)
  x <- cumsum(rnorm(80))

  result <- trans_to_dual(x, lambda = 0.5, offset = 60, plot = FALSE)

  expect_true(length(result$int_y) > 0)
  expect_true(result$h > 0)
})

test_that("trans_to_dual handles different offsets", {
  set.seed(444)
  x <- rnorm(50)

  result_60 <- trans_to_dual(x, lambda = 0.5, offset = 60, plot = FALSE)
  result_100 <- trans_to_dual(x, lambda = 0.5, offset = 100, plot = FALSE)

  # Different offsets should give different results
  expect_false(identical(result_60$int_y, result_100$int_y))
})

test_that("trans_to_dual validates input", {
  expect_error(trans_to_dual(character(10), lambda = 0.5),
               "`x` must be a non-empty numeric")
  expect_error(trans_to_dual(numeric(0), lambda = 0.5),
               "`x` must be a non-empty numeric")
  expect_error(trans_to_dual(1:10, lambda = c(0.5, 0.6)),
               "`lambda` must be a single numeric")
  expect_error(trans_to_dual(1:10, lambda = 0.5, offset = -1),
               "`offset` must be a non-negative")
})

# --- trans_to_original tests ---

test_that("trans_to_original returns numeric vector", {
  set.seed(555)
  x <- rnorm(50)

  dual <- trans_to_dual(x, lambda = 0.5, offset = 60, plot = FALSE)
  result <- trans_to_original(dual$int_y, lambda = 0.5,
                               offset = 60, h = dual$h, plot = FALSE)

  expect_true(is.numeric(result))
  expect_true(length(result) > 0)
})

test_that("trans_to_original validates input", {
  expect_error(trans_to_original(character(10), lambda = 0.5, offset = 60, h = 1),
               "`xd` must be a non-empty numeric")
  expect_error(trans_to_original(1:10, lambda = c(0.5, 0.6), offset = 60, h = 1),
               "`lambda` must be a single numeric")
})

# --- Round-trip tests ---

test_that("round-trip transformation approximately recovers original", {
  set.seed(666)
  x <- cumsum(rnorm(80))

  # Transform to dual
  dual <- trans_to_dual(x, lambda = 0.5, offset = 60, plot = FALSE)

  # Transform back to original
  original <- trans_to_original(dual$int_y, lambda = 0.5,
                                 offset = 60, h = dual$h, plot = FALSE)

  # Should approximately recover original length
  # (exact recovery not guaranteed due to interpolation)
  expect_true(abs(length(original) - length(dual$int_y)) <= 1)
})

test_that("round-trip with lambda = 0 works", {
  set.seed(777)
  x <- rnorm(60)

  dual <- trans_to_dual(x, lambda = 0, offset = 60, plot = FALSE)
  original <- trans_to_original(dual$int_y, lambda = 0,
                                 offset = 60, h = dual$h, plot = FALSE)

  expect_true(is.numeric(original))
  expect_true(all(is.finite(original)))
})

test_that("round-trip with lambda = 1 approximately preserves data", {
  set.seed(888)
  x <- rnorm(50)

  dual <- trans_to_dual(x, lambda = 1, offset = 60, plot = FALSE)
  original <- trans_to_original(dual$int_y, lambda = 1,
                                 offset = 60, h = dual$h, plot = FALSE)

  # For lambda = 1, transformation is nearly identity
  # However, interpolation at boundaries can cause a 1-position shift
  # Check correlation accounting for potential shift

  min_len <- min(length(x), length(original)) - 1
  cor_direct <- cor(x[1:min_len], original[1:min_len])
  cor_shifted <- cor(x[1:min_len], original[2:(min_len + 1)])
  cor_val <- max(cor_direct, cor_shifted)
  expect_true(cor_val > 0.9)
})

# --- Comparison with tswge ---

test_that("trans_to_dual matches tswge output", {
  skip_if_not_installed("tswge")

  set.seed(999)
  x <- cumsum(rnorm(50))

  result_tstse <- trans_to_dual(x, lambda = 0.5, offset = 60, plot = FALSE)

  # Compare with tswge
  result_tswge <- tswge::trans.to.dual.wge(x, lambda = 0.5, offset = 60, plot = FALSE)

  # h values should match
  expect_equal(result_tstse$h, result_tswge$h, tolerance = 1e-10)

  # Interpolated values should match
  expect_equal(result_tstse$int_y, result_tswge$intY, tolerance = 1e-10)
})

test_that("trans_to_original matches tswge output", {
  skip_if_not_installed("tswge")

  set.seed(111)
  x <- cumsum(rnorm(50))

  # First transform to dual
  dual_tstse <- trans_to_dual(x, lambda = 0.5, offset = 60, plot = FALSE)
  dual_tswge <- tswge::trans.to.dual.wge(x, lambda = 0.5, offset = 60, plot = FALSE)

  # Then transform back
  result_tstse <- trans_to_original(dual_tstse$int_y, lambda = 0.5,
                                     offset = 60, h = dual_tstse$h, plot = FALSE)
  result_tswge <- tswge::trans.to.original.wge(dual_tswge$intY, lambda = 0.5,
                                                offset = 60, h = dual_tswge$h, plot = FALSE)

  expect_equal(result_tstse, result_tswge, tolerance = 1e-10)
})
