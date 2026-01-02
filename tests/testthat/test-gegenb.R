# Tests for gegenb()
# Compares against tswge::gegenb.wge for equivalence

test_that("gegenb matches tswge::gegenb.wge output", {
  skip_if_not_installed("tswge")

  # Test case 1: typical GARMA parameters
  expect_equal(
    gegenb(u = 0.5, d = 0.3, n_coef = 10),
    tswge::gegenb.wge(u = 0.5, d = 0.3, n = 10)
  )

  # Test case 2: different u and d values
  expect_equal(
    gegenb(u = 0.8, d = 0.4, n_coef = 20),
    tswge::gegenb.wge(u = 0.8, d = 0.4, n = 20)
  )

  # Test case 3: u = 1 (frequency 0, fractional differencing)
  expect_equal(
    gegenb(u = 1, d = 0.25, n_coef = 15),
    tswge::gegenb.wge(u = 1, d = 0.25, n = 15)
  )

  # Test case 4: u = -1 (frequency pi)
  expect_equal(
    gegenb(u = -1, d = 0.3, n_coef = 10),
    tswge::gegenb.wge(u = -1, d = 0.3, n = 10)
  )

  # Test case 5: u = 0 (frequency pi/2)
  expect_equal(
    gegenb(u = 0, d = 0.4, n_coef = 12),
    tswge::gegenb.wge(u = 0, d = 0.4, n = 12)
  )

  # Test case 6: negative u
  expect_equal(
    gegenb(u = -0.7, d = 0.35, n_coef = 8),
    tswge::gegenb.wge(u = -0.7, d = 0.35, n = 8)
  )

  # Test case 7: larger n
  expect_equal(
    gegenb(u = 0.6, d = 0.2, n_coef = 100),
    tswge::gegenb.wge(u = 0.6, d = 0.2, n = 100)
  )

  # Test case 8: small d (weak long memory)
  expect_equal(
    gegenb(u = 0.5, d = 0.05, n_coef = 10),
    tswge::gegenb.wge(u = 0.5, d = 0.05, n = 10)
  )

  # Test case 9: d approaching 0.5 (strong long memory)
  expect_equal(
    gegenb(u = 0.5, d = 0.49, n_coef = 10),
    tswge::gegenb.wge(u = 0.5, d = 0.49, n = 10)
  )
})


test_that("gegenb returns correct first two coefficients", {
  # C_0 = 1 always
  # C_1 = 2 * d * u
  expect_equal(gegenb(u = 0.5, d = 0.3, n_coef = 2), c(1, 2 * 0.3 * 0.5))
  expect_equal(gegenb(u = 0.8, d = 0.4, n_coef = 2), c(1, 2 * 0.4 * 0.8))
  expect_equal(gegenb(u = -0.5, d = 0.25, n_coef = 2), c(1, 2 * 0.25 * -0.5))
})


test_that("gegenb handles n_coef = 1", {
  expect_equal(gegenb(u = 0.5, d = 0.3, n_coef = 1), 1)
  expect_equal(gegenb(u = 0, d = 0.5, n_coef = 1), 1)
  expect_length(gegenb(u = 0.5, d = 0.3, n_coef = 1), 1)
})


test_that("gegenb handles n_coef = 2", {
  result <- gegenb(u = 0.5, d = 0.3, n_coef = 2)
  expect_length(result, 2)
  expect_equal(result[1], 1)
  expect_equal(result[2], 2 * 0.3 * 0.5)
})


test_that("gegenb handles d = 0 (no long memory)", {
  # When d = 0, C_0 = 1 and all other coefficients are 0
  result <- gegenb(u = 0.5, d = 0, n_coef = 10)
  expect_equal(result[1], 1)
  expect_equal(result[2:10], rep(0, 9))

  result2 <- gegenb(u = 0.8, d = 0, n_coef = 5)
  expect_equal(result2, c(1, 0, 0, 0, 0))
})


test_that("gegenb returns correct length", {
  expect_length(gegenb(u = 0.5, d = 0.3, n_coef = 1), 1)
  expect_length(gegenb(u = 0.5, d = 0.3, n_coef = 5), 5)
  expect_length(gegenb(u = 0.5, d = 0.3, n_coef = 50), 50)
  expect_length(gegenb(u = 0.5, d = 0.3, n_coef = 100), 100)
})


test_that("gegenb handles u = 0 (frequency pi/2)", {
  # When u = 0, C_1 = 0, and odd-indexed coefficients should be zero
  result <- gegenb(u = 0, d = 0.3, n_coef = 10)
  expect_equal(result[1], 1)
  expect_equal(result[2], 0)  # C_1 = 2 * d * 0 = 0
  # Odd indices (C_1, C_3, C_5, ...) should be zero
  expect_equal(result[c(2, 4, 6, 8, 10)], rep(0, 5))
})


test_that("gegenb handles u = 1 (frequency 0)",
          {
            # At u = 1, this should give fractional differencing weights
            result <- gegenb(u = 1, d = 0.3, n_coef = 5)
            expect_equal(result[1], 1)
            expect_equal(result[2], 2 * 0.3 * 1)  # = 0.6
            # All coefficients should be positive for 0 < d < 0.5
            expect_true(all(result > 0))
          })


test_that("gegenb handles u = -1 (frequency pi)", {
  # At u = -1, coefficients should alternate in sign
  result <- gegenb(u = -1, d = 0.3, n_coef = 6)
  expect_equal(result[1], 1)
  expect_equal(result[2], 2 * 0.3 * -1)  # = -0.6
  # Check alternating pattern: signs should alternate
  signs <- sign(result)
  expect_equal(signs, c(1, -1, 1, -1, 1, -1))
})


test_that("gegenb validates u parameter", {
  expect_error(gegenb(u = "a", d = 0.3, n_coef = 10), "`u` must be a single numeric")
  expect_error(gegenb(u = c(0.5, 0.6), d = 0.3, n_coef = 10), "`u` must be a single numeric")
  expect_error(gegenb(u = NA, d = 0.3, n_coef = 10), "`u` must be a single numeric")
  expect_error(gegenb(u = NULL, d = 0.3, n_coef = 10), "`u` must be a single numeric")
})


test_that("gegenb validates d parameter", {
  expect_error(gegenb(u = 0.5, d = "a", n_coef = 10), "`d` must be a single numeric")
  expect_error(gegenb(u = 0.5, d = c(0.3, 0.4), n_coef = 10), "`d` must be a single numeric")
  expect_error(gegenb(u = 0.5, d = NA, n_coef = 10), "`d` must be a single numeric")
  expect_error(gegenb(u = 0.5, d = NULL, n_coef = 10), "`d` must be a single numeric")
})


test_that("gegenb validates n_coef parameter", {
  expect_error(gegenb(u = 0.5, d = 0.3, n_coef = "a"), "`n_coef` must be a single integer")
  expect_error(gegenb(u = 0.5, d = 0.3, n_coef = c(5, 10)), "`n_coef` must be a single integer")
  expect_error(gegenb(u = 0.5, d = 0.3, n_coef = NA), "`n_coef` must be a single integer")
  expect_error(gegenb(u = 0.5, d = 0.3, n_coef = NULL), "`n_coef` must be a single integer")
  expect_error(gegenb(u = 0.5, d = 0.3, n_coef = 0), "`n_coef` must be at least 1")
  expect_error(gegenb(u = 0.5, d = 0.3, n_coef = -5), "`n_coef` must be at least 1")
})


test_that("gegenb warns for u outside [-1, 1]", {
  expect_warning(gegenb(u = 1.5, d = 0.3, n_coef = 5), "outside \\[-1, 1\\]")
  expect_warning(gegenb(u = -1.5, d = 0.3, n_coef = 5), "outside \\[-1, 1\\]")
  expect_warning(gegenb(u = 2, d = 0.3, n_coef = 5), "outside \\[-1, 1\\]")
})


test_that("gegenb handles edge values of u without warning", {
  expect_silent(gegenb(u = 1, d = 0.3, n_coef = 5))
  expect_silent(gegenb(u = -1, d = 0.3, n_coef = 5))
  expect_silent(gegenb(u = 0, d = 0.3, n_coef = 5))
})


test_that("gegenb accepts numeric n_coef and coerces to integer", {
  # Should work with numeric that can be coerced
  expect_length(gegenb(u = 0.5, d = 0.3, n_coef = 5.0), 5)
  expect_length(gegenb(u = 0.5, d = 0.3, n_coef = 10L), 10)
})


test_that("gegenb handles negative d values", {
  # Negative d is mathematically valid (anti-persistent)
  result <- gegenb(u = 0.5, d = -0.3, n_coef = 5)
  expect_length(result, 5)
  expect_equal(result[1], 1)
  expect_equal(result[2], 2 * -0.3 * 0.5)  # = -0.3
})


test_that("gegenb handles d > 0.5", {
  # d > 0.5 is non-stationary but computation should still work
  result <- gegenb(u = 0.5, d = 0.8, n_coef = 5)
  expect_length(result, 5)
  expect_equal(result[1], 1)
  expect_equal(result[2], 2 * 0.8 * 0.5)  # = 0.8
})


test_that("gegenb recurrence is numerically stable for large n", {
  # Test with large n to check for numerical issues
  result <- gegenb(u = 0.5, d = 0.3, n_coef = 500)
  expect_length(result, 500)
  expect_false(any(is.na(result)))
  expect_false(any(is.infinite(result)))
})
