# Tests for macoef_geg()
# Computes MA coefficients for Gegenbauer processes

test_that("macoef_geg matches tswge::macoef.geg.wge for single factor", {
  skip_if_not_installed("tswge")

  # Note: tswge default trun=300000 returns 300001 coefficients
  # We compare first n_coef values

  # Test case 1: typical parameters
  n <- 100
  result <- macoef_geg(u = 0.5, d = 0.3, n_coef = n)
  expected <- tswge::macoef.geg.wge(u = 0.5, lambda = 0.3, trun = n - 1)
  expect_equal(result, expected[1:n])

  # Test case 2: different u
  result <- macoef_geg(u = 0.8, d = 0.4, n_coef = n)
  expected <- tswge::macoef.geg.wge(u = 0.8, lambda = 0.4, trun = n - 1)
  expect_equal(result, expected[1:n])

  # Test case 3: u = 1
  result <- macoef_geg(u = 1, d = 0.25, n_coef = n)
  expected <- tswge::macoef.geg.wge(u = 1, lambda = 0.25, trun = n - 1)
  expect_equal(result, expected[1:n])

  # Test case 4: u = -1
  result <- macoef_geg(u = -1, d = 0.3, n_coef = n)
  expected <- tswge::macoef.geg.wge(u = -1, lambda = 0.3, trun = n - 1)
  expect_equal(result, expected[1:n])

  # Test case 5: u = 0
  result <- macoef_geg(u = 0, d = 0.35, n_coef = n)
  expected <- tswge::macoef.geg.wge(u = 0, lambda = 0.35, trun = n - 1)
  expect_equal(result, expected[1:n])
})


test_that("macoef_geg single factor equals gegenb output", {
  # For k=1, macoef_geg should return exactly the Gegenbauer coefficients
  expect_equal(
    macoef_geg(u = 0.5, d = 0.3, n_coef = 50),
    gegenb(u = 0.5, d = 0.3, n_coef = 50)
  )

  expect_equal(
    macoef_geg(u = 0.8, d = 0.4, n_coef = 100),
    gegenb(u = 0.8, d = 0.4, n_coef = 100)
  )

  expect_equal(
    macoef_geg(u = -0.5, d = 0.25, n_coef = 20),
    gegenb(u = -0.5, d = 0.25, n_coef = 20)
  )
})


test_that("macoef_geg two factors: convolution is correct", {
  # For k=2, the result should be convolution of two Gegenbauer sequences
  u <- c(0.5, 0.3)
  d <- c(0.25, 0.2)
  n <- 50

  C1 <- gegenb(u = u[1], d = d[1], n_coef = n)
  C2 <- gegenb(u = u[2], d = d[2], n_coef = n)

  # Manual convolution
  expected <- numeric(n)
  for (j in 1:n) {
    expected[j] <- sum(C1[1:j] * C2[j:1])
  }

  result <- macoef_geg(u = u, d = d, n_coef = n)
  expect_equal(result, expected)
})


test_that("macoef_geg two factors: first coefficient is 1", {
  # psi_0 = C1_0 * C2_0 = 1 * 1 = 1
  result <- macoef_geg(u = c(0.5, 0.3), d = c(0.25, 0.2), n_coef = 10)
  expect_equal(result[1], 1)
})


test_that("macoef_geg two factors: second coefficient is sum of C1_1 and C2_1", {
  # psi_1 = C1_0 * C2_1 + C1_1 * C2_0 = C2_1 + C1_1 = 2*d1*u1 + 2*d2*u2
  u <- c(0.5, 0.3)
  d <- c(0.25, 0.2)
  result <- macoef_geg(u = u, d = d, n_coef = 10)

  expected_psi1 <- 2 * d[1] * u[1] + 2 * d[2] * u[2]
  expect_equal(result[2], expected_psi1)
})


# Note: The original tswge::macoef.geg.wge has a bug in the k=2 case
# where C2[2] references itself before being computed.
# We test our corrected implementation against the mathematical definition.
test_that("macoef_geg two factors: known bug in original is fixed", {
  # The original has: C2[2]=2*u[2]*((d[2]-1)/2+1)*C2[2]-...
  # This should be C2[1] not C2[2] on the right side
  # Our implementation uses gegenb() which has the correct recurrence

  u <- c(0.5, 0.8)
  d <- c(0.3, 0.4)
  n <- 20

  # Get our result
  result <- macoef_geg(u = u, d = d, n_coef = n)

  # Compute expected from correct Gegenbauer coefficients
  C1 <- gegenb(u = u[1], d = d[1], n_coef = n)
  C2 <- gegenb(u = u[2], d = d[2], n_coef = n)

  expected <- numeric(n)
  for (j in 1:n) {
    expected[j] <- sum(C1[1:j] * C2[j:1])
  }

  expect_equal(result, expected)

  # Verify C2[3] (which corresponds to C_2) is computed correctly
  # C_2 = [2*(2-1+d)*u*C_1 - (2-2+2d)*C_0] / 2
  #     = [2*(1+d)*u*C_1 - 2d*C_0] / 2
  C2_manual <- numeric(3)
  C2_manual[1] <- 1
  C2_manual[2] <- 2 * d[2] * u[2]
  C2_manual[3] <- (2 * (1 + d[2]) * u[2] * C2_manual[2] - 2 * d[2] * C2_manual[1]) / 2

  expect_equal(C2[3], C2_manual[3])
})


test_that("macoef_geg returns correct length", {
  expect_length(macoef_geg(u = 0.5, d = 0.3, n_coef = 1), 1)
  expect_length(macoef_geg(u = 0.5, d = 0.3, n_coef = 10), 10)
  expect_length(macoef_geg(u = 0.5, d = 0.3, n_coef = 100), 100)

  expect_length(macoef_geg(u = c(0.5, 0.3), d = c(0.2, 0.15), n_coef = 1), 1)
  expect_length(macoef_geg(u = c(0.5, 0.3), d = c(0.2, 0.15), n_coef = 50), 50)
})


test_that("macoef_geg handles n_coef = 1", {
  expect_equal(macoef_geg(u = 0.5, d = 0.3, n_coef = 1), 1)
  expect_equal(macoef_geg(u = c(0.5, 0.3), d = c(0.2, 0.15), n_coef = 1), 1)
})


test_that("macoef_geg handles d = 0", {
  # No long memory: all psi after psi_0 are zero for single factor
  result <- macoef_geg(u = 0.5, d = 0, n_coef = 10)
  expect_equal(result, c(1, rep(0, 9)))

  # Two factors, both d = 0
  result2 <- macoef_geg(u = c(0.5, 0.3), d = c(0, 0), n_coef = 10)
  expect_equal(result2, c(1, rep(0, 9)))

  # Two factors, one d = 0: should equal the other factor's coefficients
  result3 <- macoef_geg(u = c(0.5, 0.3), d = c(0.3, 0), n_coef = 10)
  expected3 <- gegenb(u = 0.5, d = 0.3, n_coef = 10)
  expect_equal(result3, expected3)
})


test_that("macoef_geg validates u parameter", {
  expect_error(macoef_geg(u = "a", d = 0.3), "`u` must be a numeric")
  expect_error(macoef_geg(u = NA, d = 0.3), "`u` must be a numeric")
  expect_error(macoef_geg(u = c(0.5, NA), d = c(0.3, 0.2)), "`u` must be a numeric")
})


test_that("macoef_geg validates d parameter", {
  expect_error(macoef_geg(u = 0.5, d = "a"), "`d` must be a numeric")
  expect_error(macoef_geg(u = 0.5, d = NA), "`d` must be a numeric")
  expect_error(macoef_geg(u = c(0.5, 0.3), d = c(0.2, NA)), "`d` must be a numeric")
})


test_that("macoef_geg validates u and d have same length", {
  expect_error(macoef_geg(u = c(0.5, 0.3), d = 0.2), "same length")
  expect_error(macoef_geg(u = 0.5, d = c(0.2, 0.3)), "same length")
  expect_error(macoef_geg(u = c(0.5, 0.3, 0.1), d = c(0.2, 0.15)), "same length")
})


test_that("macoef_geg rejects more than 2 factors", {
  expect_error(
    macoef_geg(u = c(0.5, 0.3, 0.1), d = c(0.2, 0.15, 0.1)),
    "1 or 2 Gegenbauer factors"
  )
})


test_that("macoef_geg rejects zero factors", {
  expect_error(macoef_geg(u = numeric(0), d = numeric(0)), "1 or 2 Gegenbauer factors")
})


test_that("macoef_geg validates n_coef parameter", {
  expect_error(macoef_geg(u = 0.5, d = 0.3, n_coef = "a"), "`n_coef` must be a single integer")
  expect_error(macoef_geg(u = 0.5, d = 0.3, n_coef = NA), "`n_coef` must be a single integer")
  expect_error(macoef_geg(u = 0.5, d = 0.3, n_coef = c(10, 20)), "`n_coef` must be a single integer")
  expect_error(macoef_geg(u = 0.5, d = 0.3, n_coef = 0), "`n_coef` must be at least 1")
  expect_error(macoef_geg(u = 0.5, d = 0.3, n_coef = -5), "`n_coef` must be at least 1")
})


test_that("macoef_geg warns for u outside [-1, 1] via gegenb", {
  # Warning comes from gegenb() which macoef_geg calls
  expect_warning(macoef_geg(u = 1.5, d = 0.3), "outside \\[-1, 1\\]")
  expect_warning(macoef_geg(u = c(0.5, 1.5), d = c(0.3, 0.2)), "outside \\[-1, 1\\]")
})


test_that("macoef_geg uses default n_coef", {
  result <- macoef_geg(u = 0.5, d = 0.3)
  expect_length(result, 1000)  # default n_coef = 1000
})


test_that("macoef_geg handles symmetric factors", {
  # Two identical factors should give squared coefficients pattern
  u <- c(0.5, 0.5)
  d <- c(0.3, 0.3)
  n <- 20

  result <- macoef_geg(u = u, d = d, n_coef = n)

  # This is self-convolution of Gegenbauer coefficients
  C <- gegenb(u = 0.5, d = 0.3, n_coef = n)
  expected <- numeric(n)
  for (j in 1:n) {
    expected[j] <- sum(C[1:j] * C[j:1])
  }

  expect_equal(result, expected)
})


test_that("macoef_geg handles negative d", {
  # Negative d (anti-persistent) should work
  result <- macoef_geg(u = 0.5, d = -0.2, n_coef = 10)
  expect_length(result, 10)
  expect_equal(result[1], 1)
  expect_false(any(is.na(result)))
})


test_that("macoef_geg is numerically stable for moderate n_coef", {
  result <- macoef_geg(u = 0.5, d = 0.3, n_coef = 5000)
  expect_length(result, 5000)
  expect_false(any(is.na(result)))
  expect_false(any(is.infinite(result)))

  result2 <- macoef_geg(u = c(0.5, 0.3), d = c(0.25, 0.2), n_coef = 1000)
  expect_length(result2, 1000)
  expect_false(any(is.na(result2)))
  expect_false(any(is.infinite(result2)))
})
