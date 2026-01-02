# Tests for unit_circle()
# Plot roots on the unit circle

test_that("unit_circle works with complex roots", {
  skip_if_not(capabilities("png"))

  roots <- c(0.5 + 0.5i, 0.5 - 0.5i)

  png(tempfile())
  on.exit(dev.off())

  result <- unit_circle(roots = roots)

  expect_equal(result$real, c(0.5, 0.5))
  expect_equal(result$imaginary, c(0.5, -0.5))
})


test_that("unit_circle works with real/imaginary parts", {
  skip_if_not(capabilities("png"))

  png(tempfile())
  on.exit(dev.off())

  result <- unit_circle(real = c(0.5, 0.8), imaginary = c(0.3, -0.2))

  expect_equal(result$real, c(0.5, 0.8))
  expect_equal(result$imaginary, c(0.3, -0.2))
})


test_that("unit_circle returns correct structure", {
  skip_if_not(capabilities("png"))

  png(tempfile())
  on.exit(dev.off())

  result <- unit_circle(roots = c(0.5 + 0.5i, 1.2 + 0i))

  expect_true(is.data.frame(result))
  expect_named(result, c("real", "imaginary", "modulus", "inside"))
  expect_equal(nrow(result), 2)
})


test_that("unit_circle correctly identifies inside/outside", {
  skip_if_not(capabilities("png"))

  # Root inside: |0.5 + 0.5i| = sqrt(0.5) ≈ 0.707 < 1
  # Root outside: |1.2 + 0i| = 1.2 > 1
  roots <- c(0.5 + 0.5i, 1.2 + 0i)

  png(tempfile())
  on.exit(dev.off())

  result <- unit_circle(roots = roots)

  expect_true(result$inside[1])   # 0.5 + 0.5i is inside
  expect_false(result$inside[2])  # 1.2 is outside
})


test_that("unit_circle calculates modulus correctly", {
  skip_if_not(capabilities("png"))

  roots <- c(3 + 4i, 0.6 + 0.8i)  # |3+4i| = 5, |0.6+0.8i| = 1

  png(tempfile())
  on.exit(dev.off())

  result <- unit_circle(roots = roots)

  expect_equal(result$modulus[1], 5)
  expect_equal(result$modulus[2], 1)
})


test_that("unit_circle handles single root", {
  skip_if_not(capabilities("png"))

  png(tempfile())
  on.exit(dev.off())

  result <- unit_circle(real = 0.8, imaginary = 0)

  expect_equal(nrow(result), 1)
  expect_equal(result$real, 0.8)
  expect_equal(result$imaginary, 0)
  expect_true(result$inside)
})


test_that("unit_circle handles many roots", {
  skip_if_not(capabilities("png"))

  # 10 random roots
  set.seed(123)
  roots <- complex(real = runif(10, -2, 2), imaginary = runif(10, -2, 2))

  png(tempfile())
  on.exit(dev.off())

  result <- unit_circle(roots = roots)

  expect_equal(nrow(result), 10)
})


test_that("unit_circle handles pure real roots", {
  skip_if_not(capabilities("png"))

  png(tempfile())
  on.exit(dev.off())

  result <- unit_circle(real = c(-0.5, 0.5, 1.5))

  expect_equal(result$imaginary, c(0, 0, 0))
  expect_true(result$inside[1])   # -0.5 inside
  expect_true(result$inside[2])   # 0.5 inside
  expect_false(result$inside[3])  # 1.5 outside
})


test_that("unit_circle handles pure imaginary roots", {
  skip_if_not(capabilities("png"))

  png(tempfile())
  on.exit(dev.off())

  result <- unit_circle(roots = c(0.5i, -0.5i, 2i))

  expect_equal(result$real, c(0, 0, 0))
  expect_true(result$inside[1])
  expect_true(result$inside[2])
  expect_false(result$inside[3])
})


test_that("unit_circle imaginary default is 0", {
  skip_if_not(capabilities("png"))

  png(tempfile())
  on.exit(dev.off())

  result <- unit_circle(real = c(0.5, 0.8))

  expect_equal(result$imaginary, c(0, 0))
})


test_that("unit_circle validates roots parameter", {
  expect_error(unit_circle(roots = "abc"), "`roots` must be")
})


test_that("unit_circle validates real parameter", {
  expect_error(unit_circle(real = "abc"), "`real` must be numeric")
})


test_that("unit_circle validates imaginary parameter", {
  expect_error(unit_circle(real = 0.5, imaginary = "abc"),
               "`imaginary` must be numeric")
})


test_that("unit_circle validates length match", {
  expect_error(unit_circle(real = c(1, 2), imaginary = c(1, 2, 3)),
               "same length")
})


test_that("unit_circle errors on empty roots", {
  expect_error(unit_circle(roots = complex(0)), "No roots")
})


test_that("unit_circle returns invisible result", {
  skip_if_not(capabilities("png"))

  png(tempfile())
  on.exit(dev.off())

  expect_invisible(unit_circle(real = 0.5))
})


test_that("unit_circle expands limits for large roots", {
  skip_if_not(capabilities("png"))

  # Root at 3+0i should expand plot limits
  roots <- c(3 + 0i, 0 + 3i)

  png(tempfile())
  on.exit(dev.off())

  # Should not error even with roots far from unit circle
  expect_no_error(unit_circle(roots = roots))
})


test_that("unit_circle show_legend = FALSE works", {
  skip_if_not(capabilities("png"))

  roots <- c(0.5 + 0i, 1.5 + 0i)

  png(tempfile())
  on.exit(dev.off())

  expect_no_error(unit_circle(roots = roots, show_legend = FALSE))
})


test_that("unit_circle custom title works", {
  skip_if_not(capabilities("png"))

  png(tempfile())
  on.exit(dev.off())

  expect_no_error(unit_circle(real = 0.5, main = "My Custom Title"))
})


test_that("unit_circle handles root on unit circle", {
  skip_if_not(capabilities("png"))

  # Root exactly on unit circle: |exp(i*pi/4)| = 1
  root <- exp(1i * pi / 4)

  png(tempfile())
  on.exit(dev.off())

  result <- unit_circle(roots = root)

  expect_equal(result$modulus, 1, tolerance = 1e-10)
  expect_false(result$inside)  # modulus < 1 is inside, so exactly 1 is not inside
})
