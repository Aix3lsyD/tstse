# Tests for hilbert()
# Computes Hilbert transform / analytic signal

test_that("hilbert matches tswge::hilbert.wge for even-length input", {
  skip_if_not_installed("tswge")

  # Test case 1: simple even-length signal
  set.seed(123)
  x <- rnorm(100)
  expect_equal(hilbert(x), tswge::hilbert.wge(x))

  # Test case 2: sinusoid
  t <- seq(0, 1, length.out = 256)
  x <- sin(2 * pi * 5 * t)
  expect_equal(hilbert(x), tswge::hilbert.wge(x))

  # Test case 3: longer signal
  x <- rnorm(1024)
  expect_equal(hilbert(x), tswge::hilbert.wge(x))

  # Test case 4: short even signal
  x <- c(1, 2, 3, 4)
  expect_equal(hilbert(x), tswge::hilbert.wge(x))

  # Test case 5: power of 2 length
  x <- rnorm(512)
  expect_equal(hilbert(x), tswge::hilbert.wge(x))
})


test_that("hilbert matches tswge::hilbert.wge for odd-length with pad_odd = FALSE", {
  skip_if_not_installed("tswge")

  # tswge drops first element for odd length
  # We match this behavior with pad_odd = FALSE
  set.seed(456)
  x <- rnorm(101)

  expect_warning(
    result <- hilbert(x, pad_odd = FALSE),
    "dropping first element"
  )
  expect_equal(result, tswge::hilbert.wge(x))
})


test_that("hilbert returns complex vector", {
  x <- rnorm(100)
  z <- hilbert(x)
  expect_true(is.complex(z))
})


test_that("hilbert returns correct length for even input", {
  x <- rnorm(100)
  expect_length(hilbert(x), 100)

  x <- rnorm(256)
  expect_length(hilbert(x), 256)
})


test_that("hilbert returns length n+1 for odd input with pad_odd = TRUE", {
  x <- rnorm(99)
  z <- hilbert(x, pad_odd = TRUE)
  expect_length(z, 100)  # padded to even

  x <- rnorm(101)
  z <- hilbert(x, pad_odd = TRUE)
  expect_length(z, 102)
})


test_that("hilbert returns length n-1 for odd input with pad_odd = FALSE", {
  x <- rnorm(101)
  expect_warning(z <- hilbert(x, pad_odd = FALSE), "dropping")
  expect_length(z, 100)
})


test_that("hilbert real part approximately equals original signal", {
  # For real input, Re(analytic) ≈ original (up to edge effects)
  set.seed(789)
  x <- rnorm(256)
  z <- hilbert(x)

  # Real part should be very close to original
  expect_equal(Re(z), x, tolerance = 1e-10)
})


test_that("hilbert envelope of pure sinusoid is approximately constant", {
  # For a pure sinusoid, the envelope should be constant (amplitude)
  t <- seq(0, 2, length.out = 512)
  amplitude <- 3
  x <- amplitude * sin(2 * pi * 10 * t)
  z <- hilbert(x)

  envelope <- Mod(z)

  # Exclude edges where edge effects occur
  middle <- 50:462
  expect_equal(envelope[middle], rep(amplitude, length(middle)), tolerance = 0.01)
})


test_that("hilbert of cosine has sine as imaginary part", {
  # H{cos(wt)} = sin(wt)
  t <- seq(0, 2, length.out = 512)
  freq <- 5
  x_cos <- cos(2 * pi * freq * t)
  z <- hilbert(x_cos)

  x_sin_expected <- sin(2 * pi * freq * t)

  # Compare in the middle to avoid edge effects
  middle <- 50:462
  expect_equal(Im(z)[middle], x_sin_expected[middle], tolerance = 0.01)
})


test_that("hilbert of sine has negative cosine as imaginary part", {
  # H{sin(wt)} = -cos(wt)
  t <- seq(0, 2, length.out = 512)
  freq <- 5
  x_sin <- sin(2 * pi * freq * t)
  z <- hilbert(x_sin)

  x_negcos_expected <- -cos(2 * pi * freq * t)

  # Compare in the middle to avoid edge effects
  middle <- 50:462
  expect_equal(Im(z)[middle], x_negcos_expected[middle], tolerance = 0.01)
})


test_that("hilbert handles length 2 input", {
  x <- c(1, 2)
  z <- hilbert(x)
  expect_length(z, 2)
  expect_true(is.complex(z))
})


test_that("hilbert handles constant signal", {
  # Hilbert transform of constant is zero
  x <- rep(5, 100)
  z <- hilbert(x)

  expect_equal(Re(z), x)
  expect_equal(Im(z), rep(0, 100), tolerance = 1e-10)
})


test_that("hilbert handles DC + sinusoid", {
  t <- seq(0, 1, length.out = 256)
  dc <- 10
  x <- dc + sin(2 * pi * 5 * t)
  z <- hilbert(x)

  # Real part should equal original
  expect_equal(Re(z), x, tolerance = 1e-10)
})


test_that("hilbert validates x parameter", {
  expect_error(hilbert("abc"), "`x` must be a numeric")
  expect_error(hilbert(list(1, 2, 3)), "`x` must be a numeric")
  expect_error(hilbert(NULL), "`x` must be a numeric")
})


test_that("hilbert rejects empty input", {
  expect_error(hilbert(numeric(0)), "length >= 1")
})


test_that("hilbert rejects NA values", {
  expect_error(hilbert(c(1, 2, NA, 4)), "contains NA")
  expect_error(hilbert(c(NA, 1, 2)), "contains NA")
})


test_that("hilbert validates pad_odd parameter", {
  x <- rnorm(100)
  expect_error(hilbert(x, pad_odd = "yes"), "`pad_odd` must be TRUE or FALSE")
  expect_error(hilbert(x, pad_odd = NA), "`pad_odd` must be TRUE or FALSE")
  expect_error(hilbert(x, pad_odd = c(TRUE, FALSE)), "`pad_odd` must be TRUE or FALSE")
})


test_that("hilbert warns for odd length with pad_odd = FALSE", {
  x <- rnorm(101)
  expect_warning(hilbert(x, pad_odd = FALSE), "Odd-length input")
})


test_that("hilbert does not warn for odd length with pad_odd = TRUE", {
  x <- rnorm(101)
  expect_silent(hilbert(x, pad_odd = TRUE))
})


test_that("hilbert does not warn for even length", {
  x <- rnorm(100)
  expect_silent(hilbert(x, pad_odd = TRUE))
  expect_silent(hilbert(x, pad_odd = FALSE))
})


test_that("hilbert handles complex input", {
  # Should work with complex input too
  x <- complex(real = rnorm(100), imaginary = rnorm(100))
  z <- hilbert(x)
  expect_length(z, 100)
  expect_true(is.complex(z))
})


test_that("hilbert is numerically stable for long signals", {
  set.seed(111)
  x <- rnorm(10000)
  z <- hilbert(x)

  expect_length(z, 10000)
  expect_false(any(is.na(z)))
  expect_false(any(is.infinite(z)))
})


test_that("hilbert works on chirp signal", {
  # For a chirp signal, just verify it produces valid output
  t <- seq(0, 1, length.out = 1024)
  f0 <- 5  # Start frequency
  f1 <- 50 # End frequency
  # Linear chirp: f(t) = f0 + (f1 - f0) * t
  phase <- 2 * pi * (f0 * t + (f1 - f0) * t^2 / 2)
  x <- sin(phase)

  z <- hilbert(x)

  expect_length(z, 1024)
  expect_false(any(is.na(z)))
  expect_false(any(is.infinite(z)))

  # Envelope should be positive and mostly near 1 (except edges)
  envelope <- Mod(z)
  middle <- 100:924
  expect_true(all(envelope[middle] > 0.8))
  expect_true(all(envelope[middle] < 1.2))
})


test_that("hilbert default pad_odd is TRUE", {
  x <- rnorm(101)
  # Should not warn with default
  expect_silent(z <- hilbert(x))
  expect_length(z, 102)  # Padded
})
