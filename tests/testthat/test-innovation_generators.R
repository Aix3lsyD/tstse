# Tests for innovation generator factory functions

# ==============================================================================
# make_gen_norm tests
# ==============================================================================

test_that("make_gen_norm returns a function", {
  gen <- make_gen_norm()
  expect_type(gen, "closure")
})

test_that("make_gen_norm generates correct length", {
  gen <- make_gen_norm()
  x <- gen(100)
  expect_length(x, 100)
})

test_that("make_gen_norm generates numeric with no NA", {
  gen <- make_gen_norm()
  x <- gen(100)
  expect_type(x, "double")
  expect_false(any(is.na(x)))
})

test_that("make_gen_norm respects sd parameter", {
  gen <- make_gen_norm(sd = 2)
  set.seed(42)
  x <- gen(10000)
  expect_true(abs(sd(x) - 2) < 0.1)
})

test_that("make_gen_norm rejects invalid sd", {
  expect_error(make_gen_norm(sd = 0), "sd must be positive")
  expect_error(make_gen_norm(sd = -1), "sd must be positive")
})


# ==============================================================================
# make_gen_t tests
# ==============================================================================

test_that("make_gen_t returns a function", {
  gen <- make_gen_t(df = 5)
  expect_type(gen, "closure")
})

test_that("make_gen_t generates correct length", {
  gen <- make_gen_t(df = 5)
  x <- gen(100)
  expect_length(x, 100)
})

test_that("make_gen_t generates numeric with no NA", {
  gen <- make_gen_t(df = 5)
  x <- gen(100)
  expect_type(x, "double")
  expect_false(any(is.na(x)))
})

test_that("make_gen_t with scale produces unit variance", {
  gen <- make_gen_t(df = 5, scale = TRUE)
  set.seed(42)
  x <- gen(10000)
  expect_true(abs(var(x) - 1) < 0.15)
})

test_that("make_gen_t warns for df <= 2 with scale", {
  expect_warning(make_gen_t(df = 2, scale = TRUE), "infinite variance")
})

test_that("make_gen_t rejects invalid df", {
  expect_error(make_gen_t(df = 0), "df must be positive")
  expect_error(make_gen_t(df = -1), "df must be positive")
})


# ==============================================================================
# make_gen_skt tests
# ==============================================================================

test_that("make_gen_skt requires sn package", {
  skip_if_not_installed("sn")
  gen <- make_gen_skt(df = 5, alpha = 2)
  expect_type(gen, "closure")
})

test_that("make_gen_skt generates correct length", {
  skip_if_not_installed("sn")
  gen <- make_gen_skt(df = 5, alpha = 0)
  x <- gen(100)
  expect_length(x, 100)
})

test_that("make_gen_skt with scale produces unit variance", {
  skip_if_not_installed("sn")
  gen <- make_gen_skt(df = 5, alpha = 2, scale = TRUE)
  set.seed(42)
  x <- gen(10000)
  expect_true(abs(var(x) - 1) < 0.15)
})

test_that("make_gen_skt validates inputs", {
  skip_if_not_installed("sn")
  expect_error(make_gen_skt(df = 0), "df must be a single positive")
  expect_error(make_gen_skt(df = 5, alpha = c(1, 2)), "alpha must be a single")
  expect_error(make_gen_skt(df = 5, scale = "yes"), "scale must be TRUE or FALSE")
})


# ==============================================================================
# make_gen_ged tests
# ==============================================================================

test_that("make_gen_ged requires fGarch package", {
  skip_if_not_installed("fGarch")
  gen <- make_gen_ged(nu = 2)
  expect_type(gen, "closure")
})

test_that("make_gen_ged generates correct length", {
  skip_if_not_installed("fGarch")
  gen <- make_gen_ged(nu = 2)
  x <- gen(100)
  expect_length(x, 100)
})

test_that("make_gen_ged with nu=2 approximates normal", {
  skip_if_not_installed("fGarch")
  gen <- make_gen_ged(nu = 2, sd = 1)
  set.seed(42)
  x <- gen(10000)
  # Should be close to standard normal variance
  expect_true(abs(var(x) - 1) < 0.15)
})

test_that("make_gen_ged validates inputs", {
  skip_if_not_installed("fGarch")
  expect_error(make_gen_ged(nu = 0), "nu must be positive")
  expect_error(make_gen_ged(nu = 2, sd = 0), "sd must be positive")
})


# ==============================================================================
# make_gen_laplace tests
# ==============================================================================

test_that("make_gen_laplace returns a function", {
  gen <- make_gen_laplace()
  expect_type(gen, "closure")
})

test_that("make_gen_laplace generates correct length", {
  gen <- make_gen_laplace()
  x <- gen(100)
  expect_length(x, 100)
})

test_that("make_gen_laplace default has unit variance", {
  gen <- make_gen_laplace()  # default scale = 1/sqrt(2)
  set.seed(42)
  x <- gen(10000)
  expect_true(abs(var(x) - 1) < 0.15)
})

test_that("make_gen_laplace rejects invalid scale", {
  expect_error(make_gen_laplace(scale = 0), "scale must be positive")
  expect_error(make_gen_laplace(scale = -1), "scale must be positive")
})


# ==============================================================================
# make_gen_unif tests
# ==============================================================================

test_that("make_gen_unif returns a function", {
  gen <- make_gen_unif()
  expect_type(gen, "closure")
})

test_that("make_gen_unif generates correct length", {
  gen <- make_gen_unif()
  x <- gen(100)
  expect_length(x, 100)
})

test_that("make_gen_unif default has unit variance", {
  gen <- make_gen_unif()  # default half_width = sqrt(3)
  set.seed(42)
  x <- gen(10000)
  expect_true(abs(var(x) - 1) < 0.1)
})

test_that("make_gen_unif respects bounds", {
  gen <- make_gen_unif(half_width = 2)
  x <- gen(1000)
  expect_true(all(x >= -2 & x <= 2))
})

test_that("make_gen_unif rejects invalid half_width", {
  expect_error(make_gen_unif(half_width = 0), "half_width must be positive")
  expect_error(make_gen_unif(half_width = -1), "half_width must be positive")
})


# ==============================================================================
# make_gen_mixnorm tests
# ==============================================================================

test_that("make_gen_mixnorm returns a function", {
  gen <- make_gen_mixnorm()
  expect_type(gen, "closure")
})

test_that("make_gen_mixnorm generates correct length", {
  gen <- make_gen_mixnorm()
  x <- gen(100)
  expect_length(x, 100)
})

test_that("make_gen_mixnorm generates numeric with no NA", {
  gen <- make_gen_mixnorm()
  x <- gen(100)
  expect_type(x, "double")
  expect_false(any(is.na(x)))
})

test_that("make_gen_mixnorm variance matches theoretical", {
  # Variance = prob1 * sd1^2 + (1-prob1) * sd2^2 = 0.9*1 + 0.1*9 = 1.8
  gen <- make_gen_mixnorm(sd1 = 1, sd2 = 3, prob1 = 0.9)
  set.seed(42)
  x <- gen(10000)
  expect_true(abs(var(x) - 1.8) < 0.2)
})

test_that("make_gen_mixnorm validates inputs", {
  expect_error(make_gen_mixnorm(sd1 = 0), "sd1 and sd2 must be positive")
  expect_error(make_gen_mixnorm(sd2 = 0), "sd1 and sd2 must be positive")
  expect_error(make_gen_mixnorm(prob1 = 0), "prob1 must be between 0 and 1")
  expect_error(make_gen_mixnorm(prob1 = 1), "prob1 must be between 0 and 1")
})


# ==============================================================================
# make_gen_garch tests
# ==============================================================================

test_that("make_gen_garch requires rugarch package", {
  skip_if_not_installed("rugarch")
  gen <- make_gen_garch(omega = 0.1, alpha = 0.15, beta = 0.8)
  expect_type(gen, "closure")
})

test_that("make_gen_garch generates correct length", {
  skip_if_not_installed("rugarch")
  gen <- make_gen_garch(omega = 0.1, alpha = 0.15, beta = 0.8)
  x <- gen(100)
  expect_length(x, 100)
})

test_that("make_gen_garch produces heteroskedastic innovations", {
  skip_if_not_installed("rugarch")
  gen <- make_gen_garch(omega = 0.1, alpha = 0.15, beta = 0.8)
  set.seed(42)
  x <- gen(500)
  # Check for volatility clustering (squared returns autocorrelated)
  acf_sq <- acf(x^2, lag.max = 5, plot = FALSE)$acf[-1]
  expect_true(any(abs(acf_sq) > 0.05))
})

test_that("make_gen_garch warns for non-stationary params", {
  skip_if_not_installed("rugarch")
  expect_warning(
    make_gen_garch(omega = 0.1, alpha = 0.5, beta = 0.6),
    "non-stationary"
  )
})

test_that("make_gen_garch validates inputs", {
  skip_if_not_installed("rugarch")
  expect_error(make_gen_garch(omega = 0, alpha = 0.1), "omega must be positive")
  expect_error(make_gen_garch(omega = 0.1, alpha = -0.1), "alpha coefficients must be non-negative")
  expect_error(make_gen_garch(omega = 0.1, alpha = 0.1, beta = -0.1), "beta coefficients must be non-negative")
})
