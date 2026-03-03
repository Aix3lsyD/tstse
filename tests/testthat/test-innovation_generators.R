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


# ==============================================================================
# make_gen_hetero tests
# ==============================================================================

test_that("make_gen_hetero returns a function", {
  gen <- make_gen_hetero()
  expect_type(gen, "closure")
})

test_that("make_gen_hetero generates correct length", {
  gen <- make_gen_hetero()
  x <- gen(100)
  expect_length(x, 100)
})

test_that("make_gen_hetero generates numeric with no NA", {
  gen <- make_gen_hetero()
  x <- gen(100)
  expect_type(x, "double")
  expect_false(any(is.na(x)))
})

test_that("make_gen_hetero default uses linear weights from 1 to 10", {
  gen <- make_gen_hetero()
  set.seed(42)
  x <- gen(100)
  # Variance should increase with time - later values have larger weights
  var_first <- var(x[1:25])
  var_last <- var(x[76:100])
  # Last quarter should have higher variance on average (weights 76-100 vs 1-25)
  expect_true(var_last > var_first)
})

test_that("make_gen_hetero accepts function for w", {
  # sqrt weight function
  gen <- make_gen_hetero(w = function(n) sqrt(seq_len(n)))
  x <- gen(100)
  expect_length(x, 100)
  expect_type(x, "double")
})

test_that("make_gen_hetero function w must return correct length", {
  gen <- make_gen_hetero(w = function(n) seq_len(n + 1))  # Returns n+1 values
  expect_error(gen(100), "Weight function must return exactly n values")
})

test_that("make_gen_hetero accepts numeric vector for w", {
  weights <- rep(2, 100)
  gen <- make_gen_hetero(w = weights)
  x <- gen(50)  # Request fewer than available
  expect_length(x, 50)
})

test_that("make_gen_hetero vector w left-pads when n > length(w)", {
  weights <- rep(2, 50)
  gen <- make_gen_hetero(w = weights)
  set.seed(42)
  x <- gen(100)
  expect_length(x, 100)
  # First 50 positions should use w[1]=2, last 50 should also use 2
  # (all weights are 2 here, so all positions scale equally)
})

test_that("make_gen_hetero respects sd parameter", {
  gen <- make_gen_hetero(w = function(n) rep(1, n), sd = 2)
  set.seed(42)
  x <- gen(10000)
  # With unit weights, should have sd close to 2
  expect_true(abs(sd(x) - 2) < 0.1)
})

test_that("make_gen_hetero validates sd", {
  expect_error(make_gen_hetero(sd = 0), "sd must be positive")
  expect_error(make_gen_hetero(sd = -1), "sd must be positive")
})

test_that("make_gen_hetero validates w type", {
  expect_error(make_gen_hetero(w = "invalid"), "w must be NULL, a function, or a numeric vector")
  expect_error(make_gen_hetero(w = list(1, 2, 3)), "w must be NULL, a function, or a numeric vector")
})

test_that("make_gen_hetero with constant weights produces homoscedastic output", {
  gen <- make_gen_hetero(w = function(n) rep(1, n))
  set.seed(42)
  x <- gen(1000)
  # Variance should be similar in all quarters
  vars <- sapply(1:4, function(i) var(x[((i - 1) * 250 + 1):(i * 250)]))
  cv <- sd(vars) / mean(vars)  # Coefficient of variation

  expect_true(cv < 0.3)  # Should be relatively stable
})

test_that("make_gen_hetero works with common weight patterns", {
  # Exponential weights
  gen_exp <- make_gen_hetero(w = function(n) exp(seq_len(n) / n))
  x_exp <- gen_exp(50)
  expect_length(x_exp, 50)

  # Periodic weights
  gen_periodic <- make_gen_hetero(w = function(n) 1 + 0.5 * sin(2 * pi * seq_len(n) / 12))
  x_periodic <- gen_periodic(48)
  expect_length(x_periodic, 48)
})


# ==============================================================================
# make_gen_hetero shape-based tests
# ==============================================================================

test_that("make_gen_hetero shape='linear' produces correct weights", {
  gen <- make_gen_hetero(shape = "linear", from = 2, to = 8)
  set.seed(42)
  x <- gen(100)
  set.seed(42)
  z <- rnorm(100)
  t_vec <- seq(0, 1, length.out = 100)
  expected_w <- 2 + (8 - 2) * t_vec
  expect_equal(x, expected_w * z)
})

test_that("make_gen_hetero shape='sqrt' produces correct weights", {
  gen <- make_gen_hetero(shape = "sqrt", from = 1, to = 5)
  set.seed(42)
  x <- gen(100)
  set.seed(42)
  z <- rnorm(100)
  t_vec <- seq(0, 1, length.out = 100)
  expected_w <- 1 + (5 - 1) * sqrt(t_vec)
  expect_equal(x, expected_w * z)
})

test_that("make_gen_hetero shape='log' produces correct weights", {
  gen <- make_gen_hetero(shape = "log", from = 1, to = 5)
  set.seed(42)
  x <- gen(100)
  set.seed(42)
  z <- rnorm(100)
  t_vec <- seq(0, 1, length.out = 100)
  expected_w <- 1 + (5 - 1) * log1p(t_vec * (exp(1) - 1))
  expect_equal(x, expected_w * z)
})

test_that("make_gen_hetero shape='power' uses power parameter", {
  gen <- make_gen_hetero(shape = "power", from = 1, to = 5, power = 3)
  set.seed(42)
  x <- gen(100)
  set.seed(42)
  z <- rnorm(100)
  t_vec <- seq(0, 1, length.out = 100)
  expected_w <- 1 + (5 - 1) * t_vec^3
  expect_equal(x, expected_w * z)
})

test_that("make_gen_hetero shape='exp' produces exponential interpolation", {
  gen <- make_gen_hetero(shape = "exp", from = 1, to = 5)
  set.seed(42)
  x <- gen(100)
  set.seed(42)
  z <- rnorm(100)
  t_vec <- seq(0, 1, length.out = 100)
  expected_w <- 1 * (5 / 1)^t_vec
  expect_equal(x, expected_w * z)
})

test_that("make_gen_hetero shape='step' produces piecewise constant weights", {
  gen <- make_gen_hetero(shape = "step", breaks = c(0.5), levels = c(1, 5))
  set.seed(42)
  x <- gen(100)
  expect_length(x, 100)
  # First ~50 should have smaller variance, last ~50 larger
  var_first <- var(x[1:40])
  var_last <- var(x[61:100])
  expect_true(var_last > var_first)
})

test_that("make_gen_hetero shape='periodic' produces sinusoidal weights", {
  gen <- make_gen_hetero(shape = "periodic", base_w = 2, amplitude = 1, period = 12)
  set.seed(42)
  x <- gen(48)
  set.seed(42)
  z <- rnorm(48)
  expected_w <- 2 + 1 * sin(2 * pi * seq_len(48) / 12)
  expect_equal(x, expected_w * z)
})

# ==============================================================================
# make_gen_hetero from/to range tests
# ==============================================================================

test_that("make_gen_hetero from == to gives constant weights", {
  gen <- make_gen_hetero(shape = "linear", from = 3, to = 3)
  set.seed(42)
  x <- gen(100)
  set.seed(42)
  z <- rnorm(100)
  expect_equal(x, 3 * z)
})

test_that("make_gen_hetero from > to gives decreasing variance", {
  gen <- make_gen_hetero(shape = "linear", from = 10, to = 1)
  set.seed(42)
  x <- gen(1000)
  var_first <- var(x[1:250])
  var_last <- var(x[751:1000])
  expect_true(var_first > var_last)
})

# ==============================================================================
# make_gen_hetero base distribution composition tests
# ==============================================================================

test_that("make_gen_hetero with base uses custom distribution", {
  t_gen <- make_gen_t(df = 5, scale = TRUE)
  gen <- make_gen_hetero(shape = "linear", from = 1, to = 5, base = t_gen)
  x <- gen(100)
  expect_length(x, 100)
  expect_type(x, "double")
  expect_false(any(is.na(x)))
})

test_that("make_gen_hetero base=make_gen_norm(sd=2) matches sd=2", {
  gen_sd <- make_gen_hetero(shape = "linear", from = 1, to = 5, sd = 2)
  gen_base <- make_gen_hetero(shape = "linear", from = 1, to = 5,
                               base = make_gen_norm(sd = 2))
  set.seed(42)
  x1 <- gen_sd(100)
  set.seed(42)
  x2 <- gen_base(100)
  expect_equal(x1, x2)
})

test_that("make_gen_hetero with base and legacy w works", {
  t_gen <- make_gen_t(df = 5, scale = TRUE)
  gen <- make_gen_hetero(w = function(n) rep(2, n), base = t_gen)
  x <- gen(50)
  expect_length(x, 50)
  expect_type(x, "double")
})

# ==============================================================================
# make_gen_hetero mutual exclusivity tests
# ==============================================================================

test_that("make_gen_hetero errors when both shape and w provided", {
  expect_error(
    make_gen_hetero(shape = "linear", w = function(n) rep(1, n)),
    "Cannot specify both"
  )
})

test_that("make_gen_hetero errors when both sd != 1 and base provided", {
  expect_error(
    make_gen_hetero(sd = 2, base = make_gen_t(df = 5)),
    "Cannot specify both"
  )
})

test_that("make_gen_hetero allows sd=1 (default) with base", {
  gen <- make_gen_hetero(shape = "linear", from = 1, to = 5,
                          base = make_gen_t(df = 5))
  expect_type(gen, "closure")
})

# ==============================================================================
# make_gen_hetero shape parameter validation tests
# ==============================================================================

test_that("make_gen_hetero shape='power' requires positive power", {
  expect_error(make_gen_hetero(shape = "power", power = 0), "'power' must be")
  expect_error(make_gen_hetero(shape = "power", power = -1), "'power' must be")
})

test_that("make_gen_hetero shape='step' requires valid breaks and levels", {
  expect_error(make_gen_hetero(shape = "step"), "'breaks' must be")
  expect_error(make_gen_hetero(shape = "step", breaks = c(0.5)),
               "'levels' must be")
  expect_error(make_gen_hetero(shape = "step", breaks = c(0.5), levels = c(1)),
               "length\\(breaks\\) \\+ 1")
  expect_error(make_gen_hetero(shape = "step", breaks = c(1.5), levels = c(1, 2)),
               "must be in \\(0, 1\\)")
})

test_that("make_gen_hetero shape='periodic' validates parameters", {
  expect_error(make_gen_hetero(shape = "periodic", base_w = 0), "'base_w' must be")
  expect_error(make_gen_hetero(shape = "periodic", period = 0), "'period' must be")
})

test_that("make_gen_hetero validates base is a function", {
  expect_error(make_gen_hetero(base = 5), "'base' must be a function")
  expect_error(make_gen_hetero(base = "norm"), "'base' must be a function")
})

test_that("make_gen_hetero rejects invalid shape name", {
  expect_error(make_gen_hetero(shape = "invalid"))
})

# ==============================================================================
# make_gen_hetero edge case tests
# ==============================================================================

test_that("make_gen_hetero shape works with n = 1", {
  gen <- make_gen_hetero(shape = "linear", from = 2, to = 10)
  x <- gen(1)
  expect_length(x, 1)
  expect_type(x, "double")
})

test_that("make_gen_hetero shape='step' with multiple breaks works", {
  gen <- make_gen_hetero(shape = "step", breaks = c(0.25, 0.5, 0.75),
                          levels = c(1, 2, 3, 4))
  x <- gen(100)
  expect_length(x, 100)
})
