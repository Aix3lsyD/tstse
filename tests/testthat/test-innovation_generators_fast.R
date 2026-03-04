# Tests for fast innovation generator factory functions

.expect_fast_gen_output <- function(gen, n = 100) {
  x <- gen(n)
  expect_length(x, n)
  expect_type(x, "double")
  expect_false(any(is.na(x)))
  expect_false(any(!is.finite(x)))
}


test_that("make_gen_norm_fast basic contract", {
  gen <- make_gen_norm_fast(sd = 2)
  expect_type(gen, "closure")
  .expect_fast_gen_output(gen, 120)
})

test_that("make_gen_norm_fast validates sd", {
  expect_error(make_gen_norm_fast(sd = 0), "sd must be positive")
  expect_error(make_gen_norm_fast(sd = -1), "sd must be positive")
})


test_that("make_gen_t_fast basic contract", {
  gen <- make_gen_t_fast(df = 5, scale = TRUE)
  expect_type(gen, "closure")
  .expect_fast_gen_output(gen, 120)
})

test_that("make_gen_t_fast validates df and warns for df <= 2 with scale", {
  expect_error(make_gen_t_fast(df = 0), "df must be positive")
  expect_warning(make_gen_t_fast(df = 2, scale = TRUE), "infinite variance")
})


test_that("make_gen_skt_fast basic contract", {
  gen <- make_gen_skt_fast(df = 5, alpha = 1.5, scale = TRUE)
  expect_type(gen, "closure")
  .expect_fast_gen_output(gen, 120)
})

test_that("make_gen_skt_fast validates df", {
  expect_error(make_gen_skt_fast(df = 0), "df must be positive")
})


test_that("make_gen_ged_fast basic contract", {
  gen <- make_gen_ged_fast(nu = 1.8, sd = 1.2)
  expect_type(gen, "closure")
  .expect_fast_gen_output(gen, 120)
})

test_that("make_gen_ged_fast validates parameters", {
  expect_error(make_gen_ged_fast(nu = 0), "nu must be positive")
  expect_error(make_gen_ged_fast(nu = 2, sd = 0), "sd must be positive")
})


test_that("make_gen_laplace_fast basic contract", {
  gen <- make_gen_laplace_fast(scale = 1 / sqrt(2))
  expect_type(gen, "closure")
  .expect_fast_gen_output(gen, 120)
})

test_that("make_gen_laplace_fast validates scale", {
  expect_error(make_gen_laplace_fast(scale = 0), "scale must be positive")
})


test_that("make_gen_unif_fast basic contract", {
  gen <- make_gen_unif_fast(half_width = sqrt(3))
  expect_type(gen, "closure")
  x <- gen(300)
  expect_true(all(x >= -sqrt(3) & x <= sqrt(3)))
  expect_length(x, 300)
  expect_type(x, "double")
  expect_false(any(is.na(x)))
  expect_false(any(!is.finite(x)))
})

test_that("make_gen_unif_fast validates half_width", {
  expect_error(make_gen_unif_fast(half_width = 0), "half_width must be positive")
})


test_that("make_gen_mixnorm_fast basic contract", {
  gen <- make_gen_mixnorm_fast(sd1 = 1, sd2 = 3, prob1 = 0.9)
  expect_type(gen, "closure")
  .expect_fast_gen_output(gen, 120)
})

test_that("make_gen_mixnorm_fast validates parameters", {
  expect_error(make_gen_mixnorm_fast(sd1 = 0), "sd1 and sd2 must be positive")
  expect_error(make_gen_mixnorm_fast(prob1 = 0), "prob1 must be between 0 and 1")
  expect_error(make_gen_mixnorm_fast(prob1 = 1), "prob1 must be between 0 and 1")
})


test_that("make_gen_garch_fast basic contract", {
  gen <- make_gen_garch_fast(omega = 0.1, alpha = 0.15, beta = 0.8)
  expect_type(gen, "closure")
  .expect_fast_gen_output(gen, 200)
})

test_that("make_gen_garch_fast validates parameters", {
  expect_error(make_gen_garch_fast(omega = 0, alpha = 0.1), "omega must be positive")
  expect_error(make_gen_garch_fast(omega = 0.1, alpha = -0.1), "alpha coefficients must be non-negative")
  expect_error(make_gen_garch_fast(omega = 0.1, alpha = 0.1, beta = -0.1), "beta coefficients must be non-negative")
})


test_that("make_gen_hetero_fast basic contract", {
  gen <- make_gen_hetero_fast(shape = "linear", from = 1, to = 5)
  expect_type(gen, "closure")
  .expect_fast_gen_output(gen, 120)
})

test_that("make_gen_hetero_fast validates sd and w", {
  expect_error(make_gen_hetero_fast(sd = 0), "sd must be positive")
  expect_error(make_gen_hetero_fast(w = "bad"), "w must be NULL, a function, or a numeric vector")
  expect_error(make_gen_hetero_fast(shape = "linear", w = function(n) rep(1, n)), "Cannot specify both")
  expect_error(make_gen_hetero_fast(sd = 2, base = make_gen_t(df = 5)), "Cannot specify both 'sd'")
  expect_error(make_gen_hetero_fast(base = 1), "'base' must be a function")
})


test_that("fast generators integrate with gen_aruma_flex", {
  gens <- list(
    make_gen_norm_fast(),
    make_gen_t_fast(df = 5, scale = TRUE),
    make_gen_skt_fast(df = 5, alpha = 1, scale = TRUE),
    make_gen_ged_fast(nu = 1.8),
    make_gen_laplace_fast(),
    make_gen_unif_fast(),
    make_gen_mixnorm_fast(),
    make_gen_garch_fast(omega = 0.1, alpha = 0.15, beta = 0.8),
    make_gen_hetero_fast(shape = "linear", from = 1, to = 5)
  )

  for (gen in gens) {
    out <- gen_aruma_flex(n = 80, phi = 0.5, theta = 0.2,
                          innov_gen = gen, plot = FALSE, seed = 42)
    expect_length(out$y, 80)
    expect_false(any(is.na(out$y)))
    expect_false(any(!is.finite(out$y)))
  }
})


# ==============================================================================
# Level 2: Distributional / structural parity vs standard generators
# ==============================================================================

.sample_stats <- function(x) {
  m <- mean(x)
  s <- sd(x)
  c(
    mean = m,
    sd = s,
    skew = if (s > 0) mean(((x - m) / s)^3) else 0,
    kurt = if (s > 0) mean(((x - m) / s)^4) else 3
  )
}

test_that("norm_fast matches norm distributional moments", {
  n <- 30000
  x_fast <- make_gen_norm_fast(sd = 1.7)(n)
  x_std <- make_gen_norm(sd = 1.7)(n)

  s_fast <- .sample_stats(x_fast)
  s_std <- .sample_stats(x_std)

  expect_true(abs(s_fast[["mean"]]) < 0.05)
  expect_true(abs(s_fast[["sd"]] - 1.7) < 0.08)
  expect_true(abs(s_fast[["skew"]]) < 0.08)
  expect_true(abs(s_fast[["kurt"]] - 3) < 0.2)

  expect_true(abs(s_fast[["mean"]] - s_std[["mean"]]) < 0.05)
  expect_true(abs(s_fast[["sd"]] - s_std[["sd"]]) < 0.08)
})

test_that("t_fast matches t distributional moments", {
  n <- 30000
  x_fast <- make_gen_t_fast(df = 5, scale = TRUE)(n)
  x_std <- make_gen_t(df = 5, scale = TRUE)(n)

  s_fast <- .sample_stats(x_fast)
  s_std <- .sample_stats(x_std)

  expect_true(abs(s_fast[["mean"]]) < 0.05)
  expect_true(abs(s_fast[["sd"]] - 1) < 0.08)
  expect_true(abs(s_fast[["skew"]]) < 0.12)

  expect_true(abs(s_fast[["mean"]] - s_std[["mean"]]) < 0.06)
  expect_true(abs(s_fast[["sd"]] - s_std[["sd"]]) < 0.1)
})

test_that("laplace_fast matches laplace distributional moments", {
  n <- 40000
  scale <- 1 / sqrt(2)
  x_fast <- make_gen_laplace_fast(scale = scale)(n)
  x_std <- make_gen_laplace(scale = scale)(n)

  s_fast <- .sample_stats(x_fast)
  s_std <- .sample_stats(x_std)

  expect_true(abs(s_fast[["mean"]]) < 0.05)
  expect_true(abs(s_fast[["sd"]] - 1) < 0.08)
  expect_true(abs(s_fast[["kurt"]] - 6) < 0.7)

  expect_true(abs(s_fast[["mean"]] - s_std[["mean"]]) < 0.06)
  expect_true(abs(s_fast[["sd"]] - s_std[["sd"]]) < 0.1)
})

test_that("laplace_fast remains finite for large n", {
  x <- make_gen_laplace_fast(scale = 1 / sqrt(2))(200000)
  expect_false(any(!is.finite(x)))
  expect_false(any(is.na(x)))
})

test_that("unif_fast matches uniform bounds and variance", {
  n <- 30000
  hw <- sqrt(3)
  x_fast <- make_gen_unif_fast(hw)(n)
  x_std <- make_gen_unif(hw)(n)

  expect_true(all(x_fast >= -hw & x_fast <= hw))

  s_fast <- .sample_stats(x_fast)
  s_std <- .sample_stats(x_std)

  expect_true(abs(s_fast[["mean"]]) < 0.05)
  expect_true(abs(s_fast[["sd"]] - 1) < 0.05)
  expect_true(abs(s_fast[["sd"]] - s_std[["sd"]]) < 0.06)
})

test_that("mixnorm_fast matches mixture variance", {
  n <- 40000
  sd1 <- 1
  sd2 <- 3
  p1 <- 0.9
  target_sd <- sqrt(p1 * sd1^2 + (1 - p1) * sd2^2)

  x_fast <- make_gen_mixnorm_fast(sd1 = sd1, sd2 = sd2, prob1 = p1)(n)
  x_std <- make_gen_mixnorm(sd1 = sd1, sd2 = sd2, prob1 = p1)(n)

  s_fast <- .sample_stats(x_fast)
  s_std <- .sample_stats(x_std)

  expect_true(abs(s_fast[["mean"]]) < 0.06)
  expect_true(abs(s_fast[["sd"]] - target_sd) < 0.12)
  expect_true(abs(s_fast[["sd"]] - s_std[["sd"]]) < 0.12)
})

test_that("skt_fast has scaled moments close to skt", {
  skip_if_not_installed("sn")

  n <- 30000
  x_fast <- make_gen_skt_fast(df = 6, alpha = 2, scale = TRUE)(n)
  x_std <- make_gen_skt(df = 6, alpha = 2, scale = TRUE)(n)

  s_fast <- .sample_stats(x_fast)
  s_std <- .sample_stats(x_std)

  expect_true(abs(s_fast[["mean"]]) < 0.08)
  expect_true(abs(s_fast[["sd"]] - 1) < 0.1)

  expect_true(abs(s_fast[["mean"]] - s_std[["mean"]]) < 0.12)
  expect_true(abs(s_fast[["sd"]] - s_std[["sd"]]) < 0.15)
})

test_that("ged_fast has target sd and similar moments to ged", {
  skip_if_not_installed("fGarch")

  n <- 30000
  x_fast <- make_gen_ged_fast(nu = 1.7, sd = 1.2)(n)
  x_std <- make_gen_ged(nu = 1.7, sd = 1.2)(n)

  s_fast <- .sample_stats(x_fast)
  s_std <- .sample_stats(x_std)

  expect_true(abs(s_fast[["mean"]]) < 0.08)
  expect_true(abs(s_fast[["sd"]] - 1.2) < 0.1)

  expect_true(abs(s_fast[["mean"]] - s_std[["mean"]]) < 0.1)
  expect_true(abs(s_fast[["sd"]] - s_std[["sd"]]) < 0.12)
})

test_that("garch_fast shows volatility clustering comparable to garch", {
  skip_if_not_installed("rugarch")

  n <- 3000
  fast <- make_gen_garch_fast(omega = 0.1, alpha = 0.15, beta = 0.8)(n)
  std <- make_gen_garch(omega = 0.1, alpha = 0.15, beta = 0.8)(n)

  acf_fast <- as.numeric(acf(fast^2, lag.max = 1, plot = FALSE)$acf[2])
  acf_std <- as.numeric(acf(std^2, lag.max = 1, plot = FALSE)$acf[2])

  expect_true(acf_fast > 0.05)
  expect_true(acf_std > 0.05)
  expect_true(abs(acf_fast - acf_std) < 0.2)
})

test_that("hetero_fast aligns with hetero semantics in gen_aruma_flex", {
  n <- 60
  gen_shape <- make_gen_hetero_fast(shape = "linear", from = 1, to = 5)
  gen_vec <- make_gen_hetero_fast(w = seq(1, 5, length.out = n))

  r_shape <- gen_aruma_flex(n = n, phi = 0, innov_gen = gen_shape, plot = FALSE, seed = 123)
  r_vec <- gen_aruma_flex(n = n, phi = 0, innov_gen = gen_vec, plot = FALSE, seed = 123)

  expect_equal(r_shape$y, r_vec$y)
})

test_that("hetero_fast supports base composition parity", {
  n <- 5000
  base_fast <- make_gen_t_fast(df = 5, scale = TRUE)
  base_std <- make_gen_t(df = 5, scale = TRUE)

  gen_fast <- make_gen_hetero_fast(shape = "linear", from = 1, to = 5, base = base_fast)
  gen_std <- make_gen_hetero(shape = "linear", from = 1, to = 5, base = base_std)

  x_fast <- gen_fast(n)
  x_std <- gen_std(n)

  expect_length(x_fast, n)
  expect_length(x_std, n)
  expect_true(abs(mean(x_fast) - mean(x_std)) < 0.08)
  expect_true(abs(sd(x_fast) - sd(x_std)) < 0.12)
})
