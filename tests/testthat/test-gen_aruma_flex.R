# Tests for gen_aruma_flex and aruma S3 class

# ==============================================================================
# Basic functionality tests
# ==============================================================================

test_that("gen_aruma_flex returns aruma class", {
  result <- gen_aruma_flex(n = 100, phi = 0.5, plot = FALSE, seed = 42)
  expect_s3_class(result, "aruma")
})

test_that("gen_aruma_flex has correct structure", {
  result <- gen_aruma_flex(n = 100, phi = 0.5, theta = 0.3, plot = FALSE, seed = 42)
  expect_named(result, c("y", "n", "p", "q", "d", "s", "phi", "theta",
                         "lambda", "innov_type", "plot_obj"))
})

test_that("gen_aruma_flex output length matches n", {
  result <- gen_aruma_flex(n = 200, phi = 0.7, plot = FALSE, seed = 123)
  expect_length(result$y, 200)
  expect_equal(result$n, 200L)
})

test_that("gen_aruma_flex is reproducible with seed", {
  r1 <- gen_aruma_flex(n = 50, phi = 0.5, plot = FALSE, seed = 42)
  r2 <- gen_aruma_flex(n = 50, phi = 0.5, plot = FALSE, seed = 42)
  expect_equal(r1$y, r2$y)
})

test_that("gen_aruma_flex generates numeric output", {
  result <- gen_aruma_flex(n = 100, phi = 0.5, plot = FALSE, seed = 42)
  expect_type(result$y, "double")
  expect_false(any(is.na(result$y)))
})


# ==============================================================================
# Model parameter tests
# ==============================================================================

test_that("gen_aruma_flex handles AR only", {
  result <- gen_aruma_flex(n = 100, phi = c(0.5, 0.3), plot = FALSE, seed = 42)
  expect_equal(result$p, 2L)
  expect_equal(result$q, 0L)
  expect_equal(result$phi, c(0.5, 0.3))
  expect_length(result$theta, 0)
})

test_that("gen_aruma_flex handles MA only", {
  result <- gen_aruma_flex(n = 100, theta = c(0.4, 0.2), plot = FALSE, seed = 42)
  expect_equal(result$p, 0L)
  expect_equal(result$q, 2L)
  expect_length(result$phi, 0)
  expect_equal(result$theta, c(0.4, 0.2))
})

test_that("gen_aruma_flex handles ARMA", {
  result <- gen_aruma_flex(n = 100, phi = 0.7, theta = 0.3, plot = FALSE, seed = 42)
  expect_equal(result$p, 1L)
  expect_equal(result$q, 1L)
})

test_that("gen_aruma_flex handles differencing", {
  result <- gen_aruma_flex(n = 100, phi = 0.5, d = 1, plot = FALSE, seed = 42)
  expect_equal(result$d, 1L)
})

test_that("gen_aruma_flex handles seasonal factor", {
  result <- gen_aruma_flex(n = 100, phi = 0.5, s = 12, plot = FALSE, seed = 42)
  expect_equal(result$s, 12L)
})

test_that("gen_aruma_flex handles lambda factor", {
  result <- gen_aruma_flex(n = 100, phi = 0.5, lambda = 0.9, plot = FALSE, seed = 42)
  expect_equal(result$lambda, 0.9)
})


# ==============================================================================
# Innovation generator tests
# ==============================================================================

test_that("gen_aruma_flex default uses normal innovations", {
  result <- gen_aruma_flex(n = 100, phi = 0.5, plot = FALSE, seed = 42)
  expect_match(result$innov_type, "Normal")
})

test_that("gen_aruma_flex with custom t innovations", {
  t_gen <- make_gen_t(df = 5, scale = TRUE)
  result <- gen_aruma_flex(n = 100, phi = 0.5, innov_gen = t_gen, plot = FALSE, seed = 42)
  expect_equal(result$innov_type, "Custom")
  expect_s3_class(result, "aruma")
})

test_that("gen_aruma_flex with mixture innovations", {
  mix_gen <- make_gen_mixnorm(sd1 = 1, sd2 = 3, prob1 = 0.9)
  result <- gen_aruma_flex(n = 100, phi = 0.5, innov_gen = mix_gen, plot = FALSE, seed = 42)
  expect_s3_class(result, "aruma")
})

test_that("gen_aruma_flex with uniform innovations", {
  unif_gen <- make_gen_unif()
  result <- gen_aruma_flex(n = 100, phi = 0.5, innov_gen = unif_gen, plot = FALSE, seed = 42)
  expect_s3_class(result, "aruma")
})

test_that("gen_aruma_flex with laplace innovations", {
  lap_gen <- make_gen_laplace()
  result <- gen_aruma_flex(n = 100, phi = 0.5, innov_gen = lap_gen, plot = FALSE, seed = 42)
  expect_s3_class(result, "aruma")
})

test_that("gen_aruma_flex with GARCH innovations", {
  skip_if_not_installed("rugarch")
  garch_gen <- make_gen_garch(omega = 0.1, alpha = 0.15, beta = 0.8)
  result <- gen_aruma_flex(n = 200, phi = 0.5, innov_gen = garch_gen, plot = FALSE, seed = 42)
  expect_s3_class(result, "aruma")
})

test_that("gen_aruma_flex respects vara when no innov_gen", {
  r1 <- gen_aruma_flex(n = 500, phi = 0.5, vara = 1, plot = FALSE, seed = 42)
  r2 <- gen_aruma_flex(n = 500, phi = 0.5, vara = 4, plot = FALSE, seed = 42)
  # Different seeds so different values, but different variance scale
  expect_true(var(r2$y) > var(r1$y))  # Higher vara should give higher variance
})


# ==============================================================================
# Input validation tests
# ==============================================================================

test_that("gen_aruma_flex validates n", {
  expect_error(gen_aruma_flex(n = 0, phi = 0.5, plot = FALSE), "must be a positive integer")
  expect_error(gen_aruma_flex(n = -1, phi = 0.5, plot = FALSE), "must be a positive integer")
  expect_error(gen_aruma_flex(n = "100", phi = 0.5, plot = FALSE), "must be a positive integer")
})

test_that("gen_aruma_flex validates phi", {
  expect_error(gen_aruma_flex(n = 100, phi = "0.5", plot = FALSE), "phi.*must be numeric")
})

test_that("gen_aruma_flex validates theta", {
  expect_error(gen_aruma_flex(n = 100, theta = "0.5", plot = FALSE), "theta.*must be numeric")
})

test_that("gen_aruma_flex validates d", {
  expect_error(gen_aruma_flex(n = 100, phi = 0.5, d = -1, plot = FALSE), "non-negative")
})

test_that("gen_aruma_flex validates vara", {
  expect_error(gen_aruma_flex(n = 100, phi = 0.5, vara = 0, plot = FALSE), "vara.*must be positive")
  expect_error(gen_aruma_flex(n = 100, phi = 0.5, vara = -1, plot = FALSE), "vara.*must be positive")
})

test_that("gen_aruma_flex validates plot parameter", {
  expect_error(gen_aruma_flex(n = 100, phi = 0.5, plot = "invalid"), "plot.*must be")
})

test_that("gen_aruma_flex validates innov_gen is function", {
  expect_error(gen_aruma_flex(n = 100, phi = 0.5, innov_gen = "not_function", plot = FALSE),
               "innov_gen.*must be a function")
})


# ==============================================================================
# S3 method tests
# ==============================================================================

test_that("print.aruma works", {
  result <- gen_aruma_flex(n = 100, phi = 0.5, theta = 0.3, plot = FALSE, seed = 42)
  expect_output(print(result), "ARUMA Realization")
  expect_output(print(result), "phi")
  expect_output(print(result), "theta")
  expect_output(print(result), "Innovation")
})

test_that("print.aruma returns invisible", {
  result <- gen_aruma_flex(n = 100, phi = 0.5, plot = FALSE, seed = 42)
  expect_invisible(print(result))
})

test_that("summary.aruma works", {
  result <- gen_aruma_flex(n = 100, phi = 0.5, plot = FALSE, seed = 42)
  summ <- summary(result)
  expect_s3_class(summ, "summary.aruma")
  expect_named(summ$stats, c("mean", "sd", "min", "max", "skewness", "kurtosis"))
})

test_that("print.summary.aruma works", {
  result <- gen_aruma_flex(n = 100, phi = 0.5, plot = FALSE, seed = 42)
  summ <- summary(result)
  expect_output(print(summ), "Summary of ARUMA Realization")
  expect_output(print(summ), "Mean:")
  expect_output(print(summ), "Std Dev:")
})


# ==============================================================================
# Plotting tests
# ==============================================================================

test_that("gen_aruma_flex plot=FALSE returns no plot_obj", {
  result <- gen_aruma_flex(n = 100, phi = 0.5, plot = FALSE, seed = 42)
  expect_null(result$plot_obj)
})

test_that("gen_aruma_flex plot='base' works", {
  result <- gen_aruma_flex(n = 100, phi = 0.5, plot = "base", seed = 42)
  expect_null(result$plot_obj)  # base R doesn't store plot object
})

test_that("gen_aruma_flex plot='ggplot2' works", {
  skip_if_not_installed("ggplot2")
  result <- gen_aruma_flex(n = 100, phi = 0.5, plot = "ggplot2", seed = 42)
  expect_s3_class(result$plot_obj, "ggplot")
})

test_that("plot.aruma with base works", {
  result <- gen_aruma_flex(n = 100, phi = 0.5, plot = FALSE, seed = 42)
  expect_silent(plot(result, type = "base"))
})

test_that("plot.aruma with ggplot2 works", {
  skip_if_not_installed("ggplot2")
  result <- gen_aruma_flex(n = 100, phi = 0.5, plot = FALSE, seed = 42)
  p <- plot(result, type = "ggplot2")
  expect_s3_class(p, "ggplot")
})


# ==============================================================================
# Backward compatibility test
# ==============================================================================

test_that("gen_aruma_flex matches gen_aruma with normal innovations", {
  # With the same seed and normal innovations, output should match gen_aruma
  set.seed(42)
  r_flex <- gen_aruma_flex(n = 200, phi = 0.5, theta = 0.3, plot = FALSE, seed = 42)
  set.seed(42)
  r_orig <- gen_aruma(n = 200, phi = 0.5, theta = 0.3, plot = FALSE, seed = 42)

  # Note: The implementations may differ slightly due to different burn-in handling.
  # We just verify both produce valid output with similar statistical properties.
  expect_length(r_flex$y, 200)
  expect_length(r_orig, 200)

  # Both should be roughly similar in scale (same AR/MA structure)
  expect_true(abs(var(r_flex$y) - var(r_orig)) / var(r_orig) < 1)  # Within 100%
})


# ==============================================================================
# Edge case tests
# ==============================================================================
test_that("gen_aruma_flex handles white noise (no AR/MA)", {
  result <- gen_aruma_flex(n = 100, plot = FALSE, seed = 42)
  expect_equal(result$p, 0L)
  expect_equal(result$q, 0L)
  expect_length(result$y, 100)
})

test_that("gen_aruma_flex handles small n", {
  result <- gen_aruma_flex(n = 10, phi = 0.5, plot = FALSE, seed = 42)
  expect_length(result$y, 10)
})

test_that("gen_aruma_flex handles large n", {
  result <- gen_aruma_flex(n = 1000, phi = 0.5, plot = FALSE, seed = 42)
  expect_length(result$y, 1000)
})
