# Tests for compare_garch function

# ==============================================================================
# Basic functionality tests
# ==============================================================================

test_that("compare_garch returns garch_comparison class", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("WeightedPortTest")

  # Generate simple GARCH data
  set.seed(42)
  n <- 500
  omega <- 0.1
  alpha <- 0.15
  beta <- 0.75
  e <- rnorm(n)
  y <- numeric(n)
  h <- numeric(n)
  h[1] <- omega / (1 - alpha - beta)
  y[1] <- sqrt(h[1]) * e[1]
  for (t in 2:n) {
    h[t] <- omega + alpha * y[t-1]^2 + beta * h[t-1]
    y[t] <- sqrt(h[t]) * e[t]
  }

  result <- compare_garch(y, arch_range = 0:1, garch_range = 0:1, parallel = FALSE)
  expect_s3_class(result, "garch_comparison")
})

test_that("compare_garch has correct structure", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("WeightedPortTest")

  set.seed(42)
  y <- rnorm(300)

  result <- compare_garch(y, arch_range = 1, garch_range = 0:1, parallel = FALSE)

  expect_named(result, c("comparison", "fits", "distribution", "n", "call"))
  expect_s3_class(result$comparison, "data.frame")
  expect_type(result$fits, "list")
  expect_equal(result$distribution, "norm")
  expect_equal(result$n, 300)
})

test_that("compare_garch excludes (0,0) model", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("WeightedPortTest")

  set.seed(42)
  y <- rnorm(300)

  result <- compare_garch(y, arch_range = 0:1, garch_range = 0:1, parallel = FALSE)

  # Should have 3 models: ARCH(1), GARCH(0,1), GARCH(1,1)
  # (0,0) is excluded
  expect_false("ARCH(0)" %in% result$comparison$Model)
  expect_true(all(result$comparison$ARCH > 0 | result$comparison$GARCH > 0))
})

test_that("compare_garch comparison has all expected columns", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("WeightedPortTest")

  set.seed(42)
  y <- rnorm(300)

  result <- compare_garch(y, arch_range = 1, garch_range = 1, parallel = FALSE)

  expected_cols <- c("Model", "ARCH", "GARCH", "AIC", "AICC", "BIC",
                     "WLB1", "WLB2", "WLB3", "Nyblom", "Nyblom_crit",
                     "SignBias", "n_sig", "n_coef")
  expect_true(all(expected_cols %in% names(result$comparison)))
})


# ==============================================================================
# Model fitting tests
# ==============================================================================

test_that("compare_garch fits correct number of models", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("WeightedPortTest")

  set.seed(42)
  y <- rnorm(300)

  # arch_range = 0:2, garch_range = 0:1
  # Grid: (0,0), (1,0), (2,0), (0,1), (1,1), (2,1)
  # Minus (0,0) = 5 models
  result <- compare_garch(y, arch_range = 0:2, garch_range = 0:1, parallel = FALSE)
  expect_lte(nrow(result$comparison), 5)  # May have convergence failures
  expect_equal(length(result$fits), nrow(result$comparison))
})

test_that("compare_garch handles different distributions", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("WeightedPortTest")

  set.seed(42)
  y <- rnorm(300)

  result <- compare_garch(y, arch_range = 1, garch_range = 1,
                          distribution = "std", parallel = FALSE)
  expect_equal(result$distribution, "std")
})

test_that("compare_garch fits list contains uGARCHfit objects", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("WeightedPortTest")

  set.seed(42)
  y <- rnorm(300)

  result <- compare_garch(y, arch_range = 1, garch_range = 1, parallel = FALSE)

  expect_true(all(sapply(result$fits, function(x) inherits(x, "uGARCHfit"))))
})


# ==============================================================================
# Parallel processing tests
# ==============================================================================

test_that("compare_garch parallel and sequential produce same results", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("WeightedPortTest")
  skip_on_cran()  # Parallel tests can be flaky on CRAN

  set.seed(42)
  y <- rnorm(300)

  result_seq <- compare_garch(y, arch_range = 1, garch_range = 0:1,
                              parallel = FALSE)
  result_par <- compare_garch(y, arch_range = 1, garch_range = 0:1,
                              parallel = TRUE, cores = 2)

  # Same number of models
  expect_equal(nrow(result_seq$comparison), nrow(result_par$comparison))

  # Same models (may be in different order)
  expect_setequal(result_seq$comparison$Model, result_par$comparison$Model)
})


# ==============================================================================
# S3 method tests
# ==============================================================================

test_that("print.garch_comparison works", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("WeightedPortTest")

  set.seed(42)
  y <- rnorm(300)

  result <- compare_garch(y, arch_range = 1, garch_range = 0:1, parallel = FALSE)

  expect_output(print(result), "GARCH Model Comparison")
  expect_output(print(result), "Distribution:")
  expect_output(print(result), "Best by AIC")
})

test_that("print.garch_comparison returns invisible", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("WeightedPortTest")

  set.seed(42)
  y <- rnorm(300)

  result <- compare_garch(y, arch_range = 1, garch_range = 1, parallel = FALSE)
  expect_invisible(print(result))
})

test_that("summary.garch_comparison works", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("WeightedPortTest")

  set.seed(42)
  y <- rnorm(300)

  result <- compare_garch(y, arch_range = 1, garch_range = 1, parallel = FALSE)
  expect_output(summary(result), "Model")
})


# ==============================================================================
# Input validation tests
# ==============================================================================

test_that("compare_garch validates data", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("WeightedPortTest")

  expect_error(compare_garch("not numeric"), "numeric vector")
  expect_error(compare_garch(1:5), "at least 10 observations")
})

test_that("compare_garch validates distribution", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("WeightedPortTest")

  set.seed(42)
  y <- rnorm(100)

  expect_error(compare_garch(y, distribution = "invalid"), "must be one of")
})

test_that("compare_garch errors when no valid orders", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("WeightedPortTest")

  set.seed(42)
  y <- rnorm(100)

  expect_error(compare_garch(y, arch_range = 0, garch_range = 0),
               "No valid model orders")
})

test_that("compare_garch warns on NA values", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("WeightedPortTest")

  y <- rnorm(300)
  y[50] <- NA

  # Should warn about NA values (may also fail to converge, which is expected)
  expect_warning(
    tryCatch(
      compare_garch(y, arch_range = 1, garch_range = 1, parallel = FALSE),
      error = function(e) NULL
    ),
    "NA values"
  )
})


# ==============================================================================
# Helper function tests
# ==============================================================================

test_that(".garch_label produces correct labels", {
  # Access internal function
  garch_label <- tstse:::.garch_label

  expect_equal(garch_label(1, 0), "ARCH(1)")
  expect_equal(garch_label(2, 0), "ARCH(2)")
  expect_equal(garch_label(1, 1), "GARCH(1,1)")
  expect_equal(garch_label(2, 1), "GARCH(2,1)")
  expect_equal(garch_label(1, 2), "GARCH(1,2)")
})


# ==============================================================================
# Diagnostics computation tests
# ==============================================================================

test_that("compare_garch computes reasonable AIC/BIC values", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("WeightedPortTest")

  set.seed(42)
  y <- rnorm(300)

  result <- compare_garch(y, arch_range = 1:2, garch_range = 1, parallel = FALSE)

  # AIC and BIC should be finite

  expect_true(all(is.finite(result$comparison$AIC)))
  expect_true(all(is.finite(result$comparison$BIC)))

  # BIC should penalize more than AIC for larger models
  # (both should be positive for these models)
  expect_true(all(result$comparison$BIC >= result$comparison$AIC))
})

test_that("compare_garch WLB p-values are in [0,1]", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("WeightedPortTest")

  set.seed(42)
  y <- rnorm(300)

  result <- compare_garch(y, arch_range = 1, garch_range = 1, parallel = FALSE)

  wlb_vals <- c(result$comparison$WLB1, result$comparison$WLB2, result$comparison$WLB3)
  wlb_vals <- wlb_vals[!is.na(wlb_vals)]

  expect_true(all(wlb_vals >= 0 & wlb_vals <= 1))
})

test_that("compare_garch n_sig <= n_coef", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("WeightedPortTest")

  set.seed(42)
  y <- rnorm(300)

  result <- compare_garch(y, arch_range = 1:2, garch_range = 0:1, parallel = FALSE)

  expect_true(all(result$comparison$n_sig <= result$comparison$n_coef))
})


# ==============================================================================
# Edge case tests
# ==============================================================================

test_that("compare_garch handles convergence failures gracefully", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("WeightedPortTest")

  # Very short series may cause some models to fail
  set.seed(42)
  y <- rnorm(50)

  # Should not error, just fit what it can
  result <- compare_garch(y, arch_range = 1:2, garch_range = 0:1, parallel = FALSE)
  expect_s3_class(result, "garch_comparison")
  expect_gte(nrow(result$comparison), 1)
})

test_that("compare_garch handles single model", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("WeightedPortTest")

  set.seed(42)
  y <- rnorm(300)

  result <- compare_garch(y, arch_range = 1, garch_range = 1, parallel = FALSE)
  expect_equal(nrow(result$comparison), 1)
  expect_equal(result$comparison$Model, "GARCH(1,1)")
})
