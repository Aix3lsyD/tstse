# test-wbg_boot_fast.R
# Tests for the fast C++ bootstrap implementation

# ==============================================================================
# Basic functionality tests
# ==============================================================================

test_that("wbg_boot_fast returns correct structure", {
  set.seed(42)
  x <- arima.sim(list(ar = 0.7), n = 100)

  result <- wbg_boot_fast(x, nb = 49, seed = 123)

  expect_s3_class(result, "wbg_boot_fast")
  expect_true("p" %in% names(result))
  expect_true("phi" %in% names(result))
  expect_true("vara" %in% names(result))
  expect_true("pvalue" %in% names(result))
  expect_true("tco_obs" %in% names(result))
  expect_true("boot_tstats" %in% names(result))
  expect_true("n" %in% names(result))
  expect_true("nb" %in% names(result))
  expect_true("maxp" %in% names(result))
  expect_true("criterion" %in% names(result))
})

test_that("wbg_boot_fast returns correct dimensions", {
  set.seed(42)
  x <- arima.sim(list(ar = 0.7), n = 100)

  result <- wbg_boot_fast(x, nb = 99, seed = 123)

  expect_equal(result$n, 100)
  expect_equal(result$nb, 99)
  expect_equal(length(result$boot_tstats), 99)
  expect_equal(length(result$phi), result$p)
})


# ==============================================================================
# Reproducibility tests
# ==============================================================================

test_that("wbg_boot_fast is reproducible with seed", {
  set.seed(42)
  x <- arima.sim(list(ar = 0.7), n = 100)

  result1 <- wbg_boot_fast(x, nb = 49, seed = 123)
  result2 <- wbg_boot_fast(x, nb = 49, seed = 123)

  expect_identical(result1$boot_tstats, result2$boot_tstats)
  expect_identical(result1$pvalue, result2$pvalue)
  expect_identical(result1$tco_obs, result2$tco_obs)
})

test_that("wbg_boot_fast produces different results with different seeds", {
  set.seed(42)
  x <- arima.sim(list(ar = 0.7), n = 100)

  result1 <- wbg_boot_fast(x, nb = 49, seed = 123)
  result2 <- wbg_boot_fast(x, nb = 49, seed = 456)

  expect_false(identical(result1$boot_tstats, result2$boot_tstats))
})


# ==============================================================================
# Statistical validity tests
# ==============================================================================

test_that("wbg_boot_fast detects no trend correctly", {
  set.seed(42)
  # Pure AR series without trend
  x <- arima.sim(list(ar = 0.7), n = 200)

  result <- wbg_boot_fast(x, nb = 199, seed = 123)

  # p-value should be non-significant (> 0.05 typically)
  # Using a very lenient threshold since we're using few bootstrap samples
  expect_gt(result$pvalue, 0.01)
})

test_that("wbg_boot_fast detects strong trend correctly", {
  set.seed(42)
  # AR series with strong linear trend
  x <- arima.sim(list(ar = 0.7), n = 200) + 0.1 * (1:200)

  result <- wbg_boot_fast(x, nb = 199, seed = 123)

  # p-value should be significant
  expect_lt(result$pvalue, 0.05)
})

test_that("wbg_boot_fast pvalue is in valid range", {
  set.seed(42)
  x <- arima.sim(list(ar = 0.7), n = 100)

  result <- wbg_boot_fast(x, nb = 99, seed = 123)

  expect_gte(result$pvalue, 0)
  expect_lte(result$pvalue, 1)
})


# ==============================================================================
# Parameter validation tests
# ==============================================================================

test_that("wbg_boot_fast validates input", {
  expect_error(wbg_boot_fast(c(1, 2, 3), nb = 99), "too short")
  expect_error(wbg_boot_fast("not numeric", nb = 99), "numeric")
  expect_error(wbg_boot_fast(numeric(0), nb = 99), "non-empty")
})

test_that("wbg_boot_fast validates nb", {
  set.seed(42)
  x <- arima.sim(list(ar = 0.7), n = 100)

  expect_error(wbg_boot_fast(x, nb = 0), "positive")
  expect_error(wbg_boot_fast(x, nb = -1), "positive")
})

test_that("wbg_boot_fast validates maxp", {
  set.seed(42)
  x <- arima.sim(list(ar = 0.7), n = 100)

  expect_error(wbg_boot_fast(x, maxp = 0), "positive")
  expect_error(wbg_boot_fast(x, maxp = -1), "positive")
})

test_that("wbg_boot_fast validates criterion", {
  set.seed(42)
  x <- arima.sim(list(ar = 0.7), n = 100)

  expect_error(wbg_boot_fast(x, criterion = "invalid"), "arg")
})


# ==============================================================================
# Different criteria tests
# ==============================================================================

test_that("wbg_boot_fast works with all criteria", {
  set.seed(42)
  x <- arima.sim(list(ar = 0.7), n = 100)

  for (crit in c("aic", "aicc", "bic")) {
    result <- wbg_boot_fast(x, nb = 49, criterion = crit, seed = 123)
    expect_s3_class(result, "wbg_boot_fast")
    expect_equal(result$criterion, crit)
  }
})


# ==============================================================================
# Print method tests
# ==============================================================================

test_that("print.wbg_boot_fast works", {
  set.seed(42)
  x <- arima.sim(list(ar = 0.7), n = 100)
  result <- wbg_boot_fast(x, nb = 49, seed = 123)

  expect_output(print(result), "WBG Bootstrap")
  expect_output(print(result), "p-value")
  expect_output(print(result), "t-stat")
})

test_that("summary.wbg_boot_fast works", {
  set.seed(42)
  x <- arima.sim(list(ar = 0.7), n = 100)
  result <- wbg_boot_fast(x, nb = 49, seed = 123)

  expect_output(summary(result), "Summary")
  expect_output(summary(result), "Null Hypothesis")
})


# ==============================================================================
# AR(0) case tests
# ==============================================================================

test_that("wbg_boot_fast handles AR(0) selection", {
  set.seed(42)
  # White noise should select AR(0)
  x <- rnorm(200)

  result <- wbg_boot_fast(x, nb = 49, maxp = 3, seed = 123)

  # Should still work even if AR(0) is selected
  expect_s3_class(result, "wbg_boot_fast")
  expect_gte(result$p, 0)
})


# ==============================================================================
# Bootstrap t-statistics tests
# ==============================================================================

test_that("wbg_boot_fast bootstrap t-stats are finite", {
  set.seed(42)
  x <- arima.sim(list(ar = 0.7), n = 150)

  result <- wbg_boot_fast(x, nb = 99, seed = 123)

  expect_true(all(is.finite(result$boot_tstats)))
  expect_false(any(is.na(result$boot_tstats)))
})

test_that("wbg_boot_fast bootstrap t-stats have reasonable distribution", {
  set.seed(42)
  x <- arima.sim(list(ar = 0.7), n = 200)

  result <- wbg_boot_fast(x, nb = 199, seed = 123)

  # Under null, t-stats should be roughly centered around 0
  expect_lt(abs(mean(result$boot_tstats)), 1)

  # Should have some variation
  expect_gt(sd(result$boot_tstats), 0.5)
})


# ==============================================================================
# Comparison with wbg_boot tests
# ==============================================================================

test_that("wbg_boot_fast t-stat matches wbg_boot", {
  set.seed(42)
  x <- arima.sim(list(ar = 0.7), n = 150)

  fast_result <- wbg_boot_fast(x, nb = 49, seed = 123)
  orig_result <- wbg_boot(x, nb = 49, method = "burg", seed = 123)

  # Observed t-statistics should match closely
  # They use same CO algorithm
  expect_equal(fast_result$tco_obs, orig_result$tco_obs, tolerance = 0.01)
})
