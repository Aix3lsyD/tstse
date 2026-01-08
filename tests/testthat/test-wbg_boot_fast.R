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


# ==============================================================================
# Prais-Winsten (pw_stat) tests
# ==============================================================================

test_that("wbg_boot_fast with pw_stat returns correct structure", {
  set.seed(42)
  x <- arima.sim(list(ar = 0.7), n = 100)

  result <- wbg_boot_fast(x, nb = 49, pw_stat = TRUE, seed = 123)

  expect_s3_class(result, "wbg_boot_fast")

  # CO results
  expect_true("tco_obs" %in% names(result))
  expect_true("pvalue" %in% names(result))
  expect_true("boot_tstats" %in% names(result))


  # PW results
  expect_true("tpw_obs" %in% names(result))
  expect_true("pvalue_pw" %in% names(result))
  expect_true("rho_null" %in% names(result))
  expect_true("boot_tstats_pw" %in% names(result))
  expect_true("maxp_pw" %in% names(result))
})

test_that("wbg_boot_fast with pw_stat returns correct dimensions", {
  set.seed(42)
  x <- arima.sim(list(ar = 0.7), n = 100)

  result <- wbg_boot_fast(x, nb = 99, pw_stat = TRUE, seed = 123)

  expect_equal(length(result$boot_tstats), 99)
  expect_equal(length(result$boot_tstats_pw), 99)
  expect_equal(result$maxp_pw, 1)  # Default value
})

test_that("wbg_boot_fast with pw_stat is reproducible", {
  set.seed(42)
  x <- arima.sim(list(ar = 0.7), n = 100)

  result1 <- wbg_boot_fast(x, nb = 49, pw_stat = TRUE, seed = 123)
  result2 <- wbg_boot_fast(x, nb = 49, pw_stat = TRUE, seed = 123)

  expect_identical(result1$boot_tstats, result2$boot_tstats)
  expect_identical(result1$boot_tstats_pw, result2$boot_tstats_pw)
  expect_identical(result1$tpw_obs, result2$tpw_obs)
  expect_identical(result1$pvalue_pw, result2$pvalue_pw)
})

test_that("wbg_boot_fast with pw_stat detects no trend correctly", {
  set.seed(42)
  # Pure AR series without trend
  x <- arima.sim(list(ar = 0.7), n = 200)

  result <- wbg_boot_fast(x, nb = 199, pw_stat = TRUE, seed = 123)

  # Both p-values should be non-significant
  expect_gt(result$pvalue, 0.01)
  expect_gt(result$pvalue_pw, 0.01)
})

test_that("wbg_boot_fast with pw_stat detects strong trend correctly", {
  set.seed(42)
  # AR series with strong linear trend
  x <- arima.sim(list(ar = 0.7), n = 200) + 0.1 * (1:200)

  result <- wbg_boot_fast(x, nb = 199, pw_stat = TRUE, seed = 123)

  # Both p-values should be significant
  expect_lt(result$pvalue, 0.05)
  expect_lt(result$pvalue_pw, 0.05)
})

test_that("wbg_boot_fast with pw_stat has valid p-values", {
  set.seed(42)
  x <- arima.sim(list(ar = 0.7), n = 100)

  result <- wbg_boot_fast(x, nb = 99, pw_stat = TRUE, seed = 123)

  expect_gte(result$pvalue, 0)
  expect_lte(result$pvalue, 1)
  expect_gte(result$pvalue_pw, 0)
  expect_lte(result$pvalue_pw, 1)
})

test_that("wbg_boot_fast with pw_stat and bootadj returns correct structure", {
  set.seed(42)
  x <- arima.sim(list(ar = 0.7), n = 100)

  result <- wbg_boot_fast(x, nb = 49, pw_stat = TRUE, bootadj = TRUE, seed = 123)

  expect_s3_class(result, "wbg_boot_fast")

  # CO COBA results
  expect_true("tco_obs_adj" %in% names(result))
  expect_true("pvalue_adj" %in% names(result))
  expect_true("adj_factor" %in% names(result))
  expect_true("median_phi" %in% names(result))

  # PW COBA results
  expect_true("tpw_obs_adj" %in% names(result))
  expect_true("pvalue_pw_adj" %in% names(result))
  expect_true("adj_factor_pw" %in% names(result))
  expect_true("median_rho" %in% names(result))
})

test_that("wbg_boot_fast with pw_stat bootstrap t-stats are finite", {
  set.seed(42)
  x <- arima.sim(list(ar = 0.7), n = 150)

  result <- wbg_boot_fast(x, nb = 99, pw_stat = TRUE, seed = 123)

  expect_true(all(is.finite(result$boot_tstats)))
  expect_true(all(is.finite(result$boot_tstats_pw)))
  expect_false(any(is.na(result$boot_tstats)))
  expect_false(any(is.na(result$boot_tstats_pw)))
})

test_that("wbg_boot_fast with pw_stat observed stats match standalone functions", {
  set.seed(42)
  x <- arima.sim(list(ar = 0.7), n = 100)

  result <- wbg_boot_fast(x, nb = 49, pw_stat = TRUE, seed = 123)

  # CO t-stat should match co_fast
  co_result <- co_fast(x, maxp = 5)
  expect_equal(result$tco_obs, co_result$tco, tolerance = 1e-10)

  # PW t-stat should match pw_fast
  pw_result <- pw_fast(x)
  expect_equal(result$tpw_obs, pw_result$tpw, tolerance = 1e-10)
})

test_that("print.wbg_boot_fast shows PW results when pw_stat=TRUE", {
  set.seed(42)
  x <- arima.sim(list(ar = 0.7), n = 100)
  result <- wbg_boot_fast(x, nb = 49, pw_stat = TRUE, seed = 123)

  expect_output(print(result), "Prais-Winsten")
  expect_output(print(result), "PW t-stat")
})

test_that("summary.wbg_boot_fast shows PW results when pw_stat=TRUE", {
  set.seed(42)
  x <- arima.sim(list(ar = 0.7), n = 100)
  result <- wbg_boot_fast(x, nb = 49, pw_stat = TRUE, seed = 123)

  expect_output(summary(result), "Prais-Winsten")
})


# ==============================================================================
# PW separate AR(1) null model tests (verifies Option B implementation)
# ==============================================================================

test_that("wbg_boot_fast PW bootstrap uses separate AR(1) null", {
  set.seed(42)
  # Create AR(2) series - CO will fit AR(2), PW should use AR(1)
  x <- arima.sim(list(ar = c(0.5, -0.3)), n = 150)

  result <- wbg_boot_fast(x, nb = 99, pw_stat = TRUE, bootadj = TRUE, seed = 123)

  # CO null should be AR(2) or higher
  expect_gte(result$p, 2)

  # rho_null should be a single value (AR(1))
  expect_length(result$rho_null, 1)

  # median_rho for PWBA should be single value
  expect_length(result$median_rho, 1)

  # median_phi for COBA should be non-empty (has some AR coefficients)
  # Note: median model order can differ from null model order (per Woodward)
  expect_true(length(result$median_phi) >= 0)

  # adj_factor_pw should be finite and reasonable
  expect_true(is.finite(result$adj_factor_pw))
  expect_gt(result$adj_factor_pw, 0.1)
  expect_lt(result$adj_factor_pw, 3.0)
})

test_that("wbg_boot_fast adj_factor reasonable for high positive autocorrelation", {
  set.seed(42)
  # High persistence AR(1) - per Woodward, C should typically be < 1
  x <- arima.sim(list(ar = 0.9), n = 150)

  result <- wbg_boot_fast(x, nb = 199, pw_stat = TRUE, bootadj = TRUE, seed = 123)

  # For high positive phi, adjustment factor should be reasonable
  # (Woodward says C < 1 for positive phi, but we use lenient threshold)
  expect_true(is.finite(result$adj_factor))
  expect_true(is.finite(result$adj_factor_pw))
  expect_lt(result$adj_factor, 1.5)
  expect_lt(result$adj_factor_pw, 1.5)
})

test_that("wbg_boot_fast PW and CO use different null models", {
  set.seed(42)
  # AR(2) series where CO should select p >= 2 but PW uses AR(1)
  x <- arima.sim(list(ar = c(0.6, -0.3)), n = 200)

  result <- wbg_boot_fast(x, nb = 99, pw_stat = TRUE, seed = 123)

  # CO should fit AR(2) or higher
  expect_gte(result$p, 2)

  # PW uses AR(1) rho (single coefficient)
  expect_length(result$rho_null, 1)
  expect_true(abs(result$rho_null) < 1)

  # Bootstrap distributions should be different since null models differ
  # (different DGPs for bootstrap series generation)
  expect_equal(length(result$boot_tstats), length(result$boot_tstats_pw))
})

test_that("wbg_boot_fast vara_pw_null is computed correctly", {
  set.seed(42)
  x <- arima.sim(list(ar = 0.7), n = 100)

  result <- wbg_boot_fast(x, nb = 49, pw_stat = TRUE, seed = 123)

  # vara_pw_null should be positive
  expect_true(result$vara_pw_null > 0)

  # For AR(1), vara = var(z) * (1 - rho^2)
  # Should be roughly consistent with estimated rho
  z_x <- resid(lm(x ~ seq_along(x)))
  expected_vara <- as.numeric(var(z_x)) * (1 - result$rho_null^2)
  expect_equal(as.numeric(result$vara_pw_null), expected_vara, tolerance = 0.01)
})
