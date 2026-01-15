# Tests for wbg_boot() - WBG Bootstrap Trend Test

test_that("wbg_boot returns expected structure", {
  set.seed(123)
  x <- rnorm(50)

  result <- wbg_boot(x, nb = 49, maxp = 3, seed = 456)

  expect_type(result, "list")
  expect_named(result, c("p", "phi", "vara", "pvalue", "pvalue_upper", "pvalue_lower",
                         "pvalue_asymp", "tco_obs", "boot_tstats", "n", "nb", "boot_seeds", "master_seed"))
  expect_true(result$p >= 0)
  expect_true(result$pvalue >= 0 && result$pvalue <= 1)
  expect_equal(result$nb, 49L)
})

test_that("wbg_boot detects no trend in white noise", {
  set.seed(111)
  x <- rnorm(100)  # Pure white noise, no trend

  result <- wbg_boot(x, nb = 99, maxp = 5, seed = 222)

  # Should not detect significant trend
  expect_true(result$pvalue > 0.05)
})
test_that("wbg_boot detects clear trend", {
  set.seed(333)
  n <- 100
  x <- rnorm(n) + 0.15 * seq_len(n)  # Strong trend

  result <- wbg_boot(x, nb = 99, maxp = 5, seed = 444)

  # Should detect significant trend
  expect_true(result$pvalue < 0.10)
})

test_that("wbg_boot handles AR(1) with trend", {
  set.seed(555)
  n <- 150
  z <- arima.sim(list(ar = 0.6), n = n)
  x <- z + 0.1 * seq_len(n)  # AR(1) noise + trend

  result <- wbg_boot(x, nb = 99, maxp = 5, seed = 666)

  # Should still detect trend despite autocorrelation
  expect_true(result$pvalue < 0.20)
  # Should detect AR structure
  expect_true(result$p >= 1)
})

test_that("wbg_boot is reproducible with seed", {
  set.seed(777)
  x <- rnorm(50) + 0.05 * seq_len(50)

  result1 <- wbg_boot(x, nb = 49, seed = 888)
  result2 <- wbg_boot(x, nb = 49, seed = 888)

  expect_equal(result1$pvalue, result2$pvalue)
  expect_equal(result1$tco_obs, result2$tco_obs)
  expect_equal(result1$p, result2$p)
})

test_that("wbg_boot works with different methods", {
  set.seed(999)
  x <- rnorm(50) + 0.05 * seq_len(50)

  result_burg <- wbg_boot(x, nb = 49, method = "burg", seed = 111)
  # MLE and YW require use_fast = FALSE since co_fast only supports Burg
  result_mle <- wbg_boot(x, nb = 49, method = "mle", use_fast = FALSE, seed = 111)
  result_yw <- wbg_boot(x, nb = 49, method = "yw", use_fast = FALSE, seed = 111)

  # All should return valid results
  expect_true(result_burg$pvalue >= 0 && result_burg$pvalue <= 1)
  expect_true(result_mle$pvalue >= 0 && result_mle$pvalue <= 1)
  expect_true(result_yw$pvalue >= 0 && result_yw$pvalue <= 1)
})

test_that("wbg_boot works with different criteria", {
  set.seed(222)
  x <- rnorm(50) + 0.05 * seq_len(50)

  result_aic <- wbg_boot(x, nb = 49, type = "aic", seed = 333)
  result_bic <- wbg_boot(x, nb = 49, type = "bic", seed = 333)

  # Both should return valid results
  expect_true(result_aic$pvalue >= 0 && result_aic$pvalue <= 1)
  expect_true(result_bic$pvalue >= 0 && result_bic$pvalue <= 1)
})

test_that("wbg_boot validates input", {
  expect_error(wbg_boot(character(10)), "`x` must be a non-empty numeric vector")
  expect_error(wbg_boot(numeric(0)), "`x` must be a non-empty numeric vector")
  expect_error(wbg_boot(1:10, nb = -1), "`nb` must be a positive integer")
  expect_error(wbg_boot(1:10, maxp = 0), "`maxp` must be a positive integer")
})

test_that("wbg_boot works with short series", {
  set.seed(444)
  x <- rnorm(20)

  # Should work with short series
  result <- wbg_boot(x, nb = 29, maxp = 2, seed = 555)
  expect_true(is.finite(result$pvalue))
})

test_that("wbg_boot handles AR(0) case", {
  set.seed(666)
  x <- rnorm(100)  # White noise

  result <- wbg_boot(x, nb = 49, maxp = 5, seed = 777)

  # Might select AR(0) or low order
  expect_true(result$p >= 0 && result$p <= 5)
})

test_that("wbg_boot parallel gives consistent results", {
  skip("Parallel RNG differs between sequential and parallel execution")
  skip_on_cran()
  skip_if(parallel::detectCores() < 2, "Not enough cores for parallel test")

  set.seed(888)
  x <- rnorm(50) + 0.05 * seq_len(50)

  result_seq <- wbg_boot(x, nb = 49, cores = 1, seed = 999)
  result_par <- wbg_boot(x, nb = 49, cores = 2, seed = 999)

  # With same seed, results should match
  expect_equal(result_seq$pvalue, result_par$pvalue)
})

test_that("wbg_boot matches tswge for basic case", {
  skip_if_not_installed("tswge")

  set.seed(111)
  x <- rnorm(100) + 0.05 * seq_len(100)

  result_tstse <- wbg_boot(x, nb = 99, maxp = 5, method = "burg", seed = 222)

  # Compare with tswge (approximate due to implementation differences)
  result_tswge <- tswge::wbg.boot.wge(x, nb = 99, sn = 222)

  # AR order should match
  expect_equal(result_tstse$p, result_tswge$p)
  # p-values should be similar (not exact due to seed handling differences)
  expect_true(abs(result_tstse$pvalue - result_tswge$pv) < 0.15)
})


# =============================================================================
# Time-transform optimization tests
# =============================================================================

test_that("C++ CO t-statistic matches R implementation for various AR orders", {
  set.seed(42)

  # Test cases covering various AR structures
  test_cases <- list(
    list(phi = 0.7, n = 100),                    # AR(1) positive
    list(phi = c(0.5, -0.3), n = 100),           # AR(2)
    list(phi = c(0.3, 0.2, -0.1), n = 150),      # AR(3)
    list(phi = c(0.95), n = 50),                 # Near unit root
    list(phi = c(-0.7), n = 100)                 # Negative phi, makes phi(1) > 1
  )

  for (tc in test_cases) {
    x <- arima.sim(list(ar = tc$phi), n = tc$n) + 0.05 * seq_len(tc$n)

    # C++ implementation (uses optimized time transform)
    t_cpp <- co_tstat_cpp(x, maxp = 5, criterion = "aic")

    # R implementation
    t_r <- co(x, maxp = 5, method = "burg", type = "aic")$tco

    # Should match within floating-point precision
    expect_equal(t_cpp, t_r, tolerance = 1e-10,
      info = sprintf("AR(%d) with phi=%s", length(tc$phi),
                     paste(round(tc$phi, 2), collapse = ",")))
  }
})

test_that("CO t-statistic is scale and shift invariant",
{
  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 100) + 0.05 * seq_len(100)

  t_original <- co_tstat_cpp(x, maxp = 5, criterion = "aic")

  # Scale invariance: 2*x should give same t-stat
  t_scaled <- co_tstat_cpp(2 * x, maxp = 5, criterion = "aic")
  expect_equal(t_original, t_scaled, tolerance = 1e-10)

  # Shift invariance: x + 100 should give same t-stat
  t_shifted <- co_tstat_cpp(x + 100, maxp = 5, criterion = "aic")
  expect_equal(t_original, t_shifted, tolerance = 1e-10)
})

test_that("wbg_boot use_fast=TRUE matches use_fast=FALSE", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.6), n = 100)

  # Same seed for both to ensure identical bootstrap samples
  # Both must use method="burg" since co_fast only supports Burg
  result_fast <- wbg_boot(x, nb = 49, maxp = 5, method = "burg", use_fast = TRUE, seed = 456)
  result_slow <- wbg_boot(x, nb = 49, maxp = 5, method = "burg", use_fast = FALSE, seed = 456)

  # Core statistics should match
  expect_equal(result_fast$tco_obs, result_slow$tco_obs, tolerance = 0.01)
  expect_equal(result_fast$p, result_slow$p)
  expect_equal(result_fast$phi, result_slow$phi, tolerance = 1e-6)
  expect_equal(result_fast$pvalue, result_slow$pvalue, tolerance = 0.1)
})

test_that("wbg_boot parallel matches sequential", {
  skip_on_cran()
  skip_if(parallel::detectCores(logical = FALSE) < 2, "Not enough cores")

  set.seed(123)
  x <- arima.sim(list(ar = 0.6), n = 100)

  result_seq <- wbg_boot(x, nb = 49, maxp = 5, cores = 1, seed = 456)
  result_par <- wbg_boot(x, nb = 49, maxp = 5, cores = 2, seed = 456)

  # Results should be identical with same seed
  expect_equal(result_seq$tco_obs, result_par$tco_obs)
  expect_equal(result_seq$pvalue, result_par$pvalue)
  expect_equal(result_seq$boot_tstats, result_par$boot_tstats)
  expect_equal(result_seq$boot_seeds, result_par$boot_seeds)
})
