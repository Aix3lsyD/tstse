# Tests for wbg_boot_garch function

# ==============================================================================
# Skip if rugarch not installed
# ==============================================================================

skip_if_no_rugarch <- function() {
  skip_if_not_installed("rugarch")
}


# ==============================================================================
# Basic structure tests
# ==============================================================================

test_that("wbg_boot_garch returns correct class", {
  skip_if_no_rugarch()

  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 200)
  stat_fn <- make_stat_ols_t()

  result <- wbg_boot_garch(x, stat_fn, nb = 19, p_max = 2,
                           seed = 456, verbose = FALSE)
  expect_s3_class(result, "wbg_boot_garch")
})

test_that("wbg_boot_garch returns correct structure", {
  skip_if_no_rugarch()

  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 200)
  stat_fn <- make_stat_ols_t()

  result <- wbg_boot_garch(x, stat_fn, nb = 19, p_max = 2,
                           seed = 456, verbose = FALSE)

  # Required elements
  expect_true("obs_stat" %in% names(result))
  expect_true("boot_dist" %in% names(result))
  expect_true("pvalue" %in% names(result))
  expect_true("pvalue_upper" %in% names(result))
  expect_true("pvalue_lower" %in% names(result))
  expect_true("ar_order" %in% names(result))
  expect_true("garch_order" %in% names(result))
  expect_true("ar_coef" %in% names(result))
  expect_true("garch_coef" %in% names(result))
  expect_true("distribution" %in% names(result))
  expect_true("n" %in% names(result))
  expect_true("nb" %in% names(result))
  expect_true("ic" %in% names(result))
})

test_that("wbg_boot_garch garch_order is always c(1,1)", {
  skip_if_no_rugarch()

  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 200)
  stat_fn <- make_stat_ols_t()

  result <- wbg_boot_garch(x, stat_fn, nb = 19, p_max = 2,
                           seed = 456, verbose = FALSE)
  expect_equal(result$garch_order, c(1, 1))
})

test_that("wbg_boot_garch garch_coef has correct names", {
  skip_if_no_rugarch()

  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 200)
  stat_fn <- make_stat_ols_t()

  result <- wbg_boot_garch(x, stat_fn, nb = 19, p_max = 2,
                           seed = 456, verbose = FALSE)

  expect_true("omega" %in% names(result$garch_coef))
  expect_true("alpha1" %in% names(result$garch_coef))
  expect_true("beta1" %in% names(result$garch_coef))
})


# ==============================================================================
# COBA adjustment tests
# ==============================================================================

test_that("wbg_boot_garch COBA elements present when bootadj=TRUE", {
  skip_if_no_rugarch()
  skip_on_cran()  # COBA is slow

  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 200)
  stat_fn <- make_stat_ols_t()

  result <- wbg_boot_garch(x, stat_fn, nb = 19, p_max = 2, bootadj = TRUE,
                           seed = 456, verbose = FALSE)

  expect_true("obs_stat_adj" %in% names(result))
  expect_true("pvalue_adj" %in% names(result))
  expect_true("adj_factor" %in% names(result))
  expect_true("median_coef" %in% names(result))
})

test_that("wbg_boot_garch COBA elements absent when bootadj=FALSE", {
  skip_if_no_rugarch()

  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 200)
  stat_fn <- make_stat_ols_t()

  result <- wbg_boot_garch(x, stat_fn, nb = 19, p_max = 2, bootadj = FALSE,
                           seed = 456, verbose = FALSE)

  expect_null(result$obs_stat_adj)
  expect_null(result$pvalue_adj)
})


# ==============================================================================
# P-value tests
# ==============================================================================

test_that("wbg_boot_garch p-values are in [0, 1]", {
  skip_if_no_rugarch()

  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 200)
  stat_fn <- make_stat_ols_t()

  result <- wbg_boot_garch(x, stat_fn, nb = 29, p_max = 2,
                           seed = 456, verbose = FALSE)

  expect_gte(result$pvalue, 0)
  expect_lte(result$pvalue, 1)
  expect_gte(result$pvalue_upper, 0)
  expect_lte(result$pvalue_upper, 1)
  expect_gte(result$pvalue_lower, 0)
  expect_lte(result$pvalue_lower, 1)
})


# ==============================================================================
# No trend detection tests
# ==============================================================================

test_that("wbg_boot_garch does not reject for white noise", {
  skip_if_no_rugarch()

  set.seed(123)
  x <- rnorm(250)
  stat_fn <- make_stat_ols_t()

  result <- wbg_boot_garch(x, stat_fn, nb = 99, p_max = 3,
                           seed = 456, verbose = FALSE)

  # Should not reject at 0.05 level for white noise
  expect_gt(result$pvalue, 0.05)
})

test_that("wbg_boot_garch does not reject for AR with no trend", {
  skip_if_no_rugarch()

  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 250)
  stat_fn <- make_stat_co()

  result <- wbg_boot_garch(x, stat_fn, nb = 99, p_max = 3,
                           seed = 456, verbose = FALSE)

  # Should typically not reject for AR with no trend
  expect_gt(result$pvalue, 0.01)
})


# ==============================================================================
# Trend detection tests
# ==============================================================================

test_that("wbg_boot_garch detects trend in trending series", {
  skip_if_no_rugarch()

  set.seed(123)
  # Strong trend
  x <- 0.15 * (1:250) + arima.sim(list(ar = 0.5), n = 250)
  stat_fn <- make_stat_ols_t()

  # Warnings expected: optim convergence and non-stationary AR when testing trending series
  result <- suppressWarnings(
    wbg_boot_garch(x, stat_fn, nb = 99, p_max = 3,
                   seed = 456, verbose = FALSE)
  )

  # Should reject for series with strong trend
  expect_lt(result$pvalue, 0.10)
})


# ==============================================================================
# Reproducibility tests
# ==============================================================================

test_that("wbg_boot_garch is reproducible with seed", {
  skip_if_no_rugarch()

  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 200)
  stat_fn <- make_stat_ols_t()

  result1 <- wbg_boot_garch(x, stat_fn, nb = 29, p_max = 2,
                            seed = 789, verbose = FALSE)
  result2 <- wbg_boot_garch(x, stat_fn, nb = 29, p_max = 2,
                            seed = 789, verbose = FALSE)

  expect_equal(result1$obs_stat, result2$obs_stat)
  expect_equal(result1$pvalue, result2$pvalue)
  expect_equal(result1$boot_dist, result2$boot_dist)
})

test_that("wbg_boot_garch differs with different seeds", {
  skip_if_no_rugarch()

  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 200)
  stat_fn <- make_stat_ols_t()

  result1 <- wbg_boot_garch(x, stat_fn, nb = 29, p_max = 2,
                            seed = 111, verbose = FALSE)
  result2 <- wbg_boot_garch(x, stat_fn, nb = 29, p_max = 2,
                            seed = 222, verbose = FALSE)

  # Same observed stat (same data)
  expect_equal(result1$obs_stat, result2$obs_stat)
  # Different bootstrap distribution
  expect_false(identical(result1$boot_dist, result2$boot_dist))
})


# ==============================================================================
# Different statistic function tests
# ==============================================================================

test_that("wbg_boot_garch works with make_stat_co", {
  skip_if_no_rugarch()

  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 200)
  stat_fn <- make_stat_co()

  result <- wbg_boot_garch(x, stat_fn, nb = 19, p_max = 2,
                           seed = 456, verbose = FALSE)
  expect_s3_class(result, "wbg_boot_garch")
  expect_false(is.na(result$pvalue))
})

test_that("wbg_boot_garch works with make_stat_ols_slope", {
  skip_if_no_rugarch()

  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 200)
  stat_fn <- make_stat_ols_slope()

  result <- wbg_boot_garch(x, stat_fn, nb = 19, p_max = 2,
                           seed = 456, verbose = FALSE)
  expect_s3_class(result, "wbg_boot_garch")
  expect_false(is.na(result$pvalue))
})

test_that("wbg_boot_garch works with make_stat_spearman", {
  skip_if_no_rugarch()

  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 200)
  stat_fn <- make_stat_spearman()

  result <- wbg_boot_garch(x, stat_fn, nb = 19, p_max = 2,
                           seed = 456, verbose = FALSE)
  expect_s3_class(result, "wbg_boot_garch")
  expect_false(is.na(result$pvalue))
})


# ==============================================================================
# Criterion tests
# ==============================================================================

test_that("wbg_boot_garch respects criterion parameter", {
  skip_if_no_rugarch()

  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 200)
  stat_fn <- make_stat_ols_t()

  result_aic <- wbg_boot_garch(x, stat_fn, nb = 19, p_max = 2, criterion = "aic",
                               seed = 456, verbose = FALSE)
  result_bic <- wbg_boot_garch(x, stat_fn, nb = 19, p_max = 2, criterion = "bic",
                               seed = 456, verbose = FALSE)

  # Both should work
  expect_false(is.na(result_aic$pvalue))
  expect_false(is.na(result_bic$pvalue))

  # Criterion should be stored
  expect_equal(result_aic$criterion, "aic")
  expect_equal(result_bic$criterion, "bic")
})


# ==============================================================================
# Distribution tests
# ==============================================================================

test_that("wbg_boot_garch works with Student's t distribution", {
  skip_if_no_rugarch()

  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 200)
  stat_fn <- make_stat_ols_t()

  result <- wbg_boot_garch(x, stat_fn, nb = 19, p_max = 2,
                           distribution = "std", seed = 456, verbose = FALSE)

  expect_s3_class(result, "wbg_boot_garch")
  expect_equal(result$distribution, "std")
  expect_false(is.na(result$pvalue))
})


# ==============================================================================
# Stationarity check tests
# ==============================================================================

test_that("wbg_boot_garch stationary_check='none' skips check", {
  skip_if_no_rugarch()

  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 200)
  stat_fn <- make_stat_ols_t()

  # Should not error or warn with stationary_check='none'
  expect_no_error(
    result <- wbg_boot_garch(x, stat_fn, nb = 19, p_max = 2,
                             stationary_check = "none",
                             seed = 456, verbose = FALSE)
  )
})


# ==============================================================================
# Input validation tests
# ==============================================================================

test_that("wbg_boot_garch validates x", {
  skip_if_no_rugarch()
  stat_fn <- make_stat_ols_t()

  expect_error(wbg_boot_garch("not numeric", stat_fn), "numeric vector")
  expect_error(wbg_boot_garch(numeric(0), stat_fn), "non-empty")
})

test_that("wbg_boot_garch validates stat_fn", {
  skip_if_no_rugarch()
  set.seed(123)
  x <- rnorm(200)

  expect_error(wbg_boot_garch(x, "not a function"), "must be a function")
  expect_error(wbg_boot_garch(x, 123), "must be a function")
})

test_that("wbg_boot_garch validates nb", {
  skip_if_no_rugarch()
  set.seed(123)
  x <- rnorm(200)
  stat_fn <- make_stat_ols_t()

  expect_error(wbg_boot_garch(x, stat_fn, nb = 0), "positive integer")
  expect_error(wbg_boot_garch(x, stat_fn, nb = -1), "positive integer")
})

test_that("wbg_boot_garch validates p_max", {
  skip_if_no_rugarch()
  set.seed(123)
  x <- rnorm(200)
  stat_fn <- make_stat_ols_t()

  expect_error(wbg_boot_garch(x, stat_fn, p_max = -1, nb = 9), "non-negative integer")
})

test_that("wbg_boot_garch requires rugarch package", {
  # Mock the package not being available
  local_mocked_bindings(
    requireNamespace = function(pkg, ...) FALSE,
    .package = "base"
  )

  set.seed(123)
  x <- rnorm(200)
  stat_fn <- make_stat_ols_t()

  expect_error(wbg_boot_garch(x, stat_fn, nb = 9), "rugarch")
})


# ==============================================================================
# Print and summary method tests
# ==============================================================================

test_that("print.wbg_boot_garch works", {
  skip_if_no_rugarch()

  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 200)
  stat_fn <- make_stat_ols_t()

  result <- wbg_boot_garch(x, stat_fn, nb = 19, p_max = 2,
                           seed = 456, verbose = FALSE)

  expect_output(print(result), "WBG Bootstrap Test")
  expect_output(print(result), "AR-GARCH")
  expect_output(print(result), "Sample size")
  expect_output(print(result), "P-values")
  expect_output(print(result), "GARCH")
})

test_that("summary.wbg_boot_garch works", {
  skip_if_no_rugarch()

  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 200)
  stat_fn <- make_stat_ols_t()

  result <- wbg_boot_garch(x, stat_fn, nb = 19, p_max = 2,
                           seed = 456, verbose = FALSE)

  expect_output(summary(result), "WBG Bootstrap Test Summary")
  expect_output(summary(result), "AR-GARCH")
  expect_output(summary(result), "Bootstrap mean")
  expect_output(summary(result), "Bootstrap SD")
})


# ==============================================================================
# Boot distribution tests
# ==============================================================================

test_that("wbg_boot_garch boot_dist has correct length", {
  skip_if_no_rugarch()

  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 200)
  stat_fn <- make_stat_ols_t()

  result <- wbg_boot_garch(x, stat_fn, nb = 37, p_max = 2,
                           seed = 456, verbose = FALSE)
  expect_equal(result$nb, length(result$boot_dist))
})


# ==============================================================================
# AR order selection tests
# ==============================================================================

test_that("wbg_boot_garch selects AR order correctly", {
  skip_if_no_rugarch()

  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 200)
  stat_fn <- make_stat_ols_t()

  result <- wbg_boot_garch(x, stat_fn, nb = 19, p_max = 3,
                           seed = 456, verbose = FALSE)

  # AR order should be between 0 and p_max
  expect_gte(result$ar_order, 0)
  expect_lte(result$ar_order, 3)
})

test_that("wbg_boot_garch works with p_max=0 (pure GARCH)", {
  skip_if_no_rugarch()

  set.seed(123)
  x <- arima.sim(list(ar = 0.3), n = 200)
  stat_fn <- make_stat_ols_t()

  result <- wbg_boot_garch(x, stat_fn, nb = 19, p_max = 0,
                           seed = 456, verbose = FALSE)

  expect_equal(result$ar_order, 0)
  expect_length(result$ar_coef, 0)
})


# ==============================================================================
# Custom statistic function test
# ==============================================================================

test_that("wbg_boot_garch works with custom stat_fn", {
  skip_if_no_rugarch()

  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 200)

  # Custom statistic: range of data
  custom_stat <- function(x) max(x) - min(x)

  result <- wbg_boot_garch(x, custom_stat, nb = 19, p_max = 2,
                           seed = 456, verbose = FALSE)
  expect_s3_class(result, "wbg_boot_garch")
  expect_false(is.na(result$pvalue))
})


# ==============================================================================
# Parallel equivalence tests
# ==============================================================================

test_that("wbg_boot_garch parallel matches sequential", {
  skip_if_no_rugarch()
  skip_on_cran()
  skip_if(parallel::detectCores(logical = FALSE) < 2, "Not enough cores")

  set.seed(123)
  x <- arima.sim(list(ar = 0.6), n = 150)
  stat_fn <- make_stat_ols_t()

  result_seq <- wbg_boot_garch(x, stat_fn, nb = 29, p_max = 2,
                                cores = 1L, seed = 456, verbose = FALSE)
  result_par <- wbg_boot_garch(x, stat_fn, nb = 29, p_max = 2,
                                cores = 2L, seed = 456, verbose = FALSE)

  # Results should be identical with same seed
  expect_equal(result_seq$obs_stat, result_par$obs_stat)
  expect_equal(result_seq$pvalue, result_par$pvalue)
  expect_equal(result_seq$boot_dist, result_par$boot_dist)
  expect_equal(result_seq$boot_seeds, result_par$boot_seeds)
})

test_that("wbg_boot_garch COBA parallel matches sequential", {
  skip_if_no_rugarch()
  skip_on_cran()
  skip_if(parallel::detectCores(logical = FALSE) < 2, "Not enough cores")

  set.seed(123)
  x <- arima.sim(list(ar = 0.6), n = 150)
  stat_fn <- make_stat_ols_t()

  result_seq <- wbg_boot_garch(x, stat_fn, nb = 19, p_max = 2, bootadj = TRUE,
                                cores = 1L, seed = 456, verbose = FALSE)
  result_par <- wbg_boot_garch(x, stat_fn, nb = 19, p_max = 2, bootadj = TRUE,
                                cores = 2L, seed = 456, verbose = FALSE)

  # Results should be identical with same seed
  expect_equal(result_seq$obs_stat, result_par$obs_stat)
  expect_equal(result_seq$pvalue, result_par$pvalue)
  expect_equal(result_seq$boot_dist, result_par$boot_dist)
  expect_equal(result_seq$pvalue_adj, result_par$pvalue_adj)
  expect_equal(result_seq$boot_seeds, result_par$boot_seeds)
  expect_equal(result_seq$boot_seeds_adj, result_par$boot_seeds_adj)
})
