# Tests for wbg_boot_flex function

# ==============================================================================
# Basic structure tests
# ==============================================================================

test_that("wbg_boot_flex returns correct class", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 100)
  stat_fn <- make_stat_ols_t()

  result <- wbg_boot_flex(x, stat_fn, nb = 19, seed = 456, verbose = FALSE)
  expect_s3_class(result, "wbg_boot_flex")
})

test_that("wbg_boot_flex returns correct structure", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 100)
  stat_fn <- make_stat_ols_t()

  result <- wbg_boot_flex(x, stat_fn, nb = 19, seed = 456, verbose = FALSE)

  # Required elements
  expect_true("obs_stat" %in% names(result))
  expect_true("boot_dist" %in% names(result))
  expect_true("pvalue" %in% names(result))
  expect_true("pvalue_upper" %in% names(result))
  expect_true("pvalue_lower" %in% names(result))
  expect_true("ar_order" %in% names(result))
  expect_true("ar_phi" %in% names(result))
  expect_true("n" %in% names(result))
  expect_true("nb" %in% names(result))
})

test_that("wbg_boot_flex COBA elements present when bootadj=TRUE", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 100)
  stat_fn <- make_stat_ols_t()

  result <- wbg_boot_flex(x, stat_fn, nb = 19, bootadj = TRUE,
                          seed = 456, verbose = FALSE)

  expect_true("obs_stat_adj" %in% names(result))
  expect_true("pvalue_adj" %in% names(result))
  expect_true("adj_factor" %in% names(result))
  expect_true("median_phi" %in% names(result))
})

test_that("wbg_boot_flex COBA elements absent when bootadj=FALSE", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 100)
  stat_fn <- make_stat_ols_t()

  result <- wbg_boot_flex(x, stat_fn, nb = 19, bootadj = FALSE,
                          seed = 456, verbose = FALSE)

  expect_null(result$obs_stat_adj)
  expect_null(result$pvalue_adj)
})


# ==============================================================================
# P-value tests
# ==============================================================================

test_that("wbg_boot_flex p-values are in [0, 1]", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 100)
  stat_fn <- make_stat_ols_t()

  result <- wbg_boot_flex(x, stat_fn, nb = 29, seed = 456, verbose = FALSE)

  expect_gte(result$pvalue, 0)
  expect_lte(result$pvalue, 1)
  expect_gte(result$pvalue_upper, 0)
  expect_lte(result$pvalue_upper, 1)
  expect_gte(result$pvalue_lower, 0)
  expect_lte(result$pvalue_lower, 1)
})

test_that("wbg_boot_flex adjusted p-value is in [0, 1]", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 100)
  stat_fn <- make_stat_ols_t()

  result <- wbg_boot_flex(x, stat_fn, nb = 29, bootadj = TRUE,
                          seed = 456, verbose = FALSE)

  expect_gte(result$pvalue_adj, 0)
  expect_lte(result$pvalue_adj, 1)
})


# ==============================================================================
# No trend detection tests
# ==============================================================================

test_that("wbg_boot_flex does not reject for white noise", {
  set.seed(123)
  x <- rnorm(150)
  stat_fn <- make_stat_ols_t()

  result <- wbg_boot_flex(x, stat_fn, nb = 99, seed = 456, verbose = FALSE)

  # Should not reject at 0.05 level for white noise
  expect_gt(result$pvalue, 0.05)
})

test_that("wbg_boot_flex does not reject for AR with no trend", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 150)
  stat_fn <- make_stat_co()

  result <- wbg_boot_flex(x, stat_fn, nb = 99, seed = 456, verbose = FALSE)

  # Should typically not reject for AR with no trend
  expect_gt(result$pvalue, 0.01)
})


# ==============================================================================
# Trend detection tests
# ==============================================================================
test_that("wbg_boot_flex detects trend in trending series", {
  set.seed(123)
  # Strong trend
  x <- 0.15 * (1:150) + arima.sim(list(ar = 0.5), n = 150)
  stat_fn <- make_stat_ols_t()

  result <- wbg_boot_flex(x, stat_fn, nb = 99, seed = 456, verbose = FALSE)

  # Should reject for series with strong trend
  expect_lt(result$pvalue, 0.10)
})


# ==============================================================================
# Reproducibility tests
# ==============================================================================

test_that("wbg_boot_flex is reproducible with seed", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 100)
  stat_fn <- make_stat_ols_t()

  result1 <- wbg_boot_flex(x, stat_fn, nb = 29, seed = 789, verbose = FALSE)
  result2 <- wbg_boot_flex(x, stat_fn, nb = 29, seed = 789, verbose = FALSE)

  expect_equal(result1$obs_stat, result2$obs_stat)
  expect_equal(result1$pvalue, result2$pvalue)
  expect_equal(result1$boot_dist, result2$boot_dist)
})

test_that("wbg_boot_flex differs with different seeds", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 100)
  stat_fn <- make_stat_ols_t()

  result1 <- wbg_boot_flex(x, stat_fn, nb = 29, seed = 111, verbose = FALSE)
  result2 <- wbg_boot_flex(x, stat_fn, nb = 29, seed = 222, verbose = FALSE)

  # Same observed stat (same data)
  expect_equal(result1$obs_stat, result2$obs_stat)
  # Different bootstrap distribution
  expect_false(identical(result1$boot_dist, result2$boot_dist))
})


# ==============================================================================
# Different statistic function tests
# ==============================================================================

test_that("wbg_boot_flex works with make_stat_co", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 100)
  stat_fn <- make_stat_co()

  result <- wbg_boot_flex(x, stat_fn, nb = 19, seed = 456, verbose = FALSE)
  expect_s3_class(result, "wbg_boot_flex")
  expect_false(is.na(result$pvalue))
})

test_that("wbg_boot_flex works with make_stat_ols_slope", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 100)
  stat_fn <- make_stat_ols_slope()

  result <- wbg_boot_flex(x, stat_fn, nb = 19, seed = 456, verbose = FALSE)
  expect_s3_class(result, "wbg_boot_flex")
  expect_false(is.na(result$pvalue))
})

test_that("wbg_boot_flex works with make_stat_spearman", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 100)
  stat_fn <- make_stat_spearman()

  result <- wbg_boot_flex(x, stat_fn, nb = 19, seed = 456, verbose = FALSE)
  expect_s3_class(result, "wbg_boot_flex")
  expect_false(is.na(result$pvalue))
})

test_that("wbg_boot_flex works with make_stat_sen", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 100)
  stat_fn <- make_stat_sen()

  result <- wbg_boot_flex(x, stat_fn, nb = 19, seed = 456, verbose = FALSE)
  expect_s3_class(result, "wbg_boot_flex")
  expect_false(is.na(result$pvalue))
})

test_that("wbg_boot_flex works with make_stat_bn", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 100)
  stat_fn <- make_stat_bn()

  result <- wbg_boot_flex(x, stat_fn, nb = 19, seed = 456, verbose = FALSE)
  expect_s3_class(result, "wbg_boot_flex")
  expect_false(is.na(result$pvalue))
})

test_that("wbg_boot_flex works with make_stat_lr", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 100)
  stat_fn <- make_stat_lr()

  result <- wbg_boot_flex(x, stat_fn, nb = 19, seed = 456, verbose = FALSE)
  expect_s3_class(result, "wbg_boot_flex")
  expect_false(is.na(result$pvalue))
})


# ==============================================================================
# AR method tests
# ==============================================================================

test_that("wbg_boot_flex works with ar_method='burg'", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 100)
  stat_fn <- make_stat_ols_t()

  result <- wbg_boot_flex(x, stat_fn, nb = 19, ar_method = "burg",
                          seed = 456, verbose = FALSE)
  expect_equal(result$ar_method, "burg")
  expect_false(is.na(result$pvalue))
})

test_that("wbg_boot_flex works with ar_method='mle'", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 100)
  stat_fn <- make_stat_ols_t()

  result <- wbg_boot_flex(x, stat_fn, nb = 19, ar_method = "mle",
                          seed = 456, verbose = FALSE)
  expect_equal(result$ar_method, "mle")
  expect_false(is.na(result$pvalue))
})


# ==============================================================================
# Criterion tests
# ==============================================================================

test_that("wbg_boot_flex respects criterion parameter", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 100)
  stat_fn <- make_stat_ols_t()

  result_aic <- wbg_boot_flex(x, stat_fn, nb = 19, criterion = "aic",
                              seed = 456, verbose = FALSE)
  result_bic <- wbg_boot_flex(x, stat_fn, nb = 19, criterion = "bic",
                              seed = 456, verbose = FALSE)

  # Both should work
  expect_false(is.na(result_aic$pvalue))
  expect_false(is.na(result_bic$pvalue))
})


# ==============================================================================
# Input validation tests
# ==============================================================================

test_that("wbg_boot_flex validates x", {
  stat_fn <- make_stat_ols_t()

  expect_error(wbg_boot_flex("not numeric", stat_fn), "numeric vector")
  expect_error(wbg_boot_flex(numeric(0), stat_fn), "non-empty")
})

test_that("wbg_boot_flex validates stat_fn", {
  set.seed(123)
  x <- rnorm(100)

  expect_error(wbg_boot_flex(x, "not a function"), "must be a function")
  expect_error(wbg_boot_flex(x, 123), "must be a function")
})

test_that("wbg_boot_flex validates nb", {
  set.seed(123)
  x <- rnorm(100)
  stat_fn <- make_stat_ols_t()

  expect_error(wbg_boot_flex(x, stat_fn, nb = 0), "positive integer")
  expect_error(wbg_boot_flex(x, stat_fn, nb = -1), "positive integer")
})

test_that("wbg_boot_flex validates p_max", {
  set.seed(123)
  x <- rnorm(100)
  stat_fn <- make_stat_ols_t()

  expect_error(wbg_boot_flex(x, stat_fn, p_max = 0, nb = 9), "positive integer")
  expect_error(wbg_boot_flex(x, stat_fn, p_max = -1, nb = 9), "positive integer")
})


# ==============================================================================
# Print and summary method tests
# ==============================================================================

test_that("print.wbg_boot_flex works", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 100)
  stat_fn <- make_stat_ols_t()

  result <- wbg_boot_flex(x, stat_fn, nb = 19, seed = 456, verbose = FALSE)

  expect_output(print(result), "WBG Bootstrap Test")
  expect_output(print(result), "Sample size")
  expect_output(print(result), "P-values")
})

test_that("summary.wbg_boot_flex works", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 100)
  stat_fn <- make_stat_ols_t()

  result <- wbg_boot_flex(x, stat_fn, nb = 19, seed = 456, verbose = FALSE)

  expect_output(summary(result), "WBG Bootstrap Test Summary")
  expect_output(summary(result), "Bootstrap mean")
  expect_output(summary(result), "Bootstrap SD")
})


# ==============================================================================
# COBA adjustment tests
# ==============================================================================

test_that("wbg_boot_flex COBA produces different p-value", {
  set.seed(123)
  x <- 0.05 * (1:120) + arima.sim(list(ar = 0.7), n = 120)
  stat_fn <- make_stat_ols_t()

  result <- wbg_boot_flex(x, stat_fn, nb = 49, bootadj = TRUE,
                          seed = 456, verbose = FALSE)

  # Adjustment factor should be finite and positive
  expect_true(is.finite(result$adj_factor))
  expect_gt(result$adj_factor, 0)

  # Adjusted statistic should differ from original
  # (unless adjustment factor is exactly 1)
  if (abs(result$adj_factor - 1) > 0.01) {
    expect_false(result$obs_stat == result$obs_stat_adj)
  }
})


# ==============================================================================
# Boot distribution tests
# ==============================================================================

test_that("wbg_boot_flex boot_dist has correct length", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 100)
  stat_fn <- make_stat_ols_t()

  result <- wbg_boot_flex(x, stat_fn, nb = 37, seed = 456, verbose = FALSE)
  expect_length(result$boot_dist, 37)
})

test_that("wbg_boot_flex boot_dist contains no NA", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 100)
  stat_fn <- make_stat_ols_t()

  result <- wbg_boot_flex(x, stat_fn, nb = 29, seed = 456, verbose = FALSE)
  expect_false(any(is.na(result$boot_dist)))
})


# ==============================================================================
# Custom statistic function test
# ==============================================================================

test_that("wbg_boot_flex works with custom stat_fn", {
  set.seed(123)
  x <- arima.sim(list(ar = 0.7), n = 100)

  # Custom statistic: range of data
  custom_stat <- function(x) max(x) - min(x)

  result <- wbg_boot_flex(x, custom_stat, nb = 19, seed = 456, verbose = FALSE)
  expect_s3_class(result, "wbg_boot_flex")
  expect_false(is.na(result$pvalue))
})


# ==============================================================================
# Parallel equivalence tests
# ==============================================================================

test_that("wbg_boot_flex parallel matches sequential", {
  skip_on_cran()
  skip_if(parallel::detectCores(logical = FALSE) < 2, "Not enough cores")

  set.seed(123)
  x <- arima.sim(list(ar = 0.6), n = 100)
  stat_fn <- make_stat_ols_t()

  result_seq <- wbg_boot_flex(x, stat_fn, nb = 49, cores = 1, seed = 456, verbose = FALSE)
  result_par <- wbg_boot_flex(x, stat_fn, nb = 49, cores = 2, seed = 456, verbose = FALSE)

  expect_equal(result_seq$obs_stat, result_par$obs_stat)
  expect_equal(result_seq$pvalue, result_par$pvalue)
  expect_equal(result_seq$boot_dist, result_par$boot_dist)
  expect_equal(result_seq$boot_seeds, result_par$boot_seeds)
})
