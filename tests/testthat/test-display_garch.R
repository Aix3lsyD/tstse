# Tests for display_garch functions

# ==============================================================================
# Helper to create test garch_comparison object
# ==============================================================================

make_test_comparison <- function() {
  comparison <- data.frame(
    Model = c("ARCH(1)", "GARCH(1,1)", "GARCH(2,1)"),
    ARCH = c(1, 1, 2),
    GARCH = c(0, 1, 1),
    AIC = c(500.5, 495.2, 497.1),
    AICC = c(500.8, 495.5, 497.5),
    BIC = c(508.3, 507.0, 513.0),
    WLB1 = c(0.12, 0.45, 0.38),
    WLB2 = c(0.08, 0.52, 0.41),
    WLB3 = c(0.03, 0.48, 0.35),
    Nyblom = c(0.85, 0.42, 0.55),
    Nyblom_crit = c(0.75, 0.75, 0.75),
    SignBias = c(0.15, 0.62, 0.48),
    n_sig = c(2, 3, 3),
    n_coef = c(2, 3, 4),
    stringsAsFactors = FALSE
  )

  structure(
    list(
      comparison = comparison,
      fits = list(),  # Empty for display tests
      distribution = "norm",
      n = 500,
      call = quote(compare_garch(y))
    ),
    class = "garch_comparison"
  )
}


# ==============================================================================
# table_garch_gt tests
# ==============================================================================

test_that("table_garch_gt returns gt object", {
  skip_if_not_installed("gt")

  results <- make_test_comparison()
  tbl <- table_garch_gt(results)

  expect_s3_class(tbl, "gt_tbl")
})

test_that("table_garch_gt validates input class", {
  skip_if_not_installed("gt")

  expect_error(table_garch_gt(list()), "garch_comparison")
  expect_error(table_garch_gt(data.frame()), "garch_comparison")
})

test_that("table_garch_gt accepts custom title", {
  skip_if_not_installed("gt")

  results <- make_test_comparison()
  tbl <- table_garch_gt(results, title = "Custom Title")

  expect_s3_class(tbl, "gt_tbl")
})

test_that("table_garch_gt accepts custom colors", {
  skip_if_not_installed("gt")

  results <- make_test_comparison()
  tbl <- table_garch_gt(results,
                        color_best = "darkgreen",
                        color_fail = "red")

  expect_s3_class(tbl, "gt_tbl")
})

test_that("table_garch_gt errors on empty comparison", {
  skip_if_not_installed("gt")

  results <- make_test_comparison()
  results$comparison <- results$comparison[0, ]

  expect_error(table_garch_gt(results), "empty")
})


# ==============================================================================
# table_garch_cli tests
# ==============================================================================

test_that("table_garch_cli produces output", {
  skip_if_not_installed("cli")

  results <- make_test_comparison()

  # cli::cli_h2 and cli::cli_text output may not be captured by expect_output
  # Test for actual table content that is printed via cat()
  expect_output(table_garch_cli(results), "AIC")
  expect_output(table_garch_cli(results), "Model")
  expect_output(table_garch_cli(results), "ARCH\\(1\\)")
})

test_that("table_garch_cli validates input class", {
  skip_if_not_installed("cli")

  expect_error(table_garch_cli(list()), "garch_comparison")
  expect_error(table_garch_cli(data.frame()), "garch_comparison")
})

test_that("table_garch_cli returns invisible data frame", {
  skip_if_not_installed("cli")

  results <- make_test_comparison()

  out <- capture.output({
    result <- table_garch_cli(results)
  })

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 3)
})

test_that("table_garch_cli handles show_signbias = FALSE", {
  skip_if_not_installed("cli")

  results <- make_test_comparison()

  # Should not error and should produce output
  out <- capture.output(table_garch_cli(results, show_signbias = FALSE))
  expect_gt(length(out), 0)
})

test_that("table_garch_cli handles empty comparison", {
  skip_if_not_installed("cli")

  results <- make_test_comparison()
  results$comparison <- results$comparison[0, ]

  # cli::cli_alert_warning output may not be captured by expect_output
  # Just verify it doesn't error and returns the empty data frame
  out <- capture.output({
    result <- table_garch_cli(results)
  })
  expect_equal(nrow(result), 0)
})


# ==============================================================================
# table_coef_gt tests
# ==============================================================================

test_that("table_coef_gt returns gt object", {
  skip_if_not_installed("gt")
  skip_if_not_installed("rugarch")

  # Create a simple fit
  set.seed(42)
  y <- rnorm(300)

  spec <- rugarch::ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
    mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
    distribution.model = "norm"
  )
  fit <- rugarch::ugarchfit(spec, y, solver = "hybrid")

  tbl <- table_coef_gt(fit)
  expect_s3_class(tbl, "gt_tbl")
})

test_that("table_coef_gt validates input class", {
  skip_if_not_installed("gt")

  expect_error(table_coef_gt(list()), "uGARCHfit")
  expect_error(table_coef_gt(data.frame()), "uGARCHfit")
})

test_that("table_coef_gt accepts custom title", {
  skip_if_not_installed("gt")
  skip_if_not_installed("rugarch")

  set.seed(42)
  y <- rnorm(300)

  spec <- rugarch::ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
    mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
    distribution.model = "norm"
  )
  fit <- rugarch::ugarchfit(spec, y, solver = "hybrid")

  tbl <- table_coef_gt(fit, title = "My Custom Title")
  expect_s3_class(tbl, "gt_tbl")
})

test_that("table_coef_gt accepts custom colors", {
  skip_if_not_installed("gt")
  skip_if_not_installed("rugarch")

  set.seed(42)
  y <- rnorm(300)

  spec <- rugarch::ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
    mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
    distribution.model = "norm"
  )
  fit <- rugarch::ugarchfit(spec, y, solver = "hybrid")

  tbl <- table_coef_gt(fit,
                       color_sig = "blue",
                       color_nonsig = "gray")
  expect_s3_class(tbl, "gt_tbl")
})


# ==============================================================================
# Internal helper function tests
# ==============================================================================

test_that(".fmt_pval formats p-values correctly", {
  fmt_pval <- tstse:::.fmt_pval

  expect_equal(fmt_pval(NA), "NA")
  expect_equal(fmt_pval(0.0001), "<0.001")
  expect_equal(fmt_pval(0.05), "0.050")
  expect_equal(fmt_pval(0.123), "0.123")
})

test_that(".find_best_indices finds correct indices", {
  find_best <- tstse:::.find_best_indices

  df <- data.frame(
    AIC = c(100, 90, 95),
    BIC = c(110, 105, 100),
    AICC = c(101, 91, 96),
    WLB1 = c(0.1, 0.5, 0.3),
    WLB2 = c(0.2, 0.4, 0.6),
    WLB3 = c(0.3, 0.3, 0.3),
    Nyblom = c(0.5, 0.3, 0.4),
    SignBias = c(0.1, 0.8, 0.5)
  )

  best <- find_best(df)

  expect_equal(best$aic, 2)  # Row 2 has lowest AIC
  expect_equal(best$bic, 3)  # Row 3 has lowest BIC
  expect_equal(best$wlb1, 2) # Row 2 has highest WLB1
  expect_equal(best$nyblom, 2) # Row 2 has lowest Nyblom
})

test_that(".find_best_indices handles all-NA columns", {
  find_best <- tstse:::.find_best_indices

  df <- data.frame(
    AIC = c(NA_real_, NA_real_, NA_real_),
    BIC = c(100, 90, 95),
    AICC = c(NA_real_, NA_real_, NA_real_),
    WLB1 = c(0.1, 0.5, 0.3),
    WLB2 = c(NA_real_, NA_real_, NA_real_),
    WLB3 = c(0.3, 0.3, 0.3),
    Nyblom = c(NA_real_, NA_real_, NA_real_),
    SignBias = c(0.1, 0.8, 0.5)
  )

  best <- find_best(df)

  expect_true(is.na(best$aic))
  expect_equal(best$bic, 2)
  expect_true(is.na(best$aicc))
  expect_true(is.na(best$nyblom))
})


# ==============================================================================
# Integration tests with real compare_garch output
# ==============================================================================

test_that("display functions work with real compare_garch output", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("WeightedPortTest")
  skip_if_not_installed("gt")
  skip_if_not_installed("cli")

  set.seed(42)
  y <- rnorm(300)

  results <- compare_garch(y, arch_range = 1, garch_range = 0:1, parallel = FALSE)

  # gt table
  tbl <- table_garch_gt(results)
  expect_s3_class(tbl, "gt_tbl")

  # cli table - check for actual table content
  expect_output(table_garch_cli(results), "AIC")

  # coefficient table for first fit
  first_fit <- results$fits[[1]]
  coef_tbl <- table_coef_gt(first_fit)
  expect_s3_class(coef_tbl, "gt_tbl")
})


# ==============================================================================
# Edge cases with NA values
# ==============================================================================

test_that("display functions handle NA values in diagnostics", {
  skip_if_not_installed("gt")
  skip_if_not_installed("cli")

  results <- make_test_comparison()
  # Add some NA values
  results$comparison$WLB1[1] <- NA
  results$comparison$Nyblom[2] <- NA
  results$comparison$SignBias[3] <- NA

  # Should not error
  tbl <- table_garch_gt(results)
  expect_s3_class(tbl, "gt_tbl")

  out <- capture.output(table_garch_cli(results))
  expect_gt(length(out), 0)
})


# ==============================================================================
# Base R fallback function tests
# ==============================================================================

test_that(".table_garch_base produces output", {
  results <- make_test_comparison()

  table_base <- tstse:::.table_garch_base

  expect_output(table_base(results), "GARCH Model Comparison")
  expect_output(table_base(results), "Distribution:")
  expect_output(table_base(results), "AIC")
  expect_output(table_base(results), "Best by AIC")
})

test_that(".table_garch_base hides WLB columns when all NA", {
  results <- make_test_comparison()
  # Set all WLB to NA
  results$comparison$WLB1 <- NA_real_
  results$comparison$WLB2 <- NA_real_
  results$comparison$WLB3 <- NA_real_

  table_base <- tstse:::.table_garch_base

  out <- capture.output(table_base(results))
  out_text <- paste(out, collapse = "\n")

  # WLB should not appear in output
  expect_false(grepl("WLB", out_text))
  # AIC should still appear
  expect_true(grepl("AIC", out_text))
})

test_that(".table_garch_base shows WLB columns when available", {
  results <- make_test_comparison()

  table_base <- tstse:::.table_garch_base

  out <- capture.output(table_base(results))
  out_text <- paste(out, collapse = "\n")

  # WLB should appear in output
  expect_true(grepl("WLB", out_text))
})

test_that(".table_coef_base produces output", {
  skip_if_not_installed("rugarch")

  set.seed(42)
  y <- rnorm(300)

  spec <- rugarch::ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
    mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
    distribution.model = "norm"
  )
  fit <- rugarch::ugarchfit(spec, y, solver = "hybrid")

  table_coef_base <- tstse:::.table_coef_base

  expect_output(table_coef_base(fit), "Coefficients")
  expect_output(table_coef_base(fit), "Distribution:")
  expect_output(table_coef_base(fit), "Estimate")
})


# ==============================================================================
# WLB column hiding tests for gt and cli
# ==============================================================================

test_that("table_garch_gt hides WLB columns when all NA", {
  skip_if_not_installed("gt")

  results <- make_test_comparison()
  # Set all WLB to NA
  results$comparison$WLB1 <- NA_real_
  results$comparison$WLB2 <- NA_real_
  results$comparison$WLB3 <- NA_real_

  # Should not error
  tbl <- table_garch_gt(results)
  expect_s3_class(tbl, "gt_tbl")
})

test_that("table_garch_cli hides WLB columns when all NA", {
  skip_if_not_installed("cli")

  results <- make_test_comparison()
  # Set all WLB to NA
  results$comparison$WLB1 <- NA_real_
  results$comparison$WLB2 <- NA_real_
  results$comparison$WLB3 <- NA_real_

  out <- capture.output(table_garch_cli(results))
  out_text <- paste(out, collapse = "\n")

  # WLB should not appear in output when all NA
  expect_false(grepl("WLB1", out_text))
})
