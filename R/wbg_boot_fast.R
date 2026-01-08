#' Fast WBG Bootstrap Test for Trend (C++ Implementation)
#'
#' High-performance bootstrap hypothesis test for trend in time series.
#' Uses the Woodward-Bottone-Gray (WBG) method with C++ parallelization
#' via TBB (Intel Threading Building Blocks).
#'
#' @param x A numeric vector containing the time series data.
#' @param nb Number of bootstrap replicates. Default is 399.
#' @param maxp Maximum AR order for Cochrane-Orcutt model selection. Default is 5.
#' @param criterion Character. Information criterion for order selection:
#'   `"aic"`, `"aicc"`, or `"bic"`. Default is `"aic"`.
#' @param bootadj Logical. If `TRUE`, performs the two-stage bootstrap
#'   adjustment (variance correction) for improved small-sample performance.
#'   Default is `FALSE` for maximum speed.
#' @param pw_stat Logical. If `TRUE`, also computes Prais-Winsten t-statistics
#'   alongside the Cochrane-Orcutt statistics. Default is `FALSE`.
#' @param maxp_pw Maximum AR order for Prais-Winsten (typically 1 for AR(1)).
#'   Default is 1. Only used when `pw_stat = TRUE`.
#' @param seed Integer, random seed for reproducibility. Default NULL.
#'
#' @return A list of class "wbg_boot_fast" containing:
#'   \item{p}{AR order selected for the series under null hypothesis.}
#'   \item{phi}{AR coefficients for the null model.}
#'   \item{vara}{Innovation variance for the null model.}
#'   \item{pvalue}{Bootstrap p-value for trend test (two-sided).}
#'   \item{tco_obs}{Observed t-statistic from Cochrane-Orcutt fit.}
#'   \item{boot_tstats}{Vector of bootstrap CO t-statistics.}
#'   \item{n}{Length of input series.}
#'   \item{nb}{Number of bootstrap replicates used.}
#'   \item{maxp}{Maximum AR order considered for CO.}
#'   \item{criterion}{Information criterion used.}
#'
#'   If `bootadj = TRUE`, also includes:
#'   \item{tco_obs_adj}{Bootstrap-adjusted observed CO t-statistic.}
#'   \item{pvalue_adj}{Bootstrap-adjusted two-sided CO p-value.}
#'   \item{adj_factor}{CO variance adjustment factor (C from paper).}
#'   \item{median_phi}{AR coefficients from the median bootstrap model.}
#'
#'   If `pw_stat = TRUE`, also includes:
#'   \item{tpw_obs}{Observed t-statistic from Prais-Winsten fit.}
#'   \item{pvalue_pw}{Bootstrap p-value for PW trend test (two-sided).}
#'   \item{rho_null}{AR(1) coefficient for PW null model.}
#'   \item{boot_tstats_pw}{Vector of bootstrap PW t-statistics.}
#'
#'   If both `pw_stat = TRUE` and `bootadj = TRUE`:
#'   \item{tpw_obs_adj}{Bootstrap-adjusted observed PW t-statistic.}
#'   \item{pvalue_pw_adj}{Bootstrap-adjusted two-sided PW p-value.}
#'   \item{adj_factor_pw}{PW variance adjustment factor.}
#'   \item{median_rho}{Median AR(1) coefficient from bootstrap.}
#'
#' @details
#' This is a high-performance C++ implementation of [wbg_boot()],
#' optimized for large-scale applications. Key differences:
#'
#' \itemize{
#'   \item Uses Burg algorithm only (most stable for bootstrap)
#'   \item Parallelization via TBB (works on all platforms)
#'   \item All bootstrap iterations run in C++ (no R overhead)
#'   \item Thread-safe RNG (dqrng with xoshiro256+)
#' }
#'
#' Expected speedup is 10-25x compared to [wbg_boot()] depending on
#' series length and number of available CPU cores.
#'
#' The bootstrap procedure:
#' \enumerate{
#'   \item Fit Cochrane-Orcutt model to observed data using C++
#'   \item Fit AR(p) model under null hypothesis (no trend)
#'   \item Generate `nb` bootstrap series from AR(p) in parallel
#'   \item Compute CO t-statistic for each bootstrap series
#'   \item p-value = proportion of |bootstrap t| >= |observed t|
#' }
#'
#' **Bootstrap Adjustment:**
#' When `bootadj = TRUE`, a second-stage bootstrap adjusts the observed
#' statistic for variance inflation caused by AR coefficient estimation
#' bias. This provides better small-sample significance levels but
#' approximately doubles computation time. The adjustment uses the median
#' AR coefficients from the first stage, where the "median" is selected
#' based on phi(1) = 1 - sum(phi) as described in the original paper.
#'
#' **Performance Tuning:**
#' For optimal parallel performance, you may want to disable BLAS multi-threading
#' to avoid nested parallelism (TBB threads + BLAS threads). Set these environment
#' variables **before loading any packages** that use BLAS:
#'
#' ```
#' Sys.setenv(
#'   VECLIB_MAXIMUM_THREADS = "1",
#'   OPENBLAS_NUM_THREADS = "1",
#'   MKL_NUM_THREADS = "1",
#'   OMP_NUM_THREADS = "1",
#'   BLAS_NUM_THREADS = "1"
#' )
#' ```
#'
#' On macOS with Apple Accelerate, `VECLIB_MAXIMUM_THREADS` is the key variable.
#' Note: These must be set early in the R session (before BLAS initializes) to
#' take effect. A session restart may be required after changing them.
#'
#' @references
#' Woodward, W. A., Bottone, S., and Gray, H. L. (1997). "Improved Tests for
#' Trend in Time Series Data." *Journal of Agricultural, Biological, and
#' Environmental Statistics*, 2(4), 403-416.
#'
#' @seealso [wbg_boot()] for the original R implementation,
#'   [wbg_boot_flex()] for flexible statistics with bootstrap adjustment,
#'   [co()] for Cochrane-Orcutt estimation,
#'   [pw_fast()] for Prais-Winsten estimation.
#'
#' @examples
#' \donttest{
#' # Test series with no trend (should be non-significant)
#' set.seed(123)
#' x <- arima.sim(list(ar = 0.7), n = 100)
#' result <- wbg_boot_fast(x, nb = 199, seed = 456)
#' print(result)
#'
#' # Test series with trend (should be significant)
#' set.seed(123)
#' x_trend <- arima.sim(list(ar = 0.7), n = 100) + 0.1 * (1:100)
#' result <- wbg_boot_fast(x_trend, nb = 199, seed = 456)
#' print(result)
#'
#' # With COBA adjustment for better small-sample performance
#' result_adj <- wbg_boot_fast(x_trend, nb = 199, bootadj = TRUE, seed = 456)
#' print(result_adj)
#'
#' # With both CO and PW statistics
#' result_both <- wbg_boot_fast(x_trend, nb = 199, pw_stat = TRUE, seed = 456)
#' print(result_both)
#'
#' # With PW and COBA adjustment
#' result_full <- wbg_boot_fast(x_trend, nb = 199, pw_stat = TRUE,
#'                               bootadj = TRUE, seed = 456)
#' print(result_full)
#' }
#'
#' \dontrun{
#' # Benchmark comparison
#' library(microbenchmark)
#' x <- arima.sim(list(ar = c(0.7, -0.2)), n = 500)
#' microbenchmark(
#'   fast = wbg_boot_fast(x, nb = 199, seed = 42),
#'   fast_coba = wbg_boot_fast(x, nb = 199, bootadj = TRUE, seed = 42),
#'   original = wbg_boot(x, nb = 199, seed = 42),
#'   times = 5
#' )
#' }
#'
#' @export
wbg_boot_fast <- function(x, nb = 399L, maxp = 5L,
                          criterion = c("aic", "aicc", "bic"),
                          bootadj = FALSE,
                          pw_stat = FALSE,
                          maxp_pw = 1L,
                          seed = NULL) {

  criterion <- match.arg(criterion)
  maxp_pw <- as.integer(maxp_pw)

 # Input validation
  if (!is.numeric(x) || length(x) == 0) {
    stop("`x` must be a non-empty numeric vector", call. = FALSE)
  }
  if (!is.numeric(nb) || length(nb) != 1 || nb < 1) {
    stop("`nb` must be a positive integer", call. = FALSE)
  }
  if (!is.numeric(maxp) || length(maxp) != 1 || maxp < 1) {
    stop("`maxp` must be a positive integer", call. = FALSE)
  }
  if (maxp > 20L) {
    stop("`maxp` must be <= 20 (C++ buffer limit)", call. = FALSE)
  }

  x <- as.numeric(x)
  nb <- as.integer(nb)
  maxp <- as.integer(maxp)
  n <- length(x)

  if (n <= maxp + 10) {
    stop("Series too short: need n > maxp + 10", call. = FALSE)
  }

  # Set seed for reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Generate bootstrap seeds upfront (R's RNG, single-threaded)
  boot_seeds <- as.numeric(sample.int(.Machine$integer.max, nb))
  if (bootadj) {
    boot_seeds_adj <- as.numeric(sample.int(.Machine$integer.max, nb))
  }

  # Step 1: Compute observed CO t-statistic (C++)
  tco_obs <- co_tstat_cpp(x, maxp, criterion)

  # Also compute observed PW t-statistic if requested
  if (pw_stat) {
    tpw_obs <- pw_tstat_cpp(x)
  }

  # Step 2: Fit null model (no trend) using C++ Burg with user's criterion
  # Using burg_aic_select_cpp ensures the same IC is used throughout
  fit <- burg_aic_select_cpp(x, maxp, criterion)
  p_null <- fit$p
  phi_null <- as.numeric(fit$phi)
  vara_null <- fit$vara

  # Handle AR(0) case
  if (p_null == 0) {
    phi_null <- numeric(0)
  }

  # For PW null model, estimate AR(1) rho and variance from OLS-detrended residuals
  if (pw_stat) {
    z_x <- ols_detrend_cpp(x)
    rho_null <- sum(z_x[-1] * z_x[-n]) / sum(z_x^2)
    rho_null <- max(-0.999, min(0.999, rho_null))
    # Innovation variance for AR(1): var(e) = var(z) * (1 - rho^2)
    vara_pw_null <- var(z_x) * (1 - rho_null^2)
  }

  # Compute optimal grain size for TBB parallelization

  # Heuristic: ~4 chunks per thread for good load balancing
  num_threads <- RcppParallel::defaultNumThreads()
  grain_size <- max(16L, as.integer(nb / (4L * num_threads)))

  # Step 3: Bootstrap
  if (pw_stat && bootadj) {
    # === SEPARATE CO + PW COBA PATHS ===
    # CO uses AR(p) null, PW uses AR(1) null

    # CO bootstrap with COBA data collection
    boot_result_co <- wbg_bootstrap_coba_kernel_grain_cpp(
      n = n,
      phi = phi_null,
      vara = vara_null,
      seeds = boot_seeds,
      maxp = maxp,
      criterion = criterion,
      grain_size = grain_size
    )
    boot_tstats <- boot_result_co$tstats
    phi1_values <- boot_result_co$phi1_values
    phi_matrix <- boot_result_co$phi_matrix
    orders <- boot_result_co$orders

    # PW bootstrap with COBA data collection (generates from AR(1))
    boot_result_pw <- wbg_bootstrap_pw_coba_kernel_cpp(
      n = n,
      rho = rho_null,
      vara = vara_pw_null,
      seeds = boot_seeds,
      grain_size = grain_size
    )
    boot_tstats_pw <- boot_result_pw$pw_tstats
    rho_values <- boot_result_pw$rho_values

    # Find median AR model for CO (per paper: model with phi(1) closest to median)
    median_idx <- which.min(abs(phi1_values - median(phi1_values)))
    median_order <- orders[median_idx]
    median_phi <- if (median_order > 0) {
      phi_matrix[median_idx, seq_len(median_order)]
    } else {
      numeric(0)
    }

    # Find median rho for PW
    median_rho <- median(rho_values)

    # Second CO bootstrap using median AR(p) model
    boot_tstats_co_adj <- wbg_bootstrap_kernel_grain_cpp(
      n = n,
      phi = median_phi,
      vara = vara_null,
      seeds = boot_seeds_adj,
      maxp = maxp,
      criterion = criterion,
      grain_size = grain_size
    )

    # Second PW bootstrap using median AR(1) rho
    boot_tstats_pw_adj <- wbg_bootstrap_pw_kernel_cpp(
      n = n,
      rho = median_rho,
      vara = vara_pw_null,
      seeds = boot_seeds_adj,
      grain_size = grain_size
    )

    # CO bootstrap adjustment: Ĉ = sd(t̃*) / sd(t*) per Woodward 1997 Section 2.2
    # Ĉ < 1 for positive autocorrelation, shrinking t_obs to reduce false positives
    adj_factor <- sd(boot_tstats_co_adj) / sd(boot_tstats)
    tco_obs_adj <- adj_factor * tco_obs
    pvalue_adj <- mean(abs(boot_tstats) >= abs(tco_obs_adj))

    # PW bootstrap adjustment: Ĉ = sd(t̃*) / sd(t*) per Woodward 1997 Section 2.2
    # Ĉ < 1 for positive autocorrelation, shrinking t_obs to reduce false positives
    adj_factor_pw <- sd(boot_tstats_pw_adj) / sd(boot_tstats_pw)
    tpw_obs_adj <- adj_factor_pw * tpw_obs
    pvalue_pw_adj <- mean(abs(boot_tstats_pw) >= abs(tpw_obs_adj))

  } else if (pw_stat && !bootadj) {
    # === SEPARATE CO + PW FAST PATHS ===
    # CO uses AR(p) null, PW uses AR(1) null

    # CO bootstrap
    boot_tstats <- wbg_bootstrap_kernel_grain_cpp(
      n = n,
      phi = phi_null,
      vara = vara_null,
      seeds = boot_seeds,
      maxp = maxp,
      criterion = criterion,
      grain_size = grain_size
    )

    # PW bootstrap (generates from AR(1))
    boot_tstats_pw <- wbg_bootstrap_pw_kernel_cpp(
      n = n,
      rho = rho_null,
      vara = vara_pw_null,
      seeds = boot_seeds,
      grain_size = grain_size
    )

  } else if (bootadj) {
    # === CO-ONLY COBA PATH ===
    boot_result <- wbg_bootstrap_coba_kernel_grain_cpp(
      n = n,
      phi = phi_null,
      vara = vara_null,
      seeds = boot_seeds,
      maxp = maxp,
      criterion = criterion,
      grain_size = grain_size
    )
    boot_tstats <- boot_result$tstats
    phi1_values <- boot_result$phi1_values
    phi_matrix <- boot_result$phi_matrix
    orders <- boot_result$orders

    # Find median AR model (per paper: model with phi(1) closest to median)
    median_idx <- which.min(abs(phi1_values - median(phi1_values)))
    median_order <- orders[median_idx]
    median_phi <- if (median_order > 0) {
      phi_matrix[median_idx, seq_len(median_order)]
    } else {
      numeric(0)
    }

    # Second bootstrap using median phi (with grain size tuning)
    boot_tstats_adj <- wbg_bootstrap_kernel_grain_cpp(
      n = n,
      phi = median_phi,
      vara = vara_null,
      seeds = boot_seeds_adj,
      maxp = maxp,
      criterion = criterion,
      grain_size = grain_size
    )

    # CO bootstrap adjustment: Ĉ = sd(t̃*) / sd(t*) per Woodward 1997 Section 2.2
    # Ĉ < 1 for positive autocorrelation, shrinking t_obs to reduce false positives
    adj_factor <- sd(boot_tstats_adj) / sd(boot_tstats)
    tco_obs_adj <- adj_factor * tco_obs
    pvalue_adj <- mean(abs(boot_tstats) >= abs(tco_obs_adj))

  } else {
    # === CO-ONLY FAST PATH (no AR fitting, no PW) ===
    boot_tstats <- wbg_bootstrap_kernel_grain_cpp(
      n = n,
      phi = phi_null,
      vara = vara_null,
      seeds = boot_seeds,
      maxp = maxp,
      criterion = criterion,
      grain_size = grain_size
    )
  }

  # Step 4: Compute p-values (two-sided)
  pvalue <- mean(abs(boot_tstats) >= abs(tco_obs))
  if (pw_stat) {
    pvalue_pw <- mean(abs(boot_tstats_pw) >= abs(tpw_obs))
  }

  # Build result
  result <- list(
    p           = p_null,
    phi         = phi_null,
    vara        = vara_null,
    pvalue      = pvalue,
    tco_obs     = tco_obs,
    boot_tstats = boot_tstats,
    n           = n,
    nb          = nb,
    maxp        = maxp,
    criterion   = criterion
  )

  # Add COBA results if computed
  if (bootadj) {
    result$tco_obs_adj <- tco_obs_adj
    result$pvalue_adj <- pvalue_adj
    result$adj_factor <- adj_factor
    result$median_phi <- median_phi
  }

  # Add PW results if computed
  if (pw_stat) {
    result$tpw_obs <- tpw_obs
    result$pvalue_pw <- pvalue_pw
    result$rho_null <- rho_null
    result$vara_pw_null <- vara_pw_null
    result$boot_tstats_pw <- boot_tstats_pw
    result$maxp_pw <- maxp_pw

    # Add PW COBA results if computed
    if (bootadj) {
      result$tpw_obs_adj <- tpw_obs_adj
      result$pvalue_pw_adj <- pvalue_pw_adj
      result$adj_factor_pw <- adj_factor_pw
      result$median_rho <- median_rho
    }
  }

  structure(result, class = "wbg_boot_fast")
}


#' @export
print.wbg_boot_fast <- function(x, ...) {
  cat("\nWBG Bootstrap Test for Trend (Fast C++ Implementation)\n")
  cat("-------------------------------------------------------\n")
  cat(sprintf("  Series length: %d\n", x$n))
  cat(sprintf("  Bootstrap replicates: %d\n", x$nb))
  cat(sprintf("  Max AR order (CO): %d\n", x$maxp))
  if (!is.null(x$maxp_pw)) {
    cat(sprintf("  Max AR order (PW): %d\n", x$maxp_pw))
  }
  cat(sprintf("  IC criterion: %s\n", x$criterion))
  cat("\nNull Model (no trend):\n")
  cat(sprintf("  Selected AR order: %d\n", x$p))
  if (x$p > 0) {
    cat(sprintf("  AR coefficients: %s\n",
                paste(round(x$phi, 4), collapse = ", ")))
  }
  cat(sprintf("  Innovation variance: %.4f\n", x$vara))
  if (!is.null(x$rho_null)) {
    cat(sprintf("  PW AR(1) rho: %.4f\n", x$rho_null))
  }

  cat("\nCochrane-Orcutt Results:\n")
  cat(sprintf("  Observed CO t-stat: %.4f\n", x$tco_obs))
  cat(sprintf("  Bootstrap p-value: %.4f\n", x$pvalue))

  # CO bootstrap adjustment results if present
  if (!is.null(x$pvalue_adj)) {
    cat("  Bootstrap Adjustment:\n")
    cat(sprintf("    Adjustment factor: %.4f\n", x$adj_factor))
    cat(sprintf("    Adjusted t-stat: %.4f\n", x$tco_obs_adj))
    cat(sprintf("    Adjusted p-value: %.4f\n", x$pvalue_adj))
    if (length(x$median_phi) > 0) {
      cat(sprintf("    Median AR(%d) phi: %s\n", length(x$median_phi),
                  paste(round(x$median_phi, 4), collapse = ", ")))
    }
  }

  # PW results if present
  if (!is.null(x$tpw_obs)) {
    cat("\nPrais-Winsten Results:\n")
    cat(sprintf("  Observed PW t-stat: %.4f\n", x$tpw_obs))
    cat(sprintf("  Bootstrap p-value: %.4f\n", x$pvalue_pw))

    # PW bootstrap adjustment results if present
    if (!is.null(x$pvalue_pw_adj)) {
      cat("  Bootstrap Adjustment:\n")
      cat(sprintf("    Adjustment factor: %.4f\n", x$adj_factor_pw))
      cat(sprintf("    Adjusted t-stat: %.4f\n", x$tpw_obs_adj))
      cat(sprintf("    Adjusted p-value: %.4f\n", x$pvalue_pw_adj))
      cat(sprintf("    Median rho: %.4f\n", x$median_rho))
    }
  }
  cat("\n")
  invisible(x)
}


#' @export
summary.wbg_boot_fast <- function(object, ...) {
  cat("\nWBG Bootstrap Test Summary\n")
  cat("==========================\n\n")

  cat("Data:\n")
  cat(sprintf("  n = %d, nb = %d, maxp (CO) = %d, criterion = %s\n",
              object$n, object$nb, object$maxp, object$criterion))
  if (!is.null(object$maxp_pw)) {
    cat(sprintf("  maxp (PW) = %d\n", object$maxp_pw))
  }
  cat("\n")

  cat("Null Hypothesis: No trend (slope = 0)\n")
  cat(sprintf("  AR(%d) model with variance = %.4f\n", object$p, object$vara))
  if (!is.null(object$rho_null)) {
    cat(sprintf("  PW AR(1) rho = %.4f\n", object$rho_null))
  }
  cat("\n")

  # CO results
  cat("Cochrane-Orcutt Test Statistic:\n")
  cat(sprintf("  Observed CO t = %.4f\n", object$tco_obs))
  cat(sprintf("  Bootstrap SE = %.4f\n", sd(object$boot_tstats)))
  cat(sprintf("  Bootstrap mean = %.4f\n\n", mean(object$boot_tstats)))

  cat("CO Result:\n")
  sig_level <- if (object$pvalue < 0.001) "***"
               else if (object$pvalue < 0.01) "**"
               else if (object$pvalue < 0.05) "*"
               else if (object$pvalue < 0.1) "."
               else ""
  cat(sprintf("  p-value = %.4f %s\n", object$pvalue, sig_level))
  cat("  Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")

  # CO bootstrap adjustment results if present
  if (!is.null(object$pvalue_adj)) {
    cat("CO Bootstrap Adjustment (variance correction):\n")
    cat(sprintf("  Adjustment factor C = %.4f\n", object$adj_factor))
    cat(sprintf("  Adjusted t-stat = %.4f\n", object$tco_obs_adj))
    sig_level_adj <- if (object$pvalue_adj < 0.001) "***"
                     else if (object$pvalue_adj < 0.01) "**"
                     else if (object$pvalue_adj < 0.05) "*"
                     else if (object$pvalue_adj < 0.1) "."
                     else ""
    cat(sprintf("  Adjusted p-value = %.4f %s\n", object$pvalue_adj, sig_level_adj))
    if (length(object$median_phi) > 0) {
      cat(sprintf("  Median model: AR(%d) with phi = %s\n",
                  length(object$median_phi),
                  paste(round(object$median_phi, 4), collapse = ", ")))
    } else {
      cat("  Median model: AR(0)\n")
    }
    cat("\n")
  }

  # PW results if present
  if (!is.null(object$tpw_obs)) {
    cat("Prais-Winsten Test Statistic:\n")
    cat(sprintf("  Observed PW t = %.4f\n", object$tpw_obs))
    cat(sprintf("  Bootstrap SE = %.4f\n", sd(object$boot_tstats_pw)))
    cat(sprintf("  Bootstrap mean = %.4f\n\n", mean(object$boot_tstats_pw)))

    cat("PW Result:\n")
    sig_level_pw <- if (object$pvalue_pw < 0.001) "***"
                    else if (object$pvalue_pw < 0.01) "**"
                    else if (object$pvalue_pw < 0.05) "*"
                    else if (object$pvalue_pw < 0.1) "."
                    else ""
    cat(sprintf("  p-value = %.4f %s\n", object$pvalue_pw, sig_level_pw))
    cat("  Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")

    # PW bootstrap adjustment results if present
    if (!is.null(object$pvalue_pw_adj)) {
      cat("PW Bootstrap Adjustment (variance correction):\n")
      cat(sprintf("  Adjustment factor C = %.4f\n", object$adj_factor_pw))
      cat(sprintf("  Adjusted t-stat = %.4f\n", object$tpw_obs_adj))
      sig_level_pw_adj <- if (object$pvalue_pw_adj < 0.001) "***"
                          else if (object$pvalue_pw_adj < 0.01) "**"
                          else if (object$pvalue_pw_adj < 0.05) "*"
                          else if (object$pvalue_pw_adj < 0.1) "."
                          else ""
      cat(sprintf("  Adjusted p-value = %.4f %s\n", object$pvalue_pw_adj, sig_level_pw_adj))
      cat(sprintf("  Median rho = %.4f\n", object$median_rho))
      cat("\n")
    }
  }

  invisible(object)
}
