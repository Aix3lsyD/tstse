#' Fast WBG Bootstrap Test for Trend (C++ Implementation)
#'
#' High-performance bootstrap hypothesis test for trend in time series.
#' Uses the Woodward-Bottone-Gray (WBG) method with C++ parallelization
#' via TBB (Intel Threading Building Blocks).
#'
#' @param x A numeric vector containing the time series data.
#' @param nb Number of bootstrap replicates. Default is 399.
#' @param maxp Maximum AR order for model selection. Default is 5.
#' @param criterion Character. Information criterion for order selection:
#'   `"aic"`, `"aicc"`, or `"bic"`. Default is `"aic"`.
#' @param bootadj Logical. If `TRUE`, performs the COBA (Cochrane-Orcutt
#'   Bootstrap Adjustment) variance adjustment for improved small-sample
#'   performance. Default is `FALSE` for maximum speed.
#' @param seed Integer, random seed for reproducibility. Default NULL.
#'
#' @return A list of class "wbg_boot_fast" containing:
#'   \item{p}{AR order selected for the series under null hypothesis.}
#'   \item{phi}{AR coefficients for the null model.}
#'   \item{vara}{Innovation variance for the null model.}
#'   \item{pvalue}{Bootstrap p-value for trend test (two-sided, abs-based).}
#'   \item{pvalue_upper}{Upper-tail p-value: P(T* >= T_obs). Use for testing
#'     positive trend or for quantile-based two-sided test.}
#'   \item{pvalue_lower}{Lower-tail p-value: P(T* <= T_obs). Use for testing
#'     negative trend or for quantile-based two-sided test.}
#'   \item{tco_obs}{Observed t-statistic from Cochrane-Orcutt fit.}
#'   \item{boot_tstats}{Vector of bootstrap t-statistics.}
#'   \item{boot_seeds}{Vector of RNG seeds used for each bootstrap replicate.
#'     Useful for auditability and exact reproduction of individual replicates.}
#'   \item{n}{Length of input series.}
#'   \item{nb}{Number of bootstrap replicates used.}
#'   \item{maxp}{Maximum AR order considered.}
#'   \item{criterion}{Information criterion used.}
#'
#'   If `bootadj = TRUE`, also includes:
#'   \item{tco_obs_adj}{COBA-adjusted observed t-statistic.}
#'   \item{pvalue_adj}{COBA-adjusted two-sided p-value.}
#'   \item{adj_factor}{Variance adjustment factor (C from paper).}
#'   \item{median_phi}{AR coefficients from the median bootstrap model.}
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
#' **COBA Adjustment:**
#' When `bootadj = TRUE`, a second-stage bootstrap adjusts the observed
#' statistic for variance inflation caused by AR coefficient estimation
#' bias. This provides better small-sample significance levels but
#' approximately doubles computation time. The adjustment uses the median
#' AR coefficients from the first stage, where the "median" is selected
#' based on phi(1) = 1 - sum(phi) as described in the original paper.
#' **Note:** COBA requires storing AR fits for all bootstrap replicates,
#' using approximately `nb * maxp * 8` bytes (e.g., 16 MB for nb=100,000,
#' maxp=20). For large-scale simulation studies, consider memory constraints.
#'
#' **P-value Methods:**
#' Three p-values are returned. The main `pvalue` uses the absolute value
#' method `P(|T*| >= |T_obs|)`, which is standard for two-sided bootstrap
#' tests. Additionally, `pvalue_upper` and `pvalue_lower` provide one-sided
#' p-values for directional hypotheses or for computing the quantile-based
#' two-sided p-value: `2 * min(pvalue_upper, pvalue_lower)`. The quantile-based
#' method can differ from the abs-based method when the bootstrap distribution
#' is asymmetric.
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
#' **Reproducibility:**
#' Results are exactly reproducible given the same `seed`, platform, compiler,
#' and math library. Minor differences in bootstrap t-statistics may occur
#' across different platforms (e.g., Windows vs Linux vs macOS) due to
#' implementation-defined behavior in C++ uniform distribution generation.
#' This does not affect statistical validity - p-values and conclusions will
#' be equivalent within sampling variability.
#'
#' @references
#' Woodward, W. A., Bottone, S., and Gray, H. L. (1997). "Improved Tests for
#' Trend in Time Series Data." *Journal of Agricultural, Biological, and
#' Environmental Statistics*, 2(4), 403-416.
#'
#' @seealso [wbg_boot()] for the original R implementation,
#'   [wbg_boot_flex()] for flexible statistics with COBA,
#'   [co()] for Cochrane-Orcutt estimation.
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
                          seed = NULL) {

  criterion <- match.arg(criterion)

 # Input validation
  if (!is.numeric(x) || length(x) == 0) {
    stop("`x` must be a non-empty numeric vector", call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop("`x` must contain only finite values (no NA, NaN, or Inf)", call. = FALSE)
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


  # Step 1: Compute observed CO t-statistic using kernel implementation

  # Uses same algorithm as bootstrap kernel for bootstrap validity
  tco_obs <- co_tstat_pure_cpp(x, maxp, criterion)

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

  # Compute optimal grain size for TBB parallelization

  # Heuristic: ~4 chunks per thread for good load balancing
  num_threads <- RcppParallel::defaultNumThreads()
  grain_size <- max(16L, as.integer(nb / (4L * num_threads)))

  # Step 3: Bootstrap
  if (bootadj) {
    # === COBA PATH: First bootstrap WITH AR fitting ===
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

    # Compute adjustment factor (C = sigma_t_tilde* / sigma_t* from paper)
    adj_factor <- sd(boot_tstats_adj) / sd(boot_tstats)
    tco_obs_adj <- adj_factor * tco_obs
    pvalue_adj <- (sum(abs(boot_tstats) >= abs(tco_obs_adj)) + 1) / (nb + 1)

  } else {
    # === STANDARD PATH: No extra AR refitting for COBA median selection ===
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

  # Step 4: Compute p-values with plus-one correction
  # Per Davison & Hinkley (1997): ensures p > 0 and treats observed as "one of" the samples

  # Two-sided (abs-based): P(|T*| >= |T_obs|)
  pvalue <- (sum(abs(boot_tstats) >= abs(tco_obs)) + 1) / (nb + 1)

  # One-sided p-values for directional tests or quantile-based two-sided
  pvalue_upper <- (sum(boot_tstats >= tco_obs) + 1) / (nb + 1)  # P(T* >= T_obs)
  pvalue_lower <- (sum(boot_tstats <= tco_obs) + 1) / (nb + 1)  # P(T* <= T_obs)

  # Build result
  result <- list(
    p           = p_null,
    phi         = phi_null,
    vara        = vara_null,
    pvalue      = pvalue,
    pvalue_upper = pvalue_upper,
    pvalue_lower = pvalue_lower,
    tco_obs     = tco_obs,
    boot_tstats = boot_tstats,
    boot_seeds  = boot_seeds,
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

  structure(result, class = "wbg_boot_fast")
}


#' @export
print.wbg_boot_fast <- function(x, ...) {
  cat("\nWBG Bootstrap Test for Trend (Fast C++ Implementation)\n")
  cat("-------------------------------------------------------\n")
  cat(sprintf("  Series length: %d\n", x$n))
  cat(sprintf("  Bootstrap replicates: %d\n", x$nb))
  cat(sprintf("  Max AR order: %d\n", x$maxp))
  cat(sprintf("  IC criterion: %s\n", x$criterion))
  cat("\nNull Model (no trend):\n")
  cat(sprintf("  Selected AR order: %d\n", x$p))
  if (x$p > 0) {
    cat(sprintf("  AR coefficients: %s\n",
                paste(round(x$phi, 4), collapse = ", ")))
  }
  cat(sprintf("  Innovation variance: %.4f\n", x$vara))
  cat("\nTest Results:\n")
  cat(sprintf("  Observed CO t-stat: %.4f\n", x$tco_obs))
  cat(sprintf("  Bootstrap p-value (two-sided): %.4f\n", x$pvalue))
  cat(sprintf("  Upper-tail p-value: %.4f\n", x$pvalue_upper))
  cat(sprintf("  Lower-tail p-value: %.4f\n", x$pvalue_lower))

  # COBA results if present
  if (!is.null(x$pvalue_adj)) {
    cat("\nCOBA Adjustment:\n")
    cat(sprintf("  Adjustment factor: %.4f\n", x$adj_factor))
    cat(sprintf("  Adjusted t-stat: %.4f\n", x$tco_obs_adj))
    cat(sprintf("  Adjusted p-value: %.4f\n", x$pvalue_adj))
    if (length(x$median_phi) > 0) {
      cat(sprintf("  Median AR(%d) phi: %s\n", length(x$median_phi),
                  paste(round(x$median_phi, 4), collapse = ", ")))
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
  cat(sprintf("  n = %d, nb = %d, maxp = %d, criterion = %s\n\n",
              object$n, object$nb, object$maxp, object$criterion))

  cat("Null Hypothesis: No trend (slope = 0)\n")
  cat(sprintf("  AR(%d) model with variance = %.4f\n\n", object$p, object$vara))

  cat("Test Statistic:\n")
  cat(sprintf("  Observed CO t = %.4f\n", object$tco_obs))
  cat(sprintf("  Bootstrap SE = %.4f\n", sd(object$boot_tstats)))
  cat(sprintf("  Bootstrap mean = %.4f\n\n", mean(object$boot_tstats)))

  cat("Result:\n")
  sig_level <- if (object$pvalue < 0.001) "***"
               else if (object$pvalue < 0.01) "**"
               else if (object$pvalue < 0.05) "*"
               else if (object$pvalue < 0.1) "."
               else ""
  cat(sprintf("  Two-sided p-value (abs-based) = %.4f %s\n", object$pvalue, sig_level))
  cat(sprintf("  Upper-tail p-value = %.4f\n", object$pvalue_upper))
  cat(sprintf("  Lower-tail p-value = %.4f\n", object$pvalue_lower))
  cat("  Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")

  # COBA adjustment results if present
  if (!is.null(object$pvalue_adj)) {
    cat("COBA Adjustment (variance correction):\n")
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

  invisible(object)
}
