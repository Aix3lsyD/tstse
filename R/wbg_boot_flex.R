#' Flexible WBG Bootstrap Test for Trend
#'
#' Performs a bootstrap hypothesis test for trend in a time series while
#' accounting for autocorrelation. Unlike [wbg_boot()], this function accepts
#' any test statistic function, enabling flexible testing with different
#' statistics (Cochrane-Orcutt, OLS, Mann-Kendall, Sen's slope, etc.).
#'
#' @importFrom stats median sd
#'
#' @param x Numeric vector. The observed time series.
#' @param stat_fn Function. A function that takes a numeric vector and returns
#'   a single numeric test statistic. Use factory functions like
#'   [make_stat_co()], [make_stat_ols_t()], etc., or create your own with
#'   signature `function(x) -> numeric(1)`.
#' @param nb Integer. Number of bootstrap replicates. Default is 399.
#' @param p_max Integer. Maximum AR order to consider when fitting the null
#'   model. Order is selected by information criterion. Default is 5.
#' @param ar_method Character. Method for AR estimation of null model:
#'   `"burg"` (default, always stationary) or `"mle"` (less biased near
#'   unit root but requires stationarity checking).
#' @param criterion Character. Information criterion for AR order selection:
#'   `"aic"` (default), `"aicc"`, or `"bic"`.
#' @param bootadj Logical. If TRUE (default), performs the COBA
#'   (Cochrane-Orcutt Bootstrap Adjustment) variance adjustment for
#'   improved small-sample performance.
#' @param cores Integer or NULL. Number of CPU cores for parallel processing.
#'   Default NULL uses `getOption("tstse.cores", 1)`. Set to 0 to use all
#'   available cores.
#' @param seed Integer or NULL. Random seed for reproducibility. If NULL
#'   (default), no seed is set.
#' @param verbose Logical. If TRUE (default), display progress messages.
#'
#' @return An object of class `"wbg_boot_flex"` containing:
#'   \item{obs_stat}{Observed test statistic.}
#'   \item{boot_dist}{Vector of bootstrap statistics.}
#'   \item{pvalue}{Two-sided p-value.}
#'   \item{pvalue_upper}{Upper-tail p-value.}
#'   \item{pvalue_lower}{Lower-tail p-value.}
#'   \item{ar_order}{Selected AR order under null.}
#'   \item{ar_phi}{AR coefficients under null.}
#'   \item{ar_vara}{Innovation variance under null (NA if MLE method used).}
#'   \item{ar_method}{AR estimation method used.}
#'   \item{n}{Sample size.}
#'   \item{nb}{Number of bootstrap replicates.}
#'   \item{boot_seeds}{Vector of RNG seeds used for each bootstrap replicate.}
#'   \item{master_seed}{The seed parameter used to generate boot_seeds (for reproducibility).}
#'
#'   If `bootadj = TRUE`, also includes:
#'   \item{obs_stat_adj}{COBA-adjusted observed statistic.}
#'   \item{pvalue_adj}{COBA-adjusted two-sided p-value.}
#'   \item{adj_factor}{Variance adjustment factor.}
#'   \item{median_phi}{Median AR coefficients from first stage.}
#'
#' @details
#' This function implements the Woodward-Bottone-Gray (WBG) bootstrap test
#' for trend in autocorrelated time series. The null hypothesis is that the
#' series is a stationary AR(p) process with no trend.
#'
#' **Procedure:**
#' 1. Compute the observed test statistic `stat_fn(x)`
#' 2. Fit AR(p) model to data under null hypothesis (no trend)
#' 3. Generate `nb` bootstrap series from the fitted AR model using
#'    [gen_ar_fast()]
#' 4. Compute test statistic for each bootstrap series
#' 5. Calculate p-value as proportion of |bootstrap stat| >= |observed stat|
#'
#' **COBA Adjustment:**
#' When `bootadj = TRUE`, a second-stage bootstrap adjusts the observed
#' statistic for variance inflation, providing better small-sample
#' performance. The adjustment uses the median AR coefficients from
#' the first stage.
#'
#' **Statistic Functions:**
#' Use the provided factory functions to create statistics:
#' - [make_stat_co()] - Cochrane-Orcutt t-statistic
#' - [make_stat_ols_t()] - OLS t-statistic
#' - [make_stat_ols_slope()] - OLS slope estimate
#' - [make_stat_mk()] - Mann-Kendall S (requires Kendall package)
#' - [make_stat_spearman()] - Spearman correlation
#' - [make_stat_sen()] - Sen's slope
#' - [make_stat_hac()] - HAC t-statistic (requires sandwich, lmtest)
#' - [make_stat_bn()] - Bloomfield-Nychka t-statistic
#' - [make_stat_lr()] - Likelihood ratio
#' - [make_stat_gls()] - GLS t-statistic (requires nlme)
#'
#' **Performance:**
#' Uses [gen_ar_fast()] which provides ~30-50x speedup over general
#' ARMA generators for bootstrap resampling.
#'
#' @references
#' Woodward, W. A., Bottone, S., and Gray, H. L. (1997). "Improved Tests for
#' Trend in Time Series Data." *Journal of Agricultural, Biological, and
#' Environmental Statistics*, 2(4), 403-416.
#'
#' @seealso [wbg_boot()] for the simpler interface with Cochrane-Orcutt only,
#'   [make_stat_co()] and other statistic factory functions,
#'   [gen_ar_fast()] for the fast AR generator used internally.
#'
#' @examples
#' \donttest{
#' # Test for trend in white noise (should not reject)
#' set.seed(123)
#' x <- rnorm(100)
#' stat_fn <- make_stat_co()
#' result <- wbg_boot_flex(x, stat_fn, nb = 99, seed = 456, verbose = FALSE)
#' print(result)
#'
#' # Test for trend in series with actual trend (should reject)
#' set.seed(123)
#' x_trend <- 0.1 * (1:100) + arima.sim(list(ar = 0.7), n = 100)
#' result <- wbg_boot_flex(x_trend, stat_fn, nb = 99, seed = 456, verbose = FALSE)
#' result$pvalue  # Should be small
#'
#' # Use different statistics
#' stat_ols <- make_stat_ols_t()
#' result_ols <- wbg_boot_flex(x_trend, stat_ols, nb = 99, seed = 456, verbose = FALSE)
#'
#' stat_spear <- make_stat_spearman()
#' result_spear <- wbg_boot_flex(x_trend, stat_spear, nb = 99, seed = 456, verbose = FALSE)
#' }
#'
#' @export
wbg_boot_flex <- function(x, stat_fn, nb = 399L, p_max = 5L,
                          ar_method = c("burg", "mle"),
                          criterion = c("aic", "aicc", "bic"),
                          bootadj = TRUE,
                          cores = 1L, seed = NULL,
                          verbose = TRUE) {

  ar_method <- match.arg(ar_method)
  criterion <- match.arg(criterion)

  # Input validation
  if (!is.numeric(x) || length(x) == 0) {
    stop("`x` must be a non-empty numeric vector", call. = FALSE)
  }
  if (!is.function(stat_fn)) {
    stop("`stat_fn` must be a function", call. = FALSE)
  }
  if (!is.numeric(nb) || length(nb) != 1 || nb < 1) {
    stop("`nb` must be a positive integer", call. = FALSE)
  }
  if (!is.numeric(p_max) || length(p_max) != 1 || p_max < 1) {
    stop("`p_max` must be a positive integer", call. = FALSE)
  }

  nb <- as.integer(nb)
  p_max <- as.integer(p_max)
  n <- length(x)
  cores <- get_cores(cores)

  # Set seed and generate boot_seeds for reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
  }
  boot_seeds <- sample.int(.Machine$integer.max, nb)
  if (bootadj) {
    boot_seeds_adj <- sample.int(.Machine$integer.max, nb)
  }

  # Compute observed statistic
  obs_stat <- stat_fn(x)

  # Fit AR model under H0: no trend, stationary AR
  # These run in main thread - safe to use defaults
  if (verbose) message("Fitting null model (", ar_method, ")...")

  if (ar_method == "burg") {
    ar_fit <- aic_burg(x, p = seq_len(p_max), type = criterion)
    ar_method_used <- "burg"
    ar_vara <- ar_fit$vara
  } else {
    ar_fit <- aic_ar_mle(x, p_max = p_max, criterion = criterion)
    ar_method_used <- "mle"
    ar_vara <- NA_real_  # MLE doesn't return vara directly
  }
  ar_phi <- ar_fit$phi
  ar_p <- ar_fit$p

  if (verbose) {
    message("  Selected AR(", ar_p, ") with phi = [",
            paste(round(ar_phi, 4), collapse = ", "), "]")
  }

  # --- First Bootstrap ---
  if (verbose) message("Running bootstrap (", nb, " replicates)...")

  # ===========================================================================
  # CRITICAL FIX: boot_fn runs inside parallel workers.
  # All AR fitting calls MUST use cores = 1L to prevent nested parallelization.
  # The pmap() function also detects if it's inside a worker and forces
  # sequential, but explicit cores = 1L is clearer and more defensive.
  # ===========================================================================
  boot_fn <- function(i) {
    # Generate AR series under null using fast generator
    xb <- gen_ar_fast(n, phi = ar_phi, seed = boot_seeds[i])
    stat <- stat_fn(xb)

    if (bootadj) {
      # Fit AR to bootstrap sample for COBA adjustment
      # NOTE: aic_burg() defaults to cores = 1L (safe for nested calls)
      #       aic_ar_mle() is already fully sequential (no cores param)
      ar_boot <- tryCatch({
        if (ar_method == "burg") {
          aic_burg(xb, p = seq_len(p_max), type = criterion)
        } else {
          aic_ar_mle(xb, p_max = p_max, criterion = criterion)
        }
      }, error = function(e) {
        # Fallback to Burg on error
        aic_burg(xb, p = seq_len(p_max), type = criterion)
      })
      list(stat = stat, phi = ar_boot$phi, phi1 = 1 - sum(ar_boot$phi))
    } else {
      list(stat = stat)
    }
  }

  # Run bootstrap (parallel or sequential)
  if (cores > 1L) {
    boot_results <- pmap(seq_len(nb), boot_fn, cores = cores)
  } else {
    boot_results <- lapply(seq_len(nb), boot_fn)
  }

  boot_stats <- sapply(boot_results, `[[`, "stat")

  # Compute p-values (with +1 correction for proper bootstrap p-values)
  pvalue_two <- (sum(abs(boot_stats) >= abs(obs_stat)) + 1) / (nb + 1)
  pvalue_upper <- (sum(boot_stats >= obs_stat) + 1) / (nb + 1)
  pvalue_lower <- (sum(boot_stats <= obs_stat) + 1) / (nb + 1)

  # Quantile-based two-sided p-value (paper-faithful, Section 2.1)
  pvalue_quantile <- 2 * min(pvalue_upper, pvalue_lower)

  # Asymptotic p-value (two-sided, standard normal)
  pvalue_asymp <- 2 * pnorm(-abs(obs_stat))

  result <- list(
    obs_stat = obs_stat,
    boot_dist = boot_stats,
    pvalue = pvalue_two,
    pvalue_upper = pvalue_upper,
    pvalue_lower = pvalue_lower,
    pvalue_quantile = pvalue_quantile,
    pvalue_asymp = pvalue_asymp,
    ar_order = ar_p,
    ar_phi = ar_phi,
    ar_vara = ar_vara,
    ar_method = ar_method_used,
    n = n,
    nb = nb,
    boot_seeds = boot_seeds,
    master_seed = seed
  )

  # --- COBA Adjustment (Second Bootstrap) ---
  if (bootadj) {
    if (verbose) message("Running COBA adjustment...")

    boot_phi1 <- sapply(boot_results, `[[`, "phi1")
    boot_phi_list <- lapply(boot_results, `[[`, "phi")

    # Find median AR model
    median_idx <- which.min(abs(boot_phi1 - median(boot_phi1)))
    median_phi <- boot_phi_list[[median_idx]]

    # Second stage bootstrap function - no AR fitting needed here
    boot_fn_adj <- function(i) {
      xb <- gen_ar_fast(n, phi = median_phi, seed = boot_seeds_adj[i])
      stat_fn(xb)
    }

    if (cores > 1L) {
      boot_stats_adj <- unlist(pmap(seq_len(nb), boot_fn_adj, cores = cores))
    } else {
      boot_stats_adj <- vapply(seq_len(nb), boot_fn_adj, numeric(1))
    }

    # Compute adjustment factor and adjusted statistic
    sd_boot <- sd(boot_stats)
    if (sd_boot < .Machine$double.eps) {
      adj_factor <- 1.0
      warning("Zero-variance bootstrap distribution; COBA adjustment set to 1.0",
              call. = FALSE)
    } else {
      adj_factor <- sd(boot_stats_adj) / sd_boot
    }
    obs_stat_adj <- adj_factor * obs_stat

    # Adjusted p-value
    pvalue_adj <- (sum(abs(boot_stats) >= abs(obs_stat_adj)) + 1) / (nb + 1)

    result$obs_stat_adj <- obs_stat_adj
    result$pvalue_adj <- pvalue_adj
    result$adj_factor <- adj_factor
    result$median_phi <- median_phi
  }

  if (verbose) message("Done.")

  class(result) <- "wbg_boot_flex"
  result
}


#' @describeIn wbg_boot_flex Print method for wbg_boot_flex objects
#' @param x A wbg_boot_flex object
#' @param ... Additional arguments (ignored)
#' @export
print.wbg_boot_flex <- function(x, ...) {
  cat("WBG Bootstrap Test (Flexible)\n")
  cat("=============================\n")
  cat("Sample size:", x$n, "\n")
  cat("Bootstrap replicates:", x$nb, "\n")
  cat("AR method:", x$ar_method, "\n")
  cat("AR order:", x$ar_order, "\n")
  if (length(x$ar_phi) > 0) {
    cat("AR coefficients:", paste(round(x$ar_phi, 4), collapse = ", "), "\n")
  }
  cat("\n")
  cat("Observed statistic:", round(x$obs_stat, 4), "\n")
  cat("\n")
  cat("P-values:\n")
  cat("  Two-sided:", format.pval(x$pvalue, digits = 4), "\n")
  cat("  Upper-tail:", format.pval(x$pvalue_upper, digits = 4), "\n")
  cat("  Lower-tail:", format.pval(x$pvalue_lower, digits = 4), "\n")

  if (!is.null(x$pvalue_adj)) {
    cat("\nCOBA Adjustment:\n")
    cat("  Adjustment factor:", round(x$adj_factor, 4), "\n")
    cat("  Adjusted statistic:", round(x$obs_stat_adj, 4), "\n")
    cat("  Adjusted p-value:", format.pval(x$pvalue_adj, digits = 4), "\n")
  }

  invisible(x)
}


#' @describeIn wbg_boot_flex Summary method for wbg_boot_flex objects
#' @param object A wbg_boot_flex object
#' @param ... Additional arguments (ignored)
#' @export
summary.wbg_boot_flex <- function(object, ...) {
  cat("WBG Bootstrap Test Summary\n")
  cat("--------------------------\n")
  cat("n =", object$n, ", nb =", object$nb, "\n")
  cat("AR(", object$ar_order, ") null model fitted via ", object$ar_method, "\n", sep = "")
  cat("\n")
  cat("Observed:", round(object$obs_stat, 4), "\n")
  cat("Bootstrap mean:", round(mean(object$boot_dist), 4), "\n")
  cat("Bootstrap SD:", round(sd(object$boot_dist), 4), "\n")
  cat("Bootstrap range: [", round(min(object$boot_dist), 4), ", ",
      round(max(object$boot_dist), 4), "]\n", sep = "")
  cat("\n")

  # Significance interpretation
  if (object$pvalue < 0.01) {
    sig <- "*** (p < 0.01)"
  } else if (object$pvalue < 0.05) {
    sig <- "** (p < 0.05)"
  } else if (object$pvalue < 0.10) {
    sig <- "* (p < 0.10)"
  } else {
    sig <- "(not significant)"
  }
  cat("p-value:", format.pval(object$pvalue, digits = 4), sig, "\n")

  if (!is.null(object$pvalue_adj)) {
    if (object$pvalue_adj < 0.01) {
      sig_adj <- "*** (p < 0.01)"
    } else if (object$pvalue_adj < 0.05) {
      sig_adj <- "** (p < 0.05)"
    } else if (object$pvalue_adj < 0.10) {
      sig_adj <- "* (p < 0.10)"
    } else {
      sig_adj <- "(not significant)"
    }
    cat("COBA-adjusted p-value:", format.pval(object$pvalue_adj, digits = 4),
        sig_adj, "\n")
  }

  invisible(object)
}
