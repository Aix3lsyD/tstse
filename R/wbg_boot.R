#' WBG Bootstrap Test for Trend
#'
#' Bootstrap hypothesis test for trend in time series with autocorrelated
#' errors. Uses the Woodward-Bottone-Gray (WBG) method to test H0: slope = 0
#' when residuals follow an AR(p) process.
#'
#' @param x A numeric vector containing the time series data.
#' @param nb Number of bootstrap replicates. Default is 399.
#' @param maxp Maximum AR order for model selection. Default is 5.
#' @param method Character. Method for AR estimation: `"mle"`, `"burg"`,
#'   or `"yw"`. Default is `"burg"` (matches original tswge).
#' @param type Character. Information criterion for order selection:
#'   `"aic"`, `"aicc"`, or `"bic"`. Default is `"aic"`.
#' @param cores Integer, number of cores for parallel processing.
#'   Default NULL uses `getOption("tstse.cores", 1)`.
#'   Set to 0 to use all available cores.
#' @param seed Integer, random seed for reproducibility. Default NULL.
#'
#' @return A list containing:
#'   \item{p}{AR order selected for the series under null hypothesis.}
#'   \item{phi}{AR coefficients for the null model.}
#'   \item{pvalue}{Bootstrap p-value for trend test (two-sided).}
#'   \item{tco_obs}{Observed t-statistic from Cochrane-Orcutt fit.}
#'   \item{boot_tstats}{Numeric vector of bootstrap t-statistics.}
#'   \item{nb}{Number of bootstrap replicates used.}
#'
#' @details
#' The WBG bootstrap test addresses the problem that standard t-tests for
#' trend are invalid when residuals are autocorrelated. The procedure:
#'
#' \enumerate{
#'   \item Fit Cochrane-Orcutt model to observed data, obtaining t-statistic
#'   \item Fit AR(p) model to data under null hypothesis (no trend)
#'   \item Generate `nb` bootstrap series from the fitted AR(p) model
#'   \item For each bootstrap series, fit Cochrane-Orcutt and compute t-statistic
#'   \item p-value = proportion of |bootstrap t| >= |observed t|
#' }
#'
#' This provides a more accurate p-value than the asymptotic approximation
#' used in [co()], especially for small to moderate sample sizes with
#' highly autocorrelated errors.
#'
#' @references
#' Woodward, W. A., Bottone, S., and Gray, H. L. (1997). "Improved Tests for
#' Trend in Time Series Data." *Journal of Agricultural, Biological, and
#' Environmental Statistics*, 2(4), 403-416.
#'
#' @seealso [co()] for Cochrane-Orcutt estimation,
#'   [aic_ar()] for AR order selection.
#'
#' @examples
#' # Test series with no trend (should be non-significant)
#' set.seed(123)
#' x <- arima.sim(list(ar = 0.7), n = 100)
#' result <- wbg_boot(x, nb = 99, seed = 456)
#' cat("p-value:", result$pvalue, "\n")
#'
#' # Test series with trend (should be significant)
#' set.seed(123)
#' x_trend <- arima.sim(list(ar = 0.7), n = 100) + 0.1 * (1:100)
#' result <- wbg_boot(x_trend, nb = 99, seed = 456)
#' cat("p-value:", result$pvalue, "\n")
#'
#' @export
wbg_boot <- function(x, nb = 399L, maxp = 5L,
                     method = c("burg", "mle", "yw"),
                     type = c("aic", "aicc", "bic"),
                     cores = 1L, seed = NULL) {

  method <- match.arg(method)
  type <- match.arg(type)

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

  nb <- as.integer(nb)
  maxp <- as.integer(maxp)
  n <- length(x)
  cores <- get_cores(cores)

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Step 1: Fit Cochrane-Orcutt on observed data (main thread - OK to use defaults)
  w <- co(x, maxp = maxp, method = method, type = type)
  tco_obs <- w$tco

  # Step 2: Fit AR model under null hypothesis (no trend)
  # These run in main thread before parallel loop - safe
  if (method == "burg") {
    x_aic <- aic_burg(x, p = seq_len(maxp))
  } else {
    x_aic <- aic_ar(x, p = seq_len(maxp), method = method, type = type)
  }
  p_null <- x_aic$p
  phi_null <- x_aic$phi

  # Step 3: Bootstrap loop
  # Generate seeds for reproducibility with parallel processing
  boot_seeds <- if (!is.null(seed)) sample.int(.Machine$integer.max, nb) else NULL

  # NOTE: boot_fn runs inside parallel workers.
  # co() internally uses cores = 1L, so no nested parallelization occurs.
  # The pmap() function also detects if it's inside a worker and forces sequential.
  boot_fn <- function(i) {
    # Set seed for this replicate if provided
    if (!is.null(boot_seeds)) {
      set.seed(boot_seeds[i])
    }

    # Generate AR(p) series under null (no trend)
    xb <- gen_arma(n = n, phi = phi_null, theta = 0, plot = FALSE)

    # Fit CO and get t-statistic (return signed value for distribution)
    # co() uses cores = 1L internally - safe for nested calls
    wb <- co(xb, maxp = maxp, method = method, type = type)
    wb$tco
  }

  # Run bootstrap (parallel or sequential)
  if (cores > 1L) {
    boot_results <- pmap(seq_len(nb), boot_fn, cores = cores)
    boot_tstats <- unlist(boot_results)
  } else {
    boot_tstats <- vapply(seq_len(nb), boot_fn, numeric(1))
  }

  # Step 4: Compute p-value
  # Proportion of bootstrap |t| >= observed |t|
  pvalue <- mean(abs(boot_tstats) >= abs(tco_obs))

  list(
    p           = p_null,
    phi         = phi_null,
    pvalue      = pvalue,
    tco_obs     = tco_obs,
    boot_tstats = boot_tstats,
    nb          = nb
  )
}
