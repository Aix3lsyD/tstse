#' Bootstrap Cochrane-Orcutt Trend Test with Turner Adjusted Sample Size
#'
#' Performs a bootstrap version of the Cochrane-Orcutt trend test with
#' Turner adjusted sample size. The bootstrap p-value is computed by
#' comparing the observed p-value to the distribution of p-values
#' under the null hypothesis (no trend).
#'
#' @param x A numeric vector containing the time series data.
#' @param nb Number of bootstrap replicates. Default is 399.
#' @param maxp Maximum AR order for model selection. Default is 5.
#' @param type Character. Information criterion for order selection:
#'   `"aic"`, `"aicc"`, or `"bic"`. Default is `"aic"`.
#' @param seed Optional integer seed for reproducibility. Default is NULL.
#'
#' @return A list with class `"co_tas_boot"` containing:
#'   \item{pvalue}{Bootstrap p-value for test of H0: no trend.}
#'   \item{tco}{t-statistic from the original data.}
#'   \item{pvalue_asymptotic}{Asymptotic p-value from [co_tas()].}
#'   \item{phi}{AR coefficients from the original fit.}
#'   \item{nb}{Number of bootstrap replicates used.}
#'
#' @details
#' The bootstrap procedure:
#' \enumerate{
#'   \item Fit [co_tas()] to the original data to obtain the AR model
#'   \item Generate `nb` bootstrap series from AR(phi) with no trend
#'   \item For each bootstrap series, compute the p-value using [co_tas()]
#'   \item Bootstrap p-value = proportion of bootstrap p-values smaller
#'     than the observed p-value
#' }
#'
#' This provides a more accurate p-value than the asymptotic version
#' when sample sizes are small or autocorrelation is strong.
#'
#' @references
#' Woodward, W. A., Gray, H. L., and Elliott, A. C. (2017).
#' *Applied Time Series Analysis with R*. CRC Press.
#'
#' @seealso [co_tas()] for the non-bootstrap version,
#'   [wbg_boot()] for WBG bootstrap trend test
#'
#' @examples
#' # Generate series with trend and AR(1) errors
#' set.seed(123)
#' n <- 100
#' t <- 1:n
#' z <- arima.sim(list(ar = 0.7), n = n)
#' x <- 5 + 0.1 * t + z
#'
#' # Bootstrap trend test (use more replicates in practice)
#' result <- co_tas_boot(x, nb = 99, seed = 42)
#' print(result)
#'
#' @export
co_tas_boot <- function(x, nb = 399L, maxp = 5L,
                         type = c("aic", "aicc", "bic"),
                         seed = NULL) {

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

  n <- length(x)
  nb <- as.integer(nb)
  maxp <- as.integer(maxp)

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Step 1: Fit co_tas to original data
  result <- co_tas(x, maxp = maxp, type = type)
  pval_obs <- result$pvalue
  phi <- result$phi

  # Step 2: Bootstrap under null hypothesis (no trend)
  boot_pvals <- numeric(nb)

  for (i in seq_len(nb)) {
    # Generate AR series with no trend
    boot_dat <- gen_arma(n = n, phi = phi, plot = FALSE)

    # Compute p-value for bootstrap series
    boot_result <- co_tas(boot_dat, maxp = maxp, type = type)
    boot_pvals[i] <- boot_result$pvalue
  }

  # Step 3: Bootstrap p-value = proportion of boot p-values < observed
  pvalue_boot <- sum(boot_pvals < pval_obs) / nb

  # Return result
  result_boot <- list(
    pvalue = pvalue_boot,
    tco = result$tco,
    pvalue_asymptotic = pval_obs,
    phi = phi,
    nb = nb
  )

  class(result_boot) <- "co_tas_boot"
  result_boot
}

#' @export
print.co_tas_boot <- function(x, digits = 4, ...) {
  cat("Bootstrap Cochrane-Orcutt Trend Test (Turner Adjusted Sample Size)\n")
  cat("===================================================================\n\n")
  cat("t-statistic:", round(x$tco, digits), "\n")
  cat("Bootstrap p-value:", format.pval(x$pvalue, digits = digits), "\n")
  cat("Asymptotic p-value:", format.pval(x$pvalue_asymptotic, digits = digits), "\n")
  cat("Bootstrap replicates:", x$nb, "\n")
  if (length(x$phi) > 0) {
    cat("AR coefficients:", paste(round(x$phi, digits), collapse = ", "), "\n")
  }
  invisible(x)
}
