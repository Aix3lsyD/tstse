#' Cochrane-Orcutt Trend Test with Turner Adjusted Sample Size
#'
#' Performs a trend test using the Cochrane-Orcutt procedure with
#' effective sample size adjustment based on Turner's method. This
#' accounts for autocorrelation when computing degrees of freedom
#' for the t-test.
#'
#' @param x A numeric vector containing the time series data.
#' @param maxp Maximum AR order for model selection. Default is 5.
#' @param type Character. Information criterion for order selection:
#'   `"aic"`, `"aicc"`, or `"bic"`. Default is `"aic"`.
#'
#' @return A list with class `"co_tas"` containing:
#'   \item{pvalue}{P-value for test of H0: no trend (adjusted for autocorrelation).}
#'   \item{tco}{t-statistic for the slope coefficient.}
#'   \item{phi}{AR coefficients for the fitted model.}
#'   \item{n_a}{Effective sample size (adjusted for autocorrelation).}
#'   \item{p}{AR order selected.}
#'
#' @details
#' The procedure:
#' \enumerate{
#'   \item Differences the series and reconstructs via cumulative sum
#'   \item Fits AR(p) model using Burg's method with AIC selection
#'   \item Transforms both the data and time index using AR filter
#'   \item Computes effective sample size using theoretical ACF:
#'     \deqn{n_a = (n-2) / A}
#'     where \eqn{A = 1 + 2 \sum_{k=1}^{n-3} (1 - k/(n-2)) \rho_k}
#'   \item Conducts t-test with adjusted degrees of freedom
#' }
#'
#' The effective sample size adjustment accounts for the reduced
#' information content when observations are correlated.
#'
#' @references
#' Turner, J. (Year). "Title of Turner's Paper."
#'
#' Woodward, W. A., Gray, H. L., and Elliott, A. C. (2017).
#' *Applied Time Series Analysis with R*. CRC Press.
#'
#' @seealso [co()] for standard Cochrane-Orcutt,
#'   [co_tas_boot()] for bootstrap version,
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
#' # Test for trend
#' result <- co_tas(x, maxp = 5)
#' print(result)
#'
#' @export
co_tas <- function(x, maxp = 5L, type = c("aic", "aicc", "bic")) {

 type <- match.arg(type)

  # Input validation
  if (!is.numeric(x) || length(x) == 0) {
    stop("`x` must be a non-empty numeric vector", call. = FALSE)
  }
  if (!is.numeric(maxp) || length(maxp) != 1 || maxp < 1) {
    stop("`maxp` must be a positive integer", call. = FALSE)
  }

  n <- length(x)
  maxp <- as.integer(maxp)
  t_idx <- seq_len(n)

 # Step 1: Difference the series (AR transform with phi = 1)
  x_diff <- artrans(x, phi = 1, plot = FALSE)

  # Step 2: Reconstruct via cumulative sum (centered)
  z_diff <- cumsum(c(x[1], x_diff - mean(x_diff)))

  # Step 3: Fit AR model to differenced series using Burg
  fit <- aic_burg(z_diff, p = seq_len(maxp), type = type)
  phi <- fit$phi
  p <- fit$p

  # Step 4: Transform data and time index
  x_trans <- artrans(x, phi = phi, plot = FALSE)
  t_trans <- artrans(t_idx, phi = phi, plot = FALSE)

  # Step 5: Regression on transformed variables
  d_co <- stats::lm(x_trans ~ t_trans)

  # Step 6: Compute effective sample size using theoretical ACF
  acf_vals <- true_acf(phi = phi, lag_max = n - 2, plot = FALSE)$acf
  # Remove lag 0 (which is 1)
  acf_vals <- acf_vals[-1]

  # Compute adjustment factor A
  k <- seq_along(acf_vals)
  weights <- 1 - k / (n - 2)
  A <- 1 + 2 * sum(weights * acf_vals)

  # Effective sample size
  n_a <- (n - 2) / A

  # Step 7: Conduct t-test with adjusted degrees of freedom
  # df = n_a - p, but ensure df > 0
  df <- if (n_a > p) n_a - p else n_a

  t_value <- summary(d_co)$coefficients[2, 3]
  p_value <- 2 * (1 - stats::pt(abs(t_value), df = df))

  # Return result
  result <- list(
    pvalue = p_value,
    tco = t_value,
    phi = phi,
    n_a = n_a,
    p = p
  )

  class(result) <- "co_tas"
  result
}

#' @export
print.co_tas <- function(x, digits = 4, ...) {
  cat("Cochrane-Orcutt Trend Test (Turner Adjusted Sample Size)\n")
  cat("=========================================================\n\n")
  cat("t-statistic:", round(x$tco, digits), "\n")
  cat("p-value:", format.pval(x$pvalue, digits = digits), "\n")
  cat("Effective sample size:", round(x$n_a, 2), "\n")
  cat("AR order:", x$p, "\n")
  if (length(x$phi) > 0) {
    cat("AR coefficients:", paste(round(x$phi, digits), collapse = ", "), "\n")
  }
  invisible(x)
}
