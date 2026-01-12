#' Fast Cochrane-Orcutt Trend Test with Turner Adjusted Sample Size
#'
#' C++ accelerated version of the Cochrane-Orcutt trend test with
#' Turner adjusted sample size. Approximately 10x faster than [co_tas()].
#'
#' @param x A numeric vector containing the time series data.
#' @param maxp Maximum AR order for model selection. Default is 5.
#' @param type Character. Information criterion for order selection:
#'   `"aic"`, `"aicc"`, or `"bic"`. Default is `"aic"`.
#'
#' @return A list with class `"co_tas"` containing:
#'   \item{pvalue}{P-value for test of H0: no trend.}
#'   \item{tco}{t-statistic for the slope coefficient.}
#'   \item{phi}{AR coefficients from Burg estimation.}
#'   \item{n_a}{Effective (adjusted) sample size.}
#'   \item{p}{AR order selected by information criterion.}
#'
#' @details
#' This is a C++ implementation of [co_tas()] that is approximately 10x faster.
#' The algorithm is identical:
#' \enumerate{
#'   \item Difference the series and reconstruct via cumulative sum (centered)
#'   \item Fit AR(p) model to the reconstructed series using Burg algorithm
#'   \item Transform both data and time index using the AR filter
#'   \item Perform OLS regression of transformed data on transformed time
#'   \item Compute effective sample size using Turner's adjustment
#'   \item Compute two-sided p-value with adjusted degrees of freedom
#' }
#'
#' The effective sample size accounts for autocorrelation in the data,
#' providing a more accurate p-value than naive OLS regression.
#'
#' @references
#' Turner, J. R. (1988). Effective sample size: A frequently neglected concept.
#' *Journal of Educational Statistics*, 13(2), 175-181.
#'
#' Woodward, W. A., Gray, H. L., and Elliott, A. C. (2017).
#' *Applied Time Series Analysis with R*. CRC Press.
#'
#' @seealso [co_tas()] for the R implementation,
#'   [co_tas_boot_fast()] for bootstrap version
#'
#' @examples
#' # Generate series with AR(1) errors and trend
#' set.seed(123)
#' n <- 100
#' t <- 1:n
#' z <- arima.sim(list(ar = 0.7), n = n)
#' x <- 5 + 0.1 * t + z
#'
#' # Fast CO-TAS trend test
#' result <- co_tas_fast(x, maxp = 5)
#' print(result)
#'
#' @export
co_tas_fast <- function(x, maxp = 5L, type = c("aic", "aicc", "bic")) {

  type <- match.arg(type)

  # Input validation
  if (!is.numeric(x) || length(x) == 0) {
    stop("`x` must be a non-empty numeric vector", call. = FALSE)
  }
  if (!is.numeric(maxp) || length(maxp) != 1 || maxp < 1) {
    stop("`maxp` must be a positive integer", call. = FALSE)
  }
  if (length(x) < 10) {
    stop("Series must have at least 10 observations", call. = FALSE)
  }

  # Coerce to proper types
  x <- as.numeric(x)
  maxp <- as.integer(maxp)

  # Call C++ implementation
  result <- co_tas_cpp(x, maxp = maxp, type = type)

  # Convert phi from matrix to vector (arma::vec returns as column matrix)
  result$phi <- as.numeric(result$phi)

  # Compute EXACT p-value using R's pt() instead of C++ approximation
  # The C++ code uses a normal approximation that can be inaccurate for small df
  df <- if (result$n_a > result$p) result$n_a - result$p else result$n_a
  if (df < 1) df <- 1
  result$pvalue <- 2 * (1 - stats::pt(abs(result$tco), df = df))

  # Return with same class as co_tas for compatibility
  class(result) <- "co_tas"
  result
}
