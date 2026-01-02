#' Simple Linear Regression for Time Series
#'
#' Fit a simple linear regression to detrend a time series.
#'
#' @param x Numeric vector, the time series.
#'
#' @return An object of class "slr" with components:
#'   \item{res}{Residuals (detrended series)}
#'   \item{b0hat}{Intercept estimate}
#'   \item{b1hat}{Slope estimate}
#'   \item{pvalue}{P-value for test of H0: slope = 0}
#'   \item{tstatistic}{t-statistic for slope coefficient}
#' @export
#'
#' @details
#' Fits the model \eqn{X_t = \beta_0 + \beta_1 t + \epsilon_t} where \eqn{t = 1, 2, \ldots, n}.
#' The residuals represent the detrended series.
#'
#' Note: The p-value assumes independent, normally distributed errors. For time series
#' with correlated errors, use \code{\link{co}} (Cochrane-Orcutt) for more accurate inference.
#'
#' @seealso \code{\link{co}} for regression with correlated errors
#'
#' @examples
#' # Detrend a series with linear trend
#' set.seed(123)
#' x <- 1:100 * 0.5 + rnorm(100)
#' result <- slr(x)
#' print(result)
#'
#' # Plot original vs detrended
#' plot(x, type = "l", main = "Original")
#' plot(result$res, type = "l", main = "Detrended")
slr <- function(x) {


  # Input validation

  if (!is.numeric(x)) {
    stop("`x` must be numeric")
  }
  if (length(x) < 3) {
    stop("`x` must have at least 3 observations")
  }

  n <- length(x)
  time <- seq_len(n)

  # Fit linear model
 fit <- stats::lm(x ~ time)
  fit_summary <- summary(fit)

  # Extract coefficients
  b0hat <- unname(fit$coefficients[1])
  b1hat <- unname(fit$coefficients[2])

  # Compute residuals (detrended series)
  res <- x - b0hat - b1hat * time

  # Extract test statistics for slope
  pvalue <- fit_summary$coefficients[2, 4]
  tstatistic <- fit_summary$coefficients[2, 3]

  structure(
    list(
      res = res,
      b0hat = b0hat,
      b1hat = b1hat,
      pvalue = pvalue,
      tstatistic = tstatistic
    ),
    class = "slr"
  )
}


#' @rdname slr
#' @param x An object of class "slr".
#' @param ... Additional arguments (ignored).
#' @export
print.slr <- function(x, ...) {

  cat("Simple Linear Regression\n\n")
  cat("Intercept (b0):", round(x$b0hat, 4), "\n")
  cat("Slope (b1):    ", round(x$b1hat, 4), "\n\n")
  cat("Test for slope = 0:\n")
  cat("  t-statistic:", round(x$tstatistic, 4), "\n")
  cat("  p-value:    ", format.pval(x$pvalue, digits = 4), "\n")

  invisible(x)
}
