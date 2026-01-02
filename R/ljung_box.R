#' Ljung-Box Test
#'
#' Perform the Ljung-Box portmanteau test for autocorrelation.
#'
#' @param x Numeric vector, the time series (or residuals).
#' @param K Integer, maximum lag to test (default 24).
#' @param p Integer, AR order of fitted model (default 0).
#' @param q Integer, MA order of fitted model (default 0).
#'
#' @return A list with components:
#'   \item{test}{Name of the test}
#'   \item{K}{Maximum lag tested}
#'   \item{chi_square}{Test statistic}
#'   \item{df}{Degrees of freedom (K - p - q)}
#'   \item{pval}{p-value}
#' @export
#'
#' @details
#' The test statistic is:
#' \deqn{Q = n(n+2) \sum_{k=1}^{K} \frac{\hat{\rho}_k^2}{n-k}}
#' which is compared to a chi-squared distribution with K - p - q degrees of freedom.
#'
#' When testing residuals from an ARMA(p,q) model, specify p and q to adjust
#' the degrees of freedom.
#'
#' @examples
#' # Test white noise
#' set.seed(123)
#' x <- rnorm(200)
#' ljung_box(x)
#'
#' # Test residuals from AR(2) fit
#' y <- arima.sim(list(ar = c(0.5, -0.3)), n = 200)
#' fit <- arima(y, order = c(2, 0, 0))
#' ljung_box(residuals(fit), K = 24, p = 2, q = 0)
ljung_box <- function(x, K = 24L, p = 0L, q = 0L) {

  # Input validation
  if (!is.numeric(x) || length(x) == 0) {
    stop("`x` must be a non-empty numeric vector")
  }
  if (!is.numeric(K) || K < 1) {
    stop("`K` must be a positive integer")
  }
  if (!is.numeric(p) || p < 0) {
    stop("`p` must be a non-negative integer")
  }
  if (!is.numeric(q) || q < 0) {
    stop("`q` must be a non-negative integer")
  }

  n <- length(x)

  if (K >= n) {
    stop("`K` must be less than the length of `x`")
  }

  df <- K - p - q
  if (df <= 0) {
    stop("Degrees of freedom (K - p - q) must be positive")
  }

  # Compute sample autocorrelations (lags 1 to K)
  acf_obj <- acf(x, lag.max = K, plot = FALSE, na.action = na.pass)
  rho <- acf_obj$acf[2:(K + 1)]  # exclude lag 0

  # Ljung-Box statistic
  weights <- 1 / seq.int(n - 1, n - K)
  chi_square <- n * (n + 2) * sum(weights * rho^2)

  # p-value
  pval <- pchisq(chi_square, df, lower.tail = FALSE)

  structure(
    list(
      test = "Ljung-Box test",
      K = K,
      chi_square = chi_square,
      df = df,
      pval = pval
    ),
    class = "ljung_box_test"
  )
}


#' @export
print.ljung_box_test <- function(x, ...) {
  cat("\n\t", x$test, "\n\n")
  cat("K =", x$K, "\n")
  cat("Chi-square =", format(x$chi_square, digits = 4), "\n")
  cat("df =", x$df, "\n")
  cat("p-value =", format(x$pval, digits = 4), "\n")
  invisible(x)
}
