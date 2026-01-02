#' Generate ARCH Process
#'
#' Generate a realization from an ARCH(q) process.
#'
#' @param n Integer, length of output series.
#' @param alpha0 Numeric, constant term (must be positive).
#' @param alpha Numeric vector, ARCH coefficients (q terms).
#' @param spin Integer, burn-in period (default 1000).
#' @param plot Logical, whether to plot the series (default TRUE).
#' @param seed Integer, random seed for reproducibility (default NULL).
#'
#' @return Numeric vector of length n.
#' @export
#'
#' @details
#' Generates from the ARCH(q) model:
#' \deqn{\sigma^2_t = \alpha_0 + \sum_{i=1}^{q} \alpha_i a_{t-i}^2}
#' \deqn{a_t = \epsilon_t \sigma_t}
#' where \eqn{\epsilon_t \sim N(0,1)}.
#'
#' ARCH (AutoRegressive Conditional Heteroskedasticity) models capture
#' volatility clustering, where large price changes tend to cluster together.
#' This is a special case of GARCH with no lagged variance terms.
#'
#' For stationarity, the sum of alpha coefficients should be less than 1.
#'
#' @seealso [gen_garch()] for the generalized GARCH model which includes
#'   lagged variance terms.
#'
#' @examples
#' # ARCH(1)
#' x <- gen_arch(n = 500, alpha0 = 0.1, alpha = 0.3)
#'
#' # ARCH(2)
#' x <- gen_arch(n = 500, alpha0 = 0.1, alpha = c(0.2, 0.15))
#'
#' # Reproducible generation
#' x <- gen_arch(n = 200, alpha0 = 0.1, alpha = 0.4, seed = 123, plot = FALSE)
gen_arch <- function(n,
                     alpha0,
                     alpha,
                     spin = 1000L,
                     plot = TRUE,
                     seed = NULL) {

  # Input validation
  if (!is.numeric(n) || length(n) != 1 || n < 1) {
    stop("`n` must be a positive integer")
  }
  if (!is.numeric(alpha0) || length(alpha0) != 1 || alpha0 <= 0) {
    stop("`alpha0` must be a positive number")
  }
  if (!is.numeric(alpha) || length(alpha) == 0) {
    stop("`alpha` must be a non-empty numeric vector")
  }
  if (any(alpha < 0)) {
    stop("all `alpha` coefficients must be non-negative")
  }
  if (!is.numeric(spin) || spin < 0) {
    stop("`spin` must be a non-negative integer")
  }

  # Call Rcpp implementation
  result <- gen_arch_cpp(
    n = as.integer(n),
    alpha0 = as.double(alpha0),
    alpha = as.double(alpha),
    spin = as.integer(spin),
    seed = seed
  )

  if (plot) {
    plot(result, type = "l", xlab = "Time", ylab = "",
         main = "ARCH Process")
    invisible(result)
  } else {
    result
  }
}
