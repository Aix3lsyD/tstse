#' Gegenbauer Polynomial Coefficients
#'
#' Computes the coefficients \eqn{C_k(d, u)} of the Gegenbauer polynomial expansion
#' arising from the power series:
#' \deqn{(1 - 2uB + B^2)^{-d} = \sum_{k=0}^{\infty} C_k(d, u) B^k}
#'
#' These coefficients are used in long-memory time series models, particularly
#' Gegenbauer ARMA (GARMA) processes where the spectral density has a peak at
#' a frequency determined by \eqn{u = \cos(\lambda)}.
#'
#' @param u Numeric scalar. The cosine of the Gegenbauer frequency, must satisfy
#'   \eqn{-1 \le u \le 1}. Related to the peak frequency \eqn{\lambda} by
#'   \eqn{u = \cos(\lambda)}.
#' @param d Numeric scalar. The long-memory (fractional differencing) parameter.
#'   Typically \eqn{0 < d < 0.5} for stationary long-memory processes.
#' @param n_coef Integer. The number of coefficients to compute (including
#'   \eqn{C_0}). Must be at least 1.
#'
#' @return Numeric vector of length `n_coef` containing \eqn{C_0, C_1, \ldots, C_{n-1}}.
#'
#' @details
#' The coefficients are computed using the recurrence relation:
#' \deqn{C_0 = 1}
#' \deqn{C_1 = 2du}
#' \deqn{C_k = \frac{2(k - 1 + d)u \cdot C_{k-1} - (k - 2 + 2d) C_{k-2}}{k}, \quad k \ge 2}
#'
#' This recurrence is numerically stable and efficient, requiring only O(n)
#' operations.
#'
#' @section Special Cases:
#' \itemize{
#'   \item When \eqn{d = 0}, all coefficients after \eqn{C_0} are zero (no long memory).
#'   \item When \eqn{u = 1} (frequency 0), the coefficients reduce to fractional
#'     differencing weights.
#'   \item When \eqn{u = -1} (frequency \eqn{\pi}), the coefficients exhibit an
#'     alternating sign pattern.
#'   \item When \eqn{u = 0} (frequency \eqn{\pi/2}), odd-indexed coefficients are zero.
#' }
#'
#' @references
#' Gray, H. L., Zhang, N.-F., & Woodward, W. A. (1989). On generalized fractional
#' processes. *Journal of Time Series Analysis*, 10(3), 233-257.
#'
#' Ferrara, L., & Guegan, D. (2001). Forecasting with k-factor Gegenbauer processes.
#' *Journal of Forecasting*, 20(8), 581-601.
#'
#' @seealso [gen_garma()] for generating GARMA processes,
#'   [est_garma()] for estimation.
#'
#' @examples
#' # Coefficients for u = 0.5 (frequency lambda = pi/3), d = 0.3
#' gegenb(u = 0.5, d = 0.3, n_coef = 10)
#'
#' # At u = 1 (frequency 0), reduces to fractional differencing weights
#' gegenb(u = 1, d = 0.4, n_coef = 5)
#'
#' # At u = -1 (frequency pi), alternating pattern
#' gegenb(u = -1, d = 0.4, n_coef = 5)
#'
#' # No long memory when d = 0
#' gegenb(u = 0.5, d = 0, n_coef = 5)
#'
#' @export
gegenb <- function(u, d, n_coef) {
  # Input validation
  if (!is.numeric(u) || length(u) != 1L || is.na(u)) {
    stop("`u` must be a single numeric value.", call. = FALSE)
  }
  if (!is.numeric(d) || length(d) != 1L || is.na(d)) {
    stop("`d` must be a single numeric value.", call. = FALSE)
  }
  if (!is.numeric(n_coef) || length(n_coef) != 1L || is.na(n_coef)) {
    stop("`n_coef` must be a single integer.", call. = FALSE)
  }

  n_coef <- as.integer(n_coef)
  if (n_coef < 1L) {
    stop("`n_coef` must be at least 1.", call. = FALSE)
  }

  if (u < -1 || u > 1) {
    warning("`u` is outside [-1, 1]; results may not correspond to a valid frequency.",
            call. = FALSE)
  }

  # Fast path for d = 0 (no long memory)
  if (d == 0) {
    coef <- numeric(n_coef)
    coef[1L] <- 1
    return(coef)
  }

  # Use Rcpp implementation for the recurrence relation
  gegenb_cpp(u, d, n_coef)
}
