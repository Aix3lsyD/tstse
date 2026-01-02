#' MA Coefficients for Gegenbauer Processes
#'
#' Computes the moving average (psi-weight) coefficients of the general linear
#' process representation for Gegenbauer processes with one or two factors.
#'
#' @param u Numeric vector of length 1 or 2. The cosine of the Gegenbauer
#'   frequency(ies), each value must satisfy \eqn{-1 \le u \le 1}. Related to
#'   peak frequency \eqn{\lambda} by \eqn{u = \cos(\lambda)}.
#' @param d Numeric vector of length 1 or 2. The long-memory (fractional
#'   differencing) parameter(s). Must have the same length as `u`. Typically
#'   \eqn{0 < d < 0.5} for stationary long-memory processes.
#' @param n_coef Integer. The number of MA coefficients to compute (including
#'   \eqn{\psi_0 = 1}). Must be at least 1. Default is 1000.
#'
#' @return Numeric vector of length `n_coef` containing \eqn{\psi_0, \psi_1,
#'   \ldots, \psi_{n-1}}, where \eqn{\psi_0 = 1}.
#'
#' @details
#' For a single Gegenbauer factor (\eqn{k = 1}), the MA coefficients are simply
#' the Gegenbauer polynomial coefficients computed by [gegenb()].
#'
#' For two Gegenbauer factors (\eqn{k = 2}), the process is:
#' \deqn{(1 - 2u_1 B + B^2)^{-d_1} (1 - 2u_2 B + B^2)^{-d_2} a_t}
#'
#' The MA coefficients are obtained by convolving the two sets of Gegenbauer
#' coefficients:
#' \deqn{\psi_j = \sum_{i=0}^{j} C_i^{(1)} C_{j-i}^{(2)}}
#'
#' This corresponds to the general linear process form described in Ferrara
#' and Guegan (2001), formula (8).
#'
#' @references
#' Ferrara, L., & Guegan, D. (2001). Forecasting with k-factor Gegenbauer
#' processes: Theory and applications. *Journal of Forecasting*, 20(8), 581-601.
#'
#' Woodward, W. A., Gray, H. L., & Elliott, A. C. (2017). *Applied Time Series
#' Analysis with R* (2nd ed.). CRC Press.
#'
#' @seealso [gegenb()] for computing Gegenbauer polynomial coefficients,
#'   [gen_garma()] for generating GARMA processes.
#'
#' @examples
#' # Single Gegenbauer factor
#' psi <- macoef_geg(u = 0.5, d = 0.3, n_coef = 100)
#' plot(0:99, psi, type = "h", xlab = "Lag", ylab = expression(psi))
#'
#' # Two Gegenbauer factors at different frequencies
#' psi2 <- macoef_geg(u = c(0.5, -0.3), d = c(0.25, 0.2), n_coef = 100)
#' plot(0:99, psi2, type = "h", xlab = "Lag", ylab = expression(psi))
#'
#' @export
macoef_geg <- function(u, d, n_coef = 1000L) {
  # Input validation

  if (!is.numeric(u) || anyNA(u)) {
    stop("`u` must be a numeric vector without NA values.", call. = FALSE)
  }
  if (!is.numeric(d) || anyNA(d)) {
    stop("`d` must be a numeric vector without NA values.", call. = FALSE)
  }
  if (length(u) != length(d)) {
    stop("`u` and `d` must have the same length.", call. = FALSE)
  }

  k <- length(d)
  if (k < 1L || k > 2L) {
    stop("Only 1 or 2 Gegenbauer factors are supported.", call. = FALSE)
  }

  if (!is.numeric(n_coef) || length(n_coef) != 1L || is.na(n_coef)) {
    stop("`n_coef` must be a single integer.", call. = FALSE)
  }
  n_coef <- as.integer(n_coef)
  if (n_coef < 1L) {
    stop("`n_coef` must be at least 1.", call. = FALSE)
  }

  # Note: gegenb() will warn if u values are outside [-1, 1]

  if (k == 1L) {
    # Single factor: MA coefficients are just Gegenbauer coefficients
    psi <- gegenb(u = u, d = d, n_coef = n_coef)
  } else {
    # Two factors: convolve the Gegenbauer coefficients
    # C1 and C2 are the coefficients for each factor
    C1 <- gegenb(u = u[1L], d = d[1L], n_coef = n_coef)
    C2 <- gegenb(u = u[2L], d = d[2L], n_coef = n_coef)

    # Convolve: psi_j = sum_{i=0}^{j} C1_i * C2_{j-i}
    # Use Rcpp implementation for O(n^2) convolution with minimal overhead
    psi <- convolve_truncated_cpp(C1, C2, n_coef)
  }

  psi
}
