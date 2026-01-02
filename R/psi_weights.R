#' Compute Psi Weights (MA Representation)
#'
#' Compute the psi weights (coefficients of the infinite MA representation)
#' for an ARMA model.
#'
#' @param phi Numeric vector, AR coefficients (default 0 = none).
#' @param theta Numeric vector, MA coefficients (default 0 = none).
#' @param lag_max Integer, maximum lag for psi weights (default 5).
#'
#' @return Numeric vector of psi weights from lag 1 to lag_max.
#' @export
#'
#' @details
#' The psi weights represent the infinite MA representation of an ARMA process:
#' \deqn{X_t = a_t + \psi_1 a_{t-1} + \psi_2 a_{t-2} + \cdots}
#'
#' This is a convenience wrapper around \code{\link[stats]{ARMAtoMA}} that
#' handles the sign convention used in ATSA.
#'
#' @seealso \code{\link{pi_weights}} for the dual (AR) representation.
#'
#' @examples
#' # AR(1) psi weights
#' psi_weights(phi = 0.8, lag_max = 10)
#'
#' # ARMA(1,1) psi weights
#' psi_weights(phi = 0.7, theta = 0.4, lag_max = 10)
psi_weights <- function(phi = 0, theta = 0, lag_max = 5L) {

  # Input validation
  if (!is.numeric(phi)) stop("`phi` must be numeric")
  if (!is.numeric(theta)) stop("`theta` must be numeric")
  if (!is.numeric(lag_max) || lag_max < 1) {
    stop("`lag_max` must be a positive integer")
  }

  # ARMAtoMA uses opposite sign convention for MA
  stats::ARMAtoMA(ar = phi, ma = -theta, lag.max = lag_max)
}
