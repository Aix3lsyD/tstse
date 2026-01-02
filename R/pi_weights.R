#' Compute Pi Weights (AR Representation)
#'
#' Compute the pi weights (coefficients of the infinite AR representation)
#' for an ARMA model.
#'
#' @param phi Numeric vector, AR coefficients (default 0 = none).
#' @param theta Numeric vector, MA coefficients (default 0 = none).
#' @param lag_max Integer, maximum lag for pi weights (default 5).
#'
#' @return Numeric vector of pi weights from lag 1 to lag_max.
#' @export
#'
#' @details
#' The pi weights represent the infinite AR representation of an ARMA process:
#' \deqn{X_t = \pi_1 X_{t-1} + \pi_2 X_{t-2} + \cdots + a_t}
#'
#' If the \code{astsa} package is installed, uses \code{astsa::ARMAtoAR}.
#' Otherwise, computes pi weights directly via polynomial inversion.
#'
#' @seealso \code{\link{psi_weights}} for the dual (MA) representation.
#'
#' @examples
#' # MA(1) pi weights
#' pi_weights(theta = 0.8, lag_max = 10)
#'
#' # ARMA(1,1) pi weights
#' pi_weights(phi = 0.7, theta = 0.4, lag_max = 10)
pi_weights <- function(phi = 0, theta = 0, lag_max = 5L) {

  # Input validation
  if (!is.numeric(phi)) stop("`phi` must be numeric")
  if (!is.numeric(theta)) stop("`theta` must be numeric")
  if (!is.numeric(lag_max) || lag_max < 1) {
    stop("`lag_max` must be a positive integer")
  }

  # Use astsa if available, otherwise compute directly
  if (requireNamespace("astsa", quietly = TRUE)) {
    astsa::ARMAtoAR(ar = phi, ma = -theta, lag.max = lag_max)
  } else {
    pi_weights_internal(phi = phi, theta = theta, lag_max = lag_max)
  }
}


#' Compute pi weights without astsa dependency
#' @noRd
pi_weights_internal <- function(phi, theta, lag_max) {

  # Pi weights are psi weights of the "inverse" model
  # If X_t = phi(B)^{-1} theta(B) a_t
  # Then pi(B) X_t = a_t where pi(B) = theta(B)^{-1} phi(B)
  # So pi weights are AR representation of ARMA(q, p) with swapped roles

  # Handle trivial cases
  if (all(phi == 0) && all(theta == 0)) {
    return(rep(0, lag_max))
  }

  # Pi weights satisfy: pi_j - theta_1*pi_{j-1} - ... - theta_q*pi_{j-q} = -phi_j
  # (with phi_j = 0 for j > p)
  # This is equivalent to ARMAtoMA with roles of AR and MA swapped

  # Swap AR and MA: pi weights of ARMA(phi, theta) = psi weights of ARMA(theta, phi)
  stats::ARMAtoMA(ar = theta, ma = -phi, lag.max = lag_max)
}
