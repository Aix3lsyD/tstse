#' Compute ARMA process variance
#'
#' Computes the variance of an ARMA process using the psi-weight representation.
#'
#' @param phi Numeric vector, AR coefficients.
#' @param theta Numeric vector, MA coefficients.
#' @param vara Numeric, white noise variance.
#' @param p Integer, AR order (optional, computed from phi if NULL).
#' @param q Integer, MA order (optional, computed from theta if NULL).
#'
#' @return Numeric, the process variance.
#' @noRd
.arma_variance <- function(phi, theta, vara, p = NULL, q = NULL) {

  if (is.null(p)) p <- if (all(phi == 0)) 0L else length(phi)
  if (is.null(q)) q <- if (all(theta == 0)) 0L else length(theta)

  if (p == 0 && q == 0) {
    return(vara)
  }

  # Get psi weights (MA infinity representation)
  max_lag <- max(100, 2 * max(p, q))
  psi <- c(1, ARMAtoMA(ar = phi, ma = -theta, lag.max = max_lag))

  # Process variance = vara * sum(psi^2)
  vara * sum(psi^2)
}


#' Get effective ARMA orders
#'
#' @param phi Numeric vector, AR coefficients.
#' @param theta Numeric vector, MA coefficients.
#'
#' @return List with p and q.
#' @noRd
.arma_orders <- function(phi, theta) {
  list(
    p = if (all(phi == 0)) 0L else length(phi),
    q = if (all(theta == 0)) 0L else length(theta)
  )
}
