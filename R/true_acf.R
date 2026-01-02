#' True ARMA Autocorrelations
#'
#' Compute and optionally plot the theoretical autocorrelations for an ARMA model.
#'
#' @param phi Numeric vector, AR coefficients (default 0 = none).
#' @param theta Numeric vector, MA coefficients (default 0 = none).
#' @param lag_max Integer, maximum lag (default 25).
#' @param vara Numeric, white noise variance (default 1).
#' @param plot Logical, whether to plot ACF (default TRUE).
#'
#' @return A list with components:
#'   \item{acf}{Theoretical autocorrelations (lags 0 to lag_max)}
#'   \item{acv}{Theoretical autocovariances (lags 0 to lag_max)}
#' @export
#'
#' @details
#' Computes the theoretical autocorrelation function for an ARMA(p,q) model
#' using \code{\link[stats]{ARMAacf}}.
#'
#' @seealso \code{\link{true_spec}} for theoretical spectral density.
#'
#' @examples
#' # AR(1): rho(k) = phi^k
#' true_acf(phi = 0.8, lag_max = 10)
#'
#' # AR(2)
#' true_acf(phi = c(1.5, -0.75), lag_max = 20)
#'
#' # ARMA(1,1)
#' true_acf(phi = 0.7, theta = 0.4, lag_max = 15)
#'
#' # Without plot
#' result <- true_acf(phi = 0.9, plot = FALSE)
true_acf <- function(phi = 0, theta = 0, lag_max = 25L, vara = 1, plot = TRUE) {

  # Input validation
  if (!is.numeric(phi)) stop("`phi` must be numeric")
  if (!is.numeric(theta)) stop("`theta` must be numeric")
  if (!is.numeric(lag_max) || lag_max < 0) {
    stop("`lag_max` must be a non-negative integer")
  }
  if (!is.numeric(vara) || vara <= 0) {
    stop("`vara` must be positive")
  }

  # Effective orders
  p <- if (all(phi == 0)) 0L else length(phi)
  q <- if (all(theta == 0)) 0L else length(theta)

  # Compute ACF using ARMAacf (note: opposite sign for theta)
  if (p == 0 && q == 0) {
    # White noise
    acf_values <- c(1, rep(0, lag_max))
  } else {
    acf_values <- as.numeric(ARMAacf(ar = phi, ma = -theta, lag.max = lag_max))
  }

  # Compute process variance to get ACV
  process_var <- .arma_variance(phi, theta, vara, p, q)
  acv_values <- acf_values * process_var

  # Plot if requested
  if (plot) {
    k <- 0:lag_max
    plot(k, acf_values, type = "h", ylim = c(-1, 1),
         xlab = "Lag", ylab = "Autocorrelation",
         main = "True Autocorrelations")
    abline(h = 0)
    invisible(list(acf = acf_values, acv = acv_values))
  } else {
    list(acf = acf_values, acv = acv_values)
  }
}
