#' Backcast Residuals for ARMA
#'
#' Compute residuals using the backcasting procedure for ARMA models.
#'
#' @param x Numeric vector, the time series.
#' @param phi Numeric vector, AR coefficients (default 0 = no AR).
#' @param theta Numeric vector, MA coefficients (default 0 = no MA).
#' @param n_back Integer, how far back to backcast (default 50).
#'
#' @return Numeric vector of length n containing residual estimates.
#' @export
#'
#' @examples
#' set.seed(42)
#' x <- arima.sim(model = list(ar = 0.7, ma = 0.3), n = 100)
#' resid <- backcast(x, phi = 0.7, theta = 0.3)
backcast <- function(x, phi = 0, theta = 0, n_back = 50L) {

  # Input validation
  if (!is.numeric(x) || length(x) == 0) {
    stop("`x` must be a non-empty numeric vector")
  }
  if (!is.numeric(phi)) {
    stop("`phi` must be numeric")
  }
  if (!is.numeric(theta)) {
    stop("`theta` must be numeric")
  }
  if (!is.numeric(n_back) || n_back < 1) {
    stop("`n_back` must be a positive integer
")
  }

  # Call Rcpp implementation
  backcast_cpp(
    x = as.double(x),
    phi = as.double(phi),
    theta = as.double(theta),
    n_back = as.integer(n_back)
  )
}
