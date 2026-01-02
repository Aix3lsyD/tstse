#' Generate ARMA Process
#'
#' Generate a realization from an ARMA(p,q) process.
#' This is a convenience wrapper around \code{\link{gen_arima}} with d=0.
#'
#' @param n Integer, length of output series.
#' @param phi Numeric vector, AR coefficients (default 0 = none).
#' @param theta Numeric vector, MA coefficients (default 0 = none).
#' @param mu Numeric, theoretical mean (default 0).
#' @param vara Numeric, white noise variance (default 1).
#' @param plot Logical, whether to plot the series (default TRUE).
#' @param seed Integer, random seed for reproducibility (default NULL).
#'
#' @return Numeric vector of length n.
#' @export
#'
#' @examples
#' # AR(1)
#' x <- gen_arma(n = 200, phi = 0.7)
#'
#' # MA(2)
#' x <- gen_arma(n = 200, theta = c(0.5, -0.3))
#'
#' # ARMA(2,1) with mean 10
#' x <- gen_arma(n = 200, phi = c(1.6, -0.9), theta = 0.4, mu = 10)
gen_arma <- function(n,
                     phi = 0,
                     theta = 0,
                     mu = 0,
                     vara = 1,
                     plot = TRUE,
                     seed = NULL) {

  gen_arima(
    n = n,
    phi = phi,
    theta = theta,
    d = 0L,
    s = 0L,
    mu = mu,
    vara = vara,
    plot = plot,
    seed = seed
  )
}
