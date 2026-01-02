#' ARIMA/SARIMA Forecasting
#'
#' Generate point forecasts and prediction intervals for ARIMA models.
#'
#' @param x Numeric vector, the time series to forecast.
#' @param phi Numeric vector, AR coefficients (default 0 = none).
#' @param theta Numeric vector, MA coefficients (default 0 = none).
#' @param d Integer, order of differencing (default 0, max 3).
#' @param s Integer, seasonal period (default 0 = none).
#' @param n_ahead Integer, forecast horizon (default 5).
#' @param lastn Logical, if TRUE forecast from position n-n_ahead (holdout mode).
#' @param plot Logical, whether to plot forecasts (default TRUE).
#' @param alpha Numeric, significance level for prediction intervals (default 0.05).
#' @param limits Logical, whether to compute/show prediction limits (default TRUE).
#'
#' @return An object of class "fore_arima" with components:
#'   \item{f}{Point forecasts (length n_ahead)}
#'   \item{ll}{Lower prediction limits}
#'   \item{ul}{Upper prediction limits}
#'   \item{resid}{Residuals (length n)}
#'   \item{wnv}{Estimated white noise variance}
#'   \item{xbar}{Sample mean of original series}
#'   \item{se}{Standard errors for each forecast horizon}
#'   \item{psi}{Psi weights used for prediction intervals}
#' @export
#'
#' @details
#' This function generates forecasts for ARIMA(p,d,q) and seasonal ARIMA models
#' using the Box-Jenkins approach. It is a convenience wrapper around
#' \code{\link{fore_aruma}} with \code{lambda = 0}.
#'
#' The model is specified via:
#' \itemize{
#'   \item \code{phi}: AR coefficients for the stationary component
#'   \item \code{theta}: MA coefficients
#'   \item \code{d}: Order of differencing (supports 0, 1, 2, or 3)
#'   \item \code{s}: Seasonal period (e.g., 12 for monthly data with yearly cycle)
#' }
#'
#' When \code{lastn = TRUE}, forecasts are generated from position \code{n - n_ahead}
#' rather than from position \code{n}, allowing for holdout validation.
#'
#' @seealso \code{\link{fore_aruma}}, \code{\link{fore_arma}}
#'
#' @examples
#' # Generate AR(2) data
#' set.seed(123)
#' x <- gen_arma(n = 100, phi = c(1.5, -0.75), plot = FALSE)
#'
#' # Forecast 10 steps ahead
#' fore_arima(x, phi = c(1.5, -0.75), n_ahead = 10)
#'
#' # ARIMA(1,1,0) forecast
#' fore_arima(AirPassengers, phi = 0.8, d = 1, n_ahead = 12)
fore_arima <- function(x, phi = 0, theta = 0, d = 0L, s = 0L,
                       n_ahead = 5L, lastn = FALSE, plot = TRUE,
                       alpha = 0.05, limits = TRUE) {

  # Call fore_aruma with lambda = 0
  result <- fore_aruma(x, phi = phi, theta = theta, d = d, s = s, lambda = 0,
                       n_ahead = n_ahead, lastn = lastn, plot = plot,
                       alpha = alpha, limits = limits)

  class(result) <- c("fore_arima", "fore_aruma")
  invisible(result)
}


#' @rdname fore_arima
#' @param x An object of class "fore_arima".
#' @param ... Additional arguments (ignored).
#' @export
print.fore_arima <- function(x, ...) {

  cat("ARIMA Forecast\n")
  cat("Horizon:", length(x$f), "\n\n")

  df <- data.frame(
    h = seq_along(x$f),
    Forecast = round(x$f, 4),
    Lower = round(x$ll, 4),
    Upper = round(x$ul, 4),
    SE = round(x$se, 4)
  )
  print(df, row.names = FALSE)

  cat("\nWhite noise variance:", round(x$wnv, 6), "\n")
  cat("Series mean:", round(x$xbar, 4), "\n")

  invisible(x)
}
