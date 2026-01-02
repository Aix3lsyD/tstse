#' ARMA Forecasting
#'
#' Generate point forecasts and prediction intervals for ARMA models.
#'
#' @param x Numeric vector, the time series to forecast.
#' @param phi Numeric vector, AR coefficients (default 0 = none).
#' @param theta Numeric vector, MA coefficients (default 0 = none).
#' @param n_ahead Integer, forecast horizon (default 5).
#' @param lastn Logical, if TRUE forecast from position n-n_ahead (holdout mode).
#' @param plot Logical, whether to plot forecasts (default TRUE).
#' @param alpha Numeric, significance level for prediction intervals (default 0.05).
#' @param limits Logical, whether to compute/show prediction limits (default TRUE).
#'
#' @return An object of class "fore_arma" with components:
#'   \item{f}{Point forecasts (length n_ahead)}
#'   \item{ll}{Lower prediction limits}
#'   \item{ul}{Upper prediction limits}
#'   \item{resid}{Residuals (length n)}
#'   \item{wnv}{Estimated white noise variance}
#'   \item{xbar}{Sample mean of original series}
#'   \item{se}{Standard errors for each forecast horizon}
#'   \item{psi}{Psi weights used for prediction intervals}
#'   \item{rmse}{Root mean squared error (only when lastn=TRUE)}
#'   \item{mad}{Mean absolute deviation (only when lastn=TRUE)}
#' @export
#'
#' @details
#' This is a convenience wrapper around \code{\link{fore_arima}} with \code{d=0}
#' and \code{s=0} for stationary ARMA models.
#'
#' When \code{lastn = TRUE}, forecasts are generated from position \code{n - n_ahead}
#' rather than from position \code{n}, allowing for holdout validation. In this mode,
#' RMSE (root mean squared error) and MAD (mean absolute deviation) are computed
#' by comparing forecasts against the held-out actual values.
#'
#' @seealso \code{\link{fore_arima}} for ARIMA/SARIMA forecasting.
#'
#' @examples
#' # Generate AR(2) data
#' set.seed(123)
#' x <- gen_arma(n = 100, phi = c(1.5, -0.75), plot = FALSE)
#'
#' # Forecast 5 steps ahead
#' fore_arma(x, phi = c(1.5, -0.75), n_ahead = 5)
#'
#' # Holdout validation (forecast last 10 values)
#' result <- fore_arma(x, phi = c(1.5, -0.75), n_ahead = 10, lastn = TRUE)
#' result$rmse
#' result$mad
fore_arma <- function(x, phi = 0, theta = 0, n_ahead = 5L, lastn = FALSE,
                      plot = TRUE, alpha = 0.05, limits = TRUE) {

  # Call fore_arima with d=0, s=0
  result <- fore_arima(x, phi = phi, theta = theta, d = 0L, s = 0L,
                       n_ahead = n_ahead, lastn = lastn, plot = plot,
                       alpha = alpha, limits = limits)

  # Add holdout metrics when lastn=TRUE
  if (lastn) {
    n <- length(x)
    actuals <- x[(n - n_ahead + 1):n]
    result$rmse <- sqrt(mean((result$f - actuals)^2))
    result$mad <- mean(abs(result$f - actuals))
  }

  class(result) <- c("fore_arma", "fore_arima")
  invisible(result)
}


#' @rdname fore_arma
#' @param x An object of class "fore_arma".
#' @param ... Additional arguments (ignored).
#' @export
print.fore_arma <- function(x, ...) {

  cat("ARMA Forecast\n")
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

  if (!is.null(x$rmse)) {
    cat("\nHoldout validation:\n")
    cat("  RMSE:", round(x$rmse, 4), "\n")
    cat("  MAD:", round(x$mad, 4), "\n")
  }

  invisible(x)
}
