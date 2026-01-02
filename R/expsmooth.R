#' Simple Exponential Smoothing
#'
#' Performs simple exponential smoothing on a time series with optional
#' forecasting. This is a convenience wrapper around [stats::HoltWinters()]
#' with `beta = FALSE` and `gamma = FALSE`.
#'
#' @param x Numeric vector or time series. The input data to smooth.
#' @param alpha Numeric scalar between 0 and 1, or `NULL`. The smoothing
#'   parameter. Larger values give more weight to recent observations. If
#'   `NULL` (default), the optimal value is estimated by minimizing the
#'   squared prediction error.
#' @param n_ahead Integer. Number of periods to forecast ahead. Default is 0
#'   (no forecasting, just smoothing).
#' @param plot Logical. If `TRUE` (default), plot the original series and
#'   smoothed values (and forecasts if `n_ahead > 0`).
#'
#' @return A list with class `"expsmooth"` containing:
#'   \item{alpha}{The smoothing parameter used (estimated if not supplied).}
#'   \item{fitted}{Numeric vector of smoothed/fitted values for the original
#'     series (length equal to `length(x)`).}
#'   \item{forecast}{Numeric vector of forecasted values (length `n_ahead`),
#'     or `NULL` if `n_ahead = 0`.}
#'   \item{hw}{The underlying [stats::HoltWinters] object for further analysis.}
#'
#' @details
#' Simple exponential smoothing computes the smoothed value at time \eqn{t} as:
#' \deqn{\hat{x}_t = \alpha x_{t-1} + (1 - \alpha) \hat{x}_{t-1}}
#'
#' This is equivalent to an ARIMA(0,1,1) model. It is appropriate for series
#' without trend or seasonality.
#'
#' The smoothing parameter \eqn{\alpha} controls the rate of decay:
#' \itemize{
#'   \item \eqn{\alpha} close to 1: More weight on recent observations (responsive)
#'   \item \eqn{\alpha} close to 0: More weight on past observations (smooth)
#' }
#'
#' @seealso [stats::HoltWinters()] for the underlying implementation,
#'   [fore_arma()] for ARMA-based forecasting.
#'
#' @examples
#' # Smooth a noisy series
#' set.seed(123)
#' x <- cumsum(rnorm(100))
#' result <- expsmooth(x)
#' result$alpha
#'
#' # Fixed alpha
#' result2 <- expsmooth(x, alpha = 0.3)
#'
#' # With forecasting
#' result3 <- expsmooth(x, n_ahead = 10)
#'
#' # No plot
#' result4 <- expsmooth(x, plot = FALSE)
#'
#' @export
expsmooth <- function(x, alpha = NULL, n_ahead = 0L, plot = TRUE) {
  # Input validation
  if (!is.numeric(x) || length(x) < 2L) {
    stop("`x` must be a numeric vector with at least 2 observations.", call. = FALSE)
  }
  if (anyNA(x)) {
    stop("`x` contains NA values.", call. = FALSE)
  }

  if (!is.null(alpha)) {
    if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha)) {
      stop("`alpha` must be a single numeric value or NULL.", call. = FALSE)
    }
    if (alpha <= 0 || alpha >= 1) {
      stop("`alpha` must be between 0 and 1 (exclusive).", call. = FALSE)
    }
  }

  if (!is.numeric(n_ahead) || length(n_ahead) != 1L || is.na(n_ahead)) {
    stop("`n_ahead` must be a single non-negative integer.", call. = FALSE)
  }
  n_ahead <- as.integer(n_ahead)
  if (n_ahead < 0L) {
    stop("`n_ahead` must be non-negative.", call. = FALSE)
  }

  if (!is.logical(plot) || length(plot) != 1L || is.na(plot)) {
    stop("`plot` must be TRUE or FALSE.", call. = FALSE)
  }

  n <- length(x)

  # Convert to ts if not already (HoltWinters requires ts object)
  if (!is.ts(x)) {
    x <- ts(x)
  }

  # Fit simple exponential smoothing model
  # beta = FALSE: no trend component
  # gamma = FALSE: no seasonal component
  hw <- HoltWinters(x, alpha = alpha, beta = FALSE, gamma = FALSE)

  # Extract fitted values
  # HoltWinters$fitted starts at observation 2, so we prepend the first observation
  # The first "fitted" value is conventionally set to x[1] (no prior info)
  fitted_vals <- c(x[1], hw$fitted[, "xhat"])

  # Compute forecasts if requested
  if (n_ahead > 0L) {
    forecast_vals <- as.numeric(predict(hw, n.ahead = n_ahead))
  } else {
    forecast_vals <- NULL
  }

  # Plot if requested
  if (plot) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))

    par(mar = c(4, 4, 2, 1))

    if (n_ahead == 0L) {
      # Plot original and fitted only
      plot(seq_len(n), as.numeric(x), type = "o", pch = 16, cex = 0.8,
           xlab = "Time", ylab = "Value",
           main = paste0("Exponential Smoothing (alpha = ",
                         round(hw$alpha, 4), ")"))
      lines(seq_len(n), fitted_vals, col = "blue", lwd = 2)
      legend("topleft", legend = c("Observed", "Smoothed"),
             col = c("black", "blue"), lty = c(1, 1), pch = c(16, NA),
             lwd = c(1, 2), bty = "n", cex = 0.8)
    } else {
      # Plot original, fitted, and forecasts
      all_vals <- c(fitted_vals, forecast_vals)
      time_all <- seq_len(n + n_ahead)

      ylim <- range(c(as.numeric(x), all_vals))

      plot(seq_len(n), as.numeric(x), type = "o", pch = 16, cex = 0.8,
           xlab = "Time", ylab = "Value",
           xlim = c(1, n + n_ahead), ylim = ylim,
           main = paste0("Exponential Smoothing (alpha = ",
                         round(hw$alpha, 4), ")"))
      lines(seq_len(n), fitted_vals, col = "blue", lwd = 2)
      lines((n + 1):(n + n_ahead), forecast_vals, col = "red", lwd = 2, lty = 2)
      points((n + 1):(n + n_ahead), forecast_vals, col = "red", pch = 17, cex = 0.8)

      legend("topleft", legend = c("Observed", "Smoothed", "Forecast"),
             col = c("black", "blue", "red"), lty = c(1, 1, 2),
             pch = c(16, NA, 17), lwd = c(1, 2, 2), bty = "n", cex = 0.8)
    }
  }

  # Build result
  result <- list(
    alpha = as.numeric(hw$alpha),
    fitted = fitted_vals,
    forecast = forecast_vals,
    hw = hw
  )

  class(result) <- "expsmooth"
  invisible(result)
}


#' @export
print.expsmooth <- function(x, ...) {
  cat("Simple Exponential Smoothing\n")
  cat("----------------------------\n")
  cat("Alpha:", round(x$alpha, 4), "\n")
  cat("Observations:", length(x$fitted), "\n")
  if (!is.null(x$forecast)) {
    cat("Forecasts:", length(x$forecast), "periods ahead\n")
  }
  invisible(x)
}
