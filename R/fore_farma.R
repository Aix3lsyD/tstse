#' Forecast Fractional ARMA (FARMA) Process
#'
#' Generate forecasts from a Fractional ARMA model, which incorporates
#' long memory through fractional differencing.
#'
#' @param x Numeric vector, the time series to forecast.
#' @param d Numeric, fractional differencing parameter. Typically 0 < d < 0.5
#'   for stationary long-memory processes.
#' @param phi Numeric vector, AR coefficients (default 0 = none).
#' @param theta Numeric vector, MA coefficients (default 0 = none).
#' @param n_ahead Integer, forecast horizon (default 10).
#' @param lastn Logical, if TRUE forecast from position n-n_ahead (holdout mode).
#' @param n_back Integer, number of backcasted values (default 500).
#' @param plot Logical, whether to plot forecasts (default TRUE).
#'
#' @return A list with components:
#'   \item{farma_fore}{FARMA point forecasts (length n_ahead)}
#'   \item{ar_fore}{AR model point forecasts for comparison}
#'   \item{ar_fit_order}{Order of the fitted AR model}
#' @export
#'
#' @details
#' FARMA (Fractional ARMA) models extend ARMA by incorporating long memory
#' through fractional differencing. The model is:
#' \deqn{(1-B)^d \phi(B) X_t = \theta(B) a_t}
#'
#' The forecasting algorithm:
#' \enumerate{
#'   \item Backcast using a high-order AR model
#'   \item Apply fractional differencing using Gegenbauer coefficients to
#'     remove long memory
#'   \item Forecast the resulting ARMA process
#'   \item Reverse the fractional differencing to obtain final forecasts
#' }
#'
#' This is a special case of GARMA forecasting with the Gegenbauer frequency
#' parameter \eqn{u = 1} (corresponding to frequency 0).
#'
#' An AR model is also fitted for comparison purposes.
#'
#' @seealso [gegenb()] for Gegenbauer coefficients,
#'   [fore_arma()] for stationary ARMA forecasting.
#'
#' @examples
#' # Generate a long-memory process
#' set.seed(123)
#' x <- cumsum(rnorm(200))  # Simple random walk as example
#'
#' # Forecast with FARMA model (d=0.3)
#' result <- fore_farma(x, d = 0.3, n_ahead = 10)
#'
#' # Holdout validation
#' result <- fore_farma(x, d = 0.3, n_ahead = 10, lastn = TRUE, plot = FALSE)
fore_farma <- function(x,
                       d,
                       phi = 0,
                       theta = 0,
                       n_ahead = 10L,
                       lastn = TRUE,
                       n_back = 500L,
                       plot = TRUE) {

 # Input validation
 if (!is.numeric(x) || length(x) < 10) {
   stop("`x` must be a numeric vector with at least 10 observations")
 }
 if (!is.numeric(d) || length(d) != 1) {
   stop("`d` must be a single numeric value")
 }
 if (!is.numeric(phi)) stop("`phi` must be numeric")
 if (!is.numeric(theta)) stop("`theta` must be numeric")
 if (!is.numeric(n_ahead) || n_ahead < 1) {
   stop("`n_ahead` must be a positive integer")
 }
 if (!is.numeric(n_back) || n_back < 1) {
   stop("`n_back` must be a positive integer")
 }

 n_ahead <- as.integer(n_ahead)
 n_back <- as.integer(n_back)

 # Determine AR/MA orders
 p <- if (all(phi == 0)) 0L else length(phi)
 q <- if (all(theta == 0)) 0L else length(theta)

 n <- length(x)
 x_mean <- mean(x)

 # Adjust d for Gegenbauer parameterization
 # FARMA uses (1-B)^d which corresponds to Gegenbauer with u=1 and d_geg = d/2
 d_geg <- d / 2

 # Set forecast step (negative for lastn mode)
 forecast_step <- if (lastn) -n_ahead else n_ahead

 # Backcast and center the series
 nc <- n + n_back
 z <- c(rep(0, n_back), x - x_mean)

 # Fit high-order AR model for backcasting
 ar_order <- 20L
 ar_back <- ar.yw(z[-(1:n_back)], order.max = ar_order, aic = FALSE, demean = FALSE)$ar

 # Backcast: z[n_back], z[n_back-1], ..., z[1]
 for (i in n_back:1) {
   z[i] <- sum(ar_back * z[(i + 1):(i + ar_order)])
 }

 # Apply fractional differencing using Gegenbauer with u=1
 # This removes the long memory component
 C1 <- gegenb(u = 1, d = -d_geg, n_coef = nc)
 w <- .apply_gegenbauer_filter(z, C1, n, n_back)
 w <- w - mean(w)
 w_forecast <- w

 # Fit AR model to original series for comparison
 ar_fit <- ar.yw(x, aic = TRUE, order.max = 8L)
 ar_order_fit <- ar_fit$order

 # Compute AR constant for mean adjustment
 if (ar_order_fit > 0) {
   ma_constant <- (1 - sum(ar_fit$ar)) * x_mean
 } else {
   ma_constant <- x_mean
 }

 x_ar_fore <- x

 # Generate forecasts
 if (forecast_step < 0) {
   # Lastn mode: update last |forecast_step| values using provided coefficients
   for (i in (n + forecast_step + 1):n) {
     # Apply AR filter (original behavior: sets to 0 when p=0)
     if (p > 0) {
       w_forecast[i] <- sum(w_forecast[(i - 1):(i - p)] * phi)
     } else {
       # Match original: when phi=0, the sum evaluates to 0
       w_forecast[i] <- 0
     }
     if (ar_order_fit > 0) {
       x_ar_fore[i] <- ma_constant + sum(x_ar_fore[(i - 1):(i - ar_order_fit)] * ar_fit$ar)
     } else {
       x_ar_fore[i] <- ma_constant
     }
   }
 } else {
   # Future forecast mode: re-estimate ARMA and predict
   # Note: Original tswge re-estimates coefficients rather than using provided ones
   arma_fit <- arima(w_forecast, order = c(p, 0L, q), include.mean = FALSE)
   w_pred <- predict(arma_fit, n.ahead = n_ahead)$pred
   w_forecast <- c(w_forecast, as.numeric(w_pred))

   if (ar_order_fit > 0) {
     ar_pred <- predict(ar_fit, n.ahead = n_ahead)$pred
   } else {
     ar_pred <- rep(x_mean, n_ahead)
   }
   x_ar_fore <- c(x_ar_fore, as.numeric(ar_pred))
 }

 # Backcast the extended w_forecast series
 narma2 <- length(w_forecast)
 nc2 <- n_back + narma2
 z2 <- c(rep(0, n_back), w_forecast)

 ar_back2 <- ar.yw(z2[-(1:n_back)], order.max = ar_order, aic = FALSE, demean = FALSE)$ar
 for (i in n_back:1) {
   z2[i] <- sum(ar_back2 * z2[(i + 1):(i + ar_order)])
 }

 # Reverse fractional differencing to get FARMA forecasts
 C2 <- gegenb(u = 1, d = d_geg, n_coef = nc2)
 x_forecast <- .apply_gegenbauer_filter(z2, C2, narma2, n_back)

 # Extract forecast values and add back mean
 if (forecast_step < 0) {
   idx_start <- narma2 + forecast_step + 1
   idx_end <- narma2
   farma_fore <- x_mean + x_forecast[idx_start:idx_end]
   ar_fore <- x_ar_fore[idx_start:idx_end]
 } else {
   idx_start <- n + 1
   idx_end <- narma2
   farma_fore <- x_mean + x_forecast[idx_start:idx_end]
   ar_fore <- x_ar_fore[idx_start:idx_end]
 }

 # Plot if requested
 if (plot) {
   .plot_farma_forecast(x, farma_fore, ar_fore, forecast_step, n, narma2, x_mean)
 }

 # Return results
 result <- list(
   farma_fore = as.numeric(farma_fore),
   ar_fore = as.numeric(ar_fore),
   ar_fit_order = ar_order_fit
 )

 invisible(result)
}


#' Apply Gegenbauer filter (convolution)
#' @noRd
.apply_gegenbauer_filter <- function(z, C, n_out, n_back) {
  result <- numeric(n_out)
  for (i in seq_len(n_out)) {
    idx <- n_back + i
    result[i] <- sum(C[1:(n_back + i - 1)] * z[idx:2])
  }
  result
}


#' Plot FARMA forecast
#' @noRd
.plot_farma_forecast <- function(x, farma_fore, ar_fore, forecast_step,
                                  n, narma2, x_mean) {
  op <- par(mar = c(4, 4, 3, 1))
  on.exit(par(op), add = TRUE)

  n_ahead <- abs(forecast_step)

  if (forecast_step < 0) {
    # Lastn mode
    t <- seq_len(n)
    ylim <- range(c(x, x_mean + farma_fore, ar_fore))

    plot(t, x, type = "o", pch = 16, cex = 0.7, lty = "solid",
         xlab = "Time", ylab = "Value", ylim = ylim,
         main = "Realization (Solid), FARMA Forecast (Dashed), AR Forecast (Dotted)")

    fore_idx <- (n - n_ahead + 1):n
    points(fore_idx, farma_fore, type = "o", lty = "dashed", pch = 1, cex = 0.6)
    points(fore_idx, ar_fore, type = "o", lty = "dotted", pch = 2, cex = 0.6)

    # Connect last observed to first forecast
    n1 <- n - n_ahead
    if (n1 > 0) {
      lines(c(n1, n1 + 1), c(x[n1], farma_fore[1]), lty = "dashed")
      lines(c(n1, n1 + 1), c(x[n1], ar_fore[1]), lty = "dotted")
    }
  } else {
    # Future forecast mode
    t <- seq_len(n)
    ylim <- range(c(x, x_mean + farma_fore, ar_fore))
    xlim <- c(1, n + n_ahead)

    plot(t, x, type = "o", pch = 16, cex = 0.7, lty = "solid",
         xlab = "Time", ylab = "Value", ylim = ylim, xlim = xlim,
         main = "Realization (Solid), FARMA Forecast (Dashed), AR Forecast (Dotted)")

    fore_idx <- (n + 1):(n + n_ahead)
    points(fore_idx, farma_fore, type = "o", lty = "dashed", pch = 1, cex = 0.7)
    points(fore_idx, ar_fore, type = "o", lty = "dotted", pch = 2, cex = 0.7)

    # Connect last observed to first forecast
    lines(c(n, n + 1), c(x[n], farma_fore[1]), lty = "dashed")
    lines(c(n, n + 1), c(x[n], ar_fore[1]), lty = "dotted")
  }
}
