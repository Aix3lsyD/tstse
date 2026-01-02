#' Forecast GARMA Model
#'
#' Generate forecasts from an estimated GARMA model using the inverse
#' Gegenbauer filter transformation.
#'
#' @param x Numeric vector. The time series to forecast.
#' @param u Numeric vector. Gegenbauer frequency parameter(s) (length 1 or 2).
#' @param lambda Numeric vector. Gegenbauer persistence parameter(s)
#'   (same length as `u`).
#' @param phi Numeric vector. AR coefficients (default 0 = no AR).
#' @param theta Numeric vector. MA coefficients (default 0 = no MA).
#' @param n_ahead Integer. Number of steps ahead to forecast. Default is 10.
#' @param lastn Logical. If `TRUE` (default), use holdout mode: replace the
#'   last `n_ahead` observations with forecasts for validation. If `FALSE`,
#'   forecast `n_ahead` steps into the future.
#' @param plot Logical. If `TRUE` (default), plot the forecasts.
#'
#' @return An object of class `"fore_garma"` with components:
#'   \item{f}{GARMA forecasts (length `n_ahead`)}
#'   \item{ar_f}{AR benchmark forecasts (length `n_ahead`)}
#'   \item{ar_order}{Order of the AR benchmark model}
#'   \item{original}{Original series}
#'   \item{n_ahead}{Number of steps ahead}
#'   \item{lastn}{Whether holdout mode was used}
#'
#' @details
#' The forecasting algorithm uses the inverse Gegenbauer filter approach:
#' \enumerate{
#'   \item **Forward Transform**: Remove Gegenbauer factor(s) by applying
#'     \eqn{(1 - 2uB + B^2)^{\lambda}} to convert the GARMA series to
#'     approximately ARMA.
#'   \item **Forecast**: Generate forecasts in the transformed (ARMA) domain.
#'   \item **Back-Transform**: Reapply the Gegenbauer factor(s) by applying
#'     \eqn{(1 - 2uB + B^2)^{-\lambda}} to return to the original scale.
#' }
#'
#' For two-factor models (k = 2), the transformations are applied sequentially.
#'
#' An AR benchmark forecast is computed for comparison, using AIC to select
#' the AR order (maximum 8).
#'
#' @section Holdout vs Future Mode:
#' \itemize{
#'   \item `lastn = TRUE`: Replaces the last `n_ahead` observations with
#'     forecasts, useful for model validation.
#'   \item `lastn = FALSE`: Extends the series with `n_ahead` future forecasts.
#' }
#'
#' @seealso [gen_garma()] for generating GARMA realizations,
#'   [est_garma()] for parameter estimation.
#'
#' @examples
#' # Generate GARMA data
#' x <- gen_garma(n = 200, u = 0.8, lambda = 0.3, phi = 0.5,
#'                plot = FALSE, seed = 123)
#'
#' # Holdout forecast (validation mode)
#' fit <- fore_garma(x, u = 0.8, lambda = 0.3, phi = 0.5, n_ahead = 10)
#'
#' # Future forecast
#' fit2 <- fore_garma(x, u = 0.8, lambda = 0.3, phi = 0.5,
#'                    n_ahead = 20, lastn = FALSE)
#'
#' @export
fore_garma <- function(x, u, lambda, phi = 0, theta = 0,
                       n_ahead = 10L, lastn = TRUE, plot = TRUE) {

  # Input validation
  if (!is.numeric(x) || length(x) < 20L) {
    stop("`x` must be a numeric vector with at least 20 observations.", call. = FALSE)
  }

  if (!is.numeric(u)) {
    stop("`u` must be numeric.", call. = FALSE)
  }
  if (!is.numeric(lambda)) {
    stop("`lambda` must be numeric.", call. = FALSE)
  }
  if (length(u) != length(lambda)) {
    stop("`u` and `lambda` must have the same length.", call. = FALSE)
  }
  k <- length(u)
  if (k > 2L) {
    stop("Only 1 or 2 Gegenbauer factors are supported.", call. = FALSE)
  }

  if (!is.numeric(phi)) {
    stop("`phi` must be numeric.", call. = FALSE)
  }
  if (!is.numeric(theta)) {
    stop("`theta` must be numeric.", call. = FALSE)
  }

  if (!is.numeric(n_ahead) || length(n_ahead) != 1L || n_ahead < 1L) {
    stop("`n_ahead` must be a positive integer.", call. = FALSE)
  }
  n_ahead <- as.integer(n_ahead)

  if (!is.logical(lastn) || length(lastn) != 1L || is.na(lastn)) {
    stop("`lastn` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.logical(plot) || length(plot) != 1L || is.na(plot)) {
    stop("`plot` must be TRUE or FALSE.", call. = FALSE)
  }

  # Determine effective AR/MA orders
  p <- if (all(phi == 0)) 0L else length(phi)
  q <- if (all(theta == 0)) 0L else length(theta)

  n <- length(x)
  n_back <- 500L
  ar_p <- 20L
  x_mean <- mean(x)

  # Center data and prepare backcasted sequence
  nc <- n + n_back
  z <- c(rep(0, n_back), x - x_mean)

  # Backcast using AR(20) Yule-Walker
  ar_fit_back <- ar.yw(z[-(1:n_back)], order.max = ar_p, aic = FALSE, demean = FALSE)
  ar_phi_back <- ar_fit_back$ar

  for (i in n_back:1) {
    z[i] <- sum(ar_phi_back * z[(i + 1):(i + ar_p)])
  }

  # Forecast step direction
  forecast_step <- if (lastn) -n_ahead else n_ahead

  # AR benchmark (select order by AIC, max 8)
  x_ar_est <- ar.yw(x, aic = TRUE, order.max = 8L)
  ar_order <- x_ar_est$order

  # Single factor (k = 1)
  if (k == 1L) {
    result <- .fore_garma_k1(z, x, x_mean, u, lambda, phi, theta, p, q,
                             n, n_back, ar_p, forecast_step, x_ar_est)
  } else {
    # Two factors (k = 2)
    result <- .fore_garma_k2(z, x, x_mean, u, lambda, phi, theta, p, q,
                             n, n_back, ar_p, forecast_step, x_ar_est)
  }

  x_forecast <- result$x_forecast
  x_ar_fore <- result$x_ar_fore
  n2 <- result$n2

  # Extract final forecasts
  if (forecast_step < 0) {
    # Holdout mode
    fore_idx <- (n2 + forecast_step + 1):n2
    garma_fore <- x_mean + x_forecast[fore_idx]
    ar_fore <- x_ar_fore[fore_idx]
  } else {
    # Future mode
    fore_idx <- (n2 - forecast_step + 1):n2
    garma_fore <- x_mean + x_forecast[fore_idx]
    ar_fore <- x_ar_fore[fore_idx]
  }

  # Plot if requested
  if (plot) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))

    par(mar = c(4, 4, 3, 1))

    if (forecast_step < 0) {
      # Holdout mode plot
      t <- seq_len(n)
      plot(t, x, type = "o", pch = 16, cex = 0.7,
           xlab = "Time", ylab = "Value",
           main = "GARMA Forecast (Holdout Mode)",
           ylim = range(c(x, garma_fore, ar_fore)))
      points((n + forecast_step + 1):n, garma_fore, type = "o",
             lty = 2, pch = 1, cex = 0.6, col = "blue")
      points((n + forecast_step + 1):n, ar_fore, type = "o",
             lty = 3, pch = 2, cex = 0.6, col = "red")
      legend("topleft", legend = c("Actual", "GARMA", paste0("AR(", ar_order, ")")),
             lty = c(1, 2, 3), pch = c(16, 1, 2), col = c("black", "blue", "red"),
             cex = 0.8, bty = "n")
    } else {
      # Future mode plot
      t_all <- seq_len(n + forecast_step)
      y_all <- c(x, rep(NA, forecast_step))

      plot(t_all, y_all, type = "o", pch = 16, cex = 0.7,
           xlab = "Time", ylab = "Value",
           main = "GARMA Forecast (Future Mode)",
           ylim = range(c(x, garma_fore, ar_fore), na.rm = TRUE))
      points((n + 1):(n + forecast_step), garma_fore, type = "o",
             lty = 2, pch = 1, cex = 0.6, col = "blue")
      points((n + 1):(n + forecast_step), ar_fore, type = "o",
             lty = 3, pch = 2, cex = 0.6, col = "red")
      # Connect last point to forecasts
      segments(n, x[n], n + 1, garma_fore[1], lty = 2, col = "blue")
      segments(n, x[n], n + 1, ar_fore[1], lty = 3, col = "red")

      legend("topleft", legend = c("Actual", "GARMA", paste0("AR(", ar_order, ")")),
             lty = c(1, 2, 3), pch = c(16, 1, 2), col = c("black", "blue", "red"),
             cex = 0.8, bty = "n")
    }
  }

  structure(
    list(
      f = garma_fore,
      ar_f = ar_fore,
      ar_order = ar_order,
      original = x,
      n_ahead = n_ahead,
      lastn = lastn
    ),
    class = "fore_garma"
  )
}


# Internal: Single factor GARMA forecast
.fore_garma_k1 <- function(z, x, x_mean, u, lambda, phi, theta, p, q,
                           n, n_back, ar_p, forecast_step, x_ar_est) {
  nc <- n + n_back

  # Forward transform: remove Gegenbauer factor
  C1 <- gegenb(u = u, d = -lambda, n_coef = nc)
  w <- numeric(n)
  for (i in seq_len(n)) {
    c_idx <- seq_len(n_back + i - 1)
    z_idx <- (n_back + i):2
    w[i] <- sum(C1[c_idx] * z[z_idx])
  }
  w <- w - mean(w)
  w_forecast <- w

  # AR benchmark forecast
  x_ar_fore <- x
  ar_order <- x_ar_est$order

  if (forecast_step < 0) {
    # Holdout mode: replace last |forecast_step| with AR forecasts
    if (ar_order > 0) {
      ma_const <- (1 - sum(x_ar_est$ar)) * x_mean
      for (i in (n + forecast_step + 1):n) {
        if (p > 0) {
          w_forecast[i] <- sum(w_forecast[(i - 1):(i - p)] * phi)
        } else {
          w_forecast[i] <- 0
        }
        x_ar_fore[i] <- ma_const + sum(x_ar_fore[(i - 1):(i - ar_order)] * x_ar_est$ar)
      }
    }
  } else {
    # Future mode
    if (p > 0 || q > 0) {
      pred_w <- predict(arima(w_forecast, order = c(p, 0, q)), n.ahead = forecast_step)
      w_forecast <- c(w_forecast, as.numeric(pred_w$pred))
    } else {
      # No AR/MA, just extend with zeros
      w_forecast <- c(w_forecast, rep(0, forecast_step))
    }
    ar_pred <- predict(x_ar_est, n.ahead = forecast_step)
    x_ar_fore <- c(x_ar_fore, as.numeric(ar_pred$pred))
  }

  # Back-transform: reapply Gegenbauer factor
  n2 <- length(w_forecast)
  nc2 <- n_back + n2
  z2 <- c(rep(0, n_back), w_forecast)

  # Backcast z2
  ar_phi2 <- ar.yw(z2[-(1:n_back)], order.max = ar_p, aic = FALSE, demean = FALSE)$ar
  for (i in n_back:1) {
    z2[i] <- sum(ar_phi2 * z2[(i + 1):(i + ar_p)])
  }

  # Apply Gegenbauer factor (positive lambda)
  C2 <- gegenb(u = u, d = lambda, n_coef = nc2)
  x_forecast <- numeric(n2)
  for (i in seq_len(n2)) {
    c_idx <- seq_len(n_back + i - 1)
    z_idx <- (n_back + i):2
    x_forecast[i] <- sum(C2[c_idx] * z2[z_idx])
  }

  list(x_forecast = x_forecast, x_ar_fore = x_ar_fore, n2 = n2)
}


# Internal: Two-factor GARMA forecast
.fore_garma_k2 <- function(z, x, x_mean, u, lambda, phi, theta, p, q,
                           n, n_back, ar_p, forecast_step, x_ar_est) {
  nc <- n + n_back

  # Remove first Gegenbauer factor
  C1 <- gegenb(u = u[1], d = -lambda[1], n_coef = nc)
  w1 <- numeric(n)
  for (i in seq_len(n)) {
    c_idx <- seq_len(n_back + i - 1)
    z_idx <- (n_back + i):2
    w1[i] <- sum(C1[c_idx] * z[z_idx])
  }
  w1 <- w1 - mean(w1)

  # Backcast w1
  w1_ext <- c(rep(0, n_back), w1)
  ar_phi2 <- ar.yw(w1_ext[-(1:n_back)], order.max = ar_p, aic = FALSE, demean = FALSE)$ar
  for (i in n_back:1) {
    w1_ext[i] <- sum(ar_phi2 * w1_ext[(i + 1):(i + ar_p)])
  }

  # Remove second Gegenbauer factor
  C2 <- gegenb(u = u[2], d = -lambda[2], n_coef = nc)
  w <- numeric(n)
  for (i in seq_len(n)) {
    c_idx <- seq_len(n_back + i - 1)
    z_idx <- (n_back + i):2
    w[i] <- sum(C2[c_idx] * w1_ext[z_idx])
  }
  w <- w - mean(w)
  w_forecast <- w

  # AR benchmark
  x_ar_fore <- x
  ar_order <- x_ar_est$order

  if (forecast_step < 0) {
    # Holdout mode
    for (i in (n + forecast_step + 1):n) {
      if (p > 0) {
        w_forecast[i] <- sum(w_forecast[(i - 1):(i - p)] * phi)
      }
      if (ar_order > 0) {
        x_ar_fore[i] <- sum(x_ar_fore[(i - 1):(i - ar_order)] * x_ar_est$ar)
      }
    }
  } else {
    # Future mode
    if (p > 0 || q > 0) {
      pred_w <- predict(arima(w_forecast, order = c(p, 0, q)), n.ahead = forecast_step)
      w_forecast <- c(w_forecast, as.numeric(pred_w$pred))
    } else {
      w_forecast <- c(w_forecast, rep(0, forecast_step))
    }
    ar_pred <- predict(x_ar_est, n.ahead = forecast_step)
    x_ar_fore <- c(x_ar_fore, as.numeric(ar_pred$pred))
  }

  # Back-transform step 1: add second factor back
  n2 <- length(w_forecast)
  nc2 <- n_back + n2
  z2 <- c(rep(0, n_back), w_forecast)

  ar_phi3 <- ar.yw(z2[-(1:n_back)], order.max = ar_p, aic = FALSE, demean = FALSE)$ar
  for (i in n_back:1) {
    z2[i] <- sum(ar_phi3 * z2[(i + 1):(i + ar_p)])
  }

  C22 <- gegenb(u = u[2], d = lambda[2], n_coef = nc2)
  x_forecast1 <- numeric(n2)
  for (i in seq_len(n2)) {
    c_idx <- seq_len(n_back + i - 1)
    z_idx <- (n_back + i):2
    x_forecast1[i] <- sum(C22[c_idx] * z2[z_idx])
  }

  # Back-transform step 2: add first factor back
  z22 <- c(rep(0, n_back), x_forecast1)

  ar_phi4 <- ar.yw(z22[-(1:n_back)], order.max = ar_p, aic = FALSE, demean = FALSE)$ar
  for (i in n_back:1) {
    z22[i] <- sum(ar_phi4 * z22[(i + 1):(i + ar_p)])
  }

  C12 <- gegenb(u = u[1], d = lambda[1], n_coef = nc2)
  x_forecast <- numeric(n2)
  for (i in seq_len(n2)) {
    c_idx <- seq_len(n_back + i - 1)
    z_idx <- (n_back + i):2
    x_forecast[i] <- sum(C12[c_idx] * z22[z_idx])
  }

  list(x_forecast = x_forecast, x_ar_fore = x_ar_fore, n2 = n2)
}


#' @export
print.fore_garma <- function(x, ...) {
  cat("\nGARMA Forecast\n")
  cat("==============\n\n")
  cat("Mode:", if (x$lastn) "Holdout" else "Future", "\n")
  cat("Steps ahead:", x$n_ahead, "\n")
  cat("AR benchmark order:", x$ar_order, "\n\n")

  cat("GARMA Forecasts:\n")
  print(round(x$f, 4))

  cat("\nAR Forecasts:\n")
  print(round(x$ar_f, 4))

  invisible(x)
}
