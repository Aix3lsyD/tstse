#' Forecast G-Lambda Model
#'
#' Generate forecasts using the G-lambda time transformation method.
#' Transforms to dual time scale, forecasts in the transformed domain,
#' then interpolates back to the original scale.
#'
#' @param x Numeric vector. The time series to forecast.
#' @param lambda Numeric. The G-lambda transformation parameter.
#'   Use `lambda = 0` for geometric spacing, positive values for power-law.
#' @param offset Numeric. Shift parameter for the transformation. Default is 60.
#' @param phi Numeric vector. AR coefficients for the model in dual domain.
#'   Default is 0 (no AR).
#' @param n_ahead Integer. Number of steps ahead to forecast. Default is 10.
#' @param lastn Logical. If `TRUE` (default), use holdout mode: replace the
#'   last `n_ahead` observations with forecasts. If `FALSE`, forecast
#'   `n_ahead` steps into the future.
#' @param plot Logical. If `TRUE` (default), plot the forecasts.
#'
#' @return An object of class `"fore_glambda"` with components:
#'   \item{f}{G-lambda forecasts (length `n_ahead`)}
#'   \item{ar_f}{AR benchmark forecasts (length `n_ahead`)}
#'   \item{ar_order}{Order of the AR benchmark model}
#'   \item{original}{Original series}
#'   \item{n_ahead}{Number of steps ahead}
#'   \item{lastn}{Whether holdout mode was used}
#'
#' @details
#' The G-lambda forecasting algorithm:
#' \enumerate{
#'   \item Transform data to dual time scale using [trans_to_dual()]
#'   \item Forecast in the dual domain using AR model
#'   \item Re-interpolate forecasts back to original time scale
#' }
#'
#' An AR benchmark forecast (using AIC selection) is computed for comparison.
#'
#' @section Lambda Values:
#' \itemize{
#'   \item `lambda = 0`: Geometric (exponential) time transformation
#'   \item `lambda = 1`: Linear transformation (near identity)
#'   \item `lambda > 1`: Power-law with acceleration
#' }
#'
#' @seealso [trans_to_dual()] for the transformation,
#'   [est_glambda()] for parameter estimation,
#'   [gen_glambda()] for generating G-lambda processes.
#'
#' @examples
#' # Generate data with time-varying structure
#' set.seed(123)
#' x <- gen_glambda(n = 200, lambda = 0.5, phi = 0.7, plot = FALSE)
#'
#' # Holdout forecast
#' fit <- fore_glambda(x, lambda = 0.5, offset = 20, phi = 0.7, n_ahead = 10)
#'
#' # Future forecast
#' fit2 <- fore_glambda(x, lambda = 0.5, offset = 20, phi = 0.7,
#'                      n_ahead = 20, lastn = FALSE)
#'
#' @export
fore_glambda <- function(x, lambda = 0, offset = 60, phi = 0,
                         n_ahead = 10L, lastn = TRUE, plot = TRUE) {

  # Input validation
  if (!is.numeric(x) || length(x) < 20L) {
    stop("`x` must be a numeric vector with at least 20 observations.", call. = FALSE)
  }

  if (!is.numeric(lambda) || length(lambda) != 1L) {
    stop("`lambda` must be a single numeric value.", call. = FALSE)
  }

  if (!is.numeric(offset) || length(offset) != 1L || offset < 0) {
    stop("`offset` must be a non-negative number.", call. = FALSE)
  }

  if (!is.numeric(phi)) {
    stop("`phi` must be numeric.", call. = FALSE)
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

  # Handle lambda = 0 (use small positive value to avoid singularity)
  lambda_use <- if (lambda == 0) 0.001 else lambda

  n <- length(x)
  shift <- offset
  p <- if (all(phi == 0)) 0L else length(phi)

  # Transform to dual scale
  dual_data <- trans_to_dual(x, lambda = lambda_use, offset = offset,
                             h = 0, plot = FALSE)
  dual_y <- dual_data$int_y
  h <- dual_data$h

  # Internal: re-interpolate from dual to original scale
  reinterpolate <- function(data, start, h, lambda, int_start, len) {
    if (lambda == 0 || abs(lambda) < 0.001) {
      # Geometric spacing
      m <- log(start) / log(h)
      ht <- h^seq(m, m + length(data) - 1)
    } else {
      # Power-law spacing
      m <- (start^lambda - 1) / (h * lambda)
      kk <- seq(m, m + length(data) - 1)
      ht <- (kk * h * lambda + 1)^(1 / lambda)
    }

    tt <- seq(int_start, int_start + len - 1)
    approx(ht, data, xout = tt, rule = 2)$y
  }

  # AR benchmark (AIC selection, max order 8 for future mode, 4 for holdout)
  ar_max <- if (lastn) 4L else 8L
  ar_fit <- aic(x, p = 0:ar_max, q = 0, type = "aic")

  if (!lastn) {
    # Future forecast mode
    # Compute number of dual steps needed
    if (lambda_use != 0) {
      len1 <- ceiling(((n + shift + n_ahead)^lambda_use - (n + shift)^lambda_use) /
                        (h * lambda_use))
    } else {
      len1 <- n_ahead
    }
    len1 <- max(len1, 1L)

    # Forecast in dual domain
    dual_fore <- fore_arma(dual_y, phi = phi, theta = 0,
                           n_ahead = len1, lastn = FALSE, plot = FALSE)

    # Re-interpolate to original scale
    int_start <- n + shift + 1
    combined_dual <- c(dual_y, dual_fore$f)

    glam_fore <- reinterpolate(combined_dual, start = shift + 1, h = h,
                               lambda = lambda_use, int_start = int_start,
                               len = n_ahead)

    # AR benchmark forecast
    ar_fore <- fore_arma(x, phi = ar_fit$phi, theta = 0,
                         n_ahead = n_ahead, lastn = FALSE, plot = FALSE)

  } else {
    # Holdout forecast mode
    forecast_step <- -n_ahead

    # Find how many dual points correspond to the holdout region
    len1 <- sum(dual_data$int_x > (n + offset + forecast_step))
    len1 <- max(len1, 1L)

    # Forecast in dual domain (holdout mode)
    dual_fore <- fore_arma(dual_y, phi = phi, theta = 0,
                           n_ahead = len1, lastn = TRUE, plot = FALSE)

    # Re-interpolate
    int_start <- n + forecast_step + shift + 1
    n_keep <- length(dual_y) - len1
    combined_dual <- c(dual_y[1:n_keep], dual_fore$f)

    glam_fore <- reinterpolate(combined_dual, start = shift + 1, h = h,
                               lambda = lambda_use, int_start = int_start,
                               len = abs(forecast_step))

    # AR benchmark forecast
    ar_fore <- fore_arma(x, phi = ar_fit$phi, theta = 0,
                         n_ahead = n_ahead, lastn = TRUE, plot = FALSE)
  }

  # Plot if requested
  if (plot) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))

    par(mar = c(4, 4, 3, 1))

    if (!lastn) {
      # Future mode plot
      t_all <- seq_len(n + n_ahead)
      y_all <- c(x, rep(NA, n_ahead))

      plot(t_all, y_all, type = "l", pch = 16, cex = 0.7,
           xlab = "Time", ylab = "Value",
           main = "G-Lambda Forecast (Future Mode)",
           ylim = range(c(x, glam_fore, ar_fore$f), na.rm = TRUE))

      # G-lambda forecasts
      points((n + 1):(n + n_ahead), glam_fore, type = "l",
             lty = 2, col = "blue", lwd = 1.2)
      # AR forecasts
      points((n + 1):(n + n_ahead), ar_fore$f, type = "l",
             lty = 3, col = "red", lwd = 1.2)

      # Connect to last observation
      segments(n, x[n], n + 1, glam_fore[1], lty = 2, col = "blue")
      segments(n, x[n], n + 1, ar_fore$f[1], lty = 3, col = "red")

      legend("topleft",
             legend = c("Actual", "G-Lambda", paste0("AR(", ar_fit$p, ")")),
             lty = c(1, 2, 3), col = c("black", "blue", "red"),
             cex = 0.8, bty = "n")

    } else {
      # Holdout mode plot
      plot_start <- max(n - 5 * n_ahead, 1)
      t <- plot_start:n

      plot(t, x[plot_start:n], type = "l", pch = 16, cex = 0.7,
           xlab = "Time", ylab = "Value",
           main = "G-Lambda Forecast (Holdout Mode)",
           ylim = range(c(x[plot_start:n], glam_fore, ar_fore$f)))

      # G-lambda forecasts
      fore_t <- (n - n_ahead + 1):n
      points(fore_t, glam_fore, type = "l", lty = 2, col = "blue", lwd = 1.2)

      # AR forecasts
      points(fore_t, ar_fore$f, type = "l", lty = 3, col = "red", lwd = 1.2)

      legend("topleft",
             legend = c("Actual", "G-Lambda", paste0("AR(", ar_fit$p, ")")),
             lty = c(1, 2, 3), col = c("black", "blue", "red"),
             cex = 0.8, bty = "n")
    }
  }

  structure(
    list(
      f = glam_fore,
      ar_f = ar_fore$f,
      ar_order = ar_fit$p,
      original = x,
      n_ahead = n_ahead,
      lastn = lastn
    ),
    class = "fore_glambda"
  )
}


#' @export
print.fore_glambda <- function(x, ...) {
  cat("\nG-Lambda Forecast\n")
  cat("=================\n\n")
  cat("Mode:", if (x$lastn) "Holdout" else "Future", "\n")
  cat("Steps ahead:", x$n_ahead, "\n")
  cat("AR benchmark order:", x$ar_order, "\n\n")

  cat("G-Lambda Forecasts:\n")
  print(round(x$f, 4))

  cat("\nAR Forecasts:\n")
  print(round(x$ar_f, 4))

  invisible(x)
}
