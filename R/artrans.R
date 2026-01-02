#' AR Transformation
#'
#' Apply an AR transformation Y(t) = phi(B) X(t) to a time series.
#'
#' @param x Numeric vector, the original time series.
#' @param phi Numeric vector of AR coefficients (without the leading 1).
#' @param lag_max Integer, maximum lag for ACF (default 25).
#' @param plot Logical, whether to plot diagnostics (default TRUE).
#'
#' @return Numeric vector, the transformed series (invisibly if plot = TRUE).
#' @export
#'
#' @examples
#' # First difference: Y(t) = (1 - B)X(t)
#' artrans(sunspot.year, phi = 1)
#'
#' # Second order: Y(t) = (1 - 1.6B + B^2)X(t)
#' artrans(sunspot.year, phi = c(1.6, -1))
artrans <- function(x, phi, lag_max = 25L, plot = TRUE) {


  # Input validation
  if (!is.numeric(x) || length(x) == 0) {
    stop("`x` must be a non-empty numeric vector")
  }
  if (!is.numeric(phi) || length(phi) == 0) {
    stop("`phi` must be a non-empty numeric vector")
  }
  if (!is.numeric(lag_max) || lag_max < 1) {
    stop("`lag_max` must be a positive integer")
  }

  x <- as.double(x)
  phi <- as.double(phi)
  p <- length(phi)
  n <- length(x)

  # Adjust lag_max if needed

  max_possible <- n - p - 1L
  if (lag_max > max_possible) lag_max <- max_possible

  # AR transformation via convolution filter (C-level speed)
  coeffs <- c(1, -phi)
  y_full <- stats::filter(x, coeffs, sides = 1, method = "convolution")
  y <- as.double(y_full[-(1:p)])

  # Compute ACFs
  acf_x <- acf(x, lag.max = lag_max, plot = FALSE)$acf[, 1, 1]
  acf_y <- acf(y, lag.max = lag_max, plot = FALSE)$acf[, 1, 1]

  if (plot) {
    .plot_artrans(x, y, acf_x, acf_y, lag_max)
    invisible(y)
  } else {
    y
  }
}


# Internal plotting helper
.plot_artrans <- function(x, y, acf_x, acf_y, lag_max) {
  op <- par(mfrow = c(2, 2), mar = c(4, 3, 2, 1))
  on.exit(par(op))

  lags <- 0:lag_max

  # Original series
  plot(seq_along(x), x, type = "o", pch = 16, cex = 0.5,
       xlab = "Time", ylab = "", main = "Original")

  # ACF of original
  plot(lags, acf_x, type = "h", ylim = c(-1, 1),
       xlab = "Lag", ylab = "", main = "ACF Original")
  abline(h = 0)

  # Transformed series
  plot(seq_along(y), y, type = "o", pch = 16, cex = 0.5,
       xlab = "Time", ylab = "", main = "Transformed")

  # ACF of transformed
  plot(lags, acf_y, type = "h", ylim = c(-1, 1),
       xlab = "Lag", ylab = "", main = "ACF Transformed")
  abline(h = 0)
}
