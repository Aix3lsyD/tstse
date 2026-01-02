#' Moving Average Prediction
#'
#' Compute a simple moving average and extend forecasts beyond the series.
#'
#' @param x Numeric vector, the time series.
#' @param order Integer, the moving average window size (default 3).
#' @param n_ahead Integer, number of steps to forecast beyond the series (default 1).
#' @param plot Logical, whether to plot original and smoothed series (default TRUE).
#'
#' @return A list with components:
#'   \item{x_original}{The original time series}
#'   \item{pred}{Smoothed values plus forecasts (length n + n_ahead, first `order` values are NA)}
#'   \item{order}{The moving average order used}
#'
#' @details
#' Computes a simple (backward-looking) moving average where each smoothed value
#' at time t is the average of the previous `order` observations.
#'
#' For forecasting beyond the series (`n_ahead > 1`), predicted values are used
#' as inputs for subsequent predictions, creating a recursive forecast.
#'
#' The first `order` values of `pred` are NA since a full window is required.
#'
#' @examples
#' # Smooth a random walk
#' set.seed(123)
#' x <- cumsum(rnorm(100))
#' result <- ma_pred(x, order = 5, n_ahead = 10, plot = FALSE)
#'
#' # Forecast values are in the last n_ahead positions
#' tail(result$pred, 10)
#'
#' @export
ma_pred <- function(x, order = 3L, n_ahead = 1L, plot = TRUE) {

  # Input validation
  if (!is.numeric(x) || length(x) == 0) {
    stop("`x` must be a non-empty numeric vector")
  }
  if (!is.numeric(order) || length(order) != 1 || order < 1) {
    stop("`order` must be a positive integer")
  }

  if (!is.numeric(n_ahead) || length(n_ahead) != 1 || n_ahead < 1) {
    stop("`n_ahead` must be a positive integer")
  }

  n <- length(x)
  k <- as.integer(order)
  n_ahead <- as.integer(n_ahead)

  if (k > n) {
    stop("`order` cannot exceed the length of `x`")
  }

  # Extended length for forecasts
  nka <- n + n_ahead

  # Working copy of x extended with NAs for forecasts
  xx <- c(x, rep(NA_real_, n_ahead))

  # Initialize smoothed values vector
  x_sm <- rep(NA_real_, nka)

  # First smoothed value at position k+1 (uses first k values)
  t_start <- k + 1L
  x_sm[t_start] <- mean(xx[1:k])

  # Rolling update for remaining in-sample values
  # x.sm[i] = x.sm[i-1] - xx[i-1-k]/k + xx[i-1]/k
  if (t_start < n + 1L) {
    for (i in (t_start + 1L):(n + 1L)) {
      x_sm[i] <- x_sm[i - 1L] - xx[i - 1L - k] / k + xx[i - 1L] / k
    }
  }

  # Forecast beyond the series
  if (n_ahead >= 1L) {
    # First forecast uses actual + smoothed
    xx[n + 1L] <- x_sm[n + 1L]

    # Additional forecasts use predicted values
    if (n_ahead > 1L) {
      for (i in (n + 2L):nka) {
        x_sm[i] <- mean(xx[(i - k):(i - 1L)])
        xx[i] <- x_sm[i]
      }
    }
  }

  # Plot if requested
  if (plot) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))

    plot(seq_len(n), x, type = "o", xlim = c(1, nka),
         xlab = "Time", ylab = "Value",
         main = paste0("MA(", order, ") Smoothing with ", n_ahead, "-step Forecast"))
    points(seq_len(nka), x_sm, pch = 1)
    lines(seq_len(nka), x_sm, lwd = 2)

    # Mark forecast region
    if (n_ahead > 0) {
      abline(v = n + 0.5, lty = 2, col = "gray")
    }
  }

  invisible(list(
    x_original = x,
    pred = x_sm,
    order = order
  ))
}
