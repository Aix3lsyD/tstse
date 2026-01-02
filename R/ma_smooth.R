#' Centered Moving Average Smoothing
#'
#' Apply a centered moving average smoother to a time series.
#'
#' @param x Numeric vector, the time series.
#' @param order Integer, the moving average window size (default 3).
#' @param plot Logical, whether to plot original and smoothed series (default TRUE).
#'
#' @return A list with components:
#'   \item{x}{The original time series}
#'   \item{smooth}{Smoothed values (same length as x, with NAs at ends)}
#'   \item{order}{The moving average order used}
#'
#' @details
#' Computes a centered moving average where each smoothed value at time t
#' is the average of observations centered around t.
#'
#' For odd orders, a simple centered average is used.
#'
#' For even orders, a weighted average is used where the endpoints receive
#' half weight. This is equivalent to a 2x`order` centered MA, commonly used
#' for seasonal adjustment (e.g., 2x12 MA for monthly data).
#'
#' The first and last `floor(order/2)` values are NA since a full centered
#' window is not available.
#'
#' @examples
#' # Smooth a noisy sine wave
#' t <- seq(0, 4 * pi, length.out = 200)
#' x <- sin(t) + rnorm(200, sd = 0.3)
#' result <- ma_smooth(x, order = 11, plot = FALSE)
#'
#' # Compare original and smoothed
#' plot(x, type = "l", col = "gray")
#' lines(result$smooth, col = "blue", lwd = 2)
#'
#' @export
ma_smooth <- function(x, order = 3L, plot = TRUE) {

 # Input validation
  if (!is.numeric(x) || length(x) == 0) {
    stop("`x` must be a non-empty numeric vector")
  }
  if (!is.numeric(order) || length(order) != 1 || order < 1) {
    stop("`order` must be a positive integer")
  }

  n <- length(x)
  k <- as.integer(order)

  if (k > n) {
    stop("`order` cannot exceed the length of `x`")
  }

  # Initialize smoothed values
  x_sm <- rep(NA_real_, n)

  k2 <- k %/% 2L  # floor(k/2)
  is_odd <- (k %% 2L) == 1L

  if (is_odd) {
    # Odd order: simple centered average
    # Window: [t - k2, t + k2] has k elements
    k2p1 <- k2 + 1L
    t_start <- k - k2  # = k2 + 1 for odd k

    # First smoothed value
    x_sm[t_start] <- mean(x[1:k])

    # Rolling update for remaining values
    ts1 <- t_start + 1L
    tsn <- n - k2

    if (ts1 <= tsn) {
      for (i in ts1:tsn) {
        # Remove oldest, add newest
        x_sm[i] <- x_sm[i - 1L] - x[i - k2p1] / k + x[i + k2] / k
      }
    }

  } else {
    # Even order: weighted endpoints (half weight at ends)
    # This is the 2xk centered MA
    k2p1 <- k2 + 1L
    k2m1 <- k2 - 1L
    kp1 <- k + 1L
    t_start <- k - k2 + 1L  # = k2 + 1 for even k

    # First smoothed value: endpoints get half weight
    x_sm[t_start] <- x[1] / (2 * k) + x[kp1] / (2 * k)
    for (i in 2:k) {
      x_sm[t_start] <- x_sm[t_start] + x[i] / k
    }

    # Rolling update
    ts1 <- t_start + 1L
    tsn <- n - k2

    if (ts1 <= tsn) {
      for (i in ts1:tsn) {
        # Weighted update for even order
        x_sm[i] <- x_sm[i - 1L] -
          x[i - k2p1] / (2 * k) - x[i - k2] / (2 * k) +
          x[i + k2m1] / (2 * k) + x[i + k2] / (2 * k)
      }
    }
  }

  # Determine valid range for plotting
  tsn <- n - k2
  t_start <- if (is_odd) k - k2 else k - k2 + 1L

  # Plot if requested
  if (plot) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))

    plot(seq_len(n), x, type = "l", col = "gray50",
         xlab = "Time", ylab = "Value",
         main = paste0("Centered MA(", order, ") Smoothing"))
    lines(t_start:tsn, x_sm[t_start:tsn], lwd = 2, col = "black")
    legend("topright", legend = c("Original", "Smoothed"),
           col = c("gray50", "black"), lwd = c(1, 2), bty = "n")
  }

  invisible(list(
    x = x,
    smooth = x_sm,
    order = order
  ))
}
