#' Butterworth Filter for Time Series
#'
#' Apply a Butterworth filter to a time series using zero-phase filtering.
#'
#' @param x Numeric vector, the time series to filter.
#' @param order Integer, filter order (default 4).
#' @param type Character, filter type: "low", "high", "stop", or "pass".
#' @param cutoff Numeric, cutoff frequency as a fraction of Nyquist frequency (0 to 0.5).
#'   For "low" and "high" filters, a scalar value.
#'   For "stop" and "pass" filters, a 2-element vector with lower and upper cutoffs.
#' @param plot Logical, whether to plot original and filtered series (default TRUE).
#'
#' @return A list with component:
#'   \item{x_filt}{The filtered time series}
#'
#' @details
#' This function wraps the \code{signal} package's \code{butter} and \code{filtfilt}
#' functions to apply a zero-phase Butterworth filter.
#'
#' The cutoff frequency is specified as a fraction of the Nyquist frequency
#' (half the sampling rate). For example, if your data is sampled at 100 Hz
#' and you want a 10 Hz lowpass filter, use \code{cutoff = 10/50 = 0.2}.
#'
#' Filter types:
#' \describe{
#'   \item{low}{Lowpass filter - passes frequencies below cutoff}
#'   \item{high}{Highpass filter - passes frequencies above cutoff}
#'   \item{stop}{Bandstop filter - removes frequencies between cutoffs}
#'   \item{pass}{Bandpass filter - passes frequencies between cutoffs}
#' }
#'
#' @seealso \code{\link[signal]{butter}}, \code{\link[signal]{filtfilt}}
#'
#' @examples
#' \dontrun{
#' # Generate noisy signal
#' t <- seq(0, 1, length.out = 200)
#' x <- sin(2 * pi * 5 * t) + 0.5 * sin(2 * pi * 50 * t)
#'
#' # Lowpass filter to remove high frequency component
#' result <- butterworth(x, order = 4, type = "low", cutoff = 0.2)
#'
#' # Highpass filter
#' result <- butterworth(x, order = 4, type = "high", cutoff = 0.3, plot = FALSE)
#'
#' # Bandpass filter
#' result <- butterworth(x, order = 4, type = "pass", cutoff = c(0.1, 0.3))
#' }
#'
#' @export
butterworth <- function(x, order = 4L, type = c("low", "high", "stop", "pass"),
                        cutoff, plot = TRUE) {

  # Check for signal package
 if (!requireNamespace("signal", quietly = TRUE)) {
    stop("Package 'signal' is required for butterworth().\n",
         "Install it with: install.packages('signal')")
  }

  # Input validation
  if (!is.numeric(x) || length(x) == 0) {
    stop("`x` must be a non-empty numeric vector")
  }
  if (!is.numeric(order) || length(order) != 1 || order < 1) {
    stop("`order` must be a positive integer")
  }
  type <- match.arg(type)

  if (!is.numeric(cutoff)) {
    stop("`cutoff` must be numeric")
  }

  # Validate cutoff based on filter type
  if (type %in% c("low", "high")) {
    if (length(cutoff) != 1) {
      stop("For '", type, "' filter, `cutoff` must be a single value")
    }
    if (cutoff <= 0 || cutoff >= 0.5) {
      stop("`cutoff` must be between 0 and 0.5 (exclusive)")
    }
  } else {
    # stop or pass (bandstop or bandpass)
    if (length(cutoff) != 2) {
      stop("For '", type, "' filter, `cutoff` must be a 2-element vector")
    }
    if (any(cutoff <= 0) || any(cutoff >= 0.5)) {
      stop("All `cutoff` values must be between 0 and 0.5 (exclusive)")
    }
    if (cutoff[1] >= cutoff[2]) {
      stop("`cutoff[1]` must be less than `cutoff[2]`")
    }
  }

  order <- as.integer(order)
  n <- length(x)

  # Nyquist normalization: signal::butter expects cutoff as fraction of Nyquist
  # where Nyquist = 1, so we multiply by 2
  cutoff_norm <- cutoff * 2

  # Design Butterworth filter
  bf <- signal::butter(order, cutoff_norm, type = type)

  # Apply zero-phase filtering
  x_filt <- signal::filtfilt(bf, x)

  # Plot if requested
  if (plot) {
    op <- par(mfrow = c(2, 1), mar = c(3, 3, 2, 0.5))
    on.exit(par(op))

    # Determine common y-axis limits
    ylim <- range(c(x, x_filt))
    t <- seq_len(n)

    # Original series
    plot(t, x, type = "l", ylim = ylim,
         xlab = "", ylab = "", main = "Original Series")

    # Filtered series
    plot(t, x_filt, type = "l", ylim = ylim,
         xlab = "", ylab = "", main = "Filtered Series")
    mtext("Time", side = 1, line = 2, cex = 0.8)
  }

  invisible(list(x_filt = x_filt))
}
