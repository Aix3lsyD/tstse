#' Wigner-Ville Distribution
#'
#' Computes and optionally plots the Wigner-Ville distribution, a time-frequency
#' representation of a signal useful for analyzing non-stationary time series.
#'
#' @param x Numeric vector. The time series to analyze.
#' @param n_freq Integer. Number of frequency points (0 to 0.5). Default is 101.
#' @param plot Logical. If `TRUE` (default), plot the time-frequency image.
#'
#' @return A list with class `"wv"` containing:
#'   \item{tfr}{Matrix of the time-frequency representation. Rows are frequencies,
#'     columns are time points.}
#'   \item{freq}{Numeric vector of frequencies (0 to 0.5).}
#'   \item{time}{Numeric vector of time indices.}
#'
#' @details
#' The Wigner-Ville distribution is a quadratic time-frequency representation
#' defined as:
#' \deqn{W_x(t, f) = \sum_{\tau} x(t + \tau) x^*(t - \tau) e^{-j 2\pi f \tau}}
#'
#' where \eqn{x^*} denotes the complex conjugate. The function first computes
#' the analytic signal using the Hilbert transform, then computes the
#' Wigner-Ville distribution.
#'
#' The Wigner-Ville distribution provides high resolution in both time and
#' frequency simultaneously, but can exhibit cross-term interference for
#' multi-component signals.
#'
#' @section Interpretation:
#' \itemize{
#'   \item Bright regions indicate high energy at that time-frequency location
#'   \item Horizontal bands suggest stationary frequency components
#'   \item Diagonal or curved patterns indicate frequency modulation (chirps)
#'   \item Cross-terms (interference) may appear between signal components
#' }
#'
#' @seealso [hilbert()] for the Hilbert transform used internally,
#'   [sample_spec()] for frequency-only spectral analysis.
#'
#' @examples
#' # Chirp signal (frequency increases over time)
#' t <- 1:256
#' x <- cos(2 * pi * (0.05 * t + 0.0002 * t^2))
#' result <- wv(x)
#'
#' # Sum of two sinusoids
#' x2 <- cos(2 * pi * 0.1 * t) + cos(2 * pi * 0.3 * t)
#' result2 <- wv(x2)
#'
#' # Just compute, no plot
#' result3 <- wv(x, plot = FALSE)
#' dim(result3$tfr)
#'
#' @export
wv <- function(x, n_freq = 101L, plot = TRUE) {
  # Input validation
  if (!is.numeric(x)) {
    stop("`x` must be a numeric vector.", call. = FALSE)
  }
  n <- length(x)
  if (n < 4L) {
    stop("`x` must have at least 4 observations.", call. = FALSE)
  }
  if (anyNA(x)) {
    stop("`x` contains NA values.", call. = FALSE)
  }

  if (!is.numeric(n_freq) || length(n_freq) != 1L || n_freq < 2L) {
    stop("`n_freq` must be an integer >= 2.", call. = FALSE)
  }
  n_freq <- as.integer(n_freq)

  if (!is.logical(plot) || length(plot) != 1L || is.na(plot)) {
    stop("`plot` must be TRUE or FALSE.", call. = FALSE)
  }

  # Frequency grid
  freq <- seq(0, 0.5, length.out = n_freq)

  # Compute analytic signal using Hilbert transform
  z <- hilbert(x)

  # Initialize time-frequency representation matrix
  # Rows = frequencies, Columns = time points
  tfr <- matrix(0, n_freq, n)
  time_idx <- seq_len(n)

  # Compute Wigner-Ville distribution
  for (i in seq_len(n)) {
    ti <- time_idx[i]

    # Maximum lag based on position in signal
    taumax <- min(ti - 1L, n - ti, floor(n_freq / 2) - 1L)

    if (taumax >= 0) {
      tau <- (-taumax):taumax

      # Circular indices for frequency dimension
      indices <- ((n_freq + tau) %% n_freq) + 1L

      # Instantaneous autocorrelation
      tfr[indices, i] <- z[ti + tau] * Conj(z[ti - tau])
    }

    # Handle the Nyquist frequency point
    tau_nyq <- floor(n_freq / 2) + 1L
    if ((ti <= n - tau_nyq) && (ti >= tau_nyq + 1L)) {
      tfr[tau_nyq + 1L, i] <- z[ti + tau_nyq] * Conj(z[ti - tau_nyq])
    }
  }

  # FFT along frequency dimension to get spectrum at each time
  tfr <- apply(tfr, 2, fft)

  # Take real part (imaginary should be ~0 for real input)
  tfr <- Re(tfr)

  # Plot if requested
  if (plot) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))

    par(mar = c(4, 4, 2, 1))

    # Grayscale colormap (dark = high, light = low)
    gcol <- seq(1, 0, length.out = 11)

    image(time_idx, freq, t(tfr),
          col = gray(gcol),
          xlab = "Time", ylab = "Frequency",
          main = "Wigner-Ville Distribution",
          xlim = c(0, n))
  }

  # Build result
  result <- list(
    tfr = tfr,
    freq = freq,
    time = time_idx
  )

  class(result) <- "wv"
  invisible(result)
}


#' @export
print.wv <- function(x, ...) {
  cat("Wigner-Ville Distribution\n")
  cat("-------------------------\n")
  cat("Time points:", length(x$time), "\n")
  cat("Frequency points:", length(x$freq), "\n")
  cat("TFR matrix:", nrow(x$tfr), "x", ncol(x$tfr), "\n")
  cat("\nUse $tfr to access the time-frequency matrix.\n")
  invisible(x)
}
