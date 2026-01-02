#' Hilbert Transform
#'
#' Computes the Hilbert transform of a real-valued signal using the FFT-based
#' method, returning the analytic signal.
#'
#' @param x Numeric vector. The real-valued input signal. For best FFT
#'   performance, length should be a power of 2, but any length is supported.
#' @param pad_odd Logical. If `TRUE` (default), odd-length signals are
#'
#'   zero-padded to even length. If `FALSE`, the last element is dropped
#'   (matching original `tswge` behavior, not recommended).
#'
#' @return Complex vector of the same length as `x` (or length + 1 if `x` has
#'   odd length and `pad_odd = TRUE`). This is the analytic signal where:
#'   \itemize{
#'     \item Real part: the original signal (possibly padded/truncated)
#'     \item Imaginary part: the Hilbert transform of the signal
#'   }
#'
#' @details
#' The Hilbert transform is computed using the FFT-based method:
#' \enumerate{
#'   \item Compute the FFT of the input signal
#'   \item Zero out negative frequency components
#'   \item Double the positive frequency components
#'   \item Compute the inverse FFT
#' }
#'
#' The result is the analytic signal \eqn{z(t) = x(t) + i \cdot H\{x(t)\}},
#' where \eqn{H\{x(t)\}} is the Hilbert transform.
#'
#' @section Derived Quantities:
#' From the analytic signal, you can compute:
#' \itemize{
#'   \item **Instantaneous amplitude (envelope)**: `Mod(z)` or `abs(z)`
#'   \item **Instantaneous phase**: `Arg(z)` or `atan2(Im(z), Re(z))`
#'   \item **Instantaneous frequency**: `diff(unwrap(Arg(z))) / (2 * pi * dt)`
#'   \item **Hilbert transform only**: `Im(z)`
#' }
#'
#' @references
#' Marple, S. L. (1999). Computing the discrete-time "analytic" signal via FFT.
#' *IEEE Transactions on Signal Processing*, 47(9), 2600-2603.
#'
#' @seealso [fft()] for the underlying FFT computation.
#'
#' @examples
#' # Simple sinusoid
#' t <- seq(0, 1, length.out = 256)
#' x <- sin(2 * pi * 5 * t)  # 5 Hz sine wave
#' z <- hilbert(x)
#'
#' # The envelope should be approximately 1
#' envelope <- Mod(z)
#' plot(t, envelope, type = "l", ylim = c(0, 1.5))
#'
#' # Amplitude-modulated signal
#' carrier <- sin(2 * pi * 50 * t)
#' modulator <- 0.5 + 0.5 * sin(2 * pi * 2 * t)  # 2 Hz modulation
#' am_signal <- carrier * modulator
#' z_am <- hilbert(am_signal)
#'
#' plot(t, am_signal, type = "l", col = "gray")
#' lines(t, Mod(z_am), col = "red", lwd = 2)
#' lines(t, modulator, col = "blue", lty = 2)  # True envelope
#'
#' @export
hilbert <- function(x, pad_odd = TRUE) {
  # Input validation
  if (!is.numeric(x) && !is.complex(x)) {
    stop("`x` must be a numeric or complex vector.", call. = FALSE)
  }
  if (length(x) == 0L) {
    stop("`x` must have length >= 1.", call. = FALSE)
  }
  if (anyNA(x)) {
    stop("`x` contains NA values.", call. = FALSE)
  }
  if (!is.logical(pad_odd) || length(pad_odd) != 1L || is.na(pad_odd)) {
    stop("`pad_odd` must be TRUE or FALSE.", call. = FALSE)
  }

  n <- length(x)

  # Handle odd-length signals
  if (n %% 2L == 1L) {
    if (pad_odd) {
      # Zero-pad to even length (preferred)
      x <- c(x, 0)
      n <- n + 1L
    } else {
      # Drop first element (legacy tswge behavior)
      warning("Odd-length input: dropping first element. Use `pad_odd = TRUE` to zero-pad instead.",
              call. = FALSE)
      x <- x[-1L]
      n <- n - 1L
    }
  }

  # Handle edge case of length 0 after truncation (original length was 1)
  if (n == 0L) {
    return(complex(0))
  }

  # Compute FFT
  X <- fft(x)

  # Create frequency-domain filter for analytic signal
  # h = [1, 2, 2, ..., 2, 1, 0, 0, ..., 0]
  # Index 1: DC component (keep as-is)
  # Index 2 to n/2: positive frequencies (double)
  # Index n/2 + 1: Nyquist frequency (keep as-is)
  # Index n/2 + 2 to n: negative frequencies (zero out)
  h <- numeric(n)
  h[1L] <- 1                          # DC
  h[n %/% 2L + 1L] <- 1               # Nyquist
  if (n > 2L) {
    h[2L:(n %/% 2L)] <- 2             # Positive frequencies
  }
  # Negative frequencies remain 0

  # Apply filter and inverse FFT
  z <- fft(X * h, inverse = TRUE) / n

  z
}
