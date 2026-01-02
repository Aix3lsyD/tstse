#' Sample Spectral Density
#'
#' Computes and optionally plots the sample spectral density (periodogram)
#' of a time series using the autocorrelation function.
#'
#' @param x Numeric vector. The time series to analyze.
#' @param db Logical. If `TRUE` (default), return the spectral density in
#'   decibels (dB), i.e., \eqn{10 \log_{10}(S(f))}. If `FALSE`, return the
#'   raw spectral density.
#' @param n_freq Integer. Number of frequency points to evaluate. Default is
#'   251, giving frequencies from 0 to 0.5 in increments of 0.002.
#' @param plot Logical. If `TRUE` (default), plot the spectral density.
#'
#' @return A list with class `"sample_spec"` containing:
#'   \item{freq}{Numeric vector of frequencies (0 to 0.5).}
#'   \item{spec}{Numeric vector of spectral density estimates at each frequency.}
#'   \item{db}{Logical indicating whether `spec` is in decibels.}
#'
#' @details
#' The sample spectral density is computed from the sample autocorrelation
#' function using the formula:
#' \deqn{S(f) = 1 + 2 \sum_{k=1}^{n-1} \hat{\rho}_k \cos(2\pi f k)}
#'
#' where \eqn{\hat{\rho}_k} is the sample autocorrelation at lag \eqn{k}.
#'
#' This is the raw (unsmoothed) periodogram. For a smoothed estimate using
#' the Parzen window, see [parzen()].
#'
#' @section Interpretation:
#' \itemize{
#'   \item Frequencies range from 0 (long-term trend) to 0.5 (Nyquist frequency)
#'   \item Peaks indicate dominant periodicities in the data
#'   \item The dB scale helps visualize differences across a wide range of values
#' }
#'
#' @seealso [parzen()] for smoothed spectral estimation,
#'   [true_spec()] for theoretical ARMA spectra,
#'   [stats::spectrum()] for R's built-in spectral analysis.
#'
#' @examples
#' # White noise has flat spectrum
#' set.seed(123)
#' wn <- rnorm(200)
#' sample_spec(wn)
#'
#' # AR(1) with positive phi has more power at low frequencies
#' ar1 <- arima.sim(model = list(ar = 0.8), n = 200)
#' sample_spec(ar1)
#'
#' # Seasonal pattern shows peak at seasonal frequency
#' seasonal <- sin(2 * pi * (1:200) / 12) + rnorm(200, sd = 0.5)
#' sample_spec(seasonal)
#'
#' # Non-dB scale
#' sample_spec(ar1, db = FALSE)
#'
#' @export
sample_spec <- function(x, db = TRUE, n_freq = 251L, plot = TRUE) {
  # Input validation
  if (!is.numeric(x)) {
    stop("`x` must be a numeric vector.", call. = FALSE)
  }
  n <- length(x)
  if (n < 2L) {
    stop("`x` must have at least 2 observations.", call. = FALSE)
  }
  if (anyNA(x)) {
    stop("`x` contains NA values.", call. = FALSE)
  }

  if (!is.logical(db) || length(db) != 1L || is.na(db)) {
    stop("`db` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.numeric(n_freq) || length(n_freq) != 1L || is.na(n_freq)) {
    stop("`n_freq` must be a single positive integer.", call. = FALSE)
  }
  n_freq <- as.integer(n_freq)
  if (n_freq < 2L) {
    stop("`n_freq` must be at least 2.", call. = FALSE)
  }

  if (!is.logical(plot) || length(plot) != 1L || is.na(plot)) {
    stop("`plot` must be TRUE or FALSE.", call. = FALSE)
  }

  # Compute sample autocorrelations (excluding lag 0)
  # acf() returns an array with lag 0 at index 1
  acf_result <- acf(x, lag.max = n - 1L, plot = FALSE)
  rho <- as.numeric(acf_result$acf)  # rho[1] is lag 0, rho[2] is lag 1, etc.

  # Frequency grid from 0 to 0.5
  freq <- seq(0, 0.5, length.out = n_freq)

  # Compute sample spectral density
  # S(f) = 1 + 2 * sum_{k=1}^{n-1} rho_k * cos(2*pi*f*k)
  # Using vectorized computation
  lags <- seq_len(n - 1L)  # 1 to n-1
  rho_lags <- rho[-1]      # rho at lags 1 to n-1 (exclude lag 0)

  # Compute for all frequencies at once using outer product
  # cos_matrix[i, k] = cos(2*pi*freq[i]*k)
  cos_matrix <- cos(2 * pi * outer(freq, lags))

  # spec[i] = 1 + 2 * sum_k(rho_k * cos(2*pi*freq[i]*k))
  spec <- 1 + 2 * (cos_matrix %*% rho_lags)
  spec <- as.numeric(spec)

  # Match original behavior: at frequency 0, set to exactly 1
  # (original initializes sample.sp=rep(1,251) and loop starts at i=2)
  spec[1] <- 1

  # Convert to dB if requested
  if (db) {
    # Handle potential negative values from estimation variance
    spec <- pmax(spec, .Machine$double.eps)
    spec <- 10 * log10(spec)
  }

  # Plot if requested
  if (plot) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))

    par(mar = c(4, 4, 2, 1))

    ylab <- if (db) "Spectral Density (dB)" else "Spectral Density"

    plot(freq, spec, type = "l", lwd = 1,
         xlab = "Frequency", ylab = ylab,
         main = "Sample Spectral Density")
  }

  # Build result
  result <- list(
    freq = freq,
    spec = spec,
    db = db
  )

  class(result) <- "sample_spec"
  invisible(result)
}


#' @export
print.sample_spec <- function(x, ...) {
  cat("Sample Spectral Density\n")
  cat("-----------------------\n")
  cat("Frequencies:", length(x$freq), "points from 0 to 0.5\n")
  cat("Scale:", if (x$db) "dB (decibels)" else "linear", "\n")

  # Find dominant frequency (peak)
  peak_idx <- which.max(x$spec)
  if (peak_idx > 1) {  # Ignore DC component
    cat("Peak frequency:", round(x$freq[peak_idx], 4),
        "(period:", round(1 / x$freq[peak_idx], 2), ")\n")
  }

  invisible(x)
}
