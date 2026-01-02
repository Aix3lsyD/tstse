#' Parzen Window Spectral Estimate
#'
#' Calculate and plot the smoothed periodogram using the Parzen window.
#'
#' @param x Numeric vector, the time series.
#' @param db Logical, if TRUE (default) return spectrum in decibels (10*log10).
#' @param trunc Integer, truncation point M. Default 0 uses M = 2*sqrt(n).
#' @param plot Logical, if TRUE (default) plot the spectrum.
#'
#' @return A list with components:
#'   \item{freq}{Frequency values (0 to 0.5)}
#'   \item{pzgram}{Parzen-smoothed spectral estimates}
#'   \item{M}{Truncation point used}
#' @export
parzen <- function(x, db = TRUE, trunc = 0L, plot = TRUE) {

  # Input validation
  if (!is.numeric(x) || length(x) == 0) {
    stop("`x` must be a non-empty numeric vector")
  }
  if (!is.numeric(trunc) || trunc < 0) {
    stop("`trunc` must be a non-negative integer")
  }

  n <- length(x)

  # Compute sample autocorrelations (lags 0 to n-1)
  aut <- acf(x, lag.max = n - 1, plot = FALSE)$acf[, 1, 1]

  # Frequency grid
  freq <- (1:floor(n / 2)) / n
  nf <- length(freq)

  # Truncation point
  M <- if (trunc == 0) floor(2 * sqrt(n)) else as.integer(trunc)

  # Parzen window weights (length n, indices correspond to lags 0 to n-1)
  weight <- numeric(n)
  k <- 0:(n - 1)

  # First part: 0 <= k <= M/2
  idx1 <- which(k <= floor(M / 2))
  u1 <- k[idx1] / M
  weight[idx1] <- 1 - 6 * u1^2 + 6 * u1^3

  # Second part: M/2 < k <= M
  idx2 <- which(k > floor(M / 2) & k <= M)
  u2 <- k[idx2] / M
  weight[idx2] <- 2 * (1 - u2)^3

  # k > M: weight stays 0

  # Compute Parzen-smoothed periodogram (vectorized)
  # Formula: S(f) = gamma(0)*w(0) + 2 * sum_{k=1}^{n-1} gamma(k)*w(k)*cos(2*pi*f*k)
  lags <- 1:(n - 1)

  # Precompute weighted autocorrelations (excluding lag 0)
  weighted_aut <- aut[-1] * weight[-1]

  # Create cos_matrix[i, k] = cos(2*pi*freq[i]*k) for all frequencies and lags
  cos_matrix <- cos(2 * pi * outer(freq, lags))

  # Compute all frequencies at once via matrix multiplication
  pzgram <- aut[1] * weight[1] + 2 * (cos_matrix %*% weighted_aut)
  pzgram <- as.numeric(pzgram)

  # Convert to dB if requested
  if (db) {
    pzgram <- 10 * log10(pzgram)
  }

  # Plot
  if (plot) {
    ylab <- if (db) "dB" else "Spectral Estimate"
    plot(freq, pzgram, type = "l", xlab = "Frequency", ylab = ylab,
         main = paste("Parzen Window (M =", M, ")"))
    invisible(list(freq = freq, pzgram = pzgram, M = M))
  } else {
    list(freq = freq, pzgram = pzgram, M = M)
  }
}
