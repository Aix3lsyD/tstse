#' Periodogram
#'
#' Compute and plot the periodogram of a time series.
#'
#' @param x Numeric vector, the time series.
#' @param db Logical, if TRUE (default) return in decibels (10*log10).
#' @param plot Logical, if TRUE (default) plot the periodogram.
#'
#' @return A list with components:
#'   \item{freq}{Frequency values (0 to 0.5)}
#'   \item{pgram}{Periodogram values}
#' @export
#'
#' @details
#' The periodogram is computed as:
#' \deqn{I(f) = \hat{\gamma}(0) + 2 \sum_{k=1}^{n-1} \hat{\gamma}(k) \cos(2\pi f k)}
#' where \eqn{\hat{\gamma}(k)} is the sample autocovariance at lag k.
#'
#' @seealso \code{\link{parzen}} for a smoothed spectral estimate.
#'
#' @examples
#' x <- gen_arma(n = 200, phi = c(1.5, -0.75), plot = FALSE, seed = 123)
#'
#' # Default (dB scale)
#' period(x)
#'
#' # Raw scale
#' period(x, db = FALSE)
period <- function(x, db = TRUE, plot = TRUE) {

  # Input validation
  if (!is.numeric(x) || length(x) == 0) {
    stop("`x` must be a non-empty numeric vector")
  }

  n <- length(x)

  # Compute sample autocorrelations
  aut <- acf(x, lag.max = n - 1, plot = FALSE)$acf[, 1, 1]

  # Frequency grid
  freq <- (1:floor(n / 2)) / n
  nf <- length(freq)

  # Compute periodogram
  lags <- 1:(n - 1)
  pgram <- numeric(nf)

  for (i in seq_len(nf)) {
    cosvec <- cos(2 * pi * freq[i] * lags)
    pgram[i] <- aut[1] + 2 * sum(aut[-1] * cosvec)
  }

  # Convert to dB if requested
  if (db) {
    pgram <- 10 * log10(pgram)
  }

  # Plot
  if (plot) {
    ylab <- if (db) "dB" else "Periodogram"
    .plot_periodogram(freq, pgram, ylab, n)
    invisible(list(freq = freq, pgram = pgram))
  } else {
    list(freq = freq, pgram = pgram)
  }
}


#' Plot periodogram as vertical segments
#' @noRd
.plot_periodogram <- function(freq, pgram, ylab, n) {

  min_pgram <- min(pgram, na.rm = TRUE)

  plot(freq, pgram, type = "n",
       xlab = "Frequency", ylab = ylab, main = "Periodogram")

  segments(freq, min_pgram, freq, pgram)
}
