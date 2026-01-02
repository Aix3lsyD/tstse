#' Plot Time Series Sample Analysis
#'
#' Plot a time series realization with sample autocorrelations, Parzen window
#' spectral estimate, and optionally the periodogram.
#'
#' @param x Numeric vector, the time series.
#' @param lag_max Integer, maximum lag for ACF (default 25).
#' @param trunc Integer, Parzen window truncation point. Default 0 uses 2*sqrt(n).
#' @param arlimits Logical, whether to show 95% ACF limits (default FALSE).
#' @param speclimits Numeric vector of length 2, y-axis limits for spectral plots.
#'   Default c(0, 0) uses data range.
#' @param periodogram Logical, whether to include periodogram plot (default FALSE).
#'
#' @return Invisibly returns a list with components:
#'   \item{xbar}{Sample mean}
#'   \item{acf}{Sample autocorrelations}
#'   \item{freq}{Frequencies for spectral estimates}
#'   \item{parzen_db}{Parzen spectral estimate in dB}
#'   \item{periodogram_db}{Periodogram in dB (if periodogram = TRUE)}
#' @export
#'
#' @details
#' Creates a multi-panel plot showing:
#' \itemize{
#'   \item Time series realization
#'   \item Sample autocorrelations
#'   \item Parzen window spectral estimate (dB)
#'   \item Periodogram (dB) - if \code{periodogram = TRUE}
#' }
#'
#' @examples
#' x <- gen_arma(n = 200, phi = c(1.5, -0.75), plot = FALSE, seed = 123)
#'
#' # Basic 3-panel plot
#' plotts_sample(x)
#'
#' # With periodogram (4-panel)
#' plotts_sample(x, periodogram = TRUE)
#'
#' # With ACF confidence limits
#' plotts_sample(x, arlimits = TRUE)
plotts_sample <- function(x,
                          lag_max = 25L,
                          trunc = 0L,
                          arlimits = FALSE,
                          speclimits = c(0, 0),
                          periodogram = FALSE) {

  # Input validation
  if (!is.numeric(x) || length(x) == 0) {
    stop("`x` must be a non-empty numeric vector")
  }
  if (!is.numeric(lag_max) || lag_max < 1) {
    stop("`lag_max` must be a positive integer")
  }

  n <- length(x)
  if (lag_max > n - 1) lag_max <- n - 1

  xbar <- mean(x)

  # Compute ACF
  acf_obj <- acf(x, lag.max = lag_max, plot = FALSE)
  acf_values <- acf_obj$acf[, 1, 1]

  # Compute Parzen spectral estimate
  pz <- parzen(x, db = TRUE, trunc = trunc, plot = FALSE)
  parzen_db <- pz$pzgram
  freq <- pz$freq

  # Compute periodogram if requested
  if (periodogram) {
    prd <- period(x, db = TRUE, plot = FALSE)
    periodogram_db <- prd$pgram
  }

  # Set up spectral limits
  if (sum(speclimits^2) == 0) {
    if (periodogram) {
      speclimits <- range(c(parzen_db, periodogram_db), finite = TRUE)
    } else {
      speclimits <- range(parzen_db, finite = TRUE)
    }
  }

  # Set up plot layout
  if (periodogram) {
    op <- par(mfrow = c(2, 2), mar = c(4.3, 3.5, 2, 1))
  } else {
    op <- par(mfrow = c(1, 3), mar = c(4.3, 3.5, 2, 1))
  }
  on.exit(par(op))

  # Plot 1: Time series
  .plot_realization_panel(x, n)

  # Plot 2: Sample ACF
  .plot_acf_panel(acf_values, lag_max, n, arlimits)

  # Plot 3: Parzen window
  .plot_spectrum_panel(freq, parzen_db, speclimits, "Parzen Window")

  # Plot 4: Periodogram (if requested)
  if (periodogram) {
    .plot_periodogram_panel(freq, periodogram_db, speclimits, n)
  }

  # Return results
  result <- list(
    xbar = xbar,
    acf = acf_values,
    freq = freq,
    parzen_db = parzen_db
  )

  if (periodogram) {
    result$periodogram_db <- periodogram_db
  }

  invisible(result)
}


#' Plot realization panel
#' @noRd
.plot_realization_panel <- function(x, n) {

  t <- seq_len(n)
  type <- if (n < 200) "o" else "l"
  pch <- if (n < 200) 16 else NA

  plot(t, x, type = type, pch = pch, cex = 0.5,
       xlab = "Time", ylab = "", main = "Realization")
}


#' Plot ACF panel
#' @noRd
.plot_acf_panel <- function(acf_values, lag_max, n, limits) {

  k <- 0:lag_max

  plot(k, acf_values, type = "h", ylim = c(-1, 1),
       xlab = "Lag", ylab = "", main = "Sample Autocorrelations")
  abline(h = 0)

  if (limits) {
    ul <- 2 / sqrt(n)
    abline(h = c(-ul, ul), lty = 2, col = "blue")
  }
}


#' Plot spectrum panel
#' @noRd
.plot_spectrum_panel <- function(freq, db, ylim, title) {

  plot(freq, db, type = "l",
       xlab = "Frequency", ylab = "dB", ylim = ylim, main = title)
}


#' Plot periodogram panel
#' @noRd
.plot_periodogram_panel <- function(freq, db, ylim, n) {

  min_db <- min(db, na.rm = TRUE)

  if (n <= 200) {
    # Draw as vertical segments
    plot(freq, db, type = "n",
         xlab = "Frequency", ylab = "dB", ylim = ylim, main = "Periodogram")
    segments(freq, min_db, freq, db)
  } else {
    plot(freq, db, type = "l",
         xlab = "Frequency", ylab = "dB", ylim = ylim, main = "Periodogram")
  }
}
