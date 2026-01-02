#' Plot Periodogram and Parzen-Smoothed Spectra
#'
#' Creates a multi-panel plot showing the periodogram and up to three
#' Parzen-windowed spectral density estimates with different truncation points.
#'
#' @param x Numeric vector. The time series to analyze.
#' @param trunc Additional truncation points for Parzen smoothing. A numeric
#'   vector of length 1 or 2 specifying additional truncation points beyond
#'   the default. Use `c(0, 0)` (default) for no additional truncation points,
#'   or e.g. `c(10, 24)` to add smoothed estimates with M=10 and M=24.
#' @param plot Logical. If `TRUE` (default), create the multi-panel plot.
#'   If `FALSE`, only compute and return the spectral estimates.
#'
#' @return A list with class `"plotts_parzen"` containing:
#'   \item{freq}{Numeric vector of Fourier frequencies.}
#'   \item{periodogram}{Periodogram values in dB at each frequency.}
#'   \item{parzen}{Parzen-smoothed spectrum in dB with default truncation
#'     \eqn{M = 2\sqrt{n}}.}
#'   \item{parzen1}{Parzen-smoothed spectrum in dB with truncation `trunc[1]`
#'     (NULL if not requested).}
#'   \item{parzen2}{Parzen-smoothed spectrum in dB with truncation `trunc[2]`
#'     (NULL if not requested).}
#'   \item{trunc_points}{Vector of truncation points actually used.}
#'
#' @details
#' The function computes:
#' \enumerate{
#'   \item The raw periodogram (unsmoothed)
#'   \item Parzen-smoothed spectrum with default \eqn{M = 2\sqrt{n}}
#'   \item Optionally, 1-2 additional Parzen-smoothed spectra with user-specified
#'     truncation points
#' }
#'
#' The Parzen window applies weights that taper smoothly to zero, reducing
#' spectral leakage compared to the raw periodogram. Smaller truncation points
#' give smoother estimates but may obscure fine spectral features.
#'
#' All spectral estimates are returned in decibels (dB).
#'
#' @seealso [parzen()] for computing a single Parzen-smoothed spectrum,
#'   [sample_spec()] for the raw sample spectrum,
#'   [plotts()] for time series plotting.
#'
#' @examples
#' # AR(2) process with complex roots
#' set.seed(123)
#' x <- arima.sim(model = list(ar = c(1.5, -0.75)), n = 200)
#'
#' # Default: periodogram + one Parzen estimate
#' plotts_parzen(x)
#'
#' # Add two more truncation points
#' plotts_parzen(x, trunc = c(10, 30))
#'
#' # Just compute, no plot
#' result <- plotts_parzen(x, plot = FALSE)
#' names(result)
#'
#' @export
plotts_parzen <- function(x, trunc = c(0, 0), plot = TRUE) {
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

  if (!is.numeric(trunc) || length(trunc) != 2L) {
    stop("`trunc` must be a numeric vector of length 2.", call. = FALSE)
  }
  if (any(trunc < 0)) {
    stop("`trunc` values must be non-negative.", call. = FALSE)
  }

  if (!is.logical(plot) || length(plot) != 1L || is.na(plot)) {
    stop("`plot` must be TRUE or FALSE.", call. = FALSE)
  }

  # Compute ACF for periodogram
  acf_result <- acf(x, lag.max = n - 1L, plot = FALSE)
  aut <- as.numeric(acf_result$acf)

  # Fourier frequencies (same as parzen() uses)
  freq <- (1:floor(n / 2)) / n
  nf <- length(freq)

  # Compute periodogram (unsmoothed)
  lags <- seq_len(n - 1L)
  cos_matrix <- cos(2 * pi * outer(freq, lags))
  pgram <- aut[1] + 2 * (cos_matrix %*% aut[-1])
  pgram <- as.numeric(pgram)
  pgram_db <- 10 * log10(pmax(pgram, .Machine$double.eps))

  # Use parzen() for smoothed estimates
  # Default truncation (trunc = 0 means use default M = 2*sqrt(n))
  pz_default <- parzen(x, db = TRUE, trunc = 0L, plot = FALSE)
  parzen_default <- pz_default$pzgram
  m_default <- pz_default$M

  # Additional truncation points
  m1 <- as.integer(trunc[1])
  m2 <- as.integer(trunc[2])

  parzen1 <- if (m1 > 0) parzen(x, db = TRUE, trunc = m1, plot = FALSE)$pzgram else NULL
  parzen2 <- if (m2 > 0) parzen(x, db = TRUE, trunc = m2, plot = FALSE)$pzgram else NULL

  # Collect truncation points used
  trunc_points <- m_default
  if (m1 > 0) trunc_points <- c(trunc_points, m1)
  if (m2 > 0) trunc_points <- c(trunc_points, m2)

  # Plot if requested
  if (plot) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))

    # Determine layout based on how many plots we have
    n_plots <- 1 + !is.null(parzen_default) + !is.null(parzen1) + !is.null(parzen2)
    if (n_plots <= 2) {
      par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
    } else {
      par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
    }

    # Plot periodogram (as vertical lines)
    min_db <- min(pgram_db)
    plot(freq, pgram_db, type = "n",
         xlab = "Frequency", ylab = "dB",
         main = "Periodogram")
    segments(freq, min_db, freq, pgram_db)

    # Plot Parzen with default M
    if (!is.null(parzen_default)) {
      plot(freq, parzen_default, type = "l",
           xlab = "Frequency", ylab = "dB",
           main = paste0("Parzen (M = ", m_default, ")"))
    }

    # Plot Parzen with m1
    if (!is.null(parzen1)) {
      plot(freq, parzen1, type = "l",
           xlab = "Frequency", ylab = "dB",
           main = paste0("Parzen (M = ", m1, ")"))
    }

    # Plot Parzen with m2
    if (!is.null(parzen2)) {
      plot(freq, parzen2, type = "l",
           xlab = "Frequency", ylab = "dB",
           main = paste0("Parzen (M = ", m2, ")"))
    }
  }

  # Build result
  result <- list(
    freq = freq,
    periodogram = pgram_db,
    parzen = parzen_default,
    parzen1 = parzen1,
    parzen2 = parzen2,
    trunc_points = trunc_points
  )

  class(result) <- "plotts_parzen"
  invisible(result)
}


#' @export
print.plotts_parzen <- function(x, ...) {
  cat("Periodogram and Parzen Spectral Estimates\n")
  cat("-----------------------------------------\n")
  cat("Frequencies:", length(x$freq), "Fourier frequencies\n")
  cat("Truncation points:", paste(x$trunc_points, collapse = ", "), "\n")

  cat("\nComponents:\n")
  cat("  $periodogram  - Raw periodogram (dB)\n")
  cat("  $parzen       - Parzen smoothed, M =", x$trunc_points[1], "\n")
  if (!is.null(x$parzen1)) {
    cat("  $parzen1      - Parzen smoothed, M =", x$trunc_points[2], "\n")
  }
  if (!is.null(x$parzen2)) {
    cat("  $parzen2      - Parzen smoothed, M =", x$trunc_points[3], "\n")
  }

  invisible(x)
}
