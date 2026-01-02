#' Instantaneous Spectrum (Sample-Based)
#'
#' Computes and plots the sample instantaneous spectrum of a time series
#' using the dual time transformation method.
#'
#' @param x Numeric vector. The time series to analyze.
#' @param lambda Numeric. The transformation parameter controlling how the
#'   dual time relates to original time. Common values:
#'   \itemize{
#'     \item `lambda = 0`: Logarithmic transformation
#'     \item `lambda = 1`: Linear transformation
#'     \item `lambda = 2`: Quadratic transformation
#'   }
#' @param offset Numeric. Time offset (shift) parameter. Must be positive.
#'   Affects the starting point of the dual transformation.
#' @param plot Logical. If `TRUE` (default), create the time-frequency plot.
#'
#' @return A list with class `"is_sample"` containing:
#'   \item{time}{Time indices (1 to n).}
#'   \item{freq}{Instantaneous frequencies at each time-frequency point.}
#'   \item{gfreq}{G-frequencies (transformed frequencies).}
#'   \item{spec}{Spectral values (in dB, normalized).}
#'   \item{lambda}{Transformation parameter used.}
#'   \item{offset}{Offset parameter used.}
#'
#' @details
#' The instantaneous spectrum provides a time-frequency representation
#' based on the Wold-Cramer spectral theory. Unlike the Wigner-Ville
#' distribution, this method uses a dual time transformation to map
#' a non-stationary process to a stationary one in the dual domain.
#'
#' The method:
#' \enumerate{
#'   \item Transforms the data to the dual time domain using [trans_to_dual()]
#'   \item Computes the periodogram in the dual domain
#'   \item Maps frequencies back to instantaneous frequencies using the
#'     transformation parameter lambda
#' }
#'
#' The parameter lambda controls the type of non-stationarity:
#' \itemize{
#'   \item `lambda = 0`: Exponentially evolving spectrum
#'   \item `lambda = 1`: Linearly evolving spectrum
#'   \item `lambda > 1`: Faster than linear evolution
#' }
#'
#' @seealso [trans_to_dual()] for the dual time transformation,
#'   [wv()] for Wigner-Ville distribution.
#'
#' @references
#' Jiang, H., Gray, H. L., & Woodward, W. A. (2006). Time-frequency analysis
#' of musical rhythm. *Behavior Research Methods*, 38(1), 56-63.
#'
#' @examples
#' # Chirp signal (frequency increases over time)
#' t <- 1:200
#' x <- cos(2 * pi * cumsum(seq(0.02, 0.2, length.out = 200)))
#' is_sample(x, lambda = 1, offset = 10)
#'
#' # AR process
#' set.seed(123)
#' x <- arima.sim(model = list(ar = 0.9), n = 200)
#' is_sample(x, lambda = 0, offset = 10)
#'
#' @export
is_sample <- function(x, lambda, offset, plot = TRUE) {
  # Input validation
  if (!is.numeric(x)) {
    stop("`x` must be a numeric vector.", call. = FALSE)
  }
  n <- length(x)
  if (n < 10L) {
    stop("`x` must have at least 10 observations.", call. = FALSE)
  }
  if (anyNA(x)) {
    stop("`x` contains NA values.", call. = FALSE)
  }

  if (!is.numeric(lambda) || length(lambda) != 1L) {
    stop("`lambda` must be a single numeric value.", call. = FALSE)
  }

  if (!is.numeric(offset) || length(offset) != 1L || offset <= 0) {
    stop("`offset` must be a positive number.", call. = FALSE)
  }

  if (!is.logical(plot) || length(plot) != 1L || is.na(plot)) {
    stop("`plot` must be TRUE or FALSE.", call. = FALSE)
  }

  # Helper function: compute h value for transformation
  h_value <- function(n, shift, lambda, m = n) {
    if (lambda == 0) {
      ((shift + n) / (shift + 1))^(1 / (m - 1))
    } else {
      ((shift + n)^lambda - (shift + 1)^lambda) / ((m - 1) * lambda)
    }
  }

  # Helper function: compute instantaneous frequency
  inst_freq <- function(lambda, time, gfreq) {
    if (lambda == 0) {
      1 / (time * (gfreq - 1))
    } else {
      1 / ((time^lambda + lambda / gfreq)^(1 / lambda) - time)
    }
  }

  # Transform to dual time domain
  dual_data <- trans_to_dual(x, lambda = lambda, offset = offset, plot = FALSE)

  # Get the dual series
  y_dual <- dual_data$int_y
  if (is.null(y_dual) || length(y_dual) == 0) {
    stop("trans_to_dual() did not return a valid dual series.", call. = FALSE)
  }

  # Compute periodogram in dual domain
  pgram <- spec.pgram(y_dual, taper = 0, plot = FALSE)
  per <- 10 * log10(pgram$spec)
  per <- per - min(per)  # Normalize to start at 0

  # Compute h value
  h <- h_value(n = n, shift = offset, lambda = lambda)

  # Transform frequencies
  if (lambda == 0) {
    gfreq <- h^(1 / pgram$freq)
  } else {
    gfreq <- pgram$freq / h
  }

  # Build plot matrix
  n_freq <- length(gfreq)
  n_total <- n * n_freq

  time_vec <- rep(seq_len(n), each = n_freq)
  gfreq_vec <- rep(gfreq, times = n)
  spec_vec <- rep(per, times = n)

  # Compute instantaneous frequencies
  time_shifted <- time_vec + offset
  inst_freq_vec <- mapply(inst_freq, lambda, time_shifted, gfreq_vec)

  # Plot if requested
  if (plot) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))

    par(mar = c(4, 4, 2, 1))

    # Color based on spectral intensity (darker = higher power)
    # Using HCL color space for perceptually uniform grayscale
    intensity <- 1 - spec_vec / max(spec_vec)
    colors <- hcl(h = 260, c = 1, l = 100 * intensity)

    plot(time_vec, inst_freq_vec,
         ylim = c(0, 0.5),
         pch = ".", cex = 2,
         col = colors,
         xlab = "Time", ylab = "Frequency",
         main = paste0("Instantaneous Spectrum (lambda=", lambda, ")"))
  }

  # Build result
  result <- list(
    time = time_vec,
    freq = inst_freq_vec,
    gfreq = gfreq_vec,
    spec = spec_vec,
    lambda = lambda,
    offset = offset
  )

  class(result) <- "is_sample"
  invisible(result)
}


#' @export
print.is_sample <- function(x, ...) {
  cat("Instantaneous Spectrum (Sample)\n")
  cat("-------------------------------\n")
  cat("Time points:", max(x$time), "\n")
  cat("Lambda:", x$lambda, "\n")
  cat("Offset:", x$offset, "\n")
  cat("Frequency range: 0 to 0.5\n")
  invisible(x)
}
