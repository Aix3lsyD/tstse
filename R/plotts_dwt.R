#' Plot Discrete Wavelet Transform
#'
#' Computes and plots the Discrete Wavelet Transform (DWT) of a time series,
#' displaying the original data, detail coefficients at each level, and the
#' final smooth coefficients.
#'
#' @param x Numeric vector. The time series to transform. Length must be
#'   divisible by \eqn{2^{n\_levels}}.
#' @param n_levels Integer. Number of decomposition levels. Default is 4.
#' @param wavelet Character. Wavelet filter to use. Options are:
#'   \itemize{
#'     \item `"S8"` (default): Symmlet 8 (least asymmetric, length 8)
#'     \item `"haar"`: Haar wavelet
#'     \item `"D4"`: Daubechies 4
#'     \item `"D6"`: Daubechies 6
#'     \item `"D8"`: Daubechies 8
#'   }
#' @param plot Logical. If `TRUE` (default), create the multi-panel DWT plot.
#'
#' @return A list with class `"plotts_dwt"` containing:
#'   \item{dwt}{The DWT object from [waveslim::dwt()].}
#'   \item{n_levels}{Number of decomposition levels.}
#'   \item{wavelet}{Wavelet filter used.}
#'   \item{n}{Length of the original series.}
#'
#' @details
#' The Discrete Wavelet Transform decomposes a time series into components
#' at different scales:
#' \itemize{
#'   \item **Detail coefficients (d1, d2, ...)**: Capture fluctuations at
#'     each scale. d1 captures the finest (highest frequency) details,
#'     with coarser details at higher levels.
#'   \item **Smooth coefficients (sJ)**: The remaining low-frequency trend
#'     after J levels of decomposition.
#' }
#'
#' The sample size must be divisible by \eqn{2^{n\_levels}} for the
#' standard DWT with periodic boundary conditions.
#'
#' @note
#' This function requires the `waveslim` package. Install it with
#' `install.packages("waveslim")` if needed.
#'
#' @seealso [plotts_mra()] for the multiresolution analysis plot,
#'   [waveslim::dwt()] for the underlying DWT computation.
#'
#' @examples
#' \dontrun{
#' # Generate test signal: trend + seasonal + noise
#' t <- 1:256
#' x <- 0.01 * t + 2 * sin(2 * pi * t / 32) + rnorm(256, sd = 0.5)
#'
#' # Plot DWT with 4 levels
#' result <- plotts_dwt(x, n_levels = 4)
#'
#' # Use Haar wavelet
#' result2 <- plotts_dwt(x, n_levels = 4, wavelet = "haar")
#'
#' # Just compute, no plot
#' result3 <- plotts_dwt(x, plot = FALSE)
#' }
#'
#' @export
plotts_dwt <- function(x, n_levels = 4L, wavelet = "S8", plot = TRUE) {
  # Check for waveslim package
  if (!requireNamespace("waveslim", quietly = TRUE)) {
    stop(
      "The 'waveslim' package is required for plotts_dwt().\n",
      "Install it with: install.packages(\"waveslim\")",
      call. = FALSE
    )
  }

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

  if (!is.numeric(n_levels) || length(n_levels) != 1L || n_levels < 1L) {
    stop("`n_levels` must be a positive integer.", call. = FALSE)
  }
  n_levels <- as.integer(n_levels)

  # Check sample size divisibility
  if (n %% (2^n_levels) != 0) {
    stop(paste0("Sample size (", n, ") must be divisible by 2^n_levels (",
                2^n_levels, ")."), call. = FALSE)
  }
  if (2^n_levels > n) {
    stop("2^n_levels exceeds sample size.", call. = FALSE)
  }

  # Wavelet filter mapping
  wavelet <- match.arg(wavelet, c("S8", "haar", "D4", "D6", "D8"))
  wf_map <- c(S8 = "la8", haar = "haar", D4 = "d4", D6 = "d6", D8 = "d8")
  wf <- wf_map[wavelet]

  if (!is.logical(plot) || length(plot) != 1L || is.na(plot)) {
    stop("`plot` must be TRUE or FALSE.", call. = FALSE)
  }

  # Compute DWT
  x_dwt <- waveslim::dwt(x, n.levels = n_levels, wf = wf, boundary = "periodic")

  # Plot if requested
  if (plot) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))

    n_panels <- n_levels + 2  # data + n_levels details + 1 smooth
    par(mfrow = c(n_panels, 1), mar = c(1, 4, 0.5, 1), oma = c(3, 0, 1, 0))

    tt <- seq_len(n)

    # Find global y-limits for detail/smooth coefficients
    ymin <- 0
    ymax <- 0
    for (j in seq_len(n_levels + 1)) {
      ymin <- min(ymin, x_dwt[[j]])
      ymax <- max(ymax, x_dwt[[j]])
    }

    # Plot original data
    if (n < 200) {
      plot(tt, x, type = "o", pch = 16, cex = 0.6,
           xlab = "", ylab = "", xaxt = "n", yaxt = "n",
           xlim = c(0, n))
    } else {
      plot(tt, x, type = "l",
           xlab = "", ylab = "", xaxt = "n", yaxt = "n",
           xlim = c(0, n))
    }
    mtext("Data", side = 2, line = 1, cex = 0.8)

    # Plot detail coefficients
    for (j in seq_len(n_levels)) {
      n_coef <- n / 2^j
      k <- seq_len(n_coef)
      # Center each coefficient in its time interval
      td <- k * 2^j - 0.5 * 2^j + 0.5

      plot(td, x_dwt[[j]], type = "h", lwd = 1,
           xlab = "", ylab = "", xaxt = "n", yaxt = "n",
           xlim = c(0, n), ylim = c(ymin, ymax))
      abline(h = 0, col = "gray")
      mtext(paste0("d", j), side = 2, line = 1, cex = 0.8)
    }

    # Plot smooth coefficients (last level)
    n_coef <- n / 2^n_levels
    k <- seq_len(n_coef)
    td <- k * 2^n_levels - 0.5 * 2^n_levels + 0.5

    plot(td, x_dwt[[n_levels + 1]], type = "h", lwd = 1,
         xlab = "", ylab = "", yaxt = "n",
         xlim = c(0, n), ylim = c(ymin, ymax))
    abline(h = 0, col = "gray")
    mtext(paste0("s", n_levels), side = 2, line = 1, cex = 0.8)
    mtext("Time", side = 1, line = 2, cex = 1)
  }

  # Rename wavelet attribute back to user-friendly name
  attr(x_dwt, "wavelet") <- wavelet

  # Build result
  result <- list(
    dwt = x_dwt,
    n_levels = n_levels,
    wavelet = wavelet,
    n = n
  )

  class(result) <- "plotts_dwt"
  invisible(result)
}


#' @export
print.plotts_dwt <- function(x, ...) {
  cat("Discrete Wavelet Transform\n")
  cat("--------------------------\n")
  cat("Sample size:", x$n, "\n")
  cat("Levels:", x$n_levels, "\n")
  cat("Wavelet:", x$wavelet, "\n")
  cat("\nComponents:\n")
  for (j in seq_len(x$n_levels)) {
    cat("  d", j, ": ", length(x$dwt[[j]]), " detail coefficients\n", sep = "")
  }
  cat("  s", x$n_levels, ": ", length(x$dwt[[x$n_levels + 1]]),
      " smooth coefficients\n", sep = "")
  invisible(x)
}
