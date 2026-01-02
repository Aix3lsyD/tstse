#' Plot Multiresolution Analysis
#'
#' Computes and plots the Multiresolution Analysis (MRA) of a time series,
#' showing cumulative smooth approximations at each decomposition level.
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
#' @param plot Logical. If `TRUE` (default), create the multi-panel MRA plot.
#'
#' @return A list with class `"plotts_mra"` containing:
#'   \item{mra}{The MRA object from [waveslim::mra()].}
#'   \item{smooth}{Matrix of cumulative smooth approximations. Row i contains
#'     the smooth approximation S(i-1), where S0 is the original data
#'     reconstructed and SJ is just the smooth coefficients.}
#'   \item{n_levels}{Number of decomposition levels.}
#'   \item{wavelet}{Wavelet filter used.}
#'   \item{n}{Length of the original series.}
#'
#' @details
#' Multiresolution Analysis decomposes a signal into smooth approximations
#' at different scales. Unlike [plotts_dwt()] which shows the wavelet
#' coefficients, MRA shows the cumulative reconstructions:
#'
#' \itemize{
#'   \item **S0**: Original data (reconstructed from all components)
#'   \item **S1**: S0 minus finest detail (d1 removed)
#'   \item **S2**: S1 minus next detail (d1 and d2 removed)
#'   \item **...**
#'   \item **SJ**: Smooth component only (all details removed)
#' }
#'
#' This provides a multi-scale view of the trend, with each successive
#' smooth showing a coarser approximation.
#'
#' @note
#' This function requires the `waveslim` package. Install it with
#' `install.packages("waveslim")` if needed.
#'
#' @seealso [plotts_dwt()] for the DWT coefficient plot,
#'   [waveslim::mra()] for the underlying MRA computation.
#'
#' @examples
#' \dontrun{
#' # Generate test signal: trend + seasonal + noise
#' t <- 1:256
#' x <- 0.01 * t + 2 * sin(2 * pi * t / 32) + rnorm(256, sd = 0.5)
#'
#' # Plot MRA with 4 levels
#' result <- plotts_mra(x, n_levels = 4)
#'
#' # Use Haar wavelet
#' result2 <- plotts_mra(x, n_levels = 4, wavelet = "haar")
#'
#' # Just compute, no plot
#' result3 <- plotts_mra(x, plot = FALSE)
#' }
#'
#' @export
plotts_mra <- function(x, n_levels = 4L, wavelet = "S8", plot = TRUE) {
  # Check for waveslim package
  if (!requireNamespace("waveslim", quietly = TRUE)) {
    stop(
      "The 'waveslim' package is required for plotts_mra().\n",
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

  # Compute MRA
  x_mra <- waveslim::mra(x, J = n_levels, wf = wf,
                         boundary = "periodic", method = "dwt")

  # Build cumulative smooth approximations
  # sm[J+1,] = smooth component (sJ)
  # sm[J,] = sJ + dJ
  # sm[J-1,] = sJ + dJ + d(J-1)
  # ...
  # sm[1,] = sJ + dJ + ... + d1 = original data
  smooth <- matrix(0, nrow = n_levels + 1, ncol = n)

  # Start with the smooth component
  smooth[n_levels + 1, ] <- x_mra[[n_levels + 1]]

  # Add details cumulatively (from coarsest to finest)
  for (i in seq_len(n_levels)) {
    j <- n_levels - i + 1  # Go from n_levels down to 1
    smooth[j, ] <- smooth[j + 1, ] + x_mra[[j]]
  }

  # Find global y-limits
  ymin <- min(x, smooth)
  ymax <- max(x, smooth)

  # Plot if requested
  if (plot) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))

    n_panels <- n_levels + 2  # data + (n_levels + 1) smooth approximations
    par(mfrow = c(n_panels, 1), mar = c(1, 4, 0.5, 1), oma = c(3, 0, 1, 0))

    tt <- seq_len(n)

    # Plot original data
    if (n < 200) {
      plot(tt, x, type = "o", pch = 16, cex = 0.6,
           xlab = "", ylab = "", xaxt = "n", yaxt = "n",
           xlim = c(0, n), ylim = c(ymin, ymax))
    } else {
      plot(tt, x, type = "l",
           xlab = "", ylab = "", xaxt = "n", yaxt = "n",
           xlim = c(0, n), ylim = c(ymin, ymax))
    }
    mtext("Data", side = 2, line = 1, cex = 0.8)

    # Plot smooth approximations S0, S1, ..., SJ
    for (i in seq_len(n_levels + 1)) {
      plot(tt, smooth[i, ], type = "l",
           xlab = "", ylab = "", xaxt = "n", yaxt = "n",
           xlim = c(0, n), ylim = c(ymin, ymax))
      mtext(paste0("S", i - 1), side = 2, line = 1, cex = 0.8)
    }

    # Add x-axis label to bottom panel
    axis(1)
    mtext("Time", side = 1, line = 2, cex = 1)
  }

  # Build result
  result <- list(
    mra = x_mra,
    smooth = smooth,
    n_levels = n_levels,
    wavelet = wavelet,
    n = n
  )

  class(result) <- "plotts_mra"
  invisible(result)
}


#' @export
print.plotts_mra <- function(x, ...) {
  cat("Multiresolution Analysis\n")
  cat("------------------------\n")
  cat("Sample size:", x$n, "\n")
  cat("Levels:", x$n_levels, "\n")
  cat("Wavelet:", x$wavelet, "\n")
  cat("\nSmooth approximations:\n")
  cat("  S0: Full reconstruction (", x$n, " points)\n", sep = "")
  for (j in seq_len(x$n_levels)) {
    cat("  S", j, ": d1-d", j, " removed\n", sep = "")
  }
  invisible(x)
}
