#' Theoretical Instantaneous Spectrum (G-Lambda)
#'
#' Computes and plots the theoretical instantaneous spectrum based on an AR
#' model using the G-lambda transformation.
#'
#' @param n Integer. Number of time points.
#' @param phi Numeric vector. AR coefficients. Default is 0 (white noise).
#' @param sigma2 Numeric. White noise variance. Default is 1.
#' @param lambda Numeric. The transformation parameter controlling how the
#'   dual time relates to original time. Common values:
#'   \itemize{
#'     \item `lambda = 0`: Logarithmic transformation
#'     \item `lambda = 1`: Linear transformation
#'     \item `lambda = 2`: Quadratic transformation
#'   }
#' @param offset Numeric. Time offset (shift) parameter. Must be positive.
#' @param n_freq Integer. Number of frequency points. Default is 500.
#' @param plot Logical. If `TRUE` (default), create the time-frequency plot.
#'
#' @return A list with class `"is_glambda"` containing:
#'   \item{time}{Time indices (1 to n).}
#'   \item{freq}{Instantaneous frequencies at each time-frequency point.}
#'   \item{gfreq}{G-frequencies (transformed frequencies).}
#'   \item{spec}{Theoretical spectral values (in dB, normalized).}
#'   \item{phi}{AR coefficients used.}
#'   \item{lambda}{Transformation parameter used.}
#'   \item{offset}{Offset parameter used.}
#'
#' @details
#' This function computes the **theoretical** instantaneous spectrum based on
#' the AR model parameters, rather than estimating it from data (see
#' [is_sample()] for the sample-based version).
#'
#' The theoretical AR spectrum at frequency f is:
#' \deqn{S(f) = \frac{\sigma^2}{|\phi(e^{2\pi i f})|^2}}
#'
#' where \eqn{\phi(z) = 1 - \phi_1 z - \phi_2 z^2 - \ldots}
#'
#' This spectrum is then mapped to instantaneous frequencies using the
#' G-lambda transformation.
#'
#' @seealso [is_sample()] for sample-based instantaneous spectrum,
#'   [true_spec()] for AR/ARMA theoretical spectrum.
#'
#' @examples
#' # AR(1) with phi = 0.9 (low frequency dominance)
#' is_glambda(n = 200, phi = 0.9, lambda = 1, offset = 10)
#'
#' # AR(2) with complex roots
#' is_glambda(n = 200, phi = c(1.5, -0.75), lambda = 0, offset = 10)
#'
#' # White noise
#' is_glambda(n = 200, phi = 0, lambda = 1, offset = 10)
#'
#' @export
is_glambda <- function(n, phi = 0, sigma2 = 1, lambda, offset,
                       n_freq = 500L, plot = TRUE) {
  # Input validation
  if (!is.numeric(n) || length(n) != 1L || n < 1L) {
    stop("`n` must be a positive integer.", call. = FALSE)
  }
  n <- as.integer(n)

  if (!is.numeric(phi)) {
    stop("`phi` must be numeric.", call. = FALSE)
  }
  # Handle phi = 0 as no AR component
  if (length(phi) == 1L && phi == 0) {
    phi <- numeric(0)
  }

  if (!is.numeric(sigma2) || length(sigma2) != 1L || sigma2 <= 0) {
    stop("`sigma2` must be a positive number.", call. = FALSE)
  }

  if (!is.numeric(lambda) || length(lambda) != 1L) {
    stop("`lambda` must be a single numeric value.", call. = FALSE)
  }

  if (!is.numeric(offset) || length(offset) != 1L || offset <= 0) {
    stop("`offset` must be a positive number.", call. = FALSE)
  }

  if (!is.numeric(n_freq) || length(n_freq) != 1L || n_freq < 10L) {
    stop("`n_freq` must be an integer >= 10.", call. = FALSE)
  }
  n_freq <- as.integer(n_freq)

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

  # Helper function: theoretical AR spectrum at frequency f
  spec_ar_f <- function(f, phi, sigma2) {
    p <- length(phi)
    if (p == 0) {
      return(sigma2)
    }
    # phi(exp(i*2*pi*f)) = 1 - phi_1*z - phi_2*z^2 - ...
    z <- exp(2i * pi * f)
    powers <- z^(0:p)
    phi_extended <- c(1, -phi)
    phi_func <- sum(phi_extended * powers)
    phi_norm <- Mod(phi_func)^2
    sigma2 / phi_norm
  }

  # Frequency grid (avoid exactly 0)
  freq <- seq(0.001, 0.5, length.out = n_freq)

  # Compute theoretical spectrum
  spec <- sapply(freq, spec_ar_f, phi = phi, sigma2 = sigma2)
  spec_db <- 10 * log10(spec)
  spec_db <- spec_db - min(spec_db)  # Normalize to start at 0

  # Compute h value
  h <- h_value(n = n, shift = offset, lambda = lambda)

  # Transform frequencies
  if (lambda == 0) {
    gfreq <- h^(1 / freq)
  } else {
    gfreq <- freq / h
  }

  # Build plot matrix
  n_total <- n * n_freq

  time_vec <- rep(seq_len(n), each = n_freq)
  gfreq_vec <- rep(gfreq, times = n)
  spec_vec <- rep(spec_db, times = n)

  # Compute instantaneous frequencies
  time_shifted <- time_vec + offset
  inst_freq_vec <- mapply(inst_freq, lambda, time_shifted, gfreq_vec)

  # Plot if requested
  if (plot) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))

    par(mar = c(4, 4, 2, 1))

    # Color based on spectral intensity (darker = higher power)
    intensity <- 1 - spec_vec / max(spec_vec)
    colors <- hcl(h = 260, c = 1, l = 100 * intensity)

    plot(time_vec, inst_freq_vec,
         ylim = c(0, 0.5),
         pch = ".", cex = 2,
         col = colors,
         xlab = "Time", ylab = "Frequency",
         main = paste0("Theoretical IS (AR, lambda=", lambda, ")"))
  }

  # Build result
  result <- list(
    time = time_vec,
    freq = inst_freq_vec,
    gfreq = gfreq_vec,
    spec = spec_vec,
    phi = phi,
    lambda = lambda,
    offset = offset
  )

  class(result) <- "is_glambda"
  invisible(result)
}


#' @export
print.is_glambda <- function(x, ...) {
  cat("Theoretical Instantaneous Spectrum (G-Lambda)\n")
  cat("---------------------------------------------\n")
  cat("Time points:", max(x$time), "\n")
  p <- length(x$phi)
  if (p == 0) {
    cat("Model: White noise\n")
  } else {
    cat("Model: AR(", p, ")\n", sep = "")
    cat("Phi:", round(x$phi, 4), "\n")
  }
  cat("Lambda:", x$lambda, "\n")
  cat("Offset:", x$offset, "\n")
  invisible(x)
}
