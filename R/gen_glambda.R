#' Generate G-Lambda Process Realization
#'
#' Generate a realization from a G-lambda process - an AR process on a non-linear
#' time scale with interpolation back to a regular grid.
#'
#' @param n Integer. Length of the realization to generate.
#' @param lambda Numeric. The transformation parameter controlling how the
#'   dual time relates to original time:
#'   \itemize{
#'     \item `lambda = 0`: Logarithmic (exponential) transformation
#'     \item `lambda = 1`: Linear transformation
#'     \item `lambda = 2`: Quadratic transformation
#'   }
#' @param phi Numeric vector. AR coefficients (default 0 = white noise).
#' @param offset Numeric. Time offset (shift) parameter. Must be positive.
#'   Affects the starting point of the dual transformation. Default 20.
#' @param sigma2 Numeric. Innovation variance. Default 1.
#' @param plot Logical. If `TRUE` (default), plot the realization.
#' @param seed Integer or NULL. Random seed for reproducibility.
#'   If NULL (default), no seed is set.
#'
#' @return Numeric vector of length `n` containing the G-lambda realization.
#'
#' @details
#' The G-lambda process models non-stationary time series by warping a
#' stationary AR process from a dual time domain back to regular time.
#' This is useful for modeling processes with time-varying spectral content.
#'
#' The method:
#' \enumerate{
#'   \item Computes the h parameter for the given lambda and offset
#'   \item Generates an AR series in the dual time domain
#'   \item Creates a non-uniform time grid in the dual domain
#'   \item Interpolates to equally-spaced original time points
#' }
#'
#' Different lambda values produce different non-stationarity patterns:
#' \itemize{
#'   \item `lambda = 0`: Exponentially evolving spectrum
#'   \item `lambda = 1`: Linearly evolving spectrum
#'   \item `lambda > 1`: Faster than linear evolution
#' }
#'
#' @seealso [trans_to_dual()], [trans_to_original()] for the underlying
#'   transformations, [is_sample()] for instantaneous spectrum estimation.
#'
#' @references
#' Jiang, H., Gray, H. L., & Woodward, W. A. (2006). Time-frequency analysis
#' of musical rhythm. *Behavior Research Methods*, 38(1), 56-63.
#'
#' @examples
#' # White noise on G-lambda scale (lambda = 0)
#' x <- gen_glambda(n = 200, lambda = 0, seed = 123)
#'
#' # AR(1) on G-lambda scale (lambda = 1)
#' x <- gen_glambda(n = 200, lambda = 1, phi = 0.9, seed = 123)
#'
#' # AR(2) with quadratic time warp
#' x <- gen_glambda(n = 200, lambda = 2, phi = c(1.5, -0.75), seed = 123)
#'
#' # Without plot
#' x <- gen_glambda(n = 100, lambda = 0.5, phi = 0.7, plot = FALSE, seed = 42)
#'
#' @export
gen_glambda <- function(n, lambda, phi = 0, offset = 20,
                        sigma2 = 1, plot = TRUE, seed = NULL) {

  # Input validation
  if (!is.numeric(n) || length(n) != 1L || n < 1L) {
    stop("`n` must be a positive integer.", call. = FALSE)
  }
  n <- as.integer(n)

  if (!is.numeric(lambda) || length(lambda) != 1L) {
    stop("`lambda` must be a single numeric value.", call. = FALSE)
  }

  if (!is.numeric(phi)) {
    stop("`phi` must be a numeric vector.", call. = FALSE)
  }

  if (!is.numeric(offset) || length(offset) != 1L || offset <= 0) {
    stop("`offset` must be a positive number.", call. = FALSE)
  }

  if (!is.numeric(sigma2) || length(sigma2) != 1L || sigma2 <= 0) {
    stop("`sigma2` must be a positive number.", call. = FALSE)
  }

  if (!is.logical(plot) || length(plot) != 1L || is.na(plot)) {
    stop("`plot` must be TRUE or FALSE.", call. = FALSE)
  }

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Burnin for arima.sim
  burnin <- 600L

  # Equally-spaced time points in original domain
  eq_x <- seq(offset + 1, offset + n)

  if (lambda == 0) {
    # Logarithmic transformation case
    h <- .h_value_glambda(n = n, shift = offset, lambda = 0)
    k <- log(offset + 1) / log(h)
    eu_x <- h^seq(k, k + n - 1)

    # Generate AR series in dual domain
    ar_coefs <- if (all(phi == 0)) numeric(0) else phi
    yt <- arima.sim(n = n + burnin, n.start = burnin,
                    model = list(ar = ar_coefs, ma = numeric(0)))
    yt <- yt[-seq_len(burnin)] * sqrt(sigma2)

  } else {
    # Power-law transformation case
    h <- .h_value_glambda(n = n, shift = offset, lambda = lambda)
    k <- ((offset + 1)^lambda - 1) / (h * lambda)
    eu_x <- (h * lambda * seq(k, k + n - 1) + 1)^(1 / lambda)

    # Generate AR series in dual domain
    ar_coefs <- if (all(phi == 0)) numeric(0) else phi
    yt <- arima.sim(n = n + burnin, n.start = burnin,
                    model = list(ar = ar_coefs, ma = numeric(0)))
    yt <- yt[-seq_len(burnin)] * sqrt(sigma2)
  }

  # Interpolate to equally-spaced grid
  tvfdata <- approx(eu_x, yt, xout = eq_x, rule = 2)$y

  # Plot if requested
  if (plot) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))

    par(mar = c(4, 4, 2, 1))
    t <- seq_len(n)
    plot(t, tvfdata, type = "o", cex = 0.5, pch = 16, lwd = 0.75,
         xlab = "Time", ylab = "",
         main = paste0("G-Lambda Realization (lambda=", lambda, ")"))
  }

  invisible(tvfdata)
}


# Internal helper: compute h value for G-lambda transformation
.h_value_glambda <- function(n, shift, lambda, m = n) {
  if (lambda == 0) {
    ((shift + n) / (shift + 1))^(1 / (m - 1))
  } else {
    ((shift + n)^lambda - (shift + 1)^lambda) / ((m - 1) * lambda)
  }
}
