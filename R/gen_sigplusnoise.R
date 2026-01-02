#' Generate Signal Plus Noise Process
#'
#' Generate a realization from a signal plus colored noise model.
#'
#' @param n Integer, length of output series.
#' @param b0 Numeric, intercept (default 0).
#' @param b1 Numeric, linear trend slope (default 0).
#' @param coef Numeric vector, amplitudes for cosine terms (default c(0, 0)).
#' @param freq Numeric vector, frequencies for cosine terms (default c(0, 0)).
#' @param psi Numeric vector, phase shifts for cosine terms (default c(0, 0)).
#' @param phi Numeric vector, AR coefficients for noise (default 0 = white noise).
#' @param vara Numeric, noise variance (default 1).
#' @param plot Logical, whether to plot the series (default TRUE).
#' @param seed Integer, random seed for reproducibility (default NULL).
#'
#' @return Numeric vector of length n.
#' @export
#'
#' @details
#' Generates from the model:
#' \deqn{X_t = b_0 + b_1 t + \sum_{j=1}^{k} A_j \cos(2\pi f_j t + \psi_j) + Z_t}
#' where \eqn{Z_t} is an AR process with coefficients \code{phi} and
#' innovation variance \code{vara}.
#'
#' @examples
#' # Linear trend plus white noise
#' x <- gen_sigplusnoise(n = 200, b0 = 10, b1 = 0.05)
#'
#' # Seasonal signal plus AR(1) noise
#' x <- gen_sigplusnoise(n = 200, coef = c(5, 0), freq = c(1/12, 0),
#'                       phi = 0.7, seed = 123)
#'
#' # Two cosine components plus noise
#' x <- gen_sigplusnoise(n = 200, coef = c(3, 2), freq = c(0.1, 0.25),
#'                       psi = c(0, pi/4), seed = 456)
gen_sigplusnoise <- function(n,
                             b0 = 0,
                             b1 = 0,
                             coef = c(0, 0),
                             freq = c(0, 0),
                             psi = c(0, 0),
                             phi = 0,
                             vara = 1,
                             plot = TRUE,
                             seed = NULL) {

  # Input validation
  if (!is.numeric(n) || length(n) != 1 || n < 1) {
    stop("`n` must be a positive integer")
  }
  if (!is.numeric(b0) || length(b0) != 1) {
    stop("`b0` must be a single number")
  }
  if (!is.numeric(b1) || length(b1) != 1) {
    stop("`b1` must be a single number")
  }
  if (!is.numeric(coef)) stop("`coef` must be numeric")
  if (!is.numeric(freq)) stop("`freq` must be numeric")
  if (!is.numeric(psi)) stop("`psi` must be numeric")
  if (length(coef) != length(freq) || length(freq) != length(psi)) {
    stop("`coef`, `freq`, and `psi` must have the same length")
  }
  if (!is.numeric(vara) || vara <= 0) {
    stop("`vara` must be positive")
  }

  t <- seq_len(n)

  # Generate AR noise component
  zt <- gen_arma(n = n, phi = phi, theta = 0, vara = vara,
                 seed = seed, plot = FALSE)

  # Build signal: intercept + trend + cosine components
  signal <- b0 + b1 * t

  for (j in seq_along(coef)) {
    if (coef[j] != 0) {
      signal <- signal + coef[j] * cos(2 * pi * freq[j] * t + psi[j])
    }
  }

  # Signal plus noise
  x <- signal + zt

  if (plot) {
    plot(t, x, type = "l", xlab = "Time", ylab = "",
         main = "Signal Plus Noise")
    invisible(x)
  } else {
    x
  }
}
