#' True ARMA Spectral Density
#'
#' Compute and optionally plot the theoretical spectral density for an ARMA model.
#'
#' @param phi Numeric vector, AR coefficients (default 0 = none).
#' @param theta Numeric vector, MA coefficients (default 0 = none).
#' @param vara Numeric, white noise variance (default 1).
#' @param n_freq Integer, number of frequency points (default 251).
#' @param plot Logical, whether to plot spectrum (default TRUE).
#'
#' @return A list with components:
#'   \item{freq}{Frequency values (0 to 0.5)}
#'   \item{spec}{Spectral density values (dB)}
#' @export
#'
#' @details
#' Computes the theoretical spectral density for an ARMA(p,q) model:
#' \deqn{S(f) = \frac{\sigma_a^2 |\Theta(e^{-2\pi i f})|^2}{\gamma_0 |\Phi(e^{-2\pi i f})|^2}}
#' where \eqn{\Phi} and \eqn{\Theta} are the AR and MA polynomials,
#' \eqn{\sigma_a^2} is the white noise variance, and \eqn{\gamma_0} is the
#' process variance.
#'
#' Results are returned in decibels (dB): \eqn{10 \log_{10}(S(f))}.
#'
#' @seealso \code{\link{true_acf}} for theoretical autocorrelations.
#'
#' @examples
#' # AR(1) spectrum - peak at f=0
#' true_spec(phi = 0.9)
#'
#' # AR(2) spectrum - peak near f=0.1
#' true_spec(phi = c(1.5, -0.75))
#'
#' # ARMA(1,1) spectrum
#' true_spec(phi = 0.7, theta = 0.4)
#'
#' # Without plot
#' result <- true_spec(phi = 0.9, plot = FALSE)
true_spec <- function(phi = 0, theta = 0, vara = 1, n_freq = 251L, plot = TRUE) {

  # Input validation
  if (!is.numeric(phi)) stop("`phi` must be numeric")
  if (!is.numeric(theta)) stop("`theta` must be numeric")
  if (!is.numeric(vara) || vara <= 0) {
    stop("`vara` must be positive")
  }
  if (!is.numeric(n_freq) || n_freq < 2) {
    stop("`n_freq` must be at least 2")
  }

  # Effective orders
  p <- if (all(phi == 0)) 0L else length(phi)
  q <- if (all(theta == 0)) 0L else length(theta)

  # Compute process variance
  gvar <- .arma_variance(phi, theta, vara, p, q)

  # Frequency grid: 0, 1/500, 2/500, ..., (n_freq-1)/500
  # With n_freq=251 (default), this gives 0 to 0.5 in steps of 0.002
  freq <- (seq_len(n_freq) - 1) / 500

  # Compute spectral density at each frequency
  spec <- numeric(n_freq)

  for (fi in seq_len(n_freq)) {
    f <- freq[fi]

    # MA polynomial: 1 - theta_1*z - theta_2*z^2 - ...
    num <- 1 + 0i
    if (q > 0) {
      for (k in seq_len(q)) {
        num <- num - theta[k] * exp(-2i * pi * k * f)
      }
    }

    # AR polynomial: 1 - phi_1*z - phi_2*z^2 - ...
    den <- 1 + 0i
    if (p > 0) {
      for (k in seq_len(p)) {
        den <- den - phi[k] * exp(-2i * pi * k * f)
      }
    }

    # Spectral density in dB
    spec[fi] <- 10 * log10((vara * Mod(num)^2) / (gvar * Mod(den)^2))
  }

  # Plot if requested
  if (plot) {
    plot(freq, spec, type = "l",
         xlab = "Frequency", ylab = "dB",
         main = "True Spectral Density")
    invisible(list(freq = freq, spec = spec))
  } else {
    list(freq = freq, spec = spec)
  }
}
