#' Generate Gegenbauer Process
#'
#' Generates a realization of a Gegenbauer (long-memory) process using
#' the general linear process representation.
#'
#' @param n Integer. Length of the realization to generate.
#' @param u Numeric vector. Gegenbauer frequency parameters. Each value should
#'   be in \[-1, 1\], determining the spectral peak frequency as
#'   \eqn{f = \arccos(u) / (2\pi)}.
#' @param d Numeric vector. Gegenbauer persistence parameters (same length as
#'   `u`). Values in (0, 0.5) give stationary long memory. Also known as
#'   `lambda` in some references.
#' @param n_trunc Integer. Truncation point for the infinite MA representation.
#'   Default is 1000. Larger values give more accurate approximations but
#'   slower computation.
#' @param var_a Numeric. Variance of the white noise innovations. Default is 1.
#' @param seed Integer or NULL. Random seed for reproducibility. If NULL
#'   (default), no seed is set.
#'
#' @return Numeric vector of length `n` containing the generated realization.
#'
#' @details
#' A Gegenbauer process is a generalization of ARFIMA that allows long memory
#' at frequencies other than zero. The spectral density has a pole (peak) at
#' frequency \eqn{f = \arccos(u) / (2\pi)}.
#'
#' The process is generated using the general linear process (GLP) form:
#' \deqn{X_t = \sum_{j=0}^{\infty} \psi_j a_{t-j}}
#' where \eqn{\psi_j} are the MA (Gegenbauer) coefficients computed by
#' [macoef_geg()] and \eqn{a_t} is white noise.
#'
#' Special cases:
#' \itemize{
#'   \item `u = 1`: Equivalent to ARFIMA(0, d, 0) - long memory at frequency 0
#'   \item `u = -1`: Long memory at frequency 0.5 (Nyquist)
#'   \item `u = 0`: Long memory at frequency 0.25
#' }
#'
#' @section Stationarity:
#' For stationarity with a single Gegenbauer factor:
#' \itemize{
#'   \item If \eqn{|u| < 1}: need \eqn{d < 0.5}
#'   \item If \eqn{u = \pm 1}: need \eqn{d < 0.25}
#' }
#'
#' @seealso [macoef_geg()] for computing the MA coefficients,
#'   [gegenb()] for Gegenbauer polynomial coefficients,
#'   [gen_arma()] for ARMA process generation.
#'
#' @examples
#' # Gegenbauer process with spectral peak at frequency 0.1
#' # u = cos(2*pi*0.1) ≈ 0.809
#' set.seed(123)
#' x <- gen_geg(n = 200, u = 0.809, d = 0.4)
#' plot(x, type = "l", main = "Gegenbauer Process")
#'
#' # Check spectrum - should peak near f = 0.1
#' parzen(x)
#'
#' # ARFIMA-like process (u = 1 gives long memory at f = 0)
#' x2 <- gen_geg(n = 200, u = 1, d = 0.3)
#'
#' @export
gen_geg <- function(n, u, d, n_trunc = 1000L, var_a = 1, seed = NULL) {
  # Input validation
  if (!is.numeric(n) || length(n) != 1L || n < 1L) {
    stop("`n` must be a positive integer.", call. = FALSE)
  }
  n <- as.integer(n)

  if (!is.numeric(u)) {
    stop("`u` must be numeric.", call. = FALSE)
  }
  if (!is.numeric(d)) {
    stop("`d` must be numeric.", call. = FALSE)
  }
  if (length(u) != length(d)) {
    stop("`u` and `d` must have the same length.", call. = FALSE)
  }
  if (length(u) > 2L) {
    stop("Only 1 or 2 Gegenbauer factors are supported.", call. = FALSE)
  }

  # Check for reasonable d values
  if (any(d >= 0.5)) {
    warning("d >= 0.5 may produce non-stationary process.", call. = FALSE)
  }
  if (any(d < 0)) {
    warning("d < 0 gives anti-persistent (short memory) behavior.", call. = FALSE)
  }

  if (!is.numeric(n_trunc) || length(n_trunc) != 1L || n_trunc < 1L) {
    stop("`n_trunc` must be a positive integer.", call. = FALSE)
  }
  n_trunc <- as.integer(n_trunc)

  if (!is.numeric(var_a) || length(var_a) != 1L || var_a <= 0) {
    stop("`var_a` must be a positive number.", call. = FALSE)
  }

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L) {
      stop("`seed` must be a single integer or NULL.", call. = FALSE)
    }
    set.seed(as.integer(seed))
  }

  # Compute MA (psi) coefficients using macoef_geg
  # For single factor, just pass scalar u and d
  # macoef_geg determines k from length of u
  psi <- macoef_geg(u = u, d = d, n_coef = n_trunc)

  # Generate white noise
  # Need n_trunc + n values total for the convolution
  sd_a <- sqrt(var_a)
  white_noise <- rnorm(n_trunc + n, mean = 0, sd = sd_a)

  # Generate realization via convolution
  # X_t = sum_{j=0}^{n_trunc-1} psi[j+1] * a[t-j]
  # For t = 1, ..., n, we need a[t], a[t-1], ..., a[t-n_trunc+1]
  # With white_noise indexed from 1 to (n_trunc + n),
  # we use indices (n_trunc + t) down to (t + 1)
  realization <- numeric(n)
  for (t in seq_len(n)) {
    # Indices: from (n_trunc + t) down to (t + 1), length n_trunc
    idx <- (n_trunc + t):(t + 1)
    realization[t] <- sum(psi * white_noise[idx])
  }

  realization
}
