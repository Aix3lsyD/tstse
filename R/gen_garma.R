#' Generate GARMA Process Realization
#'
#' Generate a realization from a k-factor Generalized ARMA (GARMA) model,
#' which combines Gegenbauer long-memory factors with ARMA components.
#'
#' @param n Integer. Length of the realization to generate.
#' @param u Numeric vector. Gegenbauer frequency parameters (length 1 or 2).
#'   Each value should be in \[-1, 1\], determining the spectral peak frequency
#'   as \eqn{f = \arccos(u) / (2\pi)}. Special cases:
#'   \itemize{
#'     \item `u = 1`: Long memory at frequency 0 (like ARFIMA)
#'     \item `u = 0`: Long memory at frequency 0.25
#'     \item `u = -1`: Long memory at frequency 0.5 (Nyquist)
#'   }
#' @param lambda Numeric vector. Gegenbauer persistence parameters (same length
#'   as `u`). Values in (0, 0.5) give stationary long memory.
#' @param phi Numeric vector. AR coefficients (default 0 = no AR component).
#' @param theta Numeric vector. MA coefficients (default 0 = no MA component).
#'   Uses negated sign convention (matches textbook notation).
#' @param n_trunc Integer. Truncation point for the infinite GLP approximation.
#'   Default is 1000.
#' @param burn_in Integer. Burn-in period to remove transient effects.
#'   Default is 600.
#' @param var_a Numeric. Variance of the white noise innovations. Default is 1.
#' @param plot Logical. If `TRUE` (default), plot the realization.
#' @param seed Integer or NULL. Random seed for reproducibility.
#'
#' @return Numeric vector of length `n` containing the GARMA realization.
#'   Returned invisibly if `plot = TRUE`.
#'
#' @details
#' A k-factor GARMA model has the form:
#' \deqn{\prod_{j=1}^{k}(1 - 2u_jB + B^2)^{-\lambda_j} \phi(B) X_t = \theta(B) a_t}
#'
#' where:
#' \itemize{
#'   \item \eqn{(1 - 2u_jB + B^2)^{-\lambda_j}} are Gegenbauer factors
#'   \item \eqn{\phi(B)} is the AR polynomial
#'   \item \eqn{\theta(B)} is the MA polynomial
#'   \item \eqn{a_t} is white noise
#' }
#'
#' The function generates the process using the General Linear Process (GLP)
#' representation. First, a Gegenbauer process is generated via [gen_geg()],
#' then filtered through the AR/MA components.
#'
#' @section Stationarity:
#' For stationarity with a single Gegenbauer factor:
#' \itemize{
#'   \item If \eqn{|u| < 1}: need \eqn{\lambda < 0.5}
#'   \item If \eqn{u = \pm 1}: need \eqn{\lambda < 0.25}
#' }
#'
#' @seealso [gen_geg()] for pure Gegenbauer process,
#'   [est_garma()] for parameter estimation,
#'   [fore_garma()] for forecasting,
#'   [gen_arma()] for ARMA process generation.
#'
#' @references
#' Ferrara, L., & Guegan, D. (2001). Forecasting with k-factor Gegenbauer
#' processes: Theory and applications. *Journal of Forecasting*, 20(8), 581-601.
#'
#' @examples
#' # Pure Gegenbauer process (no AR/MA)
#' x <- gen_garma(n = 200, u = 0.8, lambda = 0.4, seed = 123)
#'
#' # GARMA with AR(1) component
#' x <- gen_garma(n = 200, u = 0.5, lambda = 0.3, phi = 0.7, seed = 123)
#'
#' # GARMA with AR(2) and MA(1) components
#' x <- gen_garma(n = 200, u = 0.8, lambda = 0.4,
#'                phi = c(1.2, -0.5), theta = 0.3, seed = 123)
#'
#' # Two-factor GARMA (long memory at two frequencies)
#' x <- gen_garma(n = 200, u = c(0.8, 0.3), lambda = c(0.3, 0.2), seed = 123)
#'
#' # Without plot
#' x <- gen_garma(n = 100, u = 0.5, lambda = 0.3, plot = FALSE, seed = 42)
#'
#' @export
gen_garma <- function(n, u, lambda, phi = 0, theta = 0, n_trunc = 1000L,
                      burn_in = 600L, var_a = 1, plot = TRUE, seed = NULL) {

  # Input validation
  if (!is.numeric(n) || length(n) != 1L || n < 1L) {
    stop("`n` must be a positive integer.", call. = FALSE)
  }
  n <- as.integer(n)

  if (!is.numeric(u)) {
    stop("`u` must be numeric.", call. = FALSE)
  }
  if (!is.numeric(lambda)) {
    stop("`lambda` must be numeric.", call. = FALSE)
  }
  if (length(u) != length(lambda)) {
    stop("`u` and `lambda` must have the same length.", call. = FALSE)
  }
  k <- length(u)
  if (k > 2L) {
    stop("Only 1 or 2 Gegenbauer factors are supported.", call. = FALSE)
  }
  if (any(abs(u) > 1)) {
    stop("`u` values must be in [-1, 1].", call. = FALSE)
  }

  if (!is.numeric(phi)) {
    stop("`phi` must be numeric.", call. = FALSE)
  }
  if (!is.numeric(theta)) {
    stop("`theta` must be numeric.", call. = FALSE)
  }

  if (!is.numeric(n_trunc) || length(n_trunc) != 1L || n_trunc < 1L) {
    stop("`n_trunc` must be a positive integer.", call. = FALSE)
  }
  n_trunc <- as.integer(n_trunc)

  if (!is.numeric(burn_in) || length(burn_in) != 1L || burn_in < 0L) {
    stop("`burn_in` must be a non-negative integer.", call. = FALSE)
  }
  burn_in <- as.integer(burn_in)

  if (!is.numeric(var_a) || length(var_a) != 1L || var_a <= 0) {
    stop("`var_a` must be a positive number.", call. = FALSE)
  }

  if (!is.logical(plot) || length(plot) != 1L || is.na(plot)) {
    stop("`plot` must be TRUE or FALSE.", call. = FALSE)
  }

  # Determine effective AR/MA orders
  # If coefficients are all zero, treat as no AR/MA
  p <- if (all(phi == 0)) 0L else length(phi)
  q <- if (all(theta == 0)) 0L else length(theta)
  r <- max(p, q)

  # Generate base Gegenbauer process
  # Need extra observations for burn-in and AR/MA ramp-up
  # Note: seed is passed to gen_geg (not set here) for consistency with original
  total_n <- n + burn_in + r
  x <- gen_geg(n = total_n, u = u, d = lambda, n_trunc = n_trunc,
               var_a = var_a, seed = seed)

  # Apply AR/MA filter if needed
  if (p == 0L && q == 0L) {
    y <- x
  } else {
    y <- numeric(total_n)

    # MA coefficients: c(1, -theta) for negated sign convention
    ma_coefs <- if (q == 0L) 1 else c(1, -theta)

    for (i in (r + 1):total_n) {
      # AR part: phi[1]*y[i-1] + phi[2]*y[i-2] + ...
      ar_part <- if (p > 0L) {
        sum(phi * y[(i - 1):(i - p)])
      } else {
        0
      }

      # MA part: x[i] - theta[1]*x[i-1] - theta[2]*x[i-2] - ...
      ma_part <- if (q > 0L) {
        sum(ma_coefs * x[i:(i - q)])
      } else {
        x[i]
      }

      y[i] <- ar_part + ma_part
    }
  }

  # Extract the final n observations (discard burn-in and ramp-up)
  start_idx <- 1L + burn_in + r
  end_idx <- n + burn_in + r
  realization <- y[start_idx:end_idx]

  # Plot if requested
  if (plot) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))

    par(mar = c(4, 4, 2, 1))
    t <- seq_len(n)
    plot(t, realization, type = "o", cex = 0.5, pch = 16, lwd = 0.75,
         xlab = "Time", ylab = "",
         main = paste0("GARMA Realization (k=", k, ")"))
  }

  if (plot) invisible(realization) else realization
}
