#' True GARMA Autocorrelations
#'
#' Compute theoretical autocorrelations for a GARMA (Gegenbauer ARMA) process.
#'
#' @param u Numeric, Gegenbauer frequency parameter where u = cos(2*pi*f0).
#'   Must be in range (-1, 1). Common values: u = 1 for standard long memory,
#'   u = 0 for seasonal long memory at frequency 0.25.
#' @param lambda Numeric, Gegenbauer (long memory) parameter. Should be > 0.
#' @param phi Numeric vector, AR coefficients (default 0 = none).
#' @param theta Numeric vector, MA coefficients (default 0 = none).
#' @param lag_max Integer, maximum lag (default 50).
#' @param vara Numeric, white noise variance (default 1).
#' @param plot Logical, whether to plot ACF (default TRUE).
#'
#' @return A list with components:
#'   \item{acf}{Theoretical autocorrelations (lags 0 to lag_max)}
#'   \item{acv}{Theoretical autocovariances (lags 0 to lag_max)}
#'
#' @details
#' A GARMA(p, u, lambda, q) process has spectral density that includes a
#' Gegenbauer component causing long memory behavior centered at frequency f0.
#'
#' The autocovariances are computed by numerical integration of the spectral
#' density using adaptive quadrature. The integration is split around the
#' spectral pole at f = acos(u)/(2*pi) for numerical stability.
#'
#' Special cases:
#' \describe{
#'   \item{u near 1}{Standard fractional (ARFIMA-like) long memory at frequency 0}
#'   \item{u = 0}{Seasonal long memory at frequency 0.25}
#'   \item{u near -1}{Long memory at frequency 0.5}
#' }
#'
#' @seealso \code{\link{true_acf}} for standard ARMA autocorrelations
#'
#' @examples
#' \dontrun{
#' # Pure Gegenbauer process
#' result <- true_garma_acf(u = 0.5, lambda = 0.3, lag_max = 30)
#'
#' # GARMA with AR component
#' result <- true_garma_acf(u = 0, lambda = 0.4, phi = 0.5, lag_max = 40)
#' }
#'
#' @export
true_garma_acf <- function(u, lambda, phi = 0, theta = 0,
                           lag_max = 50L, vara = 1, plot = TRUE) {

  # Check MASS package
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required.\n",
         "Install with: install.packages('MASS')")
  }

  # Input validation
  if (!is.numeric(u) || length(u) != 1 || u <= -1 || u >= 1) {
    stop("`u` must be a single numeric value in (-1, 1)")
  }
  if (!is.numeric(lambda) || length(lambda) != 1 || lambda <= 0) {
    stop("`lambda` must be a positive number")
  }
  if (!is.numeric(phi)) stop("`phi` must be numeric")
  if (!is.numeric(theta)) stop("`theta` must be numeric")
  if (!is.numeric(lag_max) || lag_max < 0) {
    stop("`lag_max` must be a non-negative integer")
  }
  if (!is.numeric(vara) || vara <= 0) {
    stop("`vara` must be positive")
  }

  lag_max <- as.integer(lag_max)

  # Effective orders
  p <- if (all(phi == 0)) 0L else length(phi)
  q <- if (all(theta == 0)) 0L else length(theta)

  # Precompute pole location (spectral singularity)
  f_pole <- acos(u) / (2 * pi)

  # Compute autocovariances by integrating spectral density
  acv <- numeric(lag_max + 1L)

  for (k in seq_len(lag_max + 1L)) {
    lag <- k - 1L

    # Define integrand: S(f) * cos(2*pi*f*lag)
    integrand <- function(f) {
      # Handle singularity at f = acos(u)/(2*pi)
      if (f == f_pole) {
        return(0)
      }

      # MA component: |1 - theta_1*z - theta_2*z^2 - ...|^2 at z = e^{2*pi*i*f}
      if (q == 0L) {
        ma_comp <- 1
      } else {
        z_powers <- exp(2i * pi * f * (0:q))
        ma_comp <- Mod(sum(c(1, -theta) * z_powers))^2
      }

      # AR component: |1 - phi_1*z - phi_2*z^2 - ...|^2 at z = e^{2*pi*i*f}
      if (p == 0L) {
        ar_comp <- 1
      } else {
        z_powers <- exp(2i * pi * f * (0:p))
        ar_comp <- Mod(sum(c(1, -phi) * z_powers))^2
      }

      # Gegenbauer component: [4(cos(2*pi*f) - u)^2]^{-lambda}
      geg_comp <- (4 * (cos(2 * pi * f) - u)^2)^(-lambda)

      # Spectral density times cosine
      vara * ma_comp / ar_comp * geg_comp * cos(2 * pi * f * lag)
    }

    # Use MASS::area for adaptive quadrature (same as tswge)
    # Suppress warnings from numerical difficulties near the pole
    acv[k] <- 2 * suppressWarnings(MASS::area(integrand, 0, 0.5))
  }

  # Compute ACF from ACV
  acf <- acv / acv[1]

  # Plot if requested
  if (plot) {
    op <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(op))

    k <- 0:lag_max
    graphics::par(mar = c(4, 4, 2, 1))
    plot(k, acf, type = "h", ylim = c(-1, 1),
         xlab = "Lag", ylab = "Autocorrelation",
         main = paste0("True GARMA ACF (u=", round(u, 2),
                       ", lambda=", round(lambda, 2), ")"))
    graphics::abline(h = 0)
  }

  if (plot) {
    invisible(list(acf = acf, acv = acv))
  } else {
    list(acf = acf, acv = acv)
  }
}
