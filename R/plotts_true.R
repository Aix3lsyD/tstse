#' Plot True ARMA Properties
#'
#' Generate an ARMA realization and plot with theoretical ACF and spectrum.
#'
#' @param n Integer, series length (default 100).
#' @param phi Numeric vector, AR coefficients (default 0 = none).
#' @param theta Numeric vector, MA coefficients (default 0 = none).
#' @param lag_max Integer, maximum lag for ACF (default 25).
#' @param mu Numeric, process mean (default 0).
#' @param vara Numeric, white noise variance (default 1).
#' @param seed Integer or NULL, random seed for reproducibility.
#' @param plot Logical, whether to display plots (default TRUE).
#'
#' @return Invisibly, a list with components:
#'   \item{data}{Generated realization (length n)}
#'   \item{acf}{True autocorrelations (lags 0 to lag_max)}
#'   \item{acv}{True autocovariances (lags 0 to lag_max)}
#'   \item{freq}{Frequency values (0 to 0.5)}
#'   \item{spec}{Spectral density in dB}
#' @export
#'
#' @details
#' This function combines data generation with theoretical property visualization.
#' It generates an ARMA(p,q) realization using \code{\link{gen_arma}}, computes
#' theoretical autocorrelations using \code{\link{true_acf}}, and computes
#' theoretical spectral density using \code{\link{true_spec}}.
#'
#' When \code{plot = TRUE}, displays a 3-panel layout:
#' \itemize{
#'   \item Top: Time series realization
#'   \item Bottom-left: True autocorrelations
#'   \item Bottom-right: True spectral density (dB)
#' }
#'
#' @seealso \code{\link{gen_arma}}, \code{\link{true_acf}}, \code{\link{true_spec}}
#'
#' @examples
#' # AR(2) process
#' result <- plotts_true(n = 200, phi = c(1.6, -0.9), seed = 123)
#'
#' # ARMA(1,1) without plotting
#' result <- plotts_true(n = 100, phi = 0.7, theta = 0.4, plot = FALSE, seed = 456)
#' names(result)
#'
#' # MA(1) process
#' plotts_true(n = 150, theta = 0.8, seed = 789)
plotts_true <- function(n = 100L, phi = 0, theta = 0, lag_max = 25L,
                        mu = 0, vara = 1, seed = NULL, plot = TRUE) {

  # Input validation
  if (!is.numeric(n) || length(n) != 1 || n < 1) {
    stop("`n` must be a positive integer")
  }
  if (!is.numeric(phi)) stop("`phi` must be numeric")
  if (!is.numeric(theta)) stop("`theta` must be numeric")
  if (!is.numeric(lag_max) || length(lag_max) != 1 || lag_max < 0) {
    stop("`lag_max` must be a non-negative integer")
  }
  if (!is.numeric(mu) || length(mu) != 1) {
    stop("`mu` must be a single numeric value")
  }
  if (!is.numeric(vara) || length(vara) != 1 || vara <= 0) {
    stop("`vara` must be positive")
  }

  n <- as.integer(n)
  lag_max <- as.integer(lag_max)

  # Cap lag_max at n-1
  if (lag_max >= n) lag_max <- n - 1L

  # Generate realization (always)
  data <- gen_arma(n = n, phi = phi, theta = theta, mu = mu,
                   vara = vara, plot = FALSE, seed = seed)

  # Compute true ACF/ACV
  acf_result <- true_acf(phi = phi, theta = theta, lag_max = lag_max,
                         vara = vara, plot = FALSE)

  # Compute true spectrum
  spec_result <- true_spec(phi = phi, theta = theta, vara = vara,
                           n_freq = 251L, plot = FALSE)

  # Plot if requested
  if (plot) {
    # Save and restore par settings
    op <- par(no.readonly = TRUE)
    on.exit(par(op))

    # Layout: centered realization on top, ACF+spectrum below
    layout(matrix(c(0, 1, 1, 0,
                    2, 2, 3, 3), nrow = 2, byrow = TRUE))
    par(mar = c(5, 2.8, 1, 1.5))

    # Panel 1: Realization
    plot(seq_len(n), data, type = "l",
         xlab = "Time", ylab = "",
         main = "Realization")

    # Panel 2: True ACF
    k <- 0:lag_max
    plot(k, acf_result$acf, type = "h", ylim = c(-1, 1),
         xlab = "Lag", ylab = "Autocorrelations",
         main = "True Autocorrelations")
    abline(h = 0)

    # Panel 3: True Spectrum
    plot(spec_result$freq, spec_result$spec, type = "l",
         xlab = "Frequency", ylab = "dB",
         main = "True Spectral Density")
  }

  # Return result (always includes data)
  invisible(list(
    data = data,
    acf = acf_result$acf,
    acv = acf_result$acv,
    freq = spec_result$freq,
    spec = spec_result$spec
  ))
}
