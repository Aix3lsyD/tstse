#' True FARMA Autocorrelations
#'
#' Compute and optionally plot the theoretical autocorrelations for a
#' Fractional ARMA (FARMA) model.
#'
#' @param d Numeric. Fractional differencing parameter. Typically in (-0.5, 0.5).
#'   Positive values indicate long memory.
#' @param phi Numeric vector. AR coefficients (default 0 = none).
#' @param theta Numeric vector. MA coefficients (default 0 = none).
#' @param lag_max Integer. Maximum lag to compute (default 50).
#' @param trunc Integer. Truncation point for infinite sum approximation
#'   (default 1000). Larger values give more accurate results.
#' @param vara Numeric. White noise variance (default 1).
#' @param plot Logical. If `TRUE` (default), plot the ACF.
#'
#' @return A list with components:
#'   \item{acf}{Theoretical autocorrelations (lags 0 to lag_max)}
#'   \item{acv}{Theoretical autocovariances (lags 0 to lag_max)}
#'
#' @details
#' The FARMA(p,d,q) model is defined as:
#' \deqn{(1-B)^d \phi(B) X_t = \theta(B) a_t}
#'
#' where `d` is the fractional differencing parameter. This is a special case
#' of the GARMA model with `u = 1`.
#'
#' The theoretical ACF is computed by convolving the ARMA autocovariances
#' with Gegenbauer polynomial coefficients for the fractional differencing
#' operator.
#'
#' Long memory behavior:
#' \itemize{
#'   \item `d > 0`: Long memory (positive dependence)
#'   \item `d < 0`: Antipersistence (negative dependence)
#'   \item `d = 0`: Reduces to standard ARMA
#' }
#'
#' @seealso [true_acf()] for standard ARMA autocorrelations,
#'   [true_garma_acf()] for the general GARMA case.
#'
#' @references
#' Hosking, J. R. M. (1981). Fractional differencing.
#' *Biometrika*, 68(1), 165-176.
#'
#' @examples
#' # Pure fractional difference (FARMA(0,d,0))
#' true_farma_acf(d = 0.3)
#'
#' # FARMA(1,d,0) - AR(1) with long memory
#' true_farma_acf(d = 0.2, phi = 0.5)
#'
#' # FARMA(0,d,1) - MA(1) with long memory
#' true_farma_acf(d = 0.25, theta = 0.4)
#'
#' # When d=0, reduces to standard ARMA
#' true_farma_acf(d = 0, phi = 0.8, lag_max = 20)
#'
#' # Without plot
#' result <- true_farma_acf(d = 0.4, plot = FALSE)
#'
#' @export
true_farma_acf <- function(d = 0, phi = 0, theta = 0, lag_max = 50L,
                           trunc = 1000L, vara = 1, plot = TRUE) {

  # Input validation
  if (!is.numeric(d) || length(d) != 1L) {
    stop("`d` must be a single numeric value.", call. = FALSE)
  }
  if (d <= -0.5 || d >= 0.5) {
    warning("`d` outside typical range (-0.5, 0.5); results may be unreliable.")
  }

  if (!is.numeric(phi)) {
    stop("`phi` must be a numeric vector.", call. = FALSE)
  }
  if (!is.numeric(theta)) {
    stop("`theta` must be a numeric vector.", call. = FALSE)
  }

  if (!is.numeric(lag_max) || length(lag_max) != 1L || lag_max < 1L) {
    stop("`lag_max` must be a positive integer.", call. = FALSE)
  }
  lag_max <- as.integer(lag_max)

  if (!is.numeric(trunc) || length(trunc) != 1L || trunc < 1L) {
    stop("`trunc` must be a positive integer.", call. = FALSE)
  }
  trunc <- as.integer(trunc)

  if (!is.numeric(vara) || length(vara) != 1L || vara <= 0) {
    stop("`vara` must be a positive number.", call. = FALSE)
  }

  if (!is.logical(plot) || length(plot) != 1L || is.na(plot)) {
    stop("`plot` must be TRUE or FALSE.", call. = FALSE)
  }

  # Get ARMA autocovariances
  # Need lag_max + trunc lags for convolution
  arma_result <- true_acf(phi = phi, theta = theta,
                          lag_max = lag_max + trunc,
                          vara = vara, plot = FALSE)
  acv_V <- arma_result$acv

  # Compute Gegenbauer coefficients for fractional differencing
  # C_0 = gamma(1-2d) / (gamma(1-d))^2
  # C_i = C_{i-1} * (i-1+d) / (i-d)
  n_coef <- lag_max + trunc + 1L
  acv_X <- numeric(n_coef)
  acv_X[1] <- gamma(1 - 2 * d) / (gamma(1 - d))^2

  for (i in 2:n_coef) {
    acv_X[i] <- acv_X[i - 1] * (i - 2 + d) / (i - 1 - d)
  }

  # Convolve to get FARMA autocovariances
  acv_Y <- numeric(lag_max + 1L)

  for (k in seq_len(lag_max + 1L)) {
    xx <- numeric(trunc + 1L)
    for (j in seq_len(trunc + 1L)) {
      # Indices: acv_X uses 1-based indexing
      # k ranges from 1 to lag_max+1 (representing lags 0 to lag_max)
      # j ranges from 1 to trunc+1
      idx1 <- j + k - 1  # j-1 + k
      idx2 <- abs(k - j) + 1
      xx[j] <- acv_V[j] * (acv_X[idx1] + acv_X[idx2])
    }
    acv_Y[k] <- sum(xx) - acv_V[1] * acv_X[k]
  }

  # Normalize to get ACF
  acf_values <- acv_Y / acv_Y[1]

  # Plot if requested
  if (plot) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))

    par(mar = c(4, 4, 2, 1))
    k <- 0:lag_max
    plot(k, acf_values, type = "h", ylim = c(-1, 1),
         xlab = "Lag", ylab = "Autocorrelation",
         main = paste0("True FARMA ACF (d=", d, ")"))
    abline(h = 0)
  }

  list(acf = acf_values, acv = acv_Y)
}
