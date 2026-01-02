#' Estimate ARMA Model
#'
#' Estimate ARMA(p,q) model parameters using maximum likelihood.
#'
#' @param x Numeric vector, the time series.
#' @param p Integer, AR order (default 0).
#' @param q Integer, MA order (default 0).
#' @param factor Logical, if TRUE (default) print factor table for p > 0.
#'
#' @return A list with components:
#'   \item{phi}{AR coefficient estimates}
#'   \item{theta}{MA coefficient estimates}
#'   \item{res}{Residuals from backcasting}
#'   \item{avar}{Residual variance estimate}
#'   \item{xbar}{Sample mean of original series}
#'   \item{aic}{AIC value}
#'   \item{aicc}{Corrected AIC value}
#'   \item{bic}{BIC value}
#'   \item{se_phi}{Standard errors for AR coefficients}
#'   \item{se_theta}{Standard errors for MA coefficients}
#' @export
#'
#' @details
#' Uses \code{\link[stats]{arima}} for MLE estimation, then computes residuals
#' using the backcasting procedure. The MA coefficients are returned with
#' signs matching the ATSA textbook convention.
#'
#' @examples
#' # Simulate and estimate AR(2)
#' x <- gen_arma(n = 200, phi = c(1.5, -0.75), plot = FALSE, seed = 123)
#' est <- est_arma(x, p = 2)
#'
#' # Estimate ARMA(1,1)
#' x <- gen_arma(n = 200, phi = 0.7, theta = 0.4, plot = FALSE, seed = 456)
#' est <- est_arma(x, p = 1, q = 1)
est_arma <- function(x, p = 0L, q = 0L, factor = TRUE) {

  # Input validation
  if (!is.numeric(x) || length(x) == 0) {
    stop("`x` must be a non-empty numeric vector")
  }
  if (!is.numeric(p) || p < 0) {
    stop("`p` must be a non-negative integer")
  }
  if (!is.numeric(q) || q < 0) {
    stop("`q` must be a non-negative integer")
  }
  if (p == 0 && q == 0) {
    stop("At least one of `p` or `q` must be positive")
  }

  n <- length(x)
  xbar <- mean(x)
  x_centered <- x - xbar

  # Fit ARMA model (suppress convergence warnings)
  fit <- suppressWarnings(
    arima(x_centered, order = c(p, 0, q), include.mean = FALSE)
  )

  # Extract coefficients and standard errors
  coefs <- as.vector(coef(fit))
  se <- sqrt(as.vector(diag(vcov(fit))))

  # AR coefficients
  if (p == 0) {
    phi <- 0
    se_phi <- 0
  } else {
    phi <- coefs[1:p]
    se_phi <- se[1:p]
  }

  # MA coefficients (negate to match ATSA convention)
  if (q == 0) {
    theta <- 0
    se_theta <- 0
  } else {
    theta <- -coefs[(p + 1):(p + q)]
    se_theta <- se[(p + 1):(p + q)]
  }

  # Print factor table if requested
  if (factor && p > 0) {
    factor(phi = phi, theta = theta)
  }

  # Compute residuals via backcasting
  res <- backcast(x_centered, phi = phi, theta = theta, n_back = 50)

  # Residual variance
  avar <- mean(res^2)

  # Information criteria
  aic <- log(avar) + 2 * (p + q + 1) / n
  bic <- log(avar) + (p + q + 1) * log(n) / n
  aicc <- log(avar) + (n + p + q + 1) / (n - p - q - 3)

  structure(
    list(
      phi = phi,
      theta = theta,
      res = res,
      avar = avar,
      xbar = xbar,
      aic = aic,
      aicc = aicc,
      bic = bic,
      se_phi = se_phi,
      se_theta = se_theta
    ),
    class = "est_arma"
  )
}


#' @export
print.est_arma <- function(x, ...) {
  cat("\nARMA Model Estimation\n\n")

  if (!all(x$phi == 0)) {
    cat("AR Coefficients:\n")
    print(round(x$phi, 4))
    cat("SE:", round(x$se_phi, 4), "\n\n")
  }

  if (!all(x$theta == 0)) {
    cat("MA Coefficients:\n")
    print(round(x$theta, 4))
    cat("SE:", round(x$se_theta, 4), "\n\n")
  }

  cat("Mean:", round(x$xbar, 4), "\n")
  cat("Residual Variance:", round(x$avar, 4), "\n\n")
  cat("AIC: ", round(x$aic, 4), "\n")
  cat("AICC:", round(x$aicc, 4), "\n")
  cat("BIC: ", round(x$bic, 4), "\n")

  invisible(x)
}
