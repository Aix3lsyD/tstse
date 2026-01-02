#' Estimate AR Model
#'
#' Estimate AR(p) model parameters using various methods.
#'
#' @param x Numeric vector, the time series.
#' @param p Integer, AR order (default 2).
#' @param factor Logical, if TRUE (default) print factor table.
#' @param method Character, estimation method: "mle" (default), "burg", or "yw".
#'
#' @return A list with components:
#'   \item{method}{Estimation method used}
#'   \item{phi}{AR coefficient estimates}
#'   \item{res}{Residuals from backcasting}
#'   \item{avar}{Residual variance estimate}
#'   \item{xbar}{Sample mean of original series}
#'   \item{aic}{AIC value}
#'   \item{aicc}{Corrected AIC value}
#'   \item{bic}{BIC value}
#' @export
#'
#' @details
#' Three estimation methods are available:
#' \describe{
#'   \item{mle}{Maximum likelihood via \code{\link{est_arma}}}
#'   \item{burg}{Burg's algorithm via \code{\link[stats]{ar.burg}}}
#'   \item{yw}{Yule-Walker via \code{\link[stats]{ar.yw}}}
#' }
#'
#' @examples
#' x <- gen_arma(n = 200, phi = c(1.5, -0.75), plot = FALSE, seed = 123)
#'
#' # MLE (default)
#' est_ar(x, p = 2)
#'
#' # Burg's method
#' est_ar(x, p = 2, method = "burg")
#'
#' # Yule-Walker
#' est_ar(x, p = 2, method = "yw")
est_ar <- function(x, p = 2L, factor = TRUE, method = c("mle", "burg", "yw")) {

  # Input validation
  if (!is.numeric(x) || length(x) == 0) {
    stop("`x` must be a non-empty numeric vector")
  }
  if (!is.numeric(p) || p < 1) {
    stop("`p` must be a positive integer")
  }

  method <- match.arg(method)
  n <- length(x)
  xbar <- mean(x)
  x_centered <- x - xbar

  # Estimate AR coefficients
  phi <- switch(method,
                mle = {
                  fit <- est_arma(x_centered, p = p, q = 0, factor = FALSE)
                  fit$phi
                },
                burg = {
                  fit <- ar.burg(x_centered, order.max = p, aic = FALSE)
                  as.vector(fit$ar)
                },
                yw = {
                  fit <- ar.yw(x_centered, order.max = p, aic = FALSE)
                  as.vector(fit$ar)
                }
  )

  # Print factor table if requested
  if (factor) {
    factor_ts(phi = phi)
  }

  # Compute residuals via backcasting
  res <- backcast(x_centered, phi = phi, theta = 0, n_back = 50)

  # Residual variance
  avar <- mean(res^2)

  # Information criteria
  aic <- log(avar) + 2 * (p + 1) / n
  bic <- log(avar) + (p + 1) * log(n) / n
  aicc <- log(avar) + (n + p + 1) / (n - p - 3)

  structure(
    list(
      method = method,
      phi = phi,
      res = res,
      avar = avar,
      xbar = xbar,
      aic = aic,
      aicc = aicc,
      bic = bic
    ),
    class = "est_ar"
  )
}


#' @export
print.est_ar <- function(x, ...) {
  cat("\nAR Model Estimation (method:", x$method, ")\n\n")
  cat("AR Coefficients:\n")
  print(round(x$phi, 4))
  cat("\nMean:", round(x$xbar, 4), "\n")
  cat("Residual Variance:", round(x$avar, 4), "\n\n")
  cat("AIC: ", round(x$aic, 4), "\n")
  cat("AICC:", round(x$aicc, 4), "\n")
  cat("BIC: ", round(x$bic, 4), "\n")
  invisible(x)
}
