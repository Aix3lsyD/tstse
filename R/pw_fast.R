#' Fast Prais-Winsten Estimation Using C++ Primitives
#'
#' Optimized Prais-Winsten estimation for trend in time series with AR(1)
#' correlated errors. Uses C++ primitives for maximum performance.
#'
#' @param x A numeric vector containing the time series data.
#' @param iterate Logical. If `TRUE`, uses iterative Prais-Winsten that
#'   re-estimates rho until convergence. If `FALSE` (default), uses two-step
#'   method that estimates rho once from OLS residuals.
#' @param max_iter Integer. Maximum iterations for iterative method.
#'   Default is 50. Only used when `iterate = TRUE`.
#' @param tol Numeric. Convergence tolerance for rho in iterative method.
#'   Default is 1e-6. Iteration stops when |rho_new - rho_old| < tol.
#'
#' @return A list containing:
#'   \item{x}{The original time series.}
#'   \item{z_x}{Residuals from initial OLS trend fit.}
#'   \item{b0hat}{Prais-Winsten estimate of intercept.}
#'   \item{b1hat}{Prais-Winsten estimate of slope.}
#'   \item{rho}{Estimated AR(1) coefficient.}
#'   \item{pvalue}{P-value for test of H0: slope = 0.}
#'   \item{tpw}{t-statistic for the slope coefficient.}
#'
#' @details
#' This function implements Prais-Winsten estimation with two methods:
#'
#' **Two-step method** (`iterate = FALSE`, default):
#' \enumerate{
#'   \item Fit OLS trend and compute residuals.
#'   \item Estimate AR(1) coefficient (rho) from lag-1 autocorrelation.
#'   \item Apply Prais-Winsten quasi-differencing transformation.
#'   \item Regress transformed data on transformed time.
#' }
#'
#' **Iterative method** (`iterate = TRUE`):
#' \enumerate{
#'   \item Same as two-step steps 1-4.
#'   \item Compute residuals from PW regression.
#'   \item Re-estimate rho from PW residuals.
#'   \item Repeat steps 3-6 until rho converges or max_iter reached.
#' }
#'
#' The two-step method is faster and typically sufficient for most applications.
#' The iterative method may provide slightly better estimates when AR(1)
#' correlation is very high, but usually converges quickly (2-5 iterations).
#'
#' Unlike Cochrane-Orcutt which drops the first observation, Prais-Winsten
#' preserves it using the quasi-differencing transformation:
#' \deqn{W_1 = \sqrt{1-\rho^2} X_1}
#' \deqn{W_t = X_t - \rho X_{t-1} \quad \text{for } t > 1}
#'
#' This is more efficient than CO because no information is lost.
#'
#' @references
#' Prais, S. J. and Winsten, C. B. (1954). "Trend Estimators and Serial
#' Correlation." Cowles Commission Discussion Paper No. 383.
#'
#' @seealso [wbg_boot_fast()] for bootstrap-based trend test,
#'   [co_fast()] for Cochrane-Orcutt estimation
#'
#' @examples
#' # Generate series with trend and AR(1) errors
#' set.seed(123)
#' n <- 100
#' t <- 1:n
#' z <- arima.sim(list(ar = 0.7), n = n)
#' x <- 5 + 0.1 * t + z
#'
#' # Estimate trend with Prais-Winsten (two-step, default)
#' result <- pw_fast(x)
#' cat("Slope estimate:", result$b1hat, "\n")
#' cat("t-statistic:", result$tpw, "\n")
#' cat("AR(1) rho:", result$rho, "\n")
#'
#' # Iterative Prais-Winsten
#' result_iter <- pw_fast(x, iterate = TRUE)
#' cat("Iterative rho:", result_iter$rho, "\n")
#'
#' @export
pw_fast <- function(x, iterate = FALSE, max_iter = 50L, tol = 1e-6) {

  # Input validation
  if (!is.numeric(x) || length(x) == 0) {
    stop("`x` must be a non-empty numeric vector", call. = FALSE)
  }

  n <- length(x)
  if (n < 5) {
    stop("Series too short: need n >= 5", call. = FALSE)
  }

  # Get full results from C++ (tpw, rho, vara)
  # Use iterative or two-step method based on parameter
  if (iterate) {
    result <- pw_full_iterative_cpp(x, as.integer(max_iter), as.double(tol))
  } else {
    result <- pw_full_cpp(x)
  }

  # Also get OLS residuals for output
  z_x <- ols_detrend_cpp(x)

  # Compute coefficients via R for output
  # (C++ handles t-stat internally but we want explicit coefficients)
  rho <- result$rho
  sqrt_term <- sqrt(1 - rho^2)

  # PW transformation (intercept column must also be transformed)
  x_pw <- c(sqrt_term * x[1], x[-1] - rho * x[-n])
  intercept_pw <- c(sqrt_term, rep(1 - rho, n - 1))
  t_pw <- c(sqrt_term, (2:n) * (1 - rho) + rho)

  # Fast OLS via lm.fit on transformed design matrix
  fit <- stats::lm.fit(cbind(intercept_pw, t_pw), x_pw)
  b0hat <- unname(fit$coefficients[1])
  b1hat <- unname(fit$coefficients[2])

  # P-value from t-distribution
  df_resid <- n - 2
  pvalue <- 2 * stats::pt(abs(result$tpw), df_resid, lower.tail = FALSE)

  list(
    x       = x,
    z_x     = as.numeric(z_x),
    b0hat   = b0hat,
    b1hat   = b1hat,
    rho     = result$rho,
    pvalue  = pvalue,
    tpw     = result$tpw
  )
}
