#' Fast Cochrane-Orcutt Estimation Using C++ Primitives
#'
#' Optimized version of [co()] using C++ primitives for AR fitting and
#' transformation. Produces identical results to `co(..., method = "burg")`
#' but with significantly better performance.
#'
#' @param x A numeric vector containing the time series data.
#' @param maxp Maximum AR order for model selection. Default is 5.
#' @param type Character. Information criterion for order selection:
#'   `"aic"`, `"aicc"`, or `"bic"`. Default is `"aic"`.
#' @param iterate Logical. If `TRUE`, iterates by re-estimating AR coefficients
#'   from Cochrane-Orcutt residuals until convergence. If `FALSE` (default),
#'   performs a single-pass estimation.
#' @param max_iter Integer. Maximum iterations for iterative method.
#'   Default is 50. Only used when `iterate = TRUE`.
#' @param tol Numeric. Convergence tolerance for phi in iterative method.
#'   Default is 1e-6. Iteration stops when max|phi_new - phi_old| < tol.
#'
#' @return A list containing the same elements as [co()]:
#'   \item{x}{The original time series.}
#'   \item{z_x}{Residuals from initial OLS trend fit.}
#'   \item{b0hat}{Cochrane-Orcutt estimate of intercept.}
#'   \item{b1hat}{Cochrane-Orcutt estimate of slope.}
#'   \item{z_order}{AR order selected for the residuals.}
#'   \item{z_phi}{AR coefficients for the residual model.}
#'   \item{pvalue}{P-value for test of H0: slope = 0.}
#'   \item{tco}{t-statistic for the slope coefficient.}
#'   \item{method}{Always "burg" (C++ implementation uses Burg algorithm).}
#'
#' @details
#' This function uses C++ implementations for maximum performance:
#' \itemize{
#'   \item AR fitting via `burg_aic_select_cpp()` (Burg algorithm with IC selection)
#'   \item OLS detrending via `ols_detrend_cpp()`
#'   \item AR transformation via `ar_transform_cpp()`
#' }
#'
#' **Iterative mode** (`iterate = TRUE`): Classical Cochrane-Orcutt iterates
#' by re-estimating AR coefficients from the CO residuals until convergence.
#' This can provide more accurate estimates when the initial AR fit is poor.
#' Single-pass mode (`iterate = FALSE`, default) is typically sufficient and
#' is used by [wbg_boot_fast()] for performance.
#'
#' **Important note on AR order selection**: `co_fast()` uses C++ Burg which
#' computes variance via the recursive formula `var(p) = var(0) * prod(1 - phii^2)`,
#' while `co()` uses `aic_ar()` which computes variance from backcast residuals.
#' These methods can select different AR orders on some data, particularly on
#' white noise or near-white-noise series. On data with true AR structure, both
#' methods typically select the same order.
#'
#' The only limitation is that `co_fast()` only supports the Burg method for
#' AR estimation. For MLE or Yule-Walker, use [co()].
#'
#' @seealso [co()] for the original implementation with multiple AR methods,
#'   [wbg_boot_fast()] for bootstrap-based trend test
#'
#' @examples
#' # Generate series with trend and AR(1) errors
#' set.seed(123)
#' n <- 100
#' t <- 1:n
#' z <- arima.sim(list(ar = 0.7), n = n)
#' x <- 5 + 0.1 * t + z
#'
#' # Compare co() and co_fast()
#' r1 <- co(x, maxp = 5, method = "burg")
#' r2 <- co_fast(x, maxp = 5)
#'
#' # Results should be identical
#' all.equal(r1$tco, r2$tco)
#' all.equal(r1$pvalue, r2$pvalue)
#'
#' @export
co_fast <- function(x, maxp = 5L, type = c("aic", "aicc", "bic"),
                    iterate = FALSE, max_iter = 50L, tol = 1e-6) {

  type <- match.arg(type)

  # Input validation (same as co())
  if (!is.numeric(x) || length(x) == 0) {
    stop("`x` must be a non-empty numeric vector", call. = FALSE)
  }
  if (!is.numeric(maxp) || length(maxp) != 1 || maxp < 1) {
    stop("`maxp` must be a positive integer", call. = FALSE)
  }

  n <- length(x)
  maxp <- as.integer(maxp)

  # Step 1: OLS detrend to get residuals (for z_x and initial AR fitting)
  z_x <- ols_detrend_cpp(x)

  # Step 2: Fit AR model to residuals using C++ Burg
  # Use min_p=1 to match R's aic_ar(p = 1:maxp) which excludes AR(0)
  ar_fit <- burg_aic_select_cpp(z_x, maxp, type, min_p = 1L)
  pp <- ar_fit$p
  phi <- as.numeric(ar_fit$phi)

  # Helper function to perform CO transformation and regression
  co_transform_regress <- function(phi_current) {
    pp_cur <- length(phi_current)

    if (pp_cur > 0) {
      # Use C++ AR transform (fast, no ACF computation)
      x_trans <- ar_transform_cpp(x, phi_current)

      # Vectorized CO time transform
      phi_sum_weighted <- sum(phi_current * seq_len(pp_cur))
      phi_sum <- sum(phi_current)
      tt <- (pp_cur + 1):n
      t_co <- tt * (1 - phi_sum) + phi_sum_weighted
    } else {
      x_trans <- x
      t_co <- seq_len(n)
    }

    # Fast OLS regression
    X <- cbind(1, t_co)
    fit <- stats::lm.fit(X, x_trans)

    list(
      x_trans = x_trans,
      t_co = t_co,
      fit = fit,
      residuals = fit$residuals
    )
  }

  # Perform initial CO transformation and regression
  co_result <- co_transform_regress(phi)

  # Iterative CO: re-estimate AR from CO residuals until convergence

  if (iterate && max_iter > 1L) {
    for (iter in seq_len(max_iter - 1L)) {
      phi_old <- phi

      # Re-fit AR to CO residuals
      ar_fit_new <- burg_aic_select_cpp(co_result$residuals, maxp, type, min_p = 1L)
      pp <- ar_fit_new$p
      phi <- as.numeric(ar_fit_new$phi)

      # Check convergence (handle different lengths by padding with zeros)
      len_old <- length(phi_old)
      len_new <- length(phi)
      max_len <- max(len_old, len_new)
      phi_old_pad <- c(phi_old, rep(0, max_len - len_old))
      phi_new_pad <- c(phi, rep(0, max_len - len_new))

      if (max(abs(phi_new_pad - phi_old_pad)) < tol) {
        break
      }

      # Re-compute CO transformation with new phi
      co_result <- co_transform_regress(phi)
    }
  }

  # Extract final results
  fit <- co_result$fit
  b0hat <- unname(fit$coefficients[1])
  b1hat <- unname(fit$coefficients[2])

  # Compute t-statistic manually
  rss <- sum(fit$residuals^2)
  df_resid <- length(co_result$x_trans) - 2
  mse <- rss / df_resid

  # Standard error of slope: sqrt(MSE * (X'X)^{-1}[2,2])
  X <- cbind(1, co_result$t_co)
  XtX_inv <- solve(crossprod(X))
  se_b1 <- sqrt(mse * XtX_inv[2, 2])

  tco <- if (se_b1 > 1e-15) b1hat / se_b1 else 0

  # Step 5: Compute p-value from t-distribution
  pvalue <- 2 * stats::pt(abs(tco), df_resid, lower.tail = FALSE)

  list(
    x       = x,
    z_x     = as.numeric(z_x),
    b0hat   = b0hat,
    b1hat   = b1hat,
    z_order = pp,
    z_phi   = phi,
    pvalue  = pvalue,
    tco     = tco,
    method  = "burg"
  )
}
