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
co_fast <- function(x, maxp = 5L, type = c("aic", "aicc", "bic")) {

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

  # Step 1: OLS detrend to get residuals (for z_x and AR fitting)
  z_x <- ols_detrend_cpp(x)

  # Step 2: Fit AR model to residuals using C++ Burg
  # Use min_p=1 to match R's aic_ar(p = 1:maxp) which excludes AR(0)
  ar_fit <- burg_aic_select_cpp(z_x, maxp, type, min_p = 1L)
  pp <- ar_fit$p
  phi <- as.numeric(ar_fit$phi)

  # Step 3: AR transform the original data
  if (pp > 0) {
    # Use C++ AR transform (fast, no ACF computation)
    x_trans <- ar_transform_cpp(x, phi)

    # Vectorized CO time transform
    # Original loop: t_co[tt] = tt - sum(phi * (tt - seq_len(pp)))
    # Simplifies to: t_co = tt * (1 - sum(phi)) + sum(phi * seq_len(pp))
    phi_sum_weighted <- sum(phi * seq_len(pp))
    phi_sum <- sum(phi)
    tt <- (pp + 1):n
    t_co <- tt * (1 - phi_sum) + phi_sum_weighted
  } else {
    x_trans <- x
    t_co <- seq_len(n)
  }

  # Step 4: Fast OLS regression for coefficients and t-stat
  X <- cbind(1, t_co)
  fit <- stats::lm.fit(X, x_trans)

  b0hat <- unname(fit$coefficients[1])
  b1hat <- unname(fit$coefficients[2])

  # Compute t-statistic manually
  rss <- sum(fit$residuals^2)
  df_resid <- length(x_trans) - 2
  mse <- rss / df_resid

  # Standard error of slope: sqrt(MSE * (X'X)^{-1}[2,2])
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
