#' Fast Jacob Turner Trend Test with VIF-Adjusted Degrees of Freedom
#'
#' Optimized version of [jt()] using C++ primitives and vectorized operations.
#' Produces identical results to [jt()] but with significantly better performance
#' for large-scale simulations.
#'
#' @param x A numeric vector containing the time series data.
#' @param maxp Maximum AR order for model selection. Default is 5.
#' @param type Character. Information criterion for order selection:
#'   `"aic"`, `"aicc"`, or `"bic"`. Default is `"aic"`.
#'
#' @return A list containing the same elements as [jt()]:
#'   \item{x}{The original time series.}
#'   \item{z_x}{Transformed series (cumsum of demeaned differences).}
#'   \item{b0hat}{Estimate of intercept.}
#'   \item{b1hat}{Estimate of slope.}
#'   \item{z_order}{AR order selected for z_x (used for CO transform).}
#'   \item{z_phi}{AR coefficients for z_x.}
#'   \item{x_order}{AR order selected for original x (used for VIF).}
#'   \item{x_phi}{AR coefficients for original x.}
#'   \item{vif}{Variance inflation factor.}
#'   \item{rel_n}{Effective sample size: (n-2)/vif.}
#'   \item{est_df}{Adjusted degrees of freedom.}
#'   \item{pvalue}{P-value using VIF-adjusted degrees of freedom.}
#'   \item{tco}{t-statistic for the slope coefficient.}
#'
#' @details
#' This function is functionally equivalent to [jt()] but uses:
#' \itemize{
#'   \item C++ Burg algorithm via `burg_aic_select_cpp()` instead of `aic_burg()`
#'   \item C++ AR transformation via `ar_transform_cpp()` instead of `artrans()`
#'   \item Vectorized VIF computation instead of R loop
#'   \item `lm.fit()` instead of `lm()` + `summary()` for regression
#' }
#'
#' The key performance improvement comes from `ar_transform_cpp()` which
#' performs the AR filtering without computing ACFs (unlike `artrans()` which
#' computes ACFs even when `plot = FALSE`).
#'
#' @seealso [jt()] for the original implementation,
#'   [co()] for standard Cochrane-Orcutt,
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
#' # Compare jt() and jt_fast()
#' r1 <- jt(x, maxp = 5)
#' r2 <- jt_fast(x, maxp = 5)
#'
#' # Results should be identical
#' all.equal(r1$tco, r2$tco)
#' all.equal(r1$pvalue, r2$pvalue)
#'
#' @export
jt_fast <- function(x, maxp = 5L, type = c("aic", "aicc", "bic")) {

  type <- match.arg(type)

  # Input validation
  if (!is.numeric(x) || length(x) == 0) {
    stop("`x` must be a non-empty numeric vector", call. = FALSE)
  }
  if (!is.numeric(maxp) || length(maxp) != 1 || maxp < 1) {
    stop("`maxp` must be a positive integer", call. = FALSE)
  }

  n <- length(x)
  maxp <- as.integer(maxp)

  if (n <= maxp + 10) {
    stop("Series too short: need n > maxp + 10", call. = FALSE)
  }

  # Step 1: AR fit on original series
  # Use C++ Burg with min_p=1 to match R's aic_burg(p=1:maxp) behavior
  aic_x <- burg_aic_select_cpp(x, maxp, type, min_p = 1L)

  # Step 2: Transform series - difference, demean, cumsum
  x_diff <- diff(x)
  z_x_diff <- x_diff - mean(x_diff)
  z_x <- cumsum(c(x[1], z_x_diff))

  # Step 3: AR fit on z_x
  # Use C++ Burg with min_p=1 to match R's aic_burg(p=1:maxp) behavior
  aic_z <- burg_aic_select_cpp(z_x, maxp, type, min_p = 1L)
  pp <- aic_z$p
  phi <- as.numeric(aic_z$phi)

  # Step 4: AR transform (USE C++ - skips ACF entirely!)
  if (pp > 0) {
    x_trans <- ar_transform_cpp(x, phi)
  } else {
    x_trans <- x
  }

  # Step 5: CO time transform (VECTORIZED)
  # Original loop: t_co[tt] = tt - sum(phi * (tt - seq_len(pp)))
  # Simplifies to: t_co = tt * (1 - sum(phi)) + sum(phi * seq_len(pp))
  if (pp > 0) {
    phi_sum_weighted <- sum(phi * seq_len(pp))
    phi_sum <- sum(phi)
    tt <- (pp + 1):n
    t_co <- tt * (1 - phi_sum) + phi_sum_weighted
  } else {
    t_co <- seq_len(n)
  }

  # Step 6: Fast OLS regression using lm.fit()
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

  # Step 7: VIF computation (VECTORIZED)
  x_p <- aic_x$p
  x_phi <- if (x_p > 0) as.numeric(aic_x$phi) else numeric(0)

  if (x_p > 0) {
    acf_vals <- as.numeric(stats::ARMAacf(ar = x_phi, lag.max = n - 1))[-1]
  } else {
    acf_vals <- rep(0, n - 1)
  }

  # VIF = 1 + 2 * sum((1 - i/(n-2)) * rho(i)) for i = 1 to n-3
  # Vectorized version
  max_lag <- n - 3
  i <- seq_len(max_lag)
  vif <- 1 + 2 * sum((1 - i / (n - 2)) * acf_vals[i])
  vif <- max(vif, 1.0)  # Ensure VIF >= 1

  # Step 8: Adjusted df and p-value
  rel_n <- (n - 2) / vif
  est_df <- if (rel_n - pp > 0) rel_n - pp else rel_n
  est_df <- max(est_df, 1.0)  # Ensure df >= 1

  pvalue <- 2 * stats::pt(abs(tco), est_df, lower.tail = FALSE)

  # Return results (same structure as jt())
  list(
    x         = x,
    z_x       = z_x,
    b0hat     = b0hat,
    b1hat     = b1hat,
    z_order   = pp,
    z_phi     = phi,
    x_order   = x_p,
    x_phi     = x_phi,
    vif       = vif,
    rel_n     = rel_n,
    est_df    = est_df,
    pvalue    = pvalue,
    tco       = tco
  )
}
