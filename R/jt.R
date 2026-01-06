#' Jacob Turner Trend Test with VIF-Adjusted Degrees of Freedom
#'
#' Tests for linear trend in time series data using a modified Cochrane-Orcutt
#' procedure with variance inflation factor (VIF) adjustment for the degrees
#' of freedom. This provides a more accurate theoretical p-value that accounts
#' for autocorrelation-induced variance inflation.
#'
#' @param x A numeric vector containing the time series data.
#' @param maxp Maximum AR order for model selection. Default is 5.
#' @param type Character. Information criterion for order selection:
#'   `"aic"`, `"aicc"`, or `"bic"`. Default is `"aic"`.
#'
#' @return A list containing:
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
#' This function implements a modified Cochrane-Orcutt procedure that adjusts
#' the degrees of freedom for the t-test based on the variance inflation factor
#' (VIF) computed from the theoretical autocorrelation function.
#'
#' The procedure:
#' \enumerate{
#'   \item Fits AR(p) to the original series x (for VIF calculation)
#'   \item Transforms x: difference, demean, cumsum to get z_x
#'   \item Fits AR(p) to z_x (for Cochrane-Orcutt transform)
#'   \item Applies AR transform to x and regresses on transformed time
#'   \item Computes VIF from theoretical ACF of the original series AR model
#'   \item Adjusts degrees of freedom: df_adj = (n-2)/vif - p
#'   \item Computes p-value using t-distribution with adjusted df
#' }
#'
#' The VIF adjustment accounts for the fact that autocorrelation in the
#' residuals inflates the variance of the regression coefficient estimate,
#' making the standard t-test anti-conservative.
#'
#' @references
#' Turner, J. (unpublished). Modified Cochrane-Orcutt with degrees of freedom
#' adjustment.
#'
#' Cochrane, D. and Orcutt, G. H. (1949). "Application of Least Squares
#' Regression to Relationships Containing Auto-Correlated Error Terms."
#' *Journal of the American Statistical Association*, 44(245), 32-61.
#'
#' @seealso [co()] for standard Cochrane-Orcutt,
#'   [wbg_boot()] for bootstrap-based trend test
#'
#' @examples
#' # Generate series with trend and AR(1) errors
#' set.seed(123)
#' n <- 100
#' t <- 1:n
#' z <- arima.sim(list(ar = 0.7), n = n)
#' x <- 5 + 0.1 * t + z
#'
#' # Test for trend with JT adjustment
#' result <- jt(x, maxp = 5)
#' cat("Slope estimate:", result$b1hat, "\n")
#' cat("t-statistic:", result$tco, "\n")
#' cat("VIF:", result$vif, "\n")
#' cat("Effective df:", result$est_df, "\n")
#' cat("JT p-value:", result$pvalue, "\n")
#'
#' @export
jt <- function(x, maxp = 5L, type = c("aic", "aicc", "bic")) {

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

  # Step 1: AR fit on original series (for VIF calculation)
  aic_x <- aic_burg(x, p = seq_len(maxp), type = type)

  # Step 2: Transform series - difference, demean, cumsum
  # This creates z_x which is used for the CO procedure
  x_diff <- diff(x)  # First difference
  z_x_diff <- x_diff - mean(x_diff)  # Demean
  z_x <- cumsum(c(x[1], z_x_diff))  # Cumsum back

  # Step 3: AR fit on z_x (for CO transform)
  aic_z <- aic_burg(z_x, p = seq_len(maxp), type = type)
  pp <- aic_z$p
  phi <- aic_z$phi

  # Step 4: AR transform and CO regression
  x_trans <- artrans(x, phi = phi, plot = FALSE)

  # Compute transformed time index
  p1 <- pp + 1L
  t_co <- numeric(n)

  for (tt in p1:n) {
    t_co[tt] <- tt - sum(phi * (tt - seq_len(pp)))
  }

  # Regression on transformed data
  d_co <- stats::lm(x_trans ~ t_co[p1:n])
  d_co_sum <- summary(d_co)

  b0hat <- unname(d_co$coefficients[1])
  b1hat <- unname(d_co$coefficients[2])
  tco <- d_co_sum$coefficients[2, 3]

  # Step 5: VIF computation using AR from original series
  x_p <- aic_x$p
  x_phi <- aic_x$phi

  if (x_p > 0) {
    acf_vals <- as.numeric(stats::ARMAacf(ar = x_phi, lag.max = n - 1))[-1]
  } else {
    acf_vals <- rep(0, n - 1)  # White noise
  }

  # VIF = 1 + 2 * sum((1 - i/(n-2)) * rho(i)) for i = 1 to n-3
  max_lag <- n - 3
  vif <- 1.0
  for (i in seq_len(max_lag)) {
    vif <- vif + 2 * (1 - i / (n - 2)) * acf_vals[i]
  }
  vif <- max(vif, 1.0)  # Ensure VIF >= 1

  # Step 6: Adjusted df and p-value
  rel_n <- (n - 2) / vif
  est_df <- if (rel_n - pp > 0) rel_n - pp else rel_n
  est_df <- max(est_df, 1.0)  # Ensure df >= 1

  pvalue <- 2 * stats::pt(abs(tco), est_df, lower.tail = FALSE)

  # Return results
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
