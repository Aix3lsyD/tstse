#' Cochrane-Orcutt Estimation for Trend with Correlated Errors
#'
#' Estimates linear trend in time series data with autocorrelated residuals
#' using the Cochrane-Orcutt procedure. The residuals are modeled as an AR(p)
#' process with order selected by information criterion.
#'
#' @param x A numeric vector containing the time series data.
#' @param maxp Maximum AR order for model selection. Default is 5.
#' @param method Character. Method for AR estimation: `"mle"`, `"burg"`,
#'   or `"yw"`. Default is `"mle"`.
#' @param type Character. Information criterion for order selection:
#'   `"aic"`, `"aicc"`, or `"bic"`. Default is `"aic"`.
#'
#' @return A list containing:
#'   \item{x}{The original time series.}
#'   \item{z_x}{Residuals from initial OLS trend fit.}
#'   \item{b0hat}{Cochrane-Orcutt estimate of intercept.}
#'   \item{b1hat}{Cochrane-Orcutt estimate of slope.}
#'   \item{z_order}{AR order selected for the residuals.}
#'   \item{z_phi}{AR coefficients for the residual model.}
#'   \item{pvalue}{P-value for test of H0: slope = 0 (assumes uncorrelated
#'     transformed residuals; use [wbg_boot()] for more accurate
#'     p-values with correlated errors).}
#'   \item{tco}{t-statistic for the slope coefficient.}
#'   \item{method}{The AR estimation method used.}
#'
#' @details
#' Fits the model
#' \deqn{Y_t = a + bt + Z_t}
#' where \eqn{Z_t} is a stationary AR(p) process satisfying
#' \eqn{\phi(B)Z_t = a_t}.
#'
#' The procedure:
#' \enumerate{
#'   \item Obtains OLS estimates of \eqn{a} and \eqn{b}
#'   \item Computes residuals and fits AR(p) via the specified method with
#'     information criterion order selection
#'   \item Transforms the data using \eqn{\hat{\phi}(B)} to obtain
#'     \eqn{W_t = \hat{\phi}(B)Y_t}
#'   \item Regresses \eqn{W_t} on the transformed time index
#' }
#'
#' Note: The p-value returned assumes the transformed residuals are
#' uncorrelated, which is only approximate. For small to moderate sample sizes
#' with highly correlated errors, use [wbg_boot()] for bootstrap-based
#' inference with better significance level control.
#'
#' @references
#' Cochrane, D. and Orcutt, G. H. (1949). "Application of Least Squares
#' Regression to Relationships Containing Auto-Correlated Error Terms."
#' *Journal of the American Statistical Association*, 44(245), 32-61.
#'
#' Woodward, W. A., Bottone, S., and Gray, H. L. (1997). "Improved Tests for
#' Trend in Time Series Data." *Journal of Agricultural, Biological, and
#' Environmental Statistics*, 2(4), 403-416.
#'
#' @seealso [wbg_boot()] for bootstrap-based trend test,
#'   [aic_ar()], [artrans()]
#'
#' @examples
#' # Generate series with trend and AR(1) errors
#' set.seed(123)
#' n <- 100
#' t <- 1:n
#' z <- arima.sim(list(ar = 0.7), n = n)
#' x <- 5 + 0.1 * t + z
#'
#' # Estimate trend with Cochrane-Orcutt
#' result <- co(x, maxp = 5)
#' cat("Slope estimate:", result$b1hat, "\n")
#' cat("t-statistic:", result$tco, "\n")
#' cat("AR order:", result$z_order, "\n")
#'
#' @export
co <- function(x, maxp = 5L, method = c("burg", "mle", "yw"),
               type = c("aic", "aicc", "bic")) {

  method <- match.arg(method)
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

  # Step 1: OLS fit and residuals
 t1 <- seq_len(n)
  d <- stats::lm(x ~ t1)
  z_x <- stats::resid(d)

  # Step 2: Fit AR model to residuals
  # Use cores = 1L because co() is often called from within parallel contexts
  ar_fit <- aic_ar(z_x, p = seq_len(maxp), method = method, type = type, cores = 1L)
  pp <- ar_fit$p
  phi <- ar_fit$phi

  # Step 3: Transform the data
  x_trans <- artrans(x, phi = phi, plot = FALSE)

  # Step 4: Compute transformed time index
  # For AR(p), the transformed time is t - sum(phi_j * (t - j))
  p1 <- pp + 1L
  t_co <- numeric(n)

  for (tt in p1:n) {
    t_co[tt] <- tt - sum(phi * (tt - seq_len(pp)))
  }

  # Step 5: Regress transformed data on transformed time
  d_co <- stats::lm(x_trans ~ t_co[p1:n])
  d_co_sum <- summary(d_co)

  list(
    x         = x,
    z_x       = z_x,
    b0hat     = unname(d_co$coefficients[1]),
    b1hat     = unname(d_co$coefficients[2]),
    z_order   = pp,
    z_phi     = phi,
    pvalue    = d_co_sum$coefficients[2, 4],
    tco       = d_co_sum$coefficients[2, 3],
    method    = ar_fit$method
  )
}
