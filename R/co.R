#' Cochrane-Orcutt Estimation for Trend with Correlated Errors
#'
#' Estimates linear trend in time series data with autocorrelated residuals
#' using the Cochrane-Orcutt procedure. The residuals are modeled as an AR(p)
#' process with order selected by information criterion.
#'
#' @param x A numeric vector containing the time series data.
#' @param maxp Maximum AR order for model selection. Default is 5.
#' @param method Character. Method for AR estimation: `"mle"`, `"burg"`,
#'   or `"yw"`. Default is `"burg"`.
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
#'   \item (If `iterate = TRUE`) Re-estimates AR from CO residuals and repeats
#'     steps 3-4 until convergence
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
               type = c("aic", "aicc", "bic"),
               iterate = FALSE, max_iter = 50L, tol = 1e-6) {

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
  ar_method <- ar_fit$method

  # Helper function to perform CO transformation and regression
  co_transform_regress <- function(phi_current) {
    pp_cur <- length(phi_current)

    # Transform the data
    x_trans <- artrans(x, phi = phi_current, plot = FALSE)

    # Compute transformed time index
    p1 <- pp_cur + 1L
    t_co <- numeric(n)
    for (tt in p1:n) {
      t_co[tt] <- tt - sum(phi_current * (tt - seq_len(pp_cur)))
    }

    # Regress transformed data on transformed time
    d_co <- stats::lm(x_trans ~ t_co[p1:n])

    list(
      x_trans = x_trans,
      t_co = t_co,
      d_co = d_co,
      p1 = p1,
      residuals = stats::resid(d_co)
    )
  }

  # Perform initial CO transformation and regression
  co_result <- co_transform_regress(phi)

  # Iterative CO: re-estimate AR from CO residuals until convergence
  if (iterate && max_iter > 1L) {
    for (iter in seq_len(max_iter - 1L)) {
      phi_old <- phi

      # Re-fit AR to CO residuals
      ar_fit_new <- aic_ar(co_result$residuals, p = seq_len(maxp),
                           method = method, type = type, cores = 1L)
      pp <- ar_fit_new$p
      phi <- ar_fit_new$phi

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
  d_co <- co_result$d_co
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
    method    = ar_method
  )
}
