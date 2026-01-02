#' Create Cochrane-Orcutt Statistic Function
#'
#' Creates a function that computes the Cochrane-Orcutt t-statistic for
#' testing trend in a time series with autocorrelated errors.
#'
#' @importFrom stats cor median spec.ar
#'
#' @param maxp Integer. Maximum AR order to consider for the error model.
#'   Default is 5.
#' @param ar_method Character. Method for AR estimation: `"burg"` (default)
#'   or `"mle"`.
#' @param criterion Character. Information criterion for order selection:
#'   `"aic"` (default), `"aicc"`, or `"bic"`.
#'
#' @return A function that takes a numeric vector `x` and returns
#'   the Cochrane-Orcutt t-statistic for the trend coefficient.
#'
#' @seealso [wbg_boot_flex()], [co()]
#'
#' @examples
#' # Create statistic function
#' stat_fn <- make_stat_co(maxp = 5)
#'
#' # Use with bootstrap
#' set.seed(123)
#' x <- arima.sim(list(ar = 0.7), n = 100)
#' stat_fn(x)  # Returns t-statistic
#'
#' @export
make_stat_co <- function(maxp = 5L, ar_method = c("burg", "mle"),
                         criterion = c("aic", "aicc", "bic")) {
  ar_method <- match.arg(ar_method)
  criterion <- match.arg(criterion)

  # Force evaluation to capture in closure
  force(maxp)
  force(ar_method)
  force(criterion)

  function(x) {
    co(x, maxp = maxp, method = ar_method, type = criterion)$tco
  }
}


#' Create OLS t-Statistic Function
#'
#' Creates a function that computes the ordinary least squares t-statistic
#' for the slope in a linear trend regression. This statistic does not
#' account for autocorrelation in the standard error, but can be used
#' within the bootstrap framework which handles the dependence structure.
#'
#' @return A function that takes a numeric vector `x` and returns
#'   the OLS t-statistic for the trend coefficient.
#'
#' @seealso [wbg_boot_flex()], [make_stat_co()]
#'
#' @examples
#' stat_fn <- make_stat_ols_t()
#'
#' set.seed(123)
#' x <- 0.1 * (1:100) + rnorm(100)
#' stat_fn(x)  # Should be large (significant trend)
#'
#' @export
make_stat_ols_t <- function() {
  function(x) {
    fit <- lm(x ~ seq_along(x))
    summary(fit)$coefficients[2, 3]
  }
}


#' Create OLS Slope Statistic Function
#'
#' Creates a function that computes the ordinary least squares slope
#' estimate (not t-statistic) for a linear trend.
#'
#' @return A function that takes a numeric vector `x` and returns
#'   the OLS slope estimate.
#'
#' @seealso [wbg_boot_flex()], [make_stat_ols_t()]
#'
#' @examples
#' stat_fn <- make_stat_ols_slope()
#'
#' set.seed(123)
#' x <- 0.1 * (1:100) + rnorm(100)
#' stat_fn(x)  # Should be close to 0.1
#'
#' @export
make_stat_ols_slope <- function() {
  function(x) {
    lm(x ~ seq_along(x))$coefficients[2]
  }
}


#' Create Mann-Kendall Statistic Function
#'
#' Creates a function that computes the Mann-Kendall S statistic for
#' testing monotonic trend. This is a non-parametric, rank-based test
#' that is robust to outliers and non-normality.
#'
#' @return A function that takes a numeric vector `x` and returns
#'   the Mann-Kendall S statistic.
#'
#' @note Requires the \pkg{Kendall} package. Install with
#'   `install.packages("Kendall")`.
#'
#' @seealso [wbg_boot_flex()], [make_stat_spearman()]
#'
#' @examples
#' \dontrun{
#' stat_fn <- make_stat_mk()
#'
#' set.seed(123)
#' x <- 0.1 * (1:100) + rnorm(100)
#' stat_fn(x)  # Mann-Kendall S statistic
#' }
#'
#' @export
make_stat_mk <- function() {
  if (!requireNamespace("Kendall", quietly = TRUE)) {
    stop("Package 'Kendall' required. Install with: install.packages('Kendall')",
         call. = FALSE)
  }

  function(x) {
    Kendall::MannKendall(x)$S
  }
}


#' Create Spearman Correlation Statistic Function
#'
#' Creates a function that computes the Spearman rank correlation between
#' the series values and time. This is a non-parametric measure of
#' monotonic trend.
#'
#' @return A function that takes a numeric vector `x` and returns
#'   the Spearman correlation with time.
#'
#' @seealso [wbg_boot_flex()], [make_stat_mk()]
#'
#' @examples
#' stat_fn <- make_stat_spearman()
#'
#' set.seed(123)
#' x <- 0.1 * (1:100) + rnorm(100)
#' stat_fn(x)  # Should be positive (upward trend)
#'
#' @export
make_stat_spearman <- function() {
  function(x) {
    cor(x, seq_along(x), method = "spearman")
  }
}


#' Create Sen's Slope Statistic Function
#'
#' Creates a function that computes Sen's slope estimator, which is the
#' median of all pairwise slopes. This is a robust, non-parametric
#' estimate of trend magnitude.
#'
#' @return A function that takes a numeric vector `x` and returns
#'   Sen's slope estimate.
#'
#' @details
#' Sen's slope is defined as the median of all slopes computed between
#' pairs of points:
#' \deqn{\beta = \text{median}\left(\frac{x_j - x_i}{j - i}\right)}
#' for all \eqn{i < j}.
#'
#' This is more robust to outliers than OLS slope estimation.
#'
#' @seealso [wbg_boot_flex()], [make_stat_mk()]
#'
#' @examples
#' stat_fn <- make_stat_sen()
#'
#' set.seed(123)
#' x <- 0.1 * (1:100) + rnorm(100)
#' stat_fn(x)  # Should be close to 0.1
#'
#' @export
make_stat_sen <- function() {
  function(x) {
    sen_slope_cpp(x)
  }
}


#' Create HAC (Newey-West) t-Statistic Function
#'
#' Creates a function that computes the t-statistic for trend using
#' heteroskedasticity and autocorrelation consistent (HAC) standard errors
#' via the Newey-West estimator.
#'
#' @param lag Integer or NULL. The lag truncation parameter for the
#'   Newey-West estimator. If NULL (default), the bandwidth is selected
#'   automatically.
#'
#' @return A function that takes a numeric vector `x` and returns
#'   the HAC-corrected t-statistic for the trend coefficient.
#'
#' @note Requires the \pkg{sandwich} and \pkg{lmtest} packages. Install with
#'   `install.packages(c("sandwich", "lmtest"))`.
#'
#' @seealso [wbg_boot_flex()], [make_stat_bn()]
#'
#' @examples
#' \dontrun{
#' stat_fn <- make_stat_hac(lag = 4)
#'
#' set.seed(123)
#' x <- 0.1 * (1:100) + arima.sim(list(ar = 0.7), n = 100)
#' stat_fn(x)  # HAC-corrected t-statistic
#' }
#'
#' @export
make_stat_hac <- function(lag = NULL) {
  if (!requireNamespace("sandwich", quietly = TRUE)) {
    stop("Package 'sandwich' required. Install with: install.packages('sandwich')",
         call. = FALSE)
  }
  if (!requireNamespace("lmtest", quietly = TRUE)) {
    stop("Package 'lmtest' required. Install with: install.packages('lmtest')",
         call. = FALSE)
  }

  force(lag)

  function(x) {
    fit <- lm(x ~ seq_along(x))
    ct <- lmtest::coeftest(fit, vcov = sandwich::NeweyWest(fit, lag = lag))
    ct[2, 3]
  }
}


#' Create Bloomfield-Nychka t-Statistic Function
#'
#' Creates a function that computes the Bloomfield-Nychka t-statistic,
#' which corrects the standard error for trend using a spectral estimate
#' of the long-run variance (spectral density at frequency zero).
#'
#' @param order.max Integer or NULL. Maximum AR order for spectral
#'   estimation via [stats::spec.ar()]. If NULL (default),
#'   the order is selected automatically.
#'
#' @return A function that takes a numeric vector `x` and returns
#'   the Bloomfield-Nychka t-statistic for the trend coefficient.
#'
#' @details
#' The procedure:
#' \enumerate{
#'   \item Fit OLS trend and extract residuals.
#'   \item Estimate spectral density of residuals at frequency zero using
#'     an AR model.
#'   \item Compute corrected variance: \eqn{Var(\hat{\beta}) = 12 \cdot f(0) / n^3}.
#'   \item Return \eqn{t = \hat{\beta} / \sqrt{Var(\hat{\beta})}}.
#' }
#'
#' @references
#' Bloomfield, P. and Nychka, D. (1992). Climate spectra and detecting
#' climate change. *Climatic Change*, 21, 275-287.
#'
#' @seealso [wbg_boot_flex()], [make_stat_hac()]
#'
#' @examples
#' stat_fn <- make_stat_bn(order.max = 10)
#'
#' set.seed(123)
#' x <- 0.1 * (1:100) + arima.sim(list(ar = 0.7), n = 100)
#' stat_fn(x)  # Bloomfield-Nychka t-statistic
#'
#' @export
make_stat_bn <- function(order.max = NULL) {
  force(order.max)

  function(x) {
    n <- length(x)
    t_idx <- seq_along(x)

    fit <- lm(x ~ t_idx)
    b_hat <- coef(fit)[2]
    resid <- residuals(fit)

    sp <- spec.ar(resid, n.freq = 1000, plot = FALSE, order.max = order.max)
    f0 <- sp$spec[1]

    var_b <- 12 * f0 / n^3
    b_hat / sqrt(var_b)
  }
}


#' Create Likelihood Ratio Statistic Function
#'
#' Creates a function that computes the likelihood ratio test statistic
#' for trend by comparing AR models with and without a linear trend term.
#'
#' @param order Integer. The AR order for the models. Default is 1.
#'
#' @return A function that takes a numeric vector `x` and returns
#'   the likelihood ratio statistic: \eqn{2 \cdot (\ell_1 - \ell_0)}
#'   where \eqn{\ell_1} is the log-likelihood with trend and \eqn{\ell_0}
#'   is without.
#'
#' @details
#' Under the null hypothesis of no trend, the likelihood ratio statistic
#' is asymptotically chi-squared with 1 degree of freedom. However, when
#' used within the bootstrap framework, the empirical null distribution
#' is used instead, which handles autocorrelation appropriately.
#'
#' @seealso [wbg_boot_flex()]
#'
#' @examples
#' stat_fn <- make_stat_lr(order = 2)
#'
#' set.seed(123)
#' x <- 0.1 * (1:100) + arima.sim(list(ar = 0.7), n = 100)
#' stat_fn(x)  # Likelihood ratio statistic
#'
#' @export
make_stat_lr <- function(order = 1L) {
  force(order)

  function(x) {
    t_idx <- seq_along(x)

    # Model without trend
    fit0 <- suppressWarnings(
      arima(x, order = c(order, 0, 0))
    )

    # Model with trend
    fit1 <- suppressWarnings(
      arima(x, order = c(order, 0, 0), xreg = t_idx)
    )

    2 * (fit1$loglik - fit0$loglik)
  }
}


#' Create GLS t-Statistic Function
#'
#' Creates a function that computes the generalized least squares t-statistic
#' for trend with AR(p) correlated errors.
#'
#' @param p Integer. The AR order for the error correlation structure.
#'   Default is 1.
#'
#' @return A function that takes a numeric vector `x` and returns
#'   the GLS t-statistic for the trend coefficient.
#'
#' @note Requires the \pkg{nlme} package. Install with
#'   `install.packages("nlme")`.
#'
#' @seealso [wbg_boot_flex()], [make_stat_co()]
#'
#' @examples
#' \dontrun{
#' stat_fn <- make_stat_gls(p = 2)
#'
#' set.seed(123)
#' x <- 0.1 * (1:100) + arima.sim(list(ar = 0.7), n = 100)
#' stat_fn(x)  # GLS t-statistic
#' }
#'
#' @export
make_stat_gls <- function(p = 1L) {
  if (!requireNamespace("nlme", quietly = TRUE)) {
    stop("Package 'nlme' required. Install with: install.packages('nlme')",
         call. = FALSE)
  }

  force(p)

  function(x) {
    df <- data.frame(x = x, t = seq_along(x))
    fit <- nlme::gls(x ~ t, data = df, correlation = nlme::corARMA(p = p, q = 0))
    summary(fit)$tTable["t", "t-value"]
  }
}
