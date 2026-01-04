# cpp_wrappers.R
# Thin R wrappers for internal C++ functions used in validation/benchmarking

#' OLS Detrend (C++ Implementation)
#'
#' Remove linear trend from a time series using OLS.
#'
#' @param x Numeric vector.
#' @return Numeric vector of residuals after removing linear trend.
#' @export
ols_detrend <- function(x) {
 ols_detrend_cpp(as.numeric(x))
}

#' Burg AR Fit (C++ Implementation)
#'
#' Fit AR coefficients using Burg's algorithm.
#'
#' @param x Numeric vector, the time series.
#' @param p Integer, AR order to fit.
#' @return Numeric vector of AR coefficients.
#' @export
burg_fit <- function(x, p) {
 as.numeric(burg_fit_cpp(as.numeric(x), as.integer(p)))
}

#' Burg AR with AIC Selection (C++ Implementation)
#'
#' Select AR order using information criterion and fit coefficients.
#'
#' @param x Numeric vector, the time series (should be centered).
#' @param maxp Integer, maximum AR order to consider.
#' @param criterion Character, information criterion: "aic", "aicc", or "bic".
#' @return List with components:
#'   \item{p}{Selected AR order.}
#'   \item{phi}{AR coefficients.}
#'   \item{vara}{Innovation variance.}
#' @export
burg_aic_select <- function(x, maxp = 5, criterion = c("aic", "aicc", "bic")) {
 criterion <- match.arg(criterion)
 burg_aic_select_cpp(as.numeric(x), as.integer(maxp), criterion)
}

#' AR Transform (C++ Implementation)
#'
#' Apply AR filter to produce residuals.
#'
#' @param x Numeric vector, the time series.
#' @param phi Numeric vector, AR coefficients.
#' @return Numeric vector of filtered residuals (length n - p).
#' @export
ar_transform <- function(x, phi) {
 as.numeric(ar_transform_cpp(as.numeric(x), as.numeric(phi)))
}

#' Cochrane-Orcutt t-statistic (C++ Implementation)
#'
#' Compute the CO t-statistic for testing trend.
#'
#' @param x Numeric vector, the time series.
#' @param maxp Integer, maximum AR order for model selection.
#' @param criterion Character, information criterion: "aic", "aicc", or "bic".
#' @return Numeric, the t-statistic for trend.
#' @export
co_tstat <- function(x, maxp = 5, criterion = c("aic", "aicc", "bic")) {
 criterion <- match.arg(criterion)
 co_tstat_cpp(as.numeric(x), as.integer(maxp), criterion)
}

#' Cochrane-Orcutt Full Results (C++ Implementation)
#'
#' Compute the CO t-statistic with full model details.
#'
#' @param x Numeric vector, the time series.
#' @param maxp Integer, maximum AR order for model selection.
#' @param criterion Character, information criterion: "aic", "aicc", or "bic".
#' @return List with components:
#'   \item{tco}{t-statistic for trend.}
#'   \item{p}{Selected AR order.}
#'   \item{phi}{AR coefficients.}
#'   \item{vara}{Innovation variance.}
#' @export
co_full <- function(x, maxp = 5, criterion = c("aic", "aicc", "bic")) {
 criterion <- match.arg(criterion)
 co_full_cpp(as.numeric(x), as.integer(maxp), criterion)
}

#' Generate AR Series with Seed (C++ Implementation)
#'
#' Generate AR series using thread-safe RNG with explicit seed.
#'
#' @param n Integer, length of series to generate.
#' @param phi Numeric vector, AR coefficients.
#' @param vara Numeric, innovation variance.
#' @param seed Integer, RNG seed for reproducibility.
#' @return Numeric vector of length n.
#' @export
gen_ar_seeded <- function(n, phi, vara = 1.0, seed) {
 as.numeric(gen_ar_seeded_cpp(as.integer(n), as.numeric(phi),
                               as.numeric(vara), as.numeric(seed)))
}
