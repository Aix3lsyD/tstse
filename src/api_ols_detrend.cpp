// =============================================================================
// FILE: api_ols_detrend.cpp
// CATEGORY: INTERFACE (R-facing)
// THREAD-SAFE: YES (pure arma operations)
//
// OLS detrending: removes linear trend from time series via regression.
// Equivalent to: resid(lm(x ~ seq_along(x)))
//
// Exports:
//   - ols_detrend_cpp(): Returns detrended residuals
// =============================================================================
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

//' OLS Detrend (C++ implementation)
//'
//' Removes linear trend from a time series via OLS regression.
//' Equivalent to: resid(lm(x ~ seq_along(x)))
//'
//' @param x Numeric vector, the time series.
//' @return Numeric vector of residuals after removing linear trend.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::vec ols_detrend_cpp(const arma::vec& x) {
  const int n = x.n_elem;

  // Time index: 1, 2, ..., n
  arma::vec t = arma::linspace(1.0, static_cast<double>(n), n);

  // Means
  const double x_mean = arma::mean(x);
  const double t_mean = (n + 1.0) / 2.0;  // Mean of 1..n

  // Centered vectors
  arma::vec x_c = x - x_mean;
  arma::vec t_c = t - t_mean;

  // OLS slope: b = sum((x - x_mean) * (t - t_mean)) / sum((t - t_mean)^2)
  const double b_hat = arma::dot(x_c, t_c) / arma::dot(t_c, t_c);

  // OLS intercept: a = x_mean - b * t_mean
  const double a_hat = x_mean - b_hat * t_mean;

  // Residuals: x - a - b*t
  arma::vec resid = x - a_hat - b_hat * t;

  return resid;
}
