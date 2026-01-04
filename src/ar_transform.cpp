// ar_transform.cpp - AR transformation for time series
// Part of wbg_boot_fast optimization
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

//' AR Transform (C++ implementation)
//'
//' Applies AR transformation phi(B) to a time series.
//' W_t = x_t - phi_1 * x_{t-1} - phi_2 * x_{t-2} - ... - phi_p * x_{t-p}
//' Equivalent to: artrans(x, phi, plot = FALSE)
//'
//' @param x Numeric vector, the time series.
//' @param phi Numeric vector, AR coefficients.
//' @return Numeric vector of transformed series (length n-p).
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::vec ar_transform_cpp(const arma::vec& x, const arma::vec& phi) {
  const int n = x.n_elem;
  const int p = phi.n_elem;

  if (p == 0) {
    // No transformation needed
    return x;
  }

  if (p >= n) {
    Rcpp::stop("AR order must be less than series length");
  }

  // Result has length n-p (we lose p observations)
  arma::vec w(n - p);

  for (int t = p; t < n; ++t) {
    double wt = x[t];
    for (int j = 0; j < p; ++j) {
      wt -= phi[j] * x[t - j - 1];
    }
    w[t - p] = wt;
  }

  return w;
}


//' CO Time Transform (C++ implementation)
//'
//' Transforms time index for Cochrane-Orcutt procedure.
//' t_co[t] = t - sum_{j=1}^p phi_j * (t - j)
//' for t = p+1, ..., n
//'
//' @param n Integer, series length.
//' @param phi Numeric vector, AR coefficients.
//' @return Numeric vector of transformed time indices (length n-p).
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::vec co_time_transform_cpp(int n, const arma::vec& phi) {
  const int p = phi.n_elem;

  if (p == 0) {
    // No transformation: return 1, 2, ..., n
    return arma::linspace(1.0, static_cast<double>(n), n);
  }

  if (p >= n) {
    Rcpp::stop("AR order must be less than series length");
  }

  // Result has length n-p
  arma::vec t_co(n - p);

  for (int idx = 0; idx < n - p; ++idx) {
    int t = idx + p + 1;  // Original time index (1-based)
    double tco_val = static_cast<double>(t);

    for (int j = 0; j < p; ++j) {
      tco_val -= phi[j] * (t - (j + 1));
    }

    t_co[idx] = tco_val;
  }

  return t_co;
}
