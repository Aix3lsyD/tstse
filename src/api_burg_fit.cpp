// =============================================================================
// FILE: api_burg_fit.cpp
// CATEGORY: INTERFACE (R-facing)
// THREAD-SAFE: NO (returns Rcpp::List)
//
// Fixed-order Burg algorithm for AR coefficient estimation.
// For automatic order selection, use burg_aic_select_cpp() instead.
//
// Exports:
//   - burg_fit_cpp(): [DEPRECATED] Use burg_aic_select_cpp()
//   - burg_fit_full_cpp(): [DEPRECATED] Use burg_aic_select_cpp()
// =============================================================================
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

//' Burg AR Fit (C++ implementation)
//'
//' Estimates AR coefficients using the Burg algorithm.
//' Equivalent to: ar.burg(x, order.max = p, aic = FALSE)$ar
//'
//' @param x Numeric vector, the centered time series.
//' @param p Integer, the AR order to fit.
//' @return Numeric vector of AR coefficients (length p).
//' @keywords internal
//' @noRd
[[deprecated("Use burg_aic_select_cpp() for automatic order selection")]]
// [[Rcpp::export]]
arma::vec burg_fit_cpp(const arma::vec& x, int p) {
  const int n = x.n_elem;

  if (p <= 0) {
    return arma::vec();  // Empty vector for AR(0)
  }

  if (p >= n) {
    Rcpp::stop("Order p must be less than series length n");
  }

  // Center the series
  const double x_mean = arma::mean(x);
  arma::vec xc = x - x_mean;

  // Initialize forward and backward prediction errors
  arma::vec ef = xc;
  arma::vec eb = xc;

  // AR coefficients (Levinson recursion)
  arma::vec a(p, arma::fill::zeros);

  for (int k_order = 1; k_order <= p; ++k_order) {
    // Compute reflection coefficient
    double num = 0.0;
    double den = 0.0;

    for (int t = k_order; t < n; ++t) {
      num += 2.0 * ef[t] * eb[t - 1];
      den += ef[t] * ef[t] + eb[t - 1] * eb[t - 1];
    }

    double k = (den > 1e-15) ? num / den : 0.0;

    // Levinson recursion: update AR coefficients
    if (k_order == 1) {
      a[0] = k;
    } else {
      arma::vec a_new(k_order, arma::fill::zeros);
      a_new[k_order - 1] = k;
      for (int j = 0; j < k_order - 1; ++j) {
        a_new[j] = a[j] - k * a[k_order - 2 - j];
      }
      a.head(k_order) = a_new;
    }

    // Update prediction errors
    // IMPORTANT: Both updates must use OLD values
    arma::vec ef_new(n);
    arma::vec eb_new(n);
    for (int t = k_order; t < n; ++t) {
      ef_new[t] = ef[t] - k * eb[t - 1];
      eb_new[t] = eb[t - 1] - k * ef[t];
    }
    ef = ef_new;
    eb = eb_new;
  }

  return a;
}


//' Burg AR Fit with Variance (C++ implementation)
//'
//' Returns AR coefficients and residual variance.
//'
//' @param x Numeric vector, the time series (not necessarily centered).
//' @param p Integer, the AR order to fit.
//' @return List with phi (coefficients) and vara (residual variance).
//' @keywords internal
//' @noRd
[[deprecated("Use burg_aic_select_cpp() for automatic order selection")]]
// [[Rcpp::export]]
Rcpp::List burg_fit_full_cpp(const arma::vec& x, int p) {
  const int n = x.n_elem;

  if (p <= 0) {
    // AR(0): just return variance of centered series
    const double x_mean = arma::mean(x);
    arma::vec xc = x - x_mean;
    double vara = arma::dot(xc, xc) / n;
    return Rcpp::List::create(
      Rcpp::Named("phi") = arma::vec(),
      Rcpp::Named("vara") = vara
    );
  }

  if (p >= n) {
    Rcpp::stop("Order p must be less than series length n");
  }

  // Center the series
  const double x_mean = arma::mean(x);
  arma::vec xc = x - x_mean;

  // Initialize forward and backward prediction errors
  arma::vec ef = xc;
  arma::vec eb = xc;

  // AR coefficients
  arma::vec a(p, arma::fill::zeros);

  for (int k_order = 1; k_order <= p; ++k_order) {
    // Compute reflection coefficient
    double num = 0.0;
    double den = 0.0;

    for (int t = k_order; t < n; ++t) {
      num += 2.0 * ef[t] * eb[t - 1];
      den += ef[t] * ef[t] + eb[t - 1] * eb[t - 1];
    }

    double k = (den > 1e-15) ? num / den : 0.0;

    // Levinson recursion
    if (k_order == 1) {
      a[0] = k;
    } else {
      arma::vec a_new(k_order, arma::fill::zeros);
      a_new[k_order - 1] = k;
      for (int j = 0; j < k_order - 1; ++j) {
        a_new[j] = a[j] - k * a[k_order - 2 - j];
      }
      a.head(k_order) = a_new;
    }

    // Update prediction errors
    // IMPORTANT: Both updates must use OLD values
    arma::vec ef_new(n);
    arma::vec eb_new(n);
    for (int t = k_order; t < n; ++t) {
      ef_new[t] = ef[t] - k * eb[t - 1];
      eb_new[t] = eb[t - 1] - k * ef[t];
    }
    ef = ef_new;
    eb = eb_new;
  }

  // Compute residual variance from final prediction errors
  double ss = 0.0;
  int count = 0;
  for (int t = p; t < n; ++t) {
    ss += ef[t] * ef[t];
    count++;
  }
  double vara = (count > 0) ? ss / count : 0.0;

  return Rcpp::List::create(
    Rcpp::Named("phi") = a,
    Rcpp::Named("vara") = vara
  );
}
