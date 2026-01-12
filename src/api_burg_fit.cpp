// =============================================================================
// FILE: api_burg_fit.cpp
// CATEGORY: INTERFACE (R-facing)
// THREAD-SAFE: NO (returns Rcpp::List)
//
// Fixed-order Burg algorithm for AR coefficient estimation.
// Delegates to burg_fit_pure() in kernel_burg_aic.cpp.
// For automatic order selection, use burg_aic_select_cpp() instead.
//
// Exports:
//   - burg_fit_cpp(): Fixed-order Burg AR fit
// =============================================================================
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// Forward declaration of kernel function (defined in kernel_burg_aic.cpp)
void burg_fit_pure(const arma::vec& x, int p, arma::vec& phi_out, double& vara_out);

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
// [[Rcpp::export]]
arma::vec burg_fit_cpp(const arma::vec& x, int p) {
  const int n = x.n_elem;

  if (p <= 0) {
    return arma::vec();  // Empty vector for AR(0)
  }

  if (p >= n) {
    Rcpp::stop("Order p must be less than series length n");
  }

  // Delegate to kernel function
  arma::vec phi;
  double vara;  // Not used, but required by burg_fit_pure signature
  burg_fit_pure(x, p, phi, vara);

  return phi;
}
