// =============================================================================
// FILE: util_gegenb.cpp
// CATEGORY: UTILITY
// THREAD-SAFE: YES (pure computation)
//
// Gegenbauer polynomial coefficients for long-memory ARUMA models.
// Used with util_convolve.cpp for fractional differencing.
//
// Exports:
//   - gegenb_cpp(): Compute Gegenbauer coefficients C_k(d, u)
// =============================================================================

#include <Rcpp.h>
using namespace Rcpp;

//' Gegenbauer Polynomial Coefficients (C++ implementation)
//'
//' Computes the coefficients C_k(d, u) using the recurrence relation.
//' This is the internal C++ implementation called by the R wrapper.
//'
//' @param u Numeric scalar, cosine of the Gegenbauer frequency.
//' @param d Numeric scalar, long-memory parameter.
//' @param n_coef Integer, number of coefficients to compute.
//' @return NumericVector of length n_coef containing C_0, C_1, ..., C_{n-1}.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
NumericVector gegenb_cpp(double u, double d, int n_coef) {

  NumericVector coef(n_coef);

  // C_0 = 1

coef[0] = 1.0;

  // Handle trivial cases
  if (n_coef == 1) {
    return coef;
  }

  // C_1 = 2 * d * u
  coef[1] = 2.0 * d * u;

  if (n_coef == 2) {
    return coef;
  }

  // Recurrence relation for k >= 2:
  // C_k = [2(k-1+d)u * C_{k-1} - (k-2+2d) * C_{k-2}] / k
  for (int i = 2; i < n_coef; ++i) {
    double k = static_cast<double>(i);
    coef[i] = (2.0 * (k - 1.0 + d) * u * coef[i - 1] -
               (k - 2.0 + 2.0 * d) * coef[i - 2]) / k;
  }

  return coef;
}
