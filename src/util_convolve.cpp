// =============================================================================
// FILE: util_convolve.cpp
// CATEGORY: UTILITY
// THREAD-SAFE: YES (pure computation)
//
// Truncated polynomial convolution for Gegenbauer coefficient sequences.
// Used with util_gegenb.cpp for long-memory ARUMA fractional differencing.
//
// Exports:
//   - convolve_truncated_cpp(): First n coefficients of polynomial product
// =============================================================================

#include <Rcpp.h>
using namespace Rcpp;

//' Truncated Polynomial Convolution (C++ implementation)
//'
//' Computes the first n coefficients of the product of two power series.
//' Used for multiplying Gegenbauer polynomial coefficient sequences.
//'
//' @param C1 NumericVector, first coefficient sequence.
//' @param C2 NumericVector, second coefficient sequence.
//' @param n Integer, number of output coefficients.
//' @return NumericVector of length n containing convolution result.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
NumericVector convolve_truncated_cpp(NumericVector C1, NumericVector C2, int n) {
  NumericVector psi(n);

  for (int j = 0; j < n; ++j) {
    double sum = 0.0;
    for (int i = 0; i <= j; ++i) {
      sum += C1[i] * C2[j - i];
    }
    psi[j] = sum;
  }

  return psi;
}
