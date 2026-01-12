// =============================================================================
// FILE: api_burg_aic.cpp
// CATEGORY: INTERFACE (R-facing)
// THREAD-SAFE: NO (returns Rcpp::List)
//
// Public Rcpp export for Burg algorithm with AIC/BIC model selection.
// This is a thin wrapper around the kernel implementation to ensure
// consistent behavior across all code paths.
//
// Exports:
//   - burg_aic_select_cpp(): Primary Burg + AIC API
// =============================================================================
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "kernel_types.h"

using namespace Rcpp;

// Forward declaration of kernel function (defined in kernel_burg_aic.cpp)
BurgResult burg_aic_select_pure(const arma::vec& x, int maxp,
                                 const std::string& criterion,
                                 int min_p);

//' Burg AR Fit with AIC Selection (C++ implementation)
//'
//' Estimates AR coefficients using Burg algorithm and selects optimal
//' order via AIC. Uses Levinson-Durbin recursion (same as R's ar.burg).
//'
//' @param x Numeric vector, the time series.
//' @param maxp Integer, maximum AR order to consider.
//' @param criterion String, information criterion: "aic", "aicc", or "bic".
//' @param min_p Integer, minimum AR order to consider. Default 0 includes AR(0).
//'   Set to 1 to exclude AR(0) and match R's aic_burg(p=1:maxp).
//' @return List with:
//'   - p: selected AR order
//'   - phi: AR coefficients (length p)
//'   - vara: residual variance
//'   - aic: AIC value for selected model
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Rcpp::List burg_aic_select_cpp(const arma::vec& x, int maxp,
                                std::string criterion = "aic",
                                int min_p = 0) {
  // Call kernel implementation (single source of truth)
  BurgResult result = burg_aic_select_pure(x, maxp, criterion, min_p);

  // Convert to R list
  return Rcpp::List::create(
    Rcpp::Named("p") = result.p,
    Rcpp::Named("phi") = result.phi,
    Rcpp::Named("vara") = result.vara,
    Rcpp::Named("aic") = result.ic
  );
}
