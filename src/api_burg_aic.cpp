// =============================================================================
// FILE: api_burg_aic.cpp
// CATEGORY: INTERFACE (R-facing)
// THREAD-SAFE: NO (returns Rcpp::List)
//
// Public Rcpp export for Burg algorithm with AIC/BIC model selection.
// Matches R's ar.burg() with var.method=1 (recursive variance).
// For parallel bootstrap, use kernel_burg_aic.cpp functions instead.
//
// Exports:
//   - burg_aic_select_cpp(): Primary Burg + AIC API
// =============================================================================
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <limits>
using namespace Rcpp;

//' Burg AR Fit with AIC Selection (C++ implementation)
//'
//' Estimates AR coefficients using Burg algorithm and selects optimal
//' order via AIC. Matches R's ar.burg() with var.method=1.
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
  const int n = x.n_elem;

  if (maxp <= 0) {
    Rcpp::stop("maxp must be positive");
  }
  if (maxp >= n - 1) {
    maxp = n - 2;
  }

  // Center the series
  const double x_mean = arma::mean(x);
  arma::vec xc = x - x_mean;

  // Initialize forward and backward prediction errors
  arma::vec ef = xc;
  arma::vec eb = xc;

  // AR(0) variance
  double vara0 = arma::dot(xc, xc) / n;

  // Storage for best model
  double best_ic = std::numeric_limits<double>::infinity();
  int best_p = 0;
  arma::vec best_phi;
  double best_vara = vara0;

  // IC for AR(0): k = 1 parameter (mean)
  // R uses: n * log(var) + 2 * (p + 1) when demean=TRUE
  // For AR(0), p=0, so IC = n * log(var0) + 2 * 1
  // Only consider AR(0) if min_p == 0
  if (min_p == 0) {
    double ic0;
    if (criterion == "aic") {
      ic0 = n * std::log(vara0) + 2.0;
    } else if (criterion == "aicc") {
      ic0 = n * std::log(vara0) + 2.0 * n / (n - 2);
    } else {  // bic
      ic0 = n * std::log(vara0) + std::log(static_cast<double>(n));
    }
    best_ic = ic0;
  }
  // When min_p > 0, best_ic stays at infinity until we find a valid AR(p) model

  // Current variance (recursive formula like R's var1)
  double var_recursive = vara0;

  // Temporary storage for Levinson recursion
  arma::vec a_prev;

  // Determine starting order (always compute from p=1 for correct Levinson recursion,
  // but only consider orders >= min_p for IC comparison)
  const int start_p = std::max(min_p, 1);

  // Single pass through all orders
  for (int p = 1; p <= maxp; ++p) {
    // Compute reflection coefficient
    double num = 0.0;
    double den = 0.0;

    for (int t = p; t < n; ++t) {
      num += ef[t] * eb[t - 1];
      den += ef[t] * ef[t] + eb[t - 1] * eb[t - 1];
    }

    double phii = (den > 1e-15) ? 2.0 * num / den : 0.0;

    // Update variance using recursive formula (R's var1)
    var_recursive = var_recursive * (1.0 - phii * phii);

    // Levinson recursion for AR coefficients
    arma::vec a_curr(p);
    a_curr(p - 1) = phii;  // Last coefficient is reflection coefficient

    if (p > 1) {
      for (int j = 0; j < p - 1; ++j) {
        a_curr(j) = a_prev(j) - phii * a_prev(p - 2 - j);
      }
    }
    a_prev = a_curr;

    // Update prediction errors for next iteration
    arma::vec ef_new = ef;
    arma::vec eb_new = eb;
    for (int t = p; t < n; ++t) {
      ef_new(t) = ef(t) - phii * eb(t - 1);
      eb_new(t) = eb(t - 1) - phii * ef(t);
    }
    ef = ef_new;
    eb = eb_new;

    // Compute information criterion
    // k = p + 1 (p AR coefficients + 1 for mean)
    int k = p + 1;
    double ic;

    if (var_recursive <= 0) {
      ic = std::numeric_limits<double>::infinity();
    } else if (criterion == "aic") {
      ic = n * std::log(var_recursive) + 2.0 * k;
    } else if (criterion == "aicc") {
      if (n - k - 1 > 0) {
        ic = n * std::log(var_recursive) + 2.0 * k * n / (n - k - 1);
      } else {
        ic = std::numeric_limits<double>::infinity();
      }
    } else {  // bic
      ic = n * std::log(var_recursive) + std::log(static_cast<double>(n)) * k;
    }

    // Update best model if this is better AND p >= min_p
    if (p >= start_p && ic < best_ic) {
      best_ic = ic;
      best_p = p;
      best_phi = a_curr;
      best_vara = var_recursive;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("p") = best_p,
    Rcpp::Named("phi") = best_phi,
    Rcpp::Named("vara") = best_vara,
    Rcpp::Named("aic") = best_ic
  );
}
