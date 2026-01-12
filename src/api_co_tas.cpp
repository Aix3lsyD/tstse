// =============================================================================
// FILE: api_co_tas.cpp
// CATEGORY: INTERFACE (R-facing)
// THREAD-SAFE: NO (returns Rcpp::List)
//
// Public Rcpp export for CO-TAS (Turner Adjusted Sample Size) trend test.
// This is a thin wrapper around the kernel implementation.
//
// Exports:
//   - co_tas_cpp(): Full CO-TAS result
//   - co_tas_pvalue_cpp(): P-value only (for R-level bootstrap)
// =============================================================================
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "kernel_types.h"

using namespace Rcpp;

// Forward declaration of kernel result struct (defined in kernel_co_tas.cpp)
struct CoTasResult {
    double pvalue;
    double tco;
    double n_a;
    arma::vec phi;
    int p;

    CoTasResult() : pvalue(1.0), tco(0.0), n_a(0.0), p(0) {}
};

// Forward declaration of kernel function
CoTasResult co_tas_pure(const arma::vec& x, int maxp, CriterionType ic_type);
double co_tas_pvalue_pure(const arma::vec& x, int maxp, CriterionType ic_type);

//' CO-TAS Trend Test (C++ implementation)
//'
//' Cochrane-Orcutt trend test with Turner's effective sample size adjustment.
//' Uses C++ for speed.
//'
//' @param x Numeric vector, the time series.
//' @param maxp Integer, maximum AR order for model selection (default 5).
//' @param type String, information criterion: "aic", "aicc", or "bic".
//' @return List with:
//'   - pvalue: P-value for trend test (adjusted for autocorrelation)
//'   - tco: t-statistic for slope
//'   - n_a: Effective sample size
//'   - phi: AR coefficients
//'   - p: AR order selected
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Rcpp::List co_tas_cpp(const arma::vec& x, int maxp = 5,
                       std::string type = "aic") {
    // Input validation
    if (x.n_elem < 10) {
        Rcpp::stop("Series must have at least 10 observations");
    }
    if (maxp < 1) {
        Rcpp::stop("maxp must be at least 1");
    }
    if (maxp > 20) {
        Rcpp::stop("maxp must be <= 20");
    }

    // Convert criterion string to enum
    CriterionType ic_type = criterion_from_string(type);

    // Call kernel
    CoTasResult result = co_tas_pure(x, maxp, ic_type);

    // Convert to R list
    return Rcpp::List::create(
        Rcpp::Named("pvalue") = result.pvalue,
        Rcpp::Named("tco") = result.tco,
        Rcpp::Named("n_a") = result.n_a,
        Rcpp::Named("phi") = result.phi,
        Rcpp::Named("p") = result.p
    );
}

//' CO-TAS P-value Only (C++ implementation)
//'
//' Returns only the p-value, for use in R-level bootstrap loops.
//'
//' @param x Numeric vector, the time series.
//' @param maxp Integer, maximum AR order (default 5).
//' @param type String, information criterion: "aic", "aicc", or "bic".
//' @return Double, the p-value.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double co_tas_pvalue_cpp(const arma::vec& x, int maxp = 5,
                          std::string type = "aic") {
    if (x.n_elem < 10) {
        return 1.0;  // Return high p-value for short series
    }

    CriterionType ic_type = criterion_from_string(type);
    return co_tas_pvalue_pure(x, maxp, ic_type);
}
