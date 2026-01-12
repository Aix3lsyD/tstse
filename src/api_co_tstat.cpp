// =============================================================================
// FILE: api_co_tstat.cpp
// CATEGORY: INTERFACE (R-facing)
// THREAD-SAFE: NO (uses Rcpp types, not for parallel use)
//
// Public Rcpp exports for Cochrane-Orcutt t-statistic.
// For parallel bootstrap, use kernel_co_tstat.cpp functions instead.
//
// Exports:
//   - co_tstat_cpp(): Primary CO t-statistic API
//   - co_full_cpp(): Full results including AR coefficients
// =============================================================================
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// Forward declarations of helper functions (defined in other files)
arma::vec ols_detrend_cpp(const arma::vec& x);
Rcpp::List burg_aic_select_cpp(const arma::vec& x, int maxp, std::string criterion, int min_p = 0);
arma::vec ar_transform_cpp(const arma::vec& x, const arma::vec& phi);
arma::vec co_time_transform_cpp(int n, const arma::vec& phi);

// Forward declaration of kernel function (defined in kernel_co_tstat.cpp)
double ols_tstat_internal(const arma::vec& y, const arma::vec& t_idx);


//' Cochrane-Orcutt t-statistic (C++ implementation)
//'
//' Computes the Cochrane-Orcutt t-statistic for testing H0: slope = 0.
//' This is a fused implementation that combines all steps in C++.
//' Equivalent to: co(x, maxp, method = "burg", type = "aic")$tco
//'
//' @param x Numeric vector, the time series.
//' @param maxp Integer, maximum AR order for model selection.
//' @param criterion String, information criterion: "aic", "aicc", or "bic".
//' @return Double, Cochrane-Orcutt t-statistic.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double co_tstat_cpp(const arma::vec& x, int maxp = 5,
                    std::string criterion = "aic") {
  const int n = x.n_elem;

  if (n < maxp + 3) {
    Rcpp::stop("Series too short for specified maxp");
  }

  // Step 1: OLS detrend
  arma::vec resid = ols_detrend_cpp(x);

  // Step 2: Burg AR selection on residuals
  Rcpp::List ar_fit = burg_aic_select_cpp(resid, maxp, criterion);
  arma::vec phi = ar_fit["phi"];

  // Step 3: AR transform the original data
  arma::vec x_trans = ar_transform_cpp(x, phi);

  // Step 4: CO time transform
  arma::vec t_co = co_time_transform_cpp(n, phi);

  // Step 5: OLS on transformed data
  double tstat = ols_tstat_internal(x_trans, t_co);

  return tstat;
}


//' Cochrane-Orcutt Full Results (C++ implementation)
//'
//' Returns full CO results including AR order and coefficients.
//' For debugging and validation purposes.
//'
//' @param x Numeric vector, the time series.
//' @param maxp Integer, maximum AR order for model selection.
//' @param criterion String, information criterion: "aic", "aicc", or "bic".
//' @return List with tco, p, phi, and vara.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Rcpp::List co_full_cpp(const arma::vec& x, int maxp = 5,
                       std::string criterion = "aic") {
  const int n = x.n_elem;

  // Step 1: OLS detrend
  arma::vec resid = ols_detrend_cpp(x);

  // Step 2: Burg AR selection on residuals
  Rcpp::List ar_fit = burg_aic_select_cpp(resid, maxp, criterion);
  int p = ar_fit["p"];
  arma::vec phi = ar_fit["phi"];
  double vara = ar_fit["vara"];

  // Step 3: AR transform the original data
  arma::vec x_trans = ar_transform_cpp(x, phi);

  // Step 4: CO time transform
  arma::vec t_co = co_time_transform_cpp(n, phi);

  // Step 5: OLS on transformed data
  double tstat = ols_tstat_internal(x_trans, t_co);

  return Rcpp::List::create(
    Rcpp::Named("tco") = tstat,
    Rcpp::Named("p") = p,
    Rcpp::Named("phi") = phi,
    Rcpp::Named("vara") = vara
  );
}
