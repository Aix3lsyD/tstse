// co_tstat.cpp - Complete Cochrane-Orcutt t-statistic
// Part of wbg_boot_fast optimization
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// Forward declarations of helper functions (defined in other files)
arma::vec ols_detrend_cpp(const arma::vec& x);
Rcpp::List burg_aic_select_cpp(const arma::vec& x, int maxp, std::string criterion);
arma::vec ar_transform_cpp(const arma::vec& x, const arma::vec& phi);
arma::vec co_time_transform_cpp(int n, const arma::vec& phi);


//' OLS t-statistic for slope (C++ implementation)
//'
//' Computes t-statistic for slope coefficient from simple linear regression.
//'
//' @param y Numeric vector, response variable.
//' @param t_idx Numeric vector, predictor (time index).
//' @return Double, t-statistic for slope.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double ols_tstat_cpp(const arma::vec& y, const arma::vec& t_idx) {
  const int n = y.n_elem;

  if (n < 3) {
    return 0.0;  // Not enough data
  }

  // Means
  const double y_mean = arma::mean(y);
  const double t_mean = arma::mean(t_idx);

  // Centered vectors
  arma::vec y_c = y - y_mean;
  arma::vec t_c = t_idx - t_mean;

  // OLS slope
  const double ss_t = arma::dot(t_c, t_c);
  if (ss_t < 1e-15) {
    return 0.0;  // No variation in predictor
  }
  const double b_hat = arma::dot(y_c, t_c) / ss_t;

  // Intercept
  const double a_hat = y_mean - b_hat * t_mean;

  // Residuals and MSE
  arma::vec resid = y - a_hat - b_hat * t_idx;
  const double mse = arma::dot(resid, resid) / (n - 2);

  // Standard error of slope
  const double se_b = std::sqrt(mse / ss_t);

  if (se_b < 1e-15) {
    return 0.0;
  }

  return b_hat / se_b;
}


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
  double tstat = ols_tstat_cpp(x_trans, t_co);

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
  double tstat = ols_tstat_cpp(x_trans, t_co);

  return Rcpp::List::create(
    Rcpp::Named("tco") = tstat,
    Rcpp::Named("p") = p,
    Rcpp::Named("phi") = phi,
    Rcpp::Named("vara") = vara
  );
}
