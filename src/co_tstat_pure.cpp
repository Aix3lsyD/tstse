// co_tstat_pure.cpp - Pure C++ Cochrane-Orcutt t-statistic
// Thread-safe version using C++ structs instead of Rcpp types
// [[Rcpp::depends(RcppArmadillo)]]

#include "types_pure.h"
#include <cmath>

// Forward declaration of pure C++ functions
BurgResult burg_aic_select_pure(const arma::vec& x, int maxp,
                                 const std::string& criterion);

// Internal: OLS detrend (pure C++)
static arma::vec ols_detrend_internal(const arma::vec& x) {
    const int n = x.n_elem;

    // Create time index 1, 2, ..., n
    arma::vec t = arma::linspace(1.0, static_cast<double>(n), n);

    // Means
    const double x_mean = arma::mean(x);
    const double t_mean = arma::mean(t);

    // Centered
    arma::vec x_c = x - x_mean;
    arma::vec t_c = t - t_mean;

    // OLS slope
    const double ss_t = arma::dot(t_c, t_c);
    const double b_hat = arma::dot(x_c, t_c) / ss_t;
    const double a_hat = x_mean - b_hat * t_mean;

    // Residuals
    return x - a_hat - b_hat * t;
}


// Internal: AR transform (pure C++)
static arma::vec ar_transform_internal(const arma::vec& x, const arma::vec& phi) {
    const int n = x.n_elem;
    const int p = phi.n_elem;

    if (p == 0) {
        return x;
    }

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


// Internal: CO time transform (pure C++)
static arma::vec co_time_transform_internal(int n, const arma::vec& phi) {
    const int p = phi.n_elem;

    if (p == 0) {
        return arma::linspace(1.0, static_cast<double>(n), n);
    }

    arma::vec t_co(n - p);
    for (int idx = 0; idx < n - p; ++idx) {
        int t = idx + p + 1;  // 1-based time
        double tco_val = static_cast<double>(t);
        for (int j = 0; j < p; ++j) {
            tco_val -= phi[j] * (t - (j + 1));
        }
        t_co[idx] = tco_val;
    }

    return t_co;
}


// Internal: OLS t-statistic for slope (pure C++)
static double ols_tstat_internal(const arma::vec& y, const arma::vec& t_idx) {
    const int n = y.n_elem;

    if (n < 3) {
        return 0.0;
    }

    const double y_mean = arma::mean(y);
    const double t_mean = arma::mean(t_idx);

    arma::vec y_c = y - y_mean;
    arma::vec t_c = t_idx - t_mean;

    const double ss_t = arma::dot(t_c, t_c);
    if (ss_t < 1e-15) {
        return 0.0;
    }

    const double b_hat = arma::dot(y_c, t_c) / ss_t;
    const double a_hat = y_mean - b_hat * t_mean;

    arma::vec resid = y - a_hat - b_hat * t_idx;
    const double mse = arma::dot(resid, resid) / (n - 2);
    const double se_b = std::sqrt(mse / ss_t);

    if (se_b < 1e-15) {
        return 0.0;
    }

    return b_hat / se_b;
}


// Pure C++ CO t-statistic
// Thread-safe: no Rcpp types used internally
double co_tstat_pure(const arma::vec& x, int maxp, const std::string& criterion) {
    const int n = x.n_elem;

    if (n < maxp + 3) {
        return 0.0;  // Series too short
    }

    // Step 1: OLS detrend
    arma::vec resid = ols_detrend_internal(x);

    // Step 2: Burg AR selection on residuals (pure C++ version)
    BurgResult ar_fit = burg_aic_select_pure(resid, maxp, criterion);

    // Step 3: AR transform original data
    arma::vec x_trans = ar_transform_internal(x, ar_fit.phi);

    // Step 4: CO time transform
    arma::vec t_co = co_time_transform_internal(n, ar_fit.phi);

    // Step 5: OLS on transformed data
    return ols_tstat_internal(x_trans, t_co);
}


// Pure C++ CO with full results
COResult co_full_pure(const arma::vec& x, int maxp, const std::string& criterion) {
    const int n = x.n_elem;

    // Step 1: OLS detrend
    arma::vec resid = ols_detrend_internal(x);

    // Step 2: Burg AR selection
    BurgResult ar_fit = burg_aic_select_pure(resid, maxp, criterion);

    // Step 3: AR transform
    arma::vec x_trans = ar_transform_internal(x, ar_fit.phi);

    // Step 4: CO time transform
    arma::vec t_co = co_time_transform_internal(n, ar_fit.phi);

    // Step 5: OLS t-stat
    double tstat = ols_tstat_internal(x_trans, t_co);

    return COResult(tstat, ar_fit.p, ar_fit.phi, ar_fit.vara);
}


// Rcpp export wrappers (for R interface and validation)
// [[Rcpp::export]]
double co_tstat_pure_export(const arma::vec& x, int maxp = 5,
                             std::string criterion = "aic") {
    return co_tstat_pure(x, maxp, criterion);
}

// [[Rcpp::export]]
Rcpp::List co_full_pure_export(const arma::vec& x, int maxp = 5,
                                std::string criterion = "aic") {
    COResult result = co_full_pure(x, maxp, criterion);

    return Rcpp::List::create(
        Rcpp::Named("tco") = result.tco,
        Rcpp::Named("p") = result.p,
        Rcpp::Named("phi") = result.phi,
        Rcpp::Named("vara") = result.vara
    );
}
