// co_tstat_pure.cpp - Pure C++ Cochrane-Orcutt t-statistic
// Thread-safe version using C++ structs instead of Rcpp types
// [[Rcpp::depends(RcppArmadillo)]]

#include "types_pure.h"
#include <cmath>

// Forward declaration of pure C++ functions
BurgResult burg_aic_select_pure(const arma::vec& x, int maxp,
                                 const std::string& criterion);

// Internal: OLS detrend (pure C++) - allocation-optimized version
// Uses closed-form formulas for t = 1, 2, ..., n to avoid temporary vectors
static arma::vec ols_detrend_internal(const arma::vec& x) {
    const int n = x.n_elem;

    // Closed-form formulas for t = 1, 2, ..., n
    const double t_mean = (n + 1.0) / 2.0;
    const double ss_t = n * (static_cast<double>(n) * n - 1.0) / 12.0;

    // Single pass: compute sum(x) and sum(t*x)
    double sum_x = 0.0, sum_tx = 0.0;
    for (int i = 0; i < n; ++i) {
        sum_x += x[i];
        sum_tx += (i + 1) * x[i];  // 1-based index
    }

    const double x_mean = sum_x / n;
    const double sxy = sum_tx - n * x_mean * t_mean;  // Σ(x-x̄)(t-t̄)
    const double b_hat = sxy / ss_t;
    const double a_hat = x_mean - b_hat * t_mean;

    // Only allocation: the return vector
    arma::vec resid(n);
    for (int i = 0; i < n; ++i) {
        resid[i] = x[i] - a_hat - b_hat * (i + 1);
    }
    return resid;
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


// Internal: OLS t-statistic for slope (pure C++) - kept for co_full_pure
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


// Fused CO t-statistic: single-pass, zero-allocation for steps 3-5
// Computes AR transform, time transform, and regression all on-the-fly
// Uses algebraic identity: SS_res = SS_y - b² * SS_t (no residual storage!)
static double co_tstat_fused(const arma::vec& x, const arma::vec& phi) {
    const int n = x.n_elem;
    const int p = phi.n_elem;
    const int m = (p == 0) ? n : n - p;

    if (m < 3) return 0.0;

    // Single pass: accumulate all sums needed for regression
    double sum_y = 0.0, sum_t = 0.0;
    double sum_yt = 0.0, sum_y2 = 0.0, sum_t2 = 0.0;

    for (int idx = 0; idx < m; ++idx) {
        // Compute AR-transformed x[t] on the fly (no vector storage)
        const int t = idx + p;
        double y_val = x[t];
        for (int j = 0; j < p; ++j) {
            y_val -= phi[j] * x[t - j - 1];
        }

        // Compute CO time transform on the fly (no vector storage)
        const int t_1based = t + 1;
        double t_val = static_cast<double>(t_1based);
        for (int j = 0; j < p; ++j) {
            t_val -= phi[j] * (t_1based - (j + 1));
        }

        // Accumulate sums
        sum_y += y_val;
        sum_t += t_val;
        sum_yt += y_val * t_val;
        sum_y2 += y_val * y_val;
        sum_t2 += t_val * t_val;
    }

    // Compute regression statistics from sums
    const double y_mean = sum_y / m;
    const double t_mean = sum_t / m;
    const double ss_t = sum_t2 - m * t_mean * t_mean;

    if (ss_t < 1e-15) return 0.0;

    const double sxy = sum_yt - m * y_mean * t_mean;
    const double b_hat = sxy / ss_t;

    // Key insight: SS_res = SS_y - b² * SS_t (algebraic identity)
    // This avoids storing residuals entirely!
    const double ss_y = sum_y2 - m * y_mean * y_mean;
    double ss_res = ss_y - b_hat * b_hat * ss_t;
    if (ss_res < 0) ss_res = 0;  // Numerical safety

    const double mse = ss_res / (m - 2);
    const double se_b = std::sqrt(mse / ss_t);

    if (se_b < 1e-15) return 0.0;

    return b_hat / se_b;
}


// Pure C++ CO t-statistic
// Thread-safe: no Rcpp types used internally
// Optimized: uses fused single-pass for steps 3-5 (zero allocations)
double co_tstat_pure(const arma::vec& x, int maxp, const std::string& criterion) {
    const int n = x.n_elem;

    if (n < maxp + 3) {
        return 0.0;  // Series too short
    }

    // Step 1: OLS detrend (allocation-optimized)
    arma::vec resid = ols_detrend_internal(x);

    // Step 2: Burg AR selection on residuals (pure C++ version)
    BurgResult ar_fit = burg_aic_select_pure(resid, maxp, criterion);

    // Steps 3-5: Fused single-pass (ZERO allocations)
    // Computes AR transform, time transform, and OLS t-stat all on-the-fly
    return co_tstat_fused(x, ar_fit.phi);
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
