// =============================================================================
// FILE: kernel_co_tstat.cpp
// CATEGORY: HOT PATH - Pure C++ Cochrane-Orcutt t-statistic
// THREAD-SAFE: YES (no Rcpp types in core functions)
//
// Contains co_tstat_ws() and co_tstat_fused() - the heart of the bootstrap.
// co_tstat_fused() uses O(1) time-transform optimization (not O(p) per obs).
//
// Key internal functions:
//   - co_tstat_ws(): Workspace-aware, ZERO allocations per call
//   - co_tstat_fused(): Single-pass with O(1) time-transform
//   - ols_detrend_internal(): Thread-safe detrending
//
// Called by: kernel_wbg_boot.cpp (WBGBootstrapWorker)
// =============================================================================
// [[Rcpp::depends(RcppArmadillo)]]

#include "kernel_types.h"
#include <cmath>

// Forward declaration of pure C++ functions
BurgResult burg_aic_select_pure(const arma::vec& x, int maxp,
                                 const std::string& criterion,
                                 int min_p = 0);

// Forward declaration of workspace-aware Burg function
BurgResult burg_aic_select_ws(const arma::vec& x, int maxp,
                               CriterionType ic_type,
                               CoBootstrapWorkspace& ws,
                               int min_p = 0);

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


// Workspace-aware OLS detrend (ZERO allocations)
// Writes detrended residuals directly into workspace.resid
static void ols_detrend_ws(const arma::vec& x, CoBootstrapWorkspace& ws) {
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
    const double sxy = sum_tx - n * x_mean * t_mean;
    const double b_hat = sxy / ss_t;
    const double a_hat = x_mean - b_hat * t_mean;

    // Write directly to workspace (no allocation!)
    for (int i = 0; i < n; ++i) {
        ws.resid[i] = x[i] - a_hat - b_hat * (i + 1);
    }
}


// OLS t-statistic for slope (pure C++)
// Non-static: also called by api_co_tstat.cpp
double ols_tstat_internal(const arma::vec& y, const arma::vec& t_idx) {
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

    // Precompute affine coefficients for time transform: t* = phi_at_1 * t + time_const
    // Mathematical identity: t* = t - Σφⱼ(t-j) = t(1-Σφⱼ) + Σ(j·φⱼ) = φ(1)·t + C
    // This reduces O(p) per observation to O(1)
    double phi_at_1 = 1.0;    // φ(1) = 1 - Σφⱼ
    double time_const = 0.0;  // C = Σ(j·φⱼ)
    for (int j = 0; j < p; ++j) {
        phi_at_1 -= phi[j];
        time_const += (j + 1) * phi[j];
    }

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

        // OPTIMIZED O(1) time transform using precomputed affine coefficients
        // t* = φ(1)·t + C where φ(1) = 1-Σφⱼ, C = Σ(j·φⱼ)
        const int t_1based = t + 1;
        double t_val = phi_at_1 * t_1based + time_const;

        // // ORIGINAL O(p) time transform - kept for reference
        // const int t_1based = t + 1;
        // double t_val = static_cast<double>(t_1based);
        // for (int j = 0; j < p; ++j) {
        //     t_val -= phi[j] * (t_1based - (j + 1));
        // }

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

    // Clamp maxp first (same logic as Burg) to allow short series to proceed
    if (maxp >= n - 1) {
        maxp = n - 2;
    }
    if (maxp < 1) {
        maxp = 1;
    }

    // Need at least 3 observations after removing p lags
    if (n < maxp + 3) {
        return 0.0;  // Series too short even after clamping
    }

    // Step 1: OLS detrend (allocation-optimized)
    arma::vec resid = ols_detrend_internal(x);

    // Step 2: Burg AR selection on residuals (pure C++ version)
    BurgResult ar_fit = burg_aic_select_pure(resid, maxp, criterion);

    // Steps 3-5: Fused single-pass (ZERO allocations)
    // Computes AR transform, time transform, and OLS t-stat all on-the-fly
    return co_tstat_fused(x, ar_fit.phi);
}


// =============================================================================
// Workspace-aware CO t-statistic (ZERO allocations in hot path)
// Uses pre-allocated workspace for all temporary vectors
// =============================================================================
double co_tstat_ws(const arma::vec& x, int maxp, CriterionType ic_type,
                   CoBootstrapWorkspace& ws) {
    const int n = x.n_elem;

    // Clamp maxp first (same logic as Burg) to allow short series to proceed
    if (maxp >= n - 1) {
        maxp = n - 2;
    }
    if (maxp < 1) {
        maxp = 1;
    }

    // Need at least 3 observations after removing p lags
    if (n < maxp + 3) {
        return 0.0;  // Series too short even after clamping
    }

    // Step 1: OLS detrend into workspace.resid (ZERO allocations)
    ols_detrend_ws(x, ws);

    // Step 2: Burg AR selection using workspace (ZERO allocations except best_phi copy)
    BurgResult ar_fit = burg_aic_select_ws(ws.resid, maxp, ic_type, ws);

    // Steps 3-5: Fused single-pass (already ZERO allocations)
    return co_tstat_fused(x, ar_fit.phi);
}


// =============================================================================
// R-facing export wrapper for kernel CO t-statistic
// Uses same implementation as bootstrap kernel for bootstrap validity
// =============================================================================

//' CO t-statistic (kernel implementation)
//'
//' Computes Cochrane-Orcutt t-statistic using the kernel implementation.
//' This function uses the exact same algorithm as the bootstrap kernel,
//' ensuring bootstrap validity when used in wbg_boot_fast().
//'
//' @param x Numeric vector, the time series.
//' @param maxp Integer, maximum AR order for model selection.
//' @param criterion String, information criterion: "aic", "aicc", or "bic".
//' @return Double, Cochrane-Orcutt t-statistic.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double co_tstat_pure_cpp(const arma::vec& x, int maxp = 5,
                          std::string criterion = "aic") {
    return co_tstat_pure(x, maxp, criterion);
}
