// =============================================================================
// FILE: kernel_pw_tstat.cpp
// CATEGORY: HOT PATH - Pure C++ Prais-Winsten t-statistic
// THREAD-SAFE: YES (no Rcpp types in core functions)
//
// Contains pw_tstat_pure() and pw_tstat_ws() - the PW equivalent of CO.
// Uses single-pass fused computation for maximum performance.
//
// Key internal functions:
//   - pw_tstat_ws(): Workspace-aware, ZERO allocations per call
//   - pw_tstat_fused(): Single-pass transform + regression
//   - estimate_rho(): AR(1) coefficient from lag-1 autocorrelation
//
// Called by: kernel_wbg_boot_pw.cpp (WBGBootstrapPWWorker)
// =============================================================================
// [[Rcpp::depends(RcppArmadillo)]]

#include "kernel_types.h"
#include <cmath>
#include <algorithm>  // for std::max, std::min

// =============================================================================
// Internal helper: Estimate rho from lag-1 autocorrelation of residuals
// rho = Σ(z_t * z_{t-1}) / Σ(z_t²)
// =============================================================================
static double estimate_rho(const arma::vec& z) {
    const int n = z.n_elem;
    if (n < 2) return 0.0;

    double sum_zz = 0.0;
    double sum_z2 = 0.0;

    for (int i = 0; i < n; ++i) {
        sum_z2 += z[i] * z[i];
        if (i > 0) {
            sum_zz += z[i] * z[i-1];
        }
    }

    if (sum_z2 < 1e-15) return 0.0;

    double rho = sum_zz / sum_z2;
    // Clamp to avoid numerical issues with near-unit-root
    return std::max(-0.999, std::min(0.999, rho));
}

// =============================================================================
// Internal helper: OLS detrend into provided vector (ZERO allocations)
// Identical to kernel_co_tstat.cpp version
// =============================================================================
static void pw_ols_detrend_ws(const arma::vec& x, arma::vec& resid) {
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

    // Write directly to output vector
    for (int i = 0; i < n; ++i) {
        resid[i] = x[i] - a_hat - b_hat * (i + 1);
    }
}

// =============================================================================
// Fused PW t-statistic: single-pass, zero-allocation
// Computes PW transform and regression on-the-fly
// Uses algebraic identity: SS_res = SS_y - b² * SS_t (no residual storage!)
// =============================================================================
static double pw_tstat_fused(const arma::vec& x, double rho) {
    const int n = x.n_elem;
    if (n < 3) return 0.0;

    const double sqrt_term = std::sqrt(1.0 - rho * rho);

    // Single pass: accumulate all sums needed for regression
    double sum_y = 0.0, sum_t = 0.0;
    double sum_yt = 0.0, sum_y2 = 0.0, sum_t2 = 0.0;

    // First observation (quasi-differencing preserves it)
    // W_1 = sqrt(1-rho^2) * X_1
    // t*_1 = sqrt(1-rho^2) * 1
    double y_val = sqrt_term * x[0];
    double t_val = sqrt_term;

    sum_y += y_val;
    sum_t += t_val;
    sum_yt += y_val * t_val;
    sum_y2 += y_val * y_val;
    sum_t2 += t_val * t_val;

    // Remaining observations (standard PW differencing)
    // W_t = X_t - rho * X_{t-1}
    // t*_t = t - rho*(t-1) = t*(1-rho) + rho
    for (int i = 1; i < n; ++i) {
        y_val = x[i] - rho * x[i-1];
        t_val = (i + 1) * (1.0 - rho) + rho;

        sum_y += y_val;
        sum_t += t_val;
        sum_yt += y_val * t_val;
        sum_y2 += y_val * y_val;
        sum_t2 += t_val * t_val;
    }

    // Compute regression statistics from sums
    const double y_mean = sum_y / n;
    const double t_mean = sum_t / n;
    const double ss_t = sum_t2 - n * t_mean * t_mean;

    if (ss_t < 1e-15) return 0.0;

    const double sxy = sum_yt - n * y_mean * t_mean;
    const double b_hat = sxy / ss_t;

    // Key insight: SS_res = SS_y - b² * SS_t (algebraic identity)
    // This avoids storing residuals entirely!
    const double ss_y = sum_y2 - n * y_mean * y_mean;
    double ss_res = ss_y - b_hat * b_hat * ss_t;
    if (ss_res < 0) ss_res = 0;  // Numerical safety

    const double mse = ss_res / (n - 2);
    const double se_b = std::sqrt(mse / ss_t);

    if (se_b < 1e-15) return 0.0;

    return b_hat / se_b;
}

// =============================================================================
// Pure C++ PW t-statistic (thread-safe, allocates internally)
// =============================================================================
double pw_tstat_pure(const arma::vec& x) {
    const int n = x.n_elem;
    if (n < 5) return 0.0;

    // Step 1: OLS detrend (allocation for residuals)
    arma::vec resid(n);
    pw_ols_detrend_ws(x, resid);

    // Step 2: Estimate rho from lag-1 autocorrelation
    double rho = estimate_rho(resid);

    // Step 3: Fused PW transform + regression (ZERO allocations)
    return pw_tstat_fused(x, rho);
}

// =============================================================================
// Workspace-aware PW t-statistic (ZERO allocations in hot path)
// Uses pre-allocated workspace for all temporary vectors
// =============================================================================
double pw_tstat_ws(const arma::vec& x, PwBootstrapWorkspace& ws) {
    const int n = x.n_elem;
    if (n < 5) return 0.0;

    // Step 1: OLS detrend into workspace.resid (ZERO allocations)
    pw_ols_detrend_ws(x, ws.resid);

    // Step 2: Estimate rho from lag-1 autocorrelation
    double rho = estimate_rho(ws.resid);

    // Step 3: Fused PW transform + regression (already ZERO allocations)
    return pw_tstat_fused(x, rho);
}

// =============================================================================
// Pure C++ PW with full results
// =============================================================================
PWResult pw_full_pure(const arma::vec& x) {
    const int n = x.n_elem;
    if (n < 5) return PWResult();

    // Step 1: OLS detrend
    arma::vec resid(n);
    pw_ols_detrend_ws(x, resid);

    // Step 2: Estimate rho
    double rho = estimate_rho(resid);

    // Step 3: Compute t-statistic
    double tpw = pw_tstat_fused(x, rho);

    // Step 4: Compute residual variance for null model
    // var(X_t) = var(e_t) / (1 - rho^2) for AR(1)
    // So var(e_t) = var(resid) * (1 - rho^2)
    double vara = arma::var(resid) * (1.0 - rho * rho);

    return PWResult(tpw, rho, vara);
}

// =============================================================================
// Estimate rho and return it (for bootstrap COBA)
// Uses workspace residuals that were already computed by pw_tstat_ws
// =============================================================================
double pw_estimate_rho_from_ws(const PwBootstrapWorkspace& ws, int n) {
    return estimate_rho(ws.resid);
}

// =============================================================================
// PW regression with residual output (needed for iterative method)
// Returns coefficients and stores residuals in provided vector
// =============================================================================
static void pw_regression_with_resid(const arma::vec& x, double rho,
                                      double& b0_hat, double& b1_hat,
                                      arma::vec& resid) {
    const int n = x.n_elem;
    const double sqrt_term = std::sqrt(1.0 - rho * rho);

    // Compute PW-transformed values
    arma::vec y_pw(n);
    arma::vec t_pw(n);

    // First observation
    y_pw[0] = sqrt_term * x[0];
    t_pw[0] = sqrt_term;

    // Remaining observations
    for (int i = 1; i < n; ++i) {
        y_pw[i] = x[i] - rho * x[i-1];
        t_pw[i] = (i + 1) * (1.0 - rho) + rho;
    }

    // OLS on transformed data
    double sum_y = 0.0, sum_t = 0.0, sum_yt = 0.0, sum_t2 = 0.0;
    for (int i = 0; i < n; ++i) {
        sum_y += y_pw[i];
        sum_t += t_pw[i];
        sum_yt += y_pw[i] * t_pw[i];
        sum_t2 += t_pw[i] * t_pw[i];
    }

    const double y_mean = sum_y / n;
    const double t_mean = sum_t / n;
    const double ss_t = sum_t2 - n * t_mean * t_mean;

    if (ss_t < 1e-15) {
        b0_hat = y_mean;
        b1_hat = 0.0;
    } else {
        const double sxy = sum_yt - n * y_mean * t_mean;
        b1_hat = sxy / ss_t;
        b0_hat = y_mean - b1_hat * t_mean;
    }

    // Compute residuals from PW regression
    for (int i = 0; i < n; ++i) {
        resid[i] = y_pw[i] - b0_hat - b1_hat * t_pw[i];
    }
}

// =============================================================================
// Iterative Prais-Winsten with full results
// Iterates until rho converges or max_iter reached
// =============================================================================
PWResult pw_full_iterative_pure(const arma::vec& x, int max_iter, double tol) {
    const int n = x.n_elem;
    if (n < 5) return PWResult();

    // Step 1: Initial OLS detrend
    arma::vec resid(n);
    pw_ols_detrend_ws(x, resid);

    // Step 2: Initial rho estimate
    double rho = estimate_rho(resid);

    // Iteration variables
    double b0_hat = 0.0, b1_hat = 0.0;
    int iter = 0;

    // Iterative loop
    for (iter = 0; iter < max_iter; ++iter) {
        // PW regression with current rho
        pw_regression_with_resid(x, rho, b0_hat, b1_hat, resid);

        // Re-estimate rho from PW residuals
        double rho_new = estimate_rho(resid);

        // Check convergence
        if (std::abs(rho_new - rho) < tol) {
            rho = rho_new;
            break;
        }

        rho = rho_new;
    }

    // Final t-statistic computation
    double tpw = pw_tstat_fused(x, rho);

    // Compute residual variance for null model
    double vara = arma::var(resid) * (1.0 - rho * rho);

    return PWResult(tpw, rho, vara);
}

// =============================================================================
// Iterative PW t-statistic only (for fast computation without full results)
// =============================================================================
double pw_tstat_iterative_pure(const arma::vec& x, int max_iter, double tol) {
    PWResult result = pw_full_iterative_pure(x, max_iter, tol);
    return result.tpw;
}
