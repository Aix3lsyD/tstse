// =============================================================================
// FILE: kernel_co_tas.cpp
// CATEGORY: HOT PATH - CO-TAS (Turner Adjusted Sample Size) Kernel
// THREAD-SAFE: YES (pure C++ functions, no Rcpp types in hot path)
//
// Implements the Cochrane-Orcutt trend test with Turner's effective sample
// size adjustment. The algorithm differs from standard CO:
//   1. Difference series and reconstruct via cumsum (centered)
//   2. Fit AR(p) to differenced series via Burg
//   3. Transform original data and time index using AR filter
//   4. Compute effective sample size using theoretical ACF
//   5. Conduct t-test with adjusted degrees of freedom
//
// Key functions:
//   - ar_acf_pure(): Theoretical ACF for AR(p) via Yule-Walker recursion
//   - effective_n_from_acf(): Turner's effective sample size
//   - co_tas_pure(): Main CO-TAS computation
//   - co_tas_pvalue_ws(): Workspace-aware p-value for bootstrap (zero alloc)
//   - co_tas_tstat_ws(): Workspace-aware t-stat for bootstrap (zero alloc)
//
// Used by: kernel_co_tas_boot.cpp (parallel bootstrap), api_co_tas.cpp
// =============================================================================
// [[Rcpp::depends(RcppArmadillo, BH)]]

#include "kernel_types.h"
#include <cmath>
#include <limits>
#include <boost/math/distributions/students_t.hpp>

// =============================================================================
// Exact t-distribution p-value using Boost
// Replaces the previous approximation for accuracy in bootstrap comparisons
// =============================================================================
static double t_pvalue_exact(double t_value, double df) {
    if (df < 1.0) df = 1.0;
    double abs_t = std::fabs(t_value);
    boost::math::students_t dist(df);
    // Two-sided p-value: P(|T| > |t|) = 2 * P(T > |t|)
    return 2.0 * boost::math::cdf(boost::math::complement(dist, abs_t));
}

// Forward declaration of Burg kernel
BurgResult burg_aic_select_pure(const arma::vec& x, int maxp,
                                 const std::string& criterion, int min_p);

// Forward declaration of workspace-aware Burg kernel
BurgResult burg_aic_select_ws(const arma::vec& x, int maxp,
                               CriterionType ic_type,
                               CoBootstrapWorkspace& ws,
                               int min_p);

// =============================================================================
// Theoretical ACF for AR(p) Process
//
// Computes the theoretical autocorrelation function for an AR(p) process.
// Uses the autocovariance Yule-Walker system (matches stats::ARMAacf):
//   1. Solve (p+1)x(p+1) system for gamma(0), gamma(1), ..., gamma(p)
//   2. Normalize: rho(k) = gamma(k) / gamma(0)
//   3. Recurse for k > p: rho(k) = sum_{j=1}^{p} phi_j * rho(k-j)
//
// The autocovariance system (with sigma^2 = 1):
//   gamma(0) - phi_1*gamma(1) - ... - phi_p*gamma(p) = 1
//   gamma(k) - sum_j phi_j*gamma(|k-j|) = 0  for k = 1, ..., p
// =============================================================================
arma::vec ar_acf_pure(const arma::vec& phi, int lag_max) {
    int p = phi.n_elem;

    // Handle AR(0) case - white noise
    if (p == 0) {
        arma::vec acf(lag_max + 1, arma::fill::zeros);
        acf(0) = 1.0;
        return acf;
    }

    arma::vec acf(lag_max + 1);
    acf(0) = 1.0;  // rho(0) = 1 always

    if (lag_max == 0) {
        return acf;
    }

    // Build (p+1)x(p+1) system for autocovariances gamma(0), ..., gamma(p)
    arma::mat A(p + 1, p + 1, arma::fill::zeros);
    arma::vec b(p + 1, arma::fill::zeros);
    b(0) = 1.0;  // sigma^2 = 1

    // Row 0: gamma(0) - phi_1*gamma(1) - ... - phi_p*gamma(p) = 1
    A(0, 0) = 1.0;
    for (int j = 0; j < p; ++j) {
        A(0, j + 1) = -phi(j);
    }

    // Rows 1..p: gamma(k) = sum_j phi_j * gamma(|k-j|)
    // Rearranged: gamma(k) - sum_j phi_j * gamma(|k-j|) = 0
    for (int k = 1; k <= p; ++k) {
        for (int col = 0; col <= p; ++col) {
            double coef = (col == k) ? 1.0 : 0.0;
            // Subtract contributions from phi_j * gamma(|k-j|)
            for (int j = 0; j < p; ++j) {
                int lag = std::abs(k - (j + 1));  // |k - (j+1)| since phi is 0-indexed
                if (lag == col) {
                    coef -= phi(j);
                }
            }
            A(k, col) = coef;
        }
    }

    // Solve the system
    arma::vec gamma;
    bool solved = arma::solve(gamma, A, b);

    if (!solved || gamma(0) <= 0) {
        // Fallback for AR(1): rho(1) = phi_1
        acf(1) = phi(0);
        for (int k = 2; k <= lag_max; ++k) {
            double rho_k = 0.0;
            for (int j = 0; j < p && j < k; ++j) {
                rho_k += phi(j) * acf(k - j - 1);
            }
            acf(k) = rho_k;
        }
        return acf;
    }

    // Normalize to get ACF: rho(k) = gamma(k) / gamma(0)
    double gamma0 = gamma(0);
    for (int k = 1; k <= p && k <= lag_max; ++k) {
        acf(k) = gamma(k) / gamma0;
    }

    // Recurse for lags > p: rho(k) = sum_{j=1}^{p} phi_j * rho(k-j)
    for (int k = p + 1; k <= lag_max; ++k) {
        double rho_k = 0.0;
        for (int j = 0; j < p; ++j) {
            rho_k += phi(j) * acf(k - j - 1);
        }
        acf(k) = rho_k;
    }

    return acf;
}

// =============================================================================
// Effective Sample Size (Turner's Adjustment)
//
// n_a = (n-2) / A
// where A = 1 + 2 * sum_{k=1}^{n-3} (1 - k/(n-2)) * rho(k)
//
// This accounts for the reduced information content when observations
// are correlated.
// =============================================================================
double effective_n_from_acf(const arma::vec& acf, int n) {
    // acf should have at least n-2 lags (excluding lag 0)
    int max_k = n - 3;  // Number of lags to sum over
    if (max_k <= 0) {
        return static_cast<double>(n);
    }

    double A = 1.0;
    for (int k = 1; k <= max_k; ++k) {
        if (k < static_cast<int>(acf.n_elem)) {
            double weight = 1.0 - static_cast<double>(k) / (n - 2);
            A += 2.0 * weight * acf(k);
        }
    }

    // Ensure A is positive to avoid division issues
    if (A <= 0.0) {
        A = 1.0;
    }

    return (n - 2) / A;
}

// =============================================================================
// Streaming Effective Sample Size (Turner's Adjustment)
//
// Fuses ACF computation with Turner's sum to avoid O(n) allocation.
// Uses ring buffer of size p for ACF recursion.
// =============================================================================
double effective_n_streaming(const arma::vec& phi, int n) {
    int p = phi.n_elem;
    int max_k = n - 3;

    // Handle edge cases
    if (max_k <= 0) {
        return static_cast<double>(n);
    }

    // AR(0) case: all ρ(k) = 0 for k > 0, so A = 1
    if (p == 0) {
        return static_cast<double>(n - 2);
    }

    // Solve Yule-Walker for first p autocorrelations
    // Build (p+1)x(p+1) system (small allocation, p typically 1-5)
    arma::mat M(p + 1, p + 1, arma::fill::zeros);
    arma::vec b(p + 1, arma::fill::zeros);
    b(0) = 1.0;

    M(0, 0) = 1.0;
    for (int j = 0; j < p; ++j) {
        M(0, j + 1) = -phi(j);
    }

    for (int k = 1; k <= p; ++k) {
        for (int col = 0; col <= p; ++col) {
            double coef = (col == k) ? 1.0 : 0.0;
            for (int j = 0; j < p; ++j) {
                int lag = std::abs(k - (j + 1));
                if (lag == col) {
                    coef -= phi(j);
                }
            }
            M(k, col) = coef;
        }
    }

    arma::vec gamma;
    bool solved = arma::solve(gamma, M, b);

    // Initialize ring buffer with ρ(1), ..., ρ(p)
    // Using stack array for small p (max 20)
    double ring[20];
    double gamma0 = (solved && gamma(0) > 0) ? gamma(0) : 1.0;

    if (solved && gamma(0) > 0) {
        for (int j = 0; j < p; ++j) {
            ring[j] = gamma(j + 1) / gamma0;
        }
    } else {
        // Fallback: AR(1) approximation
        ring[0] = phi(0);
        for (int j = 1; j < p; ++j) {
            ring[j] = 0.0;
            for (int i = 0; i < j && i < p; ++i) {
                ring[j] += phi(i) * ring[j - i - 1];
            }
        }
    }

    // Accumulate Turner's sum streaming
    double A = 1.0;
    int ring_idx = 0;  // Points to oldest value in ring

    for (int k = 1; k <= max_k; ++k) {
        double rho_k;

        if (k <= p) {
            // Use precomputed values
            rho_k = ring[k - 1];
        } else {
            // Recursion: ρ(k) = Σⱼ φⱼ ρ(k-j)
            rho_k = 0.0;
            for (int j = 0; j < p; ++j) {
                // ring[(ring_idx + p - 1 - j) % p] = ρ(k-j-1)
                int idx = (ring_idx + p - 1 - j) % p;
                rho_k += phi(j) * ring[idx];
            }
            // Update ring buffer (overwrite oldest)
            ring[ring_idx] = rho_k;
            ring_idx = (ring_idx + 1) % p;
        }

        double weight = 1.0 - static_cast<double>(k) / (n - 2);
        A += 2.0 * weight * rho_k;
    }

    if (A <= 0.0) A = 1.0;
    return (n - 2) / A;
}

// =============================================================================
// AR Transform (Apply AR filter to series)
//
// y[t] = x[t] - sum_{j=1}^{p} phi_j * x[t-j]
//
// Returns transformed series of length n - p
// =============================================================================
static arma::vec ar_transform_internal(const arma::vec& x, const arma::vec& phi) {
    int n = x.n_elem;
    int p = phi.n_elem;

    if (p == 0 || p >= n) {
        return x;
    }

    int m = n - p;
    arma::vec y(m);

    for (int t = 0; t < m; ++t) {
        double val = x(t + p);
        for (int j = 0; j < p; ++j) {
            val -= phi(j) * x(t + p - j - 1);
        }
        y(t) = val;
    }

    return y;
}

// =============================================================================
// CO-TAS Result Structure
// =============================================================================
struct CoTasResult {
    double pvalue;
    double tco;
    double n_a;
    arma::vec phi;
    int p;

    CoTasResult() : pvalue(1.0), tco(0.0), n_a(0.0), p(0) {}
};

// =============================================================================
// Main CO-TAS Kernel
//
// Algorithm:
// 1. Difference series: x_diff = x[t] - x[t-1]
// 2. Reconstruct: z_diff = cumsum(x[1], x_diff - mean(x_diff))
// 3. Fit AR(p) to z_diff using Burg with IC selection
// 4. Transform data: x_trans = phi(B) * x
// 5. Transform time: t_trans = phi(B) * t
// 6. Regress x_trans on t_trans, get t-statistic
// 7. Compute effective sample size from theoretical ACF
// 8. Compute p-value with adjusted df
// =============================================================================
CoTasResult co_tas_pure(const arma::vec& x, int maxp, CriterionType ic_type) {
    CoTasResult result;
    int n = x.n_elem;

    if (n < 10) {
        result.pvalue = 1.0;
        result.tco = 0.0;
        result.n_a = n;
        return result;
    }

    // Steps 1 & 2: Difference and reconstruct in one pass
    // diff_mean = (x[n-1] - x[0]) / (n-1) via telescoping sum
    double diff_mean = (x(n - 1) - x(0)) / (n - 1);

    arma::vec z_diff(n);
    z_diff(0) = x(0);
    for (int t = 1; t < n; ++t) {
        z_diff(t) = z_diff(t - 1) + x(t) - x(t - 1) - diff_mean;
    }

    // Step 3: Fit AR model to differenced series
    std::string criterion;
    switch (ic_type) {
        case IC_AIC: criterion = "aic"; break;
        case IC_AICC: criterion = "aicc"; break;
        case IC_BIC: criterion = "bic"; break;
        default: criterion = "aic";
    }

    BurgResult ar_fit = burg_aic_select_pure(z_diff, maxp, criterion, 1);
    arma::vec phi = ar_fit.phi;
    int p = ar_fit.p;

    result.phi = phi;
    result.p = p;

    // Step 4: Transform data only (still need x_trans)
    arma::vec x_trans = ar_transform_internal(x, phi);

    int m = x_trans.n_elem;
    if (m < 3) {
        result.pvalue = 1.0;
        result.tco = 0.0;
        result.n_a = n;
        return result;
    }

    // Step 5: Compute affine coefficients for t_trans analytically
    // AR-filtering [1,2,...,n] produces affine: t_trans[t] = alpha * t + beta
    // This avoids allocating t_idx and t_trans arrays entirely
    double alpha = 1.0;
    double beta = static_cast<double>(p + 1);
    for (int j = 0; j < p; ++j) {
        alpha -= phi(j);
        beta -= phi(j) * (p - j);
    }

    // Step 6: Streaming OLS - single pass over x_trans
    // Accumulate sums where u = t_trans[t] = alpha*t + beta
    double sum_u = 0.0, sum_uu = 0.0;
    double sum_y = 0.0, sum_yy = 0.0;
    double sum_uy = 0.0;

    for (int t = 0; t < m; ++t) {
        double u = alpha * t + beta;
        double y = x_trans(t);
        sum_u += u;
        sum_uu += u * u;
        sum_y += y;
        sum_yy += y * y;
        sum_uy += u * y;
    }

    // OLS formulas using sums (no centering needed)
    // slope = (n*Σxy - Σx*Σy) / (n*Σx² - (Σx)²)
    double denom = m * sum_uu - sum_u * sum_u;  // This is m * SS_u
    if (denom < 1e-15) {
        result.pvalue = 1.0;
        result.tco = 0.0;
        result.n_a = n;
        return result;
    }

    double b_hat = (m * sum_uy - sum_u * sum_y) / denom;
    double a_hat = (sum_y - b_hat * sum_u) / m;

    // SSE = Σy² - a*Σy - b*Σxy
    double ss_res = sum_yy - a_hat * sum_y - b_hat * sum_uy;
    if (ss_res < 0) ss_res = 0;

    // Standard error of slope: se(b) = sqrt(MSE / SS_u)
    // SS_u = Σ(u-ū)² = sum_uu - sum_u²/m = denom/m
    double ss_u = denom / m;
    double mse = ss_res / (m - 2);
    double se_b = std::sqrt(mse / ss_u);

    // t-statistic
    double t_value = (se_b > 1e-15) ? b_hat / se_b : 0.0;
    result.tco = t_value;

    // Step 7: Compute effective sample size (streaming, no O(n) allocation)
    double n_a = effective_n_streaming(phi, n);
    result.n_a = n_a;

    // Step 8: Compute p-value with adjusted degrees of freedom
    double df = (n_a > p) ? n_a - p : n_a;
    if (df < 1) df = 1;

    // Exact two-sided p-value using Boost t-distribution
    result.pvalue = t_pvalue_exact(t_value, df);

    return result;
}

// =============================================================================
// Workspace-Aware CO-TAS Kernel (ZERO allocations in hot path)
//
// Same algorithm as co_tas_pure() but uses pre-allocated workspace buffers:
//   - ws.z_diff: Differenced/reconstructed series
//   - ws.x_trans: AR-transformed series
//   - ws.xc, ws.ef, ws.eb, etc: Burg algorithm buffers
//
// Returns both p-value and t-statistic via CoTasResult
// =============================================================================
static CoTasResult co_tas_ws_internal(const arma::vec& x, int maxp,
                                       CriterionType ic_type,
                                       CoBootstrapWorkspace& ws) {
    CoTasResult result;
    int n = x.n_elem;

    if (n < 10) {
        result.pvalue = 1.0;
        result.tco = 0.0;
        result.n_a = n;
        return result;
    }

    // Steps 1 & 2: Difference and reconstruct into ws.z_diff (ZERO allocations)
    double diff_mean = (x(n - 1) - x(0)) / (n - 1);

    ws.z_diff(0) = x(0);
    for (int t = 1; t < n; ++t) {
        ws.z_diff(t) = ws.z_diff(t - 1) + x(t) - x(t - 1) - diff_mean;
    }

    // Step 3: Fit AR model using workspace-aware Burg (ZERO allocations)
    BurgResult ar_fit = burg_aic_select_ws(ws.z_diff, maxp, ic_type, ws, 1);
    const arma::vec& phi = ar_fit.phi;
    int p = ar_fit.p;

    result.phi = phi;  // This copies, but it's small (p <= 20)
    result.p = p;

    // Step 4: Transform data into ws.x_trans (ZERO allocations)
    int m = n - p;
    if (m < 3) {
        result.pvalue = 1.0;
        result.tco = 0.0;
        result.n_a = n;
        return result;
    }

    for (int t = 0; t < m; ++t) {
        double val = x(t + p);
        for (int j = 0; j < p; ++j) {
            val -= phi(j) * x(t + p - j - 1);
        }
        ws.x_trans(t) = val;
    }

    // Step 5: Compute affine coefficients for t_trans analytically
    double alpha = 1.0;
    double beta = static_cast<double>(p + 1);
    for (int j = 0; j < p; ++j) {
        alpha -= phi(j);
        beta -= phi(j) * (p - j);
    }

    // Step 6: Streaming OLS - single pass over ws.x_trans
    double sum_u = 0.0, sum_uu = 0.0;
    double sum_y = 0.0, sum_yy = 0.0;
    double sum_uy = 0.0;

    for (int t = 0; t < m; ++t) {
        double u = alpha * t + beta;
        double y = ws.x_trans(t);
        sum_u += u;
        sum_uu += u * u;
        sum_y += y;
        sum_yy += y * y;
        sum_uy += u * y;
    }

    double denom = m * sum_uu - sum_u * sum_u;
    if (denom < 1e-15) {
        result.pvalue = 1.0;
        result.tco = 0.0;
        result.n_a = n;
        return result;
    }

    double b_hat = (m * sum_uy - sum_u * sum_y) / denom;
    double a_hat = (sum_y - b_hat * sum_u) / m;

    double ss_res = sum_yy - a_hat * sum_y - b_hat * sum_uy;
    if (ss_res < 0) ss_res = 0;

    double ss_u = denom / m;
    double mse = ss_res / (m - 2);
    double se_b = std::sqrt(mse / ss_u);

    double t_value = (se_b > 1e-15) ? b_hat / se_b : 0.0;
    result.tco = t_value;

    // Step 7: Compute effective sample size (streaming, minimal allocation)
    double n_a = effective_n_streaming(phi, n);
    result.n_a = n_a;

    // Step 8: Compute p-value with adjusted degrees of freedom
    double df = (n_a > p) ? n_a - p : n_a;
    if (df < 1) df = 1;

    // Exact two-sided p-value using Boost t-distribution
    result.pvalue = t_pvalue_exact(t_value, df);

    return result;
}

// =============================================================================
// Workspace-Aware CO-TAS P-value (For Bootstrap)
//
// ZERO allocations in hot path - uses pre-allocated workspace buffers
// =============================================================================
double co_tas_pvalue_ws(const arma::vec& x, int maxp, CriterionType ic_type,
                        CoBootstrapWorkspace& ws) {
    CoTasResult result = co_tas_ws_internal(x, maxp, ic_type, ws);
    return result.pvalue;
}

// =============================================================================
// Workspace-Aware CO-TAS t-statistic (For Bootstrap with btest=TRUE)
//
// ZERO allocations in hot path - uses pre-allocated workspace buffers
// =============================================================================
double co_tas_tstat_ws(const arma::vec& x, int maxp, CriterionType ic_type,
                       CoBootstrapWorkspace& ws) {
    CoTasResult result = co_tas_ws_internal(x, maxp, ic_type, ws);
    return result.tco;
}
