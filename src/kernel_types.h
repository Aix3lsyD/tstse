// =============================================================================
// FILE: kernel_types.h
// CATEGORY: HOT PATH - Core type definitions for parallel bootstrap
// THREAD-SAFE: YES (pure C++ structs, no Rcpp types)
//
// These structs are used by the TBB parallel bootstrap kernel.
// CoBootstrapWorkspace enables zero-allocation iterations.
//
// Used by: kernel_wbg_boot.cpp, kernel_co_tstat.cpp, kernel_burg_aic.cpp
// =============================================================================
#ifndef KERNEL_TYPES_H
#define KERNEL_TYPES_H

#include <RcppArmadillo.h>
#include <string>

// Result from Burg algorithm with AIC/BIC selection
struct BurgResult {
    arma::vec phi;      // AR coefficients
    double vara;        // Residual variance
    int p;              // Selected order
    double ic;          // Information criterion value

    BurgResult() : vara(0.0), p(0), ic(0.0) {}
    BurgResult(arma::vec phi_, double vara_, int p_, double ic_)
        : phi(phi_), vara(vara_), p(p_), ic(ic_) {}
};

// Result from Cochrane-Orcutt procedure
struct COResult {
    double tco;         // CO t-statistic
    int p;              // AR order selected
    arma::vec phi;      // AR coefficients
    double vara;        // Residual variance

    COResult() : tco(0.0), p(0), vara(0.0) {}
    COResult(double tco_, int p_, arma::vec phi_, double vara_)
        : tco(tco_), p(p_), phi(phi_), vara(vara_) {}
};

// Result from Prais-Winsten procedure
struct PWResult {
    double tpw;         // PW t-statistic
    double rho;         // AR(1) coefficient
    double vara;        // Residual variance

    PWResult() : tpw(0.0), rho(0.0), vara(0.0) {}
    PWResult(double tpw_, double rho_, double vara_)
        : tpw(tpw_), rho(rho_), vara(vara_) {}
};

// Result from OLS regression
struct OLSResult {
    arma::vec residuals;
    double intercept;
    double slope;
    double tstat;       // t-statistic for slope

    OLSResult() : intercept(0.0), slope(0.0), tstat(0.0) {}
};

// =============================================================================
// Information Criterion enum (avoids string comparison in hot path)
// =============================================================================
enum CriterionType {
    IC_AIC = 0,
    IC_AICC = 1,
    IC_BIC = 2
};

// Helper to convert string to enum (call once at entry point)
// NOTE: R wrapper should validate criterion before calling C++.
// Unknown values default to AIC (most common) for safety.
inline CriterionType criterion_from_string(const std::string& criterion) {
    if (criterion == "aic") return IC_AIC;
    if (criterion == "aicc") return IC_AICC;
    if (criterion == "bic") return IC_BIC;
    return IC_AIC;  // Default to AIC for unknown strings
}

// =============================================================================
// Thread-local workspace for bootstrap iterations
// Pre-allocate once per thread, reuse across all iterations
// =============================================================================
struct CoBootstrapWorkspace {
    // AR series buffer
    arma::vec x_buf;

    // OLS detrend workspace
    arma::vec resid;

    // Burg algorithm workspace
    arma::vec xc;       // Centered series
    arma::vec ef;       // Forward prediction errors
    arma::vec eb;       // Backward prediction errors
    arma::vec a_curr;   // Current AR coefficients
    arma::vec a_prev;   // Previous AR coefficients
    arma::vec best_phi; // Best AR coefficients (avoids allocation per IC improvement)
    int best_p;         // Order of best model (for reading best_phi)

    CoBootstrapWorkspace() : best_p(0) {}

    // Resize all vectors for given n and maxp
    void resize(int n, int maxp) {
        x_buf.set_size(n);
        resid.set_size(n);
        xc.set_size(n);
        ef.set_size(n);
        eb.set_size(n);
        a_curr.set_size(maxp);
        a_prev.set_size(maxp);
        best_phi.set_size(maxp);
        best_p = 0;
    }
};

// =============================================================================
// Thread-local workspace for PW bootstrap iterations
// Simpler than CO workspace (no variable-order AR, just AR(1))
// =============================================================================
struct PwBootstrapWorkspace {
    // AR series buffer
    arma::vec x_buf;

    // OLS detrend workspace
    arma::vec resid;

    PwBootstrapWorkspace() {}

    // Resize all vectors for given n
    void resize(int n) {
        x_buf.set_size(n);
        resid.set_size(n);
    }
};

#endif // KERNEL_TYPES_H
