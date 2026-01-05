// types_pure.h - Pure C++ types for thread-safe parallel execution
// These structs avoid Rcpp types which aren't thread-safe with OpenMP
#ifndef TYPES_PURE_H
#define TYPES_PURE_H

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
inline CriterionType criterion_from_string(const std::string& criterion) {
    if (criterion == "aic") return IC_AIC;
    if (criterion == "aicc") return IC_AICC;
    return IC_BIC;
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

    CoBootstrapWorkspace() = default;

    // Resize all vectors for given n and maxp
    void resize(int n, int maxp) {
        x_buf.set_size(n);
        resid.set_size(n);
        xc.set_size(n);
        ef.set_size(n);
        eb.set_size(n);
        a_curr.set_size(maxp);
        a_prev.set_size(maxp);
    }
};

#endif // TYPES_PURE_H
