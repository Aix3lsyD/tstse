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

#endif // TYPES_PURE_H
