// =============================================================================
// FILE: kernel_co_tas_boot.cpp
// CATEGORY: HOT PATH - Parallel Bootstrap Kernel for CO-TAS
// THREAD-SAFE: YES (TBB workers use pure C++ internals)
//
// Parallel bootstrap for the Cochrane-Orcutt trend test with Turner Adjusted
// Sample Size (CO-TAS). Uses RcppParallel (TBB) for parallelization.
//
// Call chain: co_tas_boot_fast.R -> co_tas_boot_kernel_cpp()
//             -> CoTasBootstrapWorker::operator()
//                 -> gen_ar_into_ws() + co_tas_pvalue_ws()
//
// Key optimization: Uses CoBootstrapWorkspace for ZERO allocations per iteration.
// This eliminates ~3n heap allocations per bootstrap replicate from Burg algorithm.
//
// Key exports:
//   - co_tas_boot_kernel_cpp(): Parallel bootstrap returning p-values
//   - co_tas_boot_tstat_kernel_cpp(): Parallel bootstrap returning t-statistics
// =============================================================================
// [[Rcpp::depends(RcppArmadillo, RcppParallel, dqrng, BH)]]

#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include "kernel_types.h"
#include <vector>
#include <cmath>

using namespace Rcpp;

// Forward declarations - thread-safe functions from other kernels
void gen_ar_into_ws(arma::vec& out, const arma::vec& phi, double sd, int burn,
                    uint64_t rng_seed, arma::vec& burn_buf);
int calc_ar_burnin_cpp(const arma::vec& phi, int n);

// Forward declarations - Workspace-aware CO-TAS computation (ZERO allocations)
double co_tas_pvalue_ws(const arma::vec& x, int maxp, CriterionType ic_type,
                        CoBootstrapWorkspace& ws);
double co_tas_tstat_ws(const arma::vec& x, int maxp, CriterionType ic_type,
                       CoBootstrapWorkspace& ws);

// =============================================================================
// TBB Worker for CO-TAS Bootstrap
//
// Each thread generates AR series under null (no trend) and computes CO-TAS
// p-values. Uses thread-local CoBootstrapWorkspace for ZERO allocations.
// =============================================================================

struct CoTasBootstrapWorker : public RcppParallel::Worker {
    // Input parameters (constant across all iterations)
    const int n;
    const arma::vec phi;
    const double vara;
    const int maxp;
    const CriterionType ic_type;
    const std::vector<uint64_t>& seeds;

    // Precomputed AR generation parameters
    const double sd;    // sqrt(vara) - computed once
    const int burn;     // burn-in length - computed once

    // Output (thread-safe write via RVector)
    RcppParallel::RVector<double> pvalues;

    // Constructor
    CoTasBootstrapWorker(int n_, const arma::vec& phi_, double vara_, int maxp_,
                         CriterionType ic_type_,
                         const std::vector<uint64_t>& seeds_,
                         Rcpp::NumericVector pvalues_)
        : n(n_), phi(phi_), vara(vara_), maxp(maxp_), ic_type(ic_type_),
          seeds(seeds_),
          sd(std::sqrt(vara_)),
          burn(phi_.n_elem == 0 ? 0 : calc_ar_burnin_cpp(phi_, n_)),
          pvalues(pvalues_) {}

    // Parallel worker function
    void operator()(std::size_t begin, std::size_t end) {
        // Thread-local workspace for ZERO allocations per iteration
        static thread_local CoBootstrapWorkspace ws;
        static thread_local int ws_n = 0;
        static thread_local int ws_maxp = 0;

        // Resize if needed (once per thread per unique n/maxp)
        if (ws_n != n || ws_maxp != maxp) {
            ws.resize(n, maxp);
            ws_n = n;
            ws_maxp = maxp;
        }

        for (std::size_t b = begin; b < end; ++b) {
            // Generate AR series under null hypothesis (no trend)
            gen_ar_into_ws(ws.x_buf, phi, sd, burn, seeds[b], ws.burn_buf);

            // Compute CO-TAS p-value using workspace (ZERO allocations)
            pvalues[b] = co_tas_pvalue_ws(ws.x_buf, maxp, ic_type, ws);
        }
    }
};

// =============================================================================
// Main CO-TAS Bootstrap Kernel
// =============================================================================

//' CO-TAS Bootstrap Kernel (C++ Implementation)
//'
//' Runs the bootstrap loop in parallel using TBB. Each iteration generates an
//' AR series under the null hypothesis (no trend) and computes the CO-TAS
//' p-value.
//'
//' @param n Integer, series length.
//' @param phi Numeric vector, AR coefficients from null model.
//' @param vara Double, innovation variance from null model.
//' @param seeds Vector of uint64 seeds, one per bootstrap iteration.
//' @param maxp Integer, maximum AR order for CO-TAS test.
//' @param type String, IC for AR selection: "aic", "aicc", "bic".
//' @param grain_size Integer, minimum iterations per thread (default 1).
//' @return Numeric vector of bootstrap p-values.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector co_tas_boot_kernel_cpp(
    int n,
    const arma::vec& phi,
    double vara,
    const std::vector<uint64_t>& seeds,
    int maxp = 5,
    std::string type = "aic",
    std::size_t grain_size = 1
) {
    // Validate constraints
    if (maxp > 20) {
        Rcpp::stop("maxp must be <= 20");
    }
    if (static_cast<int>(phi.n_elem) > 20) {
        Rcpp::stop("AR order p must be <= 20");
    }
    if (n < 10) {
        Rcpp::stop("Series length must be at least 10");
    }

    const int nb = seeds.size();
    Rcpp::NumericVector pvalues(nb);

    // Convert type string to enum once
    CriterionType ic_type = criterion_from_string(type);

    CoTasBootstrapWorker worker(n, phi, vara, maxp, ic_type, seeds, pvalues);
    RcppParallel::parallelFor(0, nb, worker, grain_size);

    return pvalues;
}

// =============================================================================
// TBB Worker for CO-TAS t-statistic Bootstrap (btest=TRUE)
//
// Each thread generates AR series under null (no trend) and computes CO-TAS
// t-statistics. Uses thread-local CoBootstrapWorkspace for ZERO allocations.
// =============================================================================

struct CoTasTstatWorker : public RcppParallel::Worker {
    // Input parameters (constant across all iterations)
    const int n;
    const arma::vec phi;
    const double vara;
    const int maxp;
    const CriterionType ic_type;
    const std::vector<uint64_t>& seeds;

    // Precomputed AR generation parameters
    const double sd;    // sqrt(vara) - computed once
    const int burn;     // burn-in length - computed once

    // Output (thread-safe write via RVector)
    RcppParallel::RVector<double> tstats;

    // Constructor
    CoTasTstatWorker(int n_, const arma::vec& phi_, double vara_, int maxp_,
                     CriterionType ic_type_,
                     const std::vector<uint64_t>& seeds_,
                     Rcpp::NumericVector tstats_)
        : n(n_), phi(phi_), vara(vara_), maxp(maxp_), ic_type(ic_type_),
          seeds(seeds_),
          sd(std::sqrt(vara_)),
          burn(phi_.n_elem == 0 ? 0 : calc_ar_burnin_cpp(phi_, n_)),
          tstats(tstats_) {}

    // Parallel worker function
    void operator()(std::size_t begin, std::size_t end) {
        // Thread-local workspace for ZERO allocations per iteration
        static thread_local CoBootstrapWorkspace ws;
        static thread_local int ws_n = 0;
        static thread_local int ws_maxp = 0;

        // Resize if needed (once per thread per unique n/maxp)
        if (ws_n != n || ws_maxp != maxp) {
            ws.resize(n, maxp);
            ws_n = n;
            ws_maxp = maxp;
        }

        for (std::size_t b = begin; b < end; ++b) {
            // Generate AR series under null hypothesis (no trend)
            gen_ar_into_ws(ws.x_buf, phi, sd, burn, seeds[b], ws.burn_buf);

            // Compute CO-TAS t-statistic using workspace (ZERO allocations)
            tstats[b] = co_tas_tstat_ws(ws.x_buf, maxp, ic_type, ws);
        }
    }
};

// =============================================================================
// CO-TAS t-statistic Bootstrap Kernel (for btest=TRUE)
// =============================================================================

//' CO-TAS t-statistic Bootstrap Kernel (C++ Implementation)
//'
//' Runs the bootstrap loop in parallel using TBB. Each iteration generates an
//' AR series under the null hypothesis (no trend) and computes the CO-TAS
//' t-statistic. Used when btest=TRUE.
//'
//' @param n Integer, series length.
//' @param phi Numeric vector, AR coefficients from null model.
//' @param vara Double, innovation variance from null model.
//' @param seeds Vector of uint64 seeds, one per bootstrap iteration.
//' @param maxp Integer, maximum AR order for CO-TAS test.
//' @param type String, IC for AR selection: "aic", "aicc", "bic".
//' @param grain_size Integer, minimum iterations per thread (default 1).
//' @return Numeric vector of bootstrap t-statistics.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector co_tas_boot_tstat_kernel_cpp(
    int n,
    const arma::vec& phi,
    double vara,
    const std::vector<uint64_t>& seeds,
    int maxp = 5,
    std::string type = "aic",
    std::size_t grain_size = 1
) {
    // Validate constraints
    if (maxp > 20) {
        Rcpp::stop("maxp must be <= 20");
    }
    if (static_cast<int>(phi.n_elem) > 20) {
        Rcpp::stop("AR order p must be <= 20");
    }
    if (n < 10) {
        Rcpp::stop("Series length must be at least 10");
    }

    const int nb = seeds.size();
    Rcpp::NumericVector tstats(nb);

    // Convert type string to enum once
    CriterionType ic_type = criterion_from_string(type);

    CoTasTstatWorker worker(n, phi, vara, maxp, ic_type, seeds, tstats);
    RcppParallel::parallelFor(0, nb, worker, grain_size);

    return tstats;
}
