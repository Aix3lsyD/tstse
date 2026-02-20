// =============================================================================
// FILE: kernel_wbg_boot.cpp
// CATEGORY: HOT PATH - Parallel Bootstrap Kernel
// THREAD-SAFE: YES (TBB workers use pure C++ internals)
//
// This is the main entry point for wbg_boot_fast(). Every microsecond matters.
// Uses RcppParallel (TBB) to run bootstrap iterations in parallel.
//
// Call chain: wbg_boot_fast.R -> wbg_bootstrap_kernel_grain_cpp()
//             -> WBGBootstrapWorker::operator()
//                 -> gen_ar_into_ws() + co_tstat_ws()
//
// ZERO allocations per bootstrap iteration (all workspace pre-allocated).
//
// Key exports:
//   - wbg_bootstrap_kernel_grain_cpp(): Main parallel bootstrap
//   - wbg_bootstrap_coba_kernel_grain_cpp(): COBA adjustment variant
// =============================================================================
// [[Rcpp::depends(RcppArmadillo, RcppParallel, dqrng, BH)]]

#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include "kernel_types.h"
#include <vector>
#include <cmath>
#include <algorithm>  // for std::fill

using namespace Rcpp;

// Forward declarations - thread-safe functions
arma::vec gen_ar_dqrng_boxmuller(int n, const arma::vec& phi, double vara, uint64_t rng_seed);
void gen_ar_into(arma::vec& out, const arma::vec& phi, double vara, uint64_t rng_seed);
void gen_ar_into_precomputed(arma::vec& out, const arma::vec& phi, double sd, int burn, uint64_t rng_seed);
void gen_ar_into_ws(arma::vec& out, const arma::vec& phi, double sd, int burn, uint64_t rng_seed, arma::vec& burn_buf);
int calc_ar_burnin_cpp(const arma::vec& phi, int n);
double co_tstat_pure(const arma::vec& x, int maxp, const std::string& criterion,
                     int min_p = 0);
BurgResult burg_aic_select_pure(const arma::vec& x, int maxp, const std::string& criterion, int min_p = 0);

// Forward declarations - workspace-aware functions (ZERO allocations)
double co_tstat_ws(const arma::vec& x, int maxp, CriterionType ic_type,
                   CoBootstrapWorkspace& ws, int min_p = 0);
BurgResult burg_aic_select_ws(const arma::vec& x, int maxp, CriterionType ic_type,
                               CoBootstrapWorkspace& ws, int min_p = 0);


// =============================================================================
// TBB Worker for WBG Bootstrap (Workspace-optimized)
// =============================================================================

struct WBGBootstrapWorker : public RcppParallel::Worker {
    // Input parameters (constant across all iterations)
    const int n;
    const arma::vec phi;
    const double vara;
    const int maxp;
    const CriterionType ic_type;  // Pre-converted enum (no string comparison in hot path)
    const int min_p;              // Minimum AR order for CO residual model
    const std::vector<uint64_t>& seeds;

    // Precomputed AR generation parameters (constant, computed once)
    const double sd;    // sqrt(vara) - precomputed to avoid per-iteration sqrt
    const int burn;     // burn-in length - precomputed to avoid per-iteration Armadillo ops

    // Output (thread-safe write via RVector)
    RcppParallel::RVector<double> results;

    // Constructor - precomputes sd and burn for efficiency
    WBGBootstrapWorker(int n_, const arma::vec& phi_, double vara_, int maxp_,
                       CriterionType ic_type_, int min_p_,
                       const std::vector<uint64_t>& seeds_,
                       Rcpp::NumericVector results_)
        : n(n_), phi(phi_), vara(vara_), maxp(maxp_), ic_type(ic_type_),
          min_p(min_p_), seeds(seeds_),
          sd(std::sqrt(vara_)),
          burn(phi_.n_elem == 0 ? 0 : calc_ar_burnin_cpp(phi_, n_)),
          results(results_) {}

    // Parallel worker function
    void operator()(std::size_t begin, std::size_t end) {
        // Thread-local workspace with lazy resize
        // Each thread allocates once, reuses across ALL chunks it processes
        // (TBB may call operator() multiple times for the same thread via work-stealing)
        static thread_local CoBootstrapWorkspace ws;
        static thread_local int ws_n = 0, ws_maxp = 0;

        if (ws_n != n || ws_maxp != maxp) {
            ws.resize(n, maxp);
            ws_n = n;
            ws_maxp = maxp;
        }

        for (std::size_t b = begin; b < end; ++b) {
            // Generate AR series using workspace buffer (ZERO allocations!)
            gen_ar_into_ws(ws.x_buf, phi, sd, burn, seeds[b], ws.burn_buf);

            // Compute CO t-statistic using workspace (ZERO allocations)
            results[b] = co_tstat_ws(ws.x_buf, maxp, ic_type, ws, min_p);
        }
    }
};


// =============================================================================
// Main Bootstrap Kernel
// =============================================================================

//' WBG Bootstrap Kernel (C++ Implementation)
//'
//' Runs the bootstrap loop in parallel using TBB with configurable grain size.
//'
//' @param n Integer, series length.
//' @param phi Numeric vector, AR coefficients from null model.
//' @param vara Double, innovation variance from null model.
//' @param seeds Vector of uint64 seeds, one per bootstrap iteration.
//' @param maxp Integer, maximum AR order for CO test.
//' @param criterion String, IC for AR selection: "aic", "aicc", "bic".
//' @param grain_size Integer, minimum iterations per thread (default 1).
//' @param min_p Integer, minimum AR order for CO residual model (default 0).
//' @return Numeric vector of bootstrap t-statistics.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector wbg_bootstrap_kernel_grain_cpp(
    int n,
    const arma::vec& phi,
    double vara,
    const std::vector<uint64_t>& seeds,
    int maxp = 5,
    std::string criterion = "aic",
    std::size_t grain_size = 1,
    int min_p = 0
) {
    // Validate AR order constraints (C++ buffer limit is MAX_P=20)
    if (maxp > 20) {
        Rcpp::stop("maxp must be <= 20 (C++ buffer limit)");
    }
    if (static_cast<int>(phi.n_elem) > 20) {
        Rcpp::stop("AR order p must be <= 20 (C++ buffer limit)");
    }

    const int nb = seeds.size();
    Rcpp::NumericVector results(nb);

    // Convert criterion string to enum ONCE
    CriterionType ic_type = criterion_from_string(criterion);

    WBGBootstrapWorker worker(n, phi, vara, maxp, ic_type, min_p, seeds, results);
    RcppParallel::parallelFor(0, nb, worker, grain_size);

    return results;
}


// =============================================================================
// TBB Worker for COBA Bootstrap (with AR fitting) - Workspace-optimized
// =============================================================================

struct WBGBootstrapCOBAWorker : public RcppParallel::Worker {
    // Input parameters (constant across all iterations)
    const int n;
    const arma::vec phi;
    const double vara;
    const int maxp;
    const CriterionType ic_type;  // Pre-converted enum
    const int min_p;              // Minimum AR order for CO residual model
    const std::vector<uint64_t>& seeds;

    // Precomputed AR generation parameters (constant, computed once)
    const double sd;    // sqrt(vara) - precomputed to avoid per-iteration sqrt
    const int burn;     // burn-in length - precomputed to avoid per-iteration Armadillo ops

    // Pre-allocated output (thread-safe writes to separate indices)
    RcppParallel::RVector<double> tstats;
    RcppParallel::RVector<double> phi1_values;
    RcppParallel::RMatrix<double> phi_matrix;  // nb × maxp
    RcppParallel::RVector<int> orders;

    // Constructor - precomputes sd and burn for efficiency
    WBGBootstrapCOBAWorker(int n_, const arma::vec& phi_, double vara_, int maxp_,
                           CriterionType ic_type_, int min_p_,
                           const std::vector<uint64_t>& seeds_,
                           Rcpp::NumericVector tstats_,
                           Rcpp::NumericVector phi1_values_,
                           Rcpp::NumericMatrix phi_matrix_,
                           Rcpp::IntegerVector orders_)
        : n(n_), phi(phi_), vara(vara_), maxp(maxp_), ic_type(ic_type_),
          min_p(min_p_), seeds(seeds_),
          sd(std::sqrt(vara_)),
          burn(phi_.n_elem == 0 ? 0 : calc_ar_burnin_cpp(phi_, n_)),
          tstats(tstats_), phi1_values(phi1_values_),
          phi_matrix(phi_matrix_), orders(orders_) {}

    // Parallel worker function
    void operator()(std::size_t begin, std::size_t end) {
        // Thread-local workspace with lazy resize
        // Each thread allocates once, reuses across ALL chunks it processes
        // (TBB may call operator() multiple times for the same thread via work-stealing)
        static thread_local CoBootstrapWorkspace ws;
        static thread_local int ws_n = 0, ws_maxp = 0;

        if (ws_n != n || ws_maxp != maxp) {
            ws.resize(n, maxp);
            ws_n = n;
            ws_maxp = maxp;
        }

        for (std::size_t b = begin; b < end; ++b) {
            // Generate AR series using workspace buffer (ZERO allocations!)
            gen_ar_into_ws(ws.x_buf, phi, sd, burn, seeds[b], ws.burn_buf);

            // Compute CO t-statistic using workspace (ZERO allocations)
            // Uses user's min_p for CO residual model selection
            tstats[b] = co_tstat_ws(ws.x_buf, maxp, ic_type, ws, min_p);

            // Fit AR model to bootstrap sample for COBA median selection
            // Uses min_p=0 here: bootstrap samples are generated under H0 (no trend),
            // so we want the unrestricted AR fit for estimating phi(1) bias
            BurgResult ar_fit = burg_aic_select_ws(ws.x_buf, maxp, ic_type, ws, 0);

            // Store results
            orders[b] = ar_fit.p;
            phi1_values[b] = std::abs(1.0 - arma::accu(ar_fit.phi));  // |φ(1)| per WBG paper Sec 2.2

            // Copy phi coefficients to matrix row (zero-padded)
            for (int j = 0; j < ar_fit.p && j < maxp; ++j) {
                phi_matrix(b, j) = ar_fit.phi(j);
            }
        }
    }
};


// =============================================================================
// COBA Bootstrap Kernel (for bootstrap adjustment)
// =============================================================================

//' WBG Bootstrap COBA Kernel (C++ Implementation)
//'
//' Runs the first-stage bootstrap for COBA adjustment with configurable grain size.
//'
//' @param n Integer, series length.
//' @param phi Numeric vector, AR coefficients from null model.
//' @param vara Double, innovation variance from null model.
//' @param seeds Vector of uint64 seeds, one per bootstrap iteration.
//' @param maxp Integer, maximum AR order for CO test and AR fitting.
//' @param criterion String, IC for AR selection: "aic", "aicc", "bic".
//' @param grain_size Integer, minimum iterations per thread (default 1).
//' @param min_p Integer, minimum AR order for CO residual model (default 0).
//' @return List with tstats, phi1_values, phi_matrix, orders.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Rcpp::List wbg_bootstrap_coba_kernel_grain_cpp(
    int n,
    const arma::vec& phi,
    double vara,
    const std::vector<uint64_t>& seeds,
    int maxp = 5,
    std::string criterion = "aic",
    std::size_t grain_size = 1,
    int min_p = 0
) {
    // Validate AR order constraints (C++ buffer limit is MAX_P=20)
    if (maxp > 20) {
        Rcpp::stop("maxp must be <= 20 (C++ buffer limit)");
    }
    if (static_cast<int>(phi.n_elem) > 20) {
        Rcpp::stop("AR order p must be <= 20 (C++ buffer limit)");
    }

    const int nb = seeds.size();

    // Pre-allocate outputs
    Rcpp::NumericVector tstats(nb);
    Rcpp::NumericVector phi1_values(nb);
    Rcpp::NumericMatrix phi_matrix(nb, maxp);
    std::fill(phi_matrix.begin(), phi_matrix.end(), 0.0);  // Explicit zero-init
    Rcpp::IntegerVector orders(nb);

    // Convert criterion string to enum ONCE
    CriterionType ic_type = criterion_from_string(criterion);

    WBGBootstrapCOBAWorker worker(n, phi, vara, maxp, ic_type, min_p, seeds,
                                   tstats, phi1_values, phi_matrix, orders);
    RcppParallel::parallelFor(0, nb, worker, grain_size);

    return Rcpp::List::create(
        Rcpp::Named("tstats") = tstats,
        Rcpp::Named("phi1_values") = phi1_values,
        Rcpp::Named("phi_matrix") = phi_matrix,
        Rcpp::Named("orders") = orders
    );
}
