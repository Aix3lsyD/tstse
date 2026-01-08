// =============================================================================
// FILE: kernel_wbg_boot.cpp
// CATEGORY: HOT PATH - Parallel Bootstrap Kernel
// THREAD-SAFE: YES (TBB workers use pure C++ internals)
//
// This is the main entry point for wbg_boot_fast(). Every microsecond matters.
// Uses RcppParallel (TBB) to run bootstrap iterations in parallel.
//
// Call chain: wbg_boot_fast.R -> wbg_bootstrap_kernel_cpp()
//             -> WBGBootstrapWorker::operator()
//                 -> gen_ar_into() + co_tstat_ws()
//
// Key exports:
//   - wbg_bootstrap_kernel_cpp(): Main parallel bootstrap
//   - wbg_bootstrap_coba_kernel_cpp(): COBA adjustment variant
// =============================================================================
// [[Rcpp::depends(RcppArmadillo, RcppParallel, dqrng, BH)]]

#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <dqrng.h>
#include <xoshiro.h>
#include "kernel_types.h"
#include <vector>
#include <cmath>
#include <random>

// M_PI fallback for portability (not defined on MSVC without _USE_MATH_DEFINES)
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace Rcpp;

// Forward declarations - thread-safe functions
arma::vec gen_ar_dqrng_boxmuller(int n, const arma::vec& phi, double vara, uint64_t rng_seed);
void gen_ar_into(arma::vec& out, const arma::vec& phi, double vara, uint64_t rng_seed);
void gen_ar_into_precomputed(arma::vec& out, const arma::vec& phi, double sd, int burn, uint64_t rng_seed);
int calc_ar_burnin_cpp(const arma::vec& phi, int n);
double co_tstat_pure(const arma::vec& x, int maxp, const std::string& criterion);
BurgResult burg_aic_select_pure(const arma::vec& x, int maxp, const std::string& criterion, int min_p = 0);

// Forward declarations - workspace-aware functions (ZERO allocations)
double co_tstat_ws(const arma::vec& x, int maxp, CriterionType ic_type,
                   CoBootstrapWorkspace& ws);
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
    const std::vector<uint64_t>& seeds;

    // Precomputed AR generation parameters (constant, computed once)
    const double sd;    // sqrt(vara) - precomputed to avoid per-iteration sqrt
    const int burn;     // burn-in length - precomputed to avoid per-iteration Armadillo ops

    // Output (thread-safe write via RVector)
    RcppParallel::RVector<double> results;

    // Constructor - precomputes sd and burn for efficiency
    WBGBootstrapWorker(int n_, const arma::vec& phi_, double vara_, int maxp_,
                       CriterionType ic_type_,
                       const std::vector<uint64_t>& seeds_,
                       Rcpp::NumericVector results_)
        : n(n_), phi(phi_), vara(vara_), maxp(maxp_), ic_type(ic_type_),
          seeds(seeds_),
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
            // Generate AR series using precomputed sd and burn (no sqrt/burnin calc per iteration)
            gen_ar_into_precomputed(ws.x_buf, phi, sd, burn, seeds[b]);

            // Compute CO t-statistic using workspace (ZERO allocations)
            results[b] = co_tstat_ws(ws.x_buf, maxp, ic_type, ws);
        }
    }
};


// =============================================================================
// Main Bootstrap Kernel
// =============================================================================

//' WBG Bootstrap Kernel (C++ Implementation)
//'
//' Runs the bootstrap loop in parallel using TBB.
//' Each iteration generates an AR series under the null hypothesis
//' (no trend) and computes the Cochrane-Orcutt t-statistic.
//'
//' @param n Integer, series length.
//' @param phi Numeric vector, AR coefficients from null model.
//' @param vara Double, innovation variance from null model.
//' @param seeds Vector of uint64 seeds, one per bootstrap iteration.
//' @param maxp Integer, maximum AR order for CO test.
//' @param criterion String, IC for AR selection: "aic", "aicc", "bic".
//' @return Numeric vector of bootstrap t-statistics.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector wbg_bootstrap_kernel_cpp(
    int n,
    const arma::vec& phi,
    double vara,
    const std::vector<uint64_t>& seeds,
    int maxp = 5,
    std::string criterion = "aic"
) {
    const int nb = seeds.size();
    Rcpp::NumericVector results(nb);

    // Convert criterion string to enum ONCE (not in hot path)
    CriterionType ic_type = criterion_from_string(criterion);

    // Create worker and run in parallel
    WBGBootstrapWorker worker(n, phi, vara, maxp, ic_type, seeds, results);
    RcppParallel::parallelFor(0, nb, worker);

    return results;
}


//' WBG Bootstrap Kernel with Grain Size Control
//'
//' Same as wbg_bootstrap_kernel_cpp but with explicit grain size
//' for performance tuning.
//'
//' @param n Integer, series length.
//' @param phi Numeric vector, AR coefficients from null model.
//' @param vara Double, innovation variance from null model.
//' @param seeds Vector of uint64 seeds, one per bootstrap iteration.
//' @param maxp Integer, maximum AR order for CO test.
//' @param criterion String, IC for AR selection: "aic", "aicc", "bic".
//' @param grain_size Integer, minimum iterations per thread (default 1).
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
    std::size_t grain_size = 1
) {
    const int nb = seeds.size();
    Rcpp::NumericVector results(nb);

    // Convert criterion string to enum ONCE
    CriterionType ic_type = criterion_from_string(criterion);

    WBGBootstrapWorker worker(n, phi, vara, maxp, ic_type, seeds, results);
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
                           CriterionType ic_type_,
                           const std::vector<uint64_t>& seeds_,
                           Rcpp::NumericVector tstats_,
                           Rcpp::NumericVector phi1_values_,
                           Rcpp::NumericMatrix phi_matrix_,
                           Rcpp::IntegerVector orders_)
        : n(n_), phi(phi_), vara(vara_), maxp(maxp_), ic_type(ic_type_),
          seeds(seeds_),
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
            // Generate AR series using precomputed sd and burn (no sqrt/burnin calc per iteration)
            gen_ar_into_precomputed(ws.x_buf, phi, sd, burn, seeds[b]);

            // Compute CO t-statistic using workspace (ZERO allocations)
            tstats[b] = co_tstat_ws(ws.x_buf, maxp, ic_type, ws);

            // Fit AR model to bootstrap sample using workspace
            // Note: This reuses ws.xc, ws.ef, ws.eb, ws.a_curr, ws.a_prev
            BurgResult ar_fit = burg_aic_select_ws(ws.x_buf, maxp, ic_type, ws);

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

//' WBG Bootstrap Kernel with COBA (C++ Implementation)
//'
//' Runs the first-stage bootstrap for COBA adjustment.
//' Each iteration generates an AR series, computes the CO t-statistic,
//' AND fits an AR model to get phi coefficients for median model selection.
//'
//' @param n Integer, series length.
//' @param phi Numeric vector, AR coefficients from null model.
//' @param vara Double, innovation variance from null model.
//' @param seeds Vector of uint64 seeds, one per bootstrap iteration.
//' @param maxp Integer, maximum AR order for CO test and AR fitting.
//' @param criterion String, IC for AR selection: "aic", "aicc", "bic".
//' @return List with:
//'   - tstats: Numeric vector of bootstrap t-statistics
//'   - phi1_values: Numeric vector of phi(1) = 1 - sum(phi) for each bootstrap
//'   - phi_matrix: Matrix (nb x maxp) of AR coefficients (zero-padded)
//'   - orders: Integer vector of selected AR orders
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Rcpp::List wbg_bootstrap_coba_kernel_cpp(
    int n,
    const arma::vec& phi,
    double vara,
    const std::vector<uint64_t>& seeds,
    int maxp = 5,
    std::string criterion = "aic"
) {
    const int nb = seeds.size();

    // Pre-allocate outputs
    Rcpp::NumericVector tstats(nb);
    Rcpp::NumericVector phi1_values(nb);
    Rcpp::NumericMatrix phi_matrix(nb, maxp);  // Zero-initialized
    Rcpp::IntegerVector orders(nb);

    // Convert criterion string to enum ONCE
    CriterionType ic_type = criterion_from_string(criterion);

    // Create worker and run in parallel
    WBGBootstrapCOBAWorker worker(n, phi, vara, maxp, ic_type, seeds,
                                   tstats, phi1_values, phi_matrix, orders);
    RcppParallel::parallelFor(0, nb, worker);

    return Rcpp::List::create(
        Rcpp::Named("tstats") = tstats,
        Rcpp::Named("phi1_values") = phi1_values,
        Rcpp::Named("phi_matrix") = phi_matrix,
        Rcpp::Named("orders") = orders
    );
}


//' WBG Bootstrap COBA Kernel with Grain Size Control
//'
//' Same as wbg_bootstrap_coba_kernel_cpp but with explicit grain size.
//'
//' @param n Integer, series length.
//' @param phi Numeric vector, AR coefficients from null model.
//' @param vara Double, innovation variance from null model.
//' @param seeds Vector of uint64 seeds, one per bootstrap iteration.
//' @param maxp Integer, maximum AR order for CO test and AR fitting.
//' @param criterion String, IC for AR selection: "aic", "aicc", "bic".
//' @param grain_size Integer, minimum iterations per thread (default 1).
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
    std::size_t grain_size = 1
) {
    const int nb = seeds.size();

    // Pre-allocate outputs
    Rcpp::NumericVector tstats(nb);
    Rcpp::NumericVector phi1_values(nb);
    Rcpp::NumericMatrix phi_matrix(nb, maxp);
    Rcpp::IntegerVector orders(nb);

    // Convert criterion string to enum ONCE
    CriterionType ic_type = criterion_from_string(criterion);

    WBGBootstrapCOBAWorker worker(n, phi, vara, maxp, ic_type, seeds,
                                   tstats, phi1_values, phi_matrix, orders);
    RcppParallel::parallelFor(0, nb, worker, grain_size);

    return Rcpp::List::create(
        Rcpp::Named("tstats") = tstats,
        Rcpp::Named("phi1_values") = phi1_values,
        Rcpp::Named("phi_matrix") = phi_matrix,
        Rcpp::Named("orders") = orders
    );
}


// =============================================================================
// Combined CO + PW Bootstrap Kernels
// Compute both CO and PW t-statistics in a single parallel loop
// =============================================================================

// Forward declarations for PW functions (from kernel_pw_tstat.cpp)
double pw_tstat_ws(const arma::vec& x, PwBootstrapWorkspace& ws);
double pw_estimate_rho_from_ws(const PwBootstrapWorkspace& ws, int n);

// TBB Worker for Combined CO + PW Bootstrap
struct WBGBootstrapCOPWWorker : public RcppParallel::Worker {
    // Input parameters
    const int n;
    const arma::vec phi;
    const double vara;
    const int maxp_co;
    const int maxp_pw;
    const CriterionType ic_type;
    const std::vector<uint64_t>& seeds;

    // Precomputed
    const double sd;
    const int burn;

    // Output
    RcppParallel::RVector<double> co_results;
    RcppParallel::RVector<double> pw_results;

    WBGBootstrapCOPWWorker(int n_, const arma::vec& phi_, double vara_,
                           int maxp_co_, int maxp_pw_, CriterionType ic_type_,
                           const std::vector<uint64_t>& seeds_,
                           Rcpp::NumericVector co_results_,
                           Rcpp::NumericVector pw_results_)
        : n(n_), phi(phi_), vara(vara_), maxp_co(maxp_co_), maxp_pw(maxp_pw_),
          ic_type(ic_type_), seeds(seeds_),
          sd(std::sqrt(vara_)),
          burn(phi_.n_elem == 0 ? 0 : calc_ar_burnin_cpp(phi_, n_)),
          co_results(co_results_), pw_results(pw_results_) {}

    void operator()(std::size_t begin, std::size_t end) {
        // Thread-local workspaces
        static thread_local CoBootstrapWorkspace co_ws;
        static thread_local PwBootstrapWorkspace pw_ws;
        static thread_local int ws_n = 0, ws_maxp = 0;

        if (ws_n != n || ws_maxp != maxp_co) {
            co_ws.resize(n, maxp_co);
            pw_ws.resize(n);
            ws_n = n;
            ws_maxp = maxp_co;
        }

        for (std::size_t b = begin; b < end; ++b) {
            // Generate AR series once, use for both statistics
            gen_ar_into_precomputed(co_ws.x_buf, phi, sd, burn, seeds[b]);

            // Copy to PW workspace (they need separate buffers for residuals)
            pw_ws.x_buf = co_ws.x_buf;

            // Compute CO t-statistic
            co_results[b] = co_tstat_ws(co_ws.x_buf, maxp_co, ic_type, co_ws);

            // Compute PW t-statistic
            pw_results[b] = pw_tstat_ws(pw_ws.x_buf, pw_ws);
        }
    }
};


//' Combined CO + PW Bootstrap Kernel
//'
//' Computes both CO and PW bootstrap t-statistics in parallel.
//'
//' @param n Integer, series length.
//' @param phi Numeric vector, AR coefficients from null model.
//' @param vara Double, innovation variance from null model.
//' @param seeds Vector of uint64 seeds, one per bootstrap iteration.
//' @param maxp_co Integer, maximum AR order for CO test.
//' @param maxp_pw Integer, maximum AR order for PW (typically 1).
//' @param criterion String, IC for AR selection.
//' @param grain_size Integer, minimum iterations per thread.
//' @return List with co_tstats and pw_tstats.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Rcpp::List wbg_bootstrap_copw_kernel_cpp(
    int n,
    const arma::vec& phi,
    double vara,
    const std::vector<uint64_t>& seeds,
    int maxp_co = 5,
    int maxp_pw = 1,
    std::string criterion = "aic",
    std::size_t grain_size = 1
) {
    const int nb = seeds.size();
    Rcpp::NumericVector co_tstats(nb);
    Rcpp::NumericVector pw_tstats(nb);

    CriterionType ic_type = criterion_from_string(criterion);

    WBGBootstrapCOPWWorker worker(n, phi, vara, maxp_co, maxp_pw, ic_type, seeds,
                                   co_tstats, pw_tstats);
    RcppParallel::parallelFor(0, nb, worker, grain_size);

    return Rcpp::List::create(
        Rcpp::Named("co_tstats") = co_tstats,
        Rcpp::Named("pw_tstats") = pw_tstats
    );
}


// =============================================================================
// Combined CO + PW COBA Bootstrap Worker
// =============================================================================

struct WBGBootstrapCOPWCOBAWorker : public RcppParallel::Worker {
    const int n;
    const arma::vec phi;
    const double vara;
    const int maxp_co;
    const int maxp_pw;
    const CriterionType ic_type;
    const std::vector<uint64_t>& seeds;

    const double sd;
    const int burn;

    // CO outputs
    RcppParallel::RVector<double> co_tstats;
    RcppParallel::RVector<double> phi1_values;
    RcppParallel::RMatrix<double> phi_matrix;
    RcppParallel::RVector<int> orders;

    // PW outputs
    RcppParallel::RVector<double> pw_tstats;
    RcppParallel::RVector<double> rho_values;

    WBGBootstrapCOPWCOBAWorker(int n_, const arma::vec& phi_, double vara_,
                                int maxp_co_, int maxp_pw_, CriterionType ic_type_,
                                const std::vector<uint64_t>& seeds_,
                                Rcpp::NumericVector co_tstats_,
                                Rcpp::NumericVector phi1_values_,
                                Rcpp::NumericMatrix phi_matrix_,
                                Rcpp::IntegerVector orders_,
                                Rcpp::NumericVector pw_tstats_,
                                Rcpp::NumericVector rho_values_)
        : n(n_), phi(phi_), vara(vara_), maxp_co(maxp_co_), maxp_pw(maxp_pw_),
          ic_type(ic_type_), seeds(seeds_),
          sd(std::sqrt(vara_)),
          burn(phi_.n_elem == 0 ? 0 : calc_ar_burnin_cpp(phi_, n_)),
          co_tstats(co_tstats_), phi1_values(phi1_values_),
          phi_matrix(phi_matrix_), orders(orders_),
          pw_tstats(pw_tstats_), rho_values(rho_values_) {}

    void operator()(std::size_t begin, std::size_t end) {
        static thread_local CoBootstrapWorkspace co_ws;
        static thread_local PwBootstrapWorkspace pw_ws;
        static thread_local int ws_n = 0, ws_maxp = 0;

        if (ws_n != n || ws_maxp != maxp_co) {
            co_ws.resize(n, maxp_co);
            pw_ws.resize(n);
            ws_n = n;
            ws_maxp = maxp_co;
        }

        for (std::size_t b = begin; b < end; ++b) {
            // Generate AR series once
            gen_ar_into_precomputed(co_ws.x_buf, phi, sd, burn, seeds[b]);
            pw_ws.x_buf = co_ws.x_buf;

            // CO t-statistic
            co_tstats[b] = co_tstat_ws(co_ws.x_buf, maxp_co, ic_type, co_ws);

            // CO AR fit for COBA
            BurgResult ar_fit = burg_aic_select_ws(co_ws.x_buf, maxp_co, ic_type, co_ws);
            orders[b] = ar_fit.p;
            phi1_values[b] = std::abs(1.0 - arma::accu(ar_fit.phi));

            for (int j = 0; j < ar_fit.p && j < maxp_co; ++j) {
                phi_matrix(b, j) = ar_fit.phi(j);
            }

            // PW t-statistic
            pw_tstats[b] = pw_tstat_ws(pw_ws.x_buf, pw_ws);

            // PW rho for COBA
            rho_values[b] = pw_estimate_rho_from_ws(pw_ws, n);
        }
    }
};


//' Combined CO + PW COBA Bootstrap Kernel
//'
//' Computes both CO and PW bootstrap with COBA data collection.
//'
//' @param n Integer, series length.
//' @param phi Numeric vector, AR coefficients from null model.
//' @param vara Double, innovation variance from null model.
//' @param seeds Vector of uint64 seeds.
//' @param maxp_co Integer, maximum AR order for CO.
//' @param maxp_pw Integer, maximum AR order for PW (typically 1).
//' @param criterion String, IC for AR selection.
//' @param grain_size Integer, minimum iterations per thread.
//' @return List with CO and PW results for COBA.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Rcpp::List wbg_bootstrap_copw_coba_kernel_cpp(
    int n,
    const arma::vec& phi,
    double vara,
    const std::vector<uint64_t>& seeds,
    int maxp_co = 5,
    int maxp_pw = 1,
    std::string criterion = "aic",
    std::size_t grain_size = 1
) {
    const int nb = seeds.size();

    // CO outputs
    Rcpp::NumericVector co_tstats(nb);
    Rcpp::NumericVector phi1_values(nb);
    Rcpp::NumericMatrix phi_matrix(nb, maxp_co);
    Rcpp::IntegerVector orders(nb);

    // PW outputs
    Rcpp::NumericVector pw_tstats(nb);
    Rcpp::NumericVector rho_values(nb);

    CriterionType ic_type = criterion_from_string(criterion);

    WBGBootstrapCOPWCOBAWorker worker(n, phi, vara, maxp_co, maxp_pw, ic_type, seeds,
                                       co_tstats, phi1_values, phi_matrix, orders,
                                       pw_tstats, rho_values);
    RcppParallel::parallelFor(0, nb, worker, grain_size);

    return Rcpp::List::create(
        // CO results
        Rcpp::Named("co_tstats") = co_tstats,
        Rcpp::Named("phi1_values") = phi1_values,
        Rcpp::Named("phi_matrix") = phi_matrix,
        Rcpp::Named("orders") = orders,
        // PW results
        Rcpp::Named("pw_tstats") = pw_tstats,
        Rcpp::Named("rho_values") = rho_values
    );
}


// =============================================================================
// PW-Only Bootstrap Kernels (generates from AR(1) with user-specified rho)
// =============================================================================

// TBB Worker for PW-only Bootstrap (generates from AR(1))
struct WBGBootstrapPWWorker : public RcppParallel::Worker {
    const int n;
    const double rho;           // AR(1) coefficient
    const double vara;
    const std::vector<uint64_t>& seeds;

    // Precomputed
    const double sd;
    const int burn;

    // Output
    RcppParallel::RVector<double> pw_results;

    WBGBootstrapPWWorker(int n_, double rho_, double vara_,
                         const std::vector<uint64_t>& seeds_,
                         Rcpp::NumericVector pw_results_)
        : n(n_), rho(rho_), vara(vara_), seeds(seeds_),
          sd(std::sqrt(vara_)),
          burn(calc_ar1_burnin(rho_, n_)),
          pw_results(pw_results_) {}

    // Simplified burn-in calculation for AR(1)
    static int calc_ar1_burnin(double rho, int n) {
        double abs_rho = std::abs(rho);
        if (abs_rho >= 0.999) {
            return std::min(2000, std::max({50, static_cast<int>(n * 0.5), 500}));
        } else if (abs_rho < 0.01) {
            return 50;
        } else {
            double log_factor = std::abs(std::log(0.001) / std::log(abs_rho));
            return std::min(2000, std::max(50, static_cast<int>(3.0 * log_factor)));
        }
    }

    void operator()(std::size_t begin, std::size_t end) {
        static thread_local PwBootstrapWorkspace pw_ws;
        static thread_local int ws_n = 0;

        if (ws_n != n) {
            pw_ws.resize(n);
            ws_n = n;
        }

        // Create RNG for this thread
        for (std::size_t b = begin; b < end; ++b) {
            dqrng::xoshiro256plus rng(seeds[b]);
            std::uniform_real_distribution<double> unif(0.0, 1.0);
            const double two_pi = 2.0 * M_PI;

            // Generate AR(1) series directly into workspace
            double prev = 0.0;

            // Burn-in phase
            for (int t = 0; t < burn; ++t) {
                // Box-Muller for single normal
                double u1 = unif(rng);
                double u2 = unif(rng);
                while (u1 <= 1e-15) u1 = unif(rng);
                double z = sd * std::sqrt(-2.0 * std::log(u1)) * std::cos(two_pi * u2);
                prev = rho * prev + z;
            }

            // Main generation phase
            for (int t = 0; t < n; ++t) {
                double u1 = unif(rng);
                double u2 = unif(rng);
                while (u1 <= 1e-15) u1 = unif(rng);
                double z = sd * std::sqrt(-2.0 * std::log(u1)) * std::cos(two_pi * u2);
                prev = rho * prev + z;
                pw_ws.x_buf[t] = prev;
            }

            // Compute PW t-statistic
            pw_results[b] = pw_tstat_ws(pw_ws.x_buf, pw_ws);
        }
    }
};


//' PW-Only Bootstrap Kernel (generates from AR(1))
//'
//' Computes PW bootstrap t-statistics by generating series from AR(1).
//'
//' @param n Integer, series length.
//' @param rho Double, AR(1) coefficient for null model.
//' @param vara Double, innovation variance.
//' @param seeds Vector of uint64 seeds.
//' @param grain_size Integer, minimum iterations per thread.
//' @return NumericVector of PW t-statistics.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector wbg_bootstrap_pw_kernel_cpp(
    int n,
    double rho,
    double vara,
    const std::vector<uint64_t>& seeds,
    std::size_t grain_size = 1
) {
    const int nb = seeds.size();
    Rcpp::NumericVector pw_tstats(nb);

    WBGBootstrapPWWorker worker(n, rho, vara, seeds, pw_tstats);
    RcppParallel::parallelFor(0, nb, worker, grain_size);

    return pw_tstats;
}


// =============================================================================
// PW-Only Bootstrap with COBA data collection
// =============================================================================

// TBB Worker for PW-only Bootstrap with rho collection (for COBA)
struct WBGBootstrapPWCOBAWorker : public RcppParallel::Worker {
    const int n;
    const double rho;
    const double vara;
    const std::vector<uint64_t>& seeds;

    const double sd;
    const int burn;

    // Outputs
    RcppParallel::RVector<double> pw_tstats;
    RcppParallel::RVector<double> rho_values;

    WBGBootstrapPWCOBAWorker(int n_, double rho_, double vara_,
                              const std::vector<uint64_t>& seeds_,
                              Rcpp::NumericVector pw_tstats_,
                              Rcpp::NumericVector rho_values_)
        : n(n_), rho(rho_), vara(vara_), seeds(seeds_),
          sd(std::sqrt(vara_)),
          burn(WBGBootstrapPWWorker::calc_ar1_burnin(rho_, n_)),
          pw_tstats(pw_tstats_), rho_values(rho_values_) {}

    void operator()(std::size_t begin, std::size_t end) {
        static thread_local PwBootstrapWorkspace pw_ws;
        static thread_local int ws_n = 0;

        if (ws_n != n) {
            pw_ws.resize(n);
            ws_n = n;
        }

        for (std::size_t b = begin; b < end; ++b) {
            dqrng::xoshiro256plus rng(seeds[b]);
            std::uniform_real_distribution<double> unif(0.0, 1.0);
            const double two_pi = 2.0 * M_PI;

            // Generate AR(1) series
            double prev = 0.0;

            // Burn-in phase
            for (int t = 0; t < burn; ++t) {
                double u1 = unif(rng);
                double u2 = unif(rng);
                while (u1 <= 1e-15) u1 = unif(rng);
                double z = sd * std::sqrt(-2.0 * std::log(u1)) * std::cos(two_pi * u2);
                prev = rho * prev + z;
            }

            // Main generation
            for (int t = 0; t < n; ++t) {
                double u1 = unif(rng);
                double u2 = unif(rng);
                while (u1 <= 1e-15) u1 = unif(rng);
                double z = sd * std::sqrt(-2.0 * std::log(u1)) * std::cos(two_pi * u2);
                prev = rho * prev + z;
                pw_ws.x_buf[t] = prev;
            }

            // Compute PW t-statistic
            pw_tstats[b] = pw_tstat_ws(pw_ws.x_buf, pw_ws);

            // Collect rho estimate from this bootstrap series
            rho_values[b] = pw_estimate_rho_from_ws(pw_ws, n);
        }
    }
};


//' PW-Only Bootstrap Kernel with COBA data collection
//'
//' Computes PW bootstrap t-statistics and collects rho estimates for COBA.
//'
//' @param n Integer, series length.
//' @param rho Double, AR(1) coefficient for null model.
//' @param vara Double, innovation variance.
//' @param seeds Vector of uint64 seeds.
//' @param grain_size Integer, minimum iterations per thread.
//' @return List with pw_tstats and rho_values.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Rcpp::List wbg_bootstrap_pw_coba_kernel_cpp(
    int n,
    double rho,
    double vara,
    const std::vector<uint64_t>& seeds,
    std::size_t grain_size = 1
) {
    const int nb = seeds.size();
    Rcpp::NumericVector pw_tstats(nb);
    Rcpp::NumericVector rho_values(nb);

    WBGBootstrapPWCOBAWorker worker(n, rho, vara, seeds, pw_tstats, rho_values);
    RcppParallel::parallelFor(0, nb, worker, grain_size);

    return Rcpp::List::create(
        Rcpp::Named("pw_tstats") = pw_tstats,
        Rcpp::Named("rho_values") = rho_values
    );
}
