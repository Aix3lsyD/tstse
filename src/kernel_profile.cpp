// =============================================================================
// FILE: kernel_profile.cpp
// CATEGORY: Diagnostic - Sub-component profiling for bootstrap kernel
//
// Runs bootstrap iterations SEQUENTIALLY with std::chrono timing around each
// sub-operation. Used by the Shiny Performance Profile tab to show a
// parent-child breakdown of the bootstrap kernel.
//
// This is a diagnostic-only function. It does NOT affect production performance.
// The production kernel (kernel_wbg_boot.cpp) is completely untouched.
//
// Key export:
//   - wbg_profile_kernel_components_cpp(): Returns timing ratios per sub-step
// =============================================================================
// [[Rcpp::depends(RcppArmadillo, dqrng, BH)]]

#include <RcppArmadillo.h>
#include <chrono>
#include "kernel_types.h"

using namespace Rcpp;
using hrclock = std::chrono::high_resolution_clock;

// Forward declarations for kernel sub-operations
void gen_ar_into_ws(arma::vec& out, const arma::vec& phi,
                    double sd, int burn, uint64_t rng_seed, arma::vec& burn_buf);
int calc_ar_burnin_cpp(const arma::vec& phi, int n);
void ols_detrend_ws(const arma::vec& x, CoBootstrapWorkspace& ws);
BurgResult burg_aic_select_ws(const arma::vec& x, int maxp, CriterionType ic_type,
                               CoBootstrapWorkspace& ws, int min_p);
double co_tstat_fused(const arma::vec& x, const arma::vec& phi);


//' Profile Bootstrap Kernel Sub-Components
//'
//' Runs bootstrap iterations sequentially with timing around each sub-operation
//' to produce a sub-component breakdown of the bootstrap kernel. Used by the
//' Shiny Performance Profile tab for parent-child visualization.
//'
//' @param n Integer, series length.
//' @param phi Numeric vector, AR coefficients from null model.
//' @param vara Double, innovation variance from null model.
//' @param seeds Vector of uint64 seeds (determines number of iterations).
//' @param maxp Integer, maximum AR order for CO test.
//' @param criterion String, IC for AR selection: "aic", "aicc", "bic".
//' @param min_p Integer, minimum AR order for CO residual model (default 0).
//' @param coba Logical, whether to also time COBA's second Burg fit.
//' @return List with mean microseconds per sub-component:
//'   gen_ar_us, ols_us, burg_us, co_fused_us, burg_coba_us (0 if !coba),
//'   total_iter_us, nreps.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Rcpp::List wbg_profile_kernel_components_cpp(
    int n,
    const arma::vec& phi,
    double vara,
    const std::vector<uint64_t>& seeds,
    int maxp,
    std::string criterion,
    int min_p,
    bool coba) {

    const int nreps = seeds.size();
    if (nreps < 1) {
        Rcpp::stop("seeds must have at least 1 element");
    }

    // Convert criterion string to enum once
    CriterionType ic_type = criterion_from_string(criterion);

    // Precompute AR generation params (same as TBB worker constructor)
    const double sd = std::sqrt(vara);
    const int burn = (phi.n_elem == 0) ? 0 : calc_ar_burnin_cpp(phi, n);

    // Set up workspace (same as TBB thread-local workspace)
    CoBootstrapWorkspace ws;
    ws.resize(n, maxp);

    // Accumulators (nanoseconds)
    int64_t acc_gen_ar = 0;
    int64_t acc_ols = 0;
    int64_t acc_burg = 0;
    int64_t acc_fused = 0;
    int64_t acc_burg_coba = 0;
    int64_t acc_total = 0;

    for (int b = 0; b < nreps; ++b) {
        auto t_iter_start = hrclock::now();

        // 1) AR series generation
        auto t0 = hrclock::now();
        gen_ar_into_ws(ws.x_buf, phi, sd, burn, seeds[b], ws.burn_buf);
        auto t1 = hrclock::now();
        acc_gen_ar += std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();

        // 2) OLS detrend → writes to ws.resid
        t0 = hrclock::now();
        ols_detrend_ws(ws.x_buf, ws);
        t1 = hrclock::now();
        acc_ols += std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();

        // 3) Burg AR fit on detrended residuals
        t0 = hrclock::now();
        BurgResult ar_fit = burg_aic_select_ws(ws.resid, maxp, ic_type, ws, min_p);
        t1 = hrclock::now();
        acc_burg += std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();

        // 4) CO fused transform + regression
        t0 = hrclock::now();
        co_tstat_fused(ws.x_buf, ar_fit.phi);
        t1 = hrclock::now();
        acc_fused += std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();

        // 5) COBA: second Burg fit on bootstrap sample (min_p=0 per WBG paper)
        if (coba) {
            t0 = hrclock::now();
            burg_aic_select_ws(ws.x_buf, maxp, ic_type, ws, 0);
            t1 = hrclock::now();
            acc_burg_coba += std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();
        }

        auto t_iter_end = hrclock::now();
        acc_total += std::chrono::duration_cast<std::chrono::nanoseconds>(t_iter_end - t_iter_start).count();
    }

    // Convert to mean microseconds
    double us_gen_ar     = static_cast<double>(acc_gen_ar) / (nreps * 1000.0);
    double us_ols        = static_cast<double>(acc_ols) / (nreps * 1000.0);
    double us_burg       = static_cast<double>(acc_burg) / (nreps * 1000.0);
    double us_fused      = static_cast<double>(acc_fused) / (nreps * 1000.0);
    double us_burg_coba  = static_cast<double>(acc_burg_coba) / (nreps * 1000.0);
    double us_total      = static_cast<double>(acc_total) / (nreps * 1000.0);

    return Rcpp::List::create(
        Rcpp::Named("gen_ar_us")     = us_gen_ar,
        Rcpp::Named("ols_us")        = us_ols,
        Rcpp::Named("burg_us")       = us_burg,
        Rcpp::Named("co_fused_us")   = us_fused,
        Rcpp::Named("burg_coba_us")  = us_burg_coba,
        Rcpp::Named("total_iter_us") = us_total,
        Rcpp::Named("nreps")         = nreps
    );
}
