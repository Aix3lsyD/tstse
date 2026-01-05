// wbg_boot_kernel.cpp - TBB parallel bootstrap kernel for WBG test
// Uses RcppParallel for parallelization
// [[Rcpp::depends(RcppArmadillo, RcppParallel, dqrng, BH)]]

#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include "types_pure.h"
#include <vector>

using namespace Rcpp;

// Forward declarations - thread-safe functions
arma::vec gen_ar_dqrng_boxmuller(int n, const arma::vec& phi, double vara, uint64_t rng_seed);
void gen_ar_into(arma::vec& out, const arma::vec& phi, double vara, uint64_t rng_seed);
double co_tstat_pure(const arma::vec& x, int maxp, const std::string& criterion);
BurgResult burg_aic_select_pure(const arma::vec& x, int maxp, const std::string& criterion);

// Forward declarations - workspace-aware functions (ZERO allocations)
double co_tstat_ws(const arma::vec& x, int maxp, CriterionType ic_type,
                   CoBootstrapWorkspace& ws);
BurgResult burg_aic_select_ws(const arma::vec& x, int maxp, CriterionType ic_type,
                               CoBootstrapWorkspace& ws);


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

    // Output (thread-safe write via RVector)
    RcppParallel::RVector<double> results;

    // Constructor
    WBGBootstrapWorker(int n_, const arma::vec& phi_, double vara_, int maxp_,
                       CriterionType ic_type_,
                       const std::vector<uint64_t>& seeds_,
                       Rcpp::NumericVector results_)
        : n(n_), phi(phi_), vara(vara_), maxp(maxp_), ic_type(ic_type_),
          seeds(seeds_), results(results_) {}

    // Parallel worker function
    void operator()(std::size_t begin, std::size_t end) {
        // Allocate workspace ONCE for this thread chunk
        // All iterations in this chunk reuse the same buffers
        CoBootstrapWorkspace ws;
        ws.resize(n, maxp);

        for (std::size_t b = begin; b < end; ++b) {
            // Generate AR series directly into workspace buffer (no copy)
            gen_ar_into(ws.x_buf, phi, vara, seeds[b]);

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

    // Pre-allocated output (thread-safe writes to separate indices)
    RcppParallel::RVector<double> tstats;
    RcppParallel::RVector<double> phi1_values;
    RcppParallel::RMatrix<double> phi_matrix;  // nb × maxp
    RcppParallel::RVector<int> orders;

    // Constructor
    WBGBootstrapCOBAWorker(int n_, const arma::vec& phi_, double vara_, int maxp_,
                           CriterionType ic_type_,
                           const std::vector<uint64_t>& seeds_,
                           Rcpp::NumericVector tstats_,
                           Rcpp::NumericVector phi1_values_,
                           Rcpp::NumericMatrix phi_matrix_,
                           Rcpp::IntegerVector orders_)
        : n(n_), phi(phi_), vara(vara_), maxp(maxp_), ic_type(ic_type_),
          seeds(seeds_), tstats(tstats_), phi1_values(phi1_values_),
          phi_matrix(phi_matrix_), orders(orders_) {}

    // Parallel worker function
    void operator()(std::size_t begin, std::size_t end) {
        // Allocate workspace ONCE for this thread chunk
        CoBootstrapWorkspace ws;
        ws.resize(n, maxp);

        for (std::size_t b = begin; b < end; ++b) {
            // Generate AR series directly into workspace buffer (no copy)
            gen_ar_into(ws.x_buf, phi, vara, seeds[b]);

            // Compute CO t-statistic using workspace (ZERO allocations)
            tstats[b] = co_tstat_ws(ws.x_buf, maxp, ic_type, ws);

            // Fit AR model to bootstrap sample using workspace
            // Note: This reuses ws.xc, ws.ef, ws.eb, ws.a_curr, ws.a_prev
            BurgResult ar_fit = burg_aic_select_ws(ws.x_buf, maxp, ic_type, ws);

            // Store results
            orders[b] = ar_fit.p;
            phi1_values[b] = 1.0 - arma::accu(ar_fit.phi);  // φ(1) from paper

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
