// parallel_test.cpp - Parallel execution tests for bootstrap
// Uses RcppParallel (TBB) for parallelization
// [[Rcpp::depends(RcppArmadillo, RcppParallel, dqrng, BH)]]

#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include "types_pure.h"
#include <vector>

using namespace Rcpp;

// Forward declarations - thread-safe versions
arma::vec gen_ar_seeded_cpp(int n, const arma::vec& phi, double vara, uint64_t rng_seed);
double co_tstat_cpp(const arma::vec& x, int maxp, std::string criterion);
double co_tstat_pure(const arma::vec& x, int maxp, const std::string& criterion);


// =============================================================================
// Sequential Baseline (Pure C++)
// =============================================================================

//' Sequential Bootstrap Test (Baseline)
//'
//' @param n Integer, series length.
//' @param phi Numeric vector, AR coefficients.
//' @param vara Double, innovation variance.
//' @param seeds Vector of seeds for each iteration.
//' @param maxp Integer, maximum AR order.
//' @return Numeric vector of bootstrap t-statistics.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector test_parallel_seq(int n, const arma::vec& phi, double vara,
                                       const std::vector<uint64_t>& seeds,
                                       int maxp = 5) {
    const int nb = seeds.size();
    Rcpp::NumericVector results(nb);

    for (int b = 0; b < nb; ++b) {
        arma::vec x = gen_ar_seeded_cpp(n, phi, vara, seeds[b]);
        results[b] = co_tstat_pure(x, maxp, "aic");
    }

    return results;
}


// =============================================================================
// RcppParallel (TBB) with Pure C++ - Primary Implementation
// =============================================================================

struct BootstrapWorkerPure : public RcppParallel::Worker {
    const int n;
    const arma::vec phi;
    const double vara;
    const int maxp;
    const std::vector<uint64_t>& seeds;
    RcppParallel::RVector<double> results;

    BootstrapWorkerPure(int n_, const arma::vec& phi_, double vara_, int maxp_,
                        const std::vector<uint64_t>& seeds_,
                        Rcpp::NumericVector results_)
        : n(n_), phi(phi_), vara(vara_), maxp(maxp_), seeds(seeds_),
          results(results_) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t b = begin; b < end; ++b) {
            arma::vec x = gen_ar_seeded_cpp(n, phi, vara, seeds[b]);
            results[b] = co_tstat_pure(x, maxp, "aic");
        }
    }
};

//' TBB Parallel Bootstrap Test (Pure C++)
//'
//' Uses RcppParallel with pure C++ internals for thread safety.
//' This is the primary parallel implementation.
//'
//' @param n Integer, series length.
//' @param phi Numeric vector, AR coefficients.
//' @param vara Double, innovation variance.
//' @param seeds Vector of seeds for each iteration.
//' @param maxp Integer, maximum AR order.
//' @return Numeric vector of bootstrap t-statistics.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector test_parallel_tbb_pure(int n, const arma::vec& phi, double vara,
                                            const std::vector<uint64_t>& seeds,
                                            int maxp = 5) {
    const int nb = seeds.size();
    Rcpp::NumericVector results(nb);

    BootstrapWorkerPure worker(n, phi, vara, maxp, seeds, results);
    RcppParallel::parallelFor(0, nb, worker);

    return results;
}


// =============================================================================
// RcppParallel (TBB) with Rcpp Types (for comparison)
// Note: Uses co_tstat_cpp which has Rcpp::List internally
// =============================================================================

struct BootstrapWorkerRcpp : public RcppParallel::Worker {
    const int n;
    const arma::vec phi;
    const double vara;
    const int maxp;
    const std::vector<uint64_t>& seeds;
    RcppParallel::RVector<double> results;

    BootstrapWorkerRcpp(int n_, const arma::vec& phi_, double vara_, int maxp_,
                        const std::vector<uint64_t>& seeds_,
                        Rcpp::NumericVector results_)
        : n(n_), phi(phi_), vara(vara_), maxp(maxp_), seeds(seeds_),
          results(results_) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t b = begin; b < end; ++b) {
            arma::vec x = gen_ar_seeded_cpp(n, phi, vara, seeds[b]);
            // This uses Rcpp::List internally - works with TBB
            results[b] = co_tstat_cpp(x, maxp, "aic");
        }
    }
};

//' TBB Parallel Bootstrap Test (Rcpp Types)
//'
//' Uses RcppParallel with Rcpp types internally.
//' For comparison with pure C++ version.
//'
//' @param n Integer, series length.
//' @param phi Numeric vector, AR coefficients.
//' @param vara Double, innovation variance.
//' @param seeds Vector of seeds for each iteration.
//' @param maxp Integer, maximum AR order.
//' @return Numeric vector of bootstrap t-statistics.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector test_parallel_tbb_rcpp(int n, const arma::vec& phi, double vara,
                                            const std::vector<uint64_t>& seeds,
                                            int maxp = 5) {
    const int nb = seeds.size();
    Rcpp::NumericVector results(nb);

    BootstrapWorkerRcpp worker(n, phi, vara, maxp, seeds, results);
    RcppParallel::parallelFor(0, nb, worker);

    return results;
}


// =============================================================================
// Default alias - uses pure C++ version
// =============================================================================

//' TBB Parallel Bootstrap Test (Default)
//'
//' Default parallel bootstrap using pure C++ version.
//'
//' @param n Integer, series length.
//' @param phi Numeric vector, AR coefficients.
//' @param vara Double, innovation variance.
//' @param seeds Vector of seeds for each iteration.
//' @param maxp Integer, maximum AR order.
//' @return Numeric vector of bootstrap t-statistics.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector test_parallel_tbb(int n, const arma::vec& phi, double vara,
                                       const std::vector<uint64_t>& seeds,
                                       int maxp = 5) {
    return test_parallel_tbb_pure(n, phi, vara, seeds, maxp);
}
