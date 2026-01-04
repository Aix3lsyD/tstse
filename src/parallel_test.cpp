// parallel_test.cpp - Parallel execution tests for bootstrap
// Tests both Rcpp and Pure C++ paths with TBB and OpenMP
// [[Rcpp::depends(RcppArmadillo, RcppParallel, dqrng, BH)]]

#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include "types_pure.h"
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// Forward declarations - Rcpp versions (from other files)
arma::vec gen_ar_seeded_cpp(int n, const arma::vec& phi, double vara, uint64_t rng_seed);
double co_tstat_cpp(const arma::vec& x, int maxp, std::string criterion);

// Forward declarations - Pure C++ versions
double co_tstat_pure(const arma::vec& x, int maxp, const std::string& criterion);


// =============================================================================
// OpenMP Status Functions
// =============================================================================

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
int get_omp_threads() {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
bool has_openmp() {
#ifdef _OPENMP
    return true;
#else
    return false;
#endif
}


// =============================================================================
// Sequential Baseline (Pure C++)
// =============================================================================

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
// OpenMP with Pure C++ (Thread-Safe)
// =============================================================================

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector test_parallel_omp_pure(int n, const arma::vec& phi, double vara,
                                            const std::vector<uint64_t>& seeds,
                                            int maxp = 5) {
    const int nb = seeds.size();
    std::vector<double> results(nb);

#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (int b = 0; b < nb; ++b) {
        // All functions here are pure C++ - no Rcpp types
        arma::vec x = gen_ar_seeded_cpp(n, phi, vara, seeds[b]);
        results[b] = co_tstat_pure(x, maxp, "aic");
    }

    return Rcpp::wrap(results);
}


// =============================================================================
// RcppParallel (TBB) with Pure C++
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
// Note: This uses co_tstat_cpp which has Rcpp::List internally
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
            // This uses Rcpp::List internally - works with TBB but not OpenMP
            results[b] = co_tstat_cpp(x, maxp, "aic");
        }
    }
};

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
// Legacy aliases for backward compatibility
// =============================================================================

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector test_parallel_omp(int n, const arma::vec& phi, double vara,
                                       const std::vector<uint64_t>& seeds,
                                       int maxp = 5) {
    // Use pure C++ version (thread-safe)
    return test_parallel_omp_pure(n, phi, vara, seeds, maxp);
}

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector test_parallel_tbb(int n, const arma::vec& phi, double vara,
                                       const std::vector<uint64_t>& seeds,
                                       int maxp = 5) {
    // Use pure C++ version by default
    return test_parallel_tbb_pure(n, phi, vara, seeds, maxp);
}
