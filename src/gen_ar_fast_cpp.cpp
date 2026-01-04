// gen_ar_fast_cpp.cpp - Fast AR generation with multiple RNG backends
// Part of wbg_boot_fast optimization
// [[Rcpp::depends(RcppArmadillo, dqrng, BH)]]

#include <RcppArmadillo.h>
#include <dqrng.h>
#include <dqrng_distribution.h>
#include <xoshiro.h>
#include <cmath>

using namespace Rcpp;

//' Calculate Adaptive Burn-in for AR Process
//'
//' @param phi Numeric vector, AR coefficients.
//' @param n Integer, target series length.
//' @return Integer, recommended burn-in length.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
int calc_ar_burnin_cpp(const arma::vec& phi, int n) {
  // Compute persistence (sum of absolute coefficients)
  double persistence = arma::sum(arma::abs(phi));

  int burn;
  if (persistence >= 0.999) {
    // Near unit root: need long burn-in
    burn = std::max({50, static_cast<int>(n * 0.5), 500});
  } else if (persistence < 0.01) {
    // Nearly white noise: minimal burn-in
    burn = 50;
  } else {
    // Normal case: decay rate based
    double log_factor = std::abs(std::log(0.001) / std::log(persistence));
    burn = std::max(50, static_cast<int>(3.0 * log_factor));
  }

  // Cap at 2000
  return std::min(burn, 2000);
}


// =============================================================================
// VARIANT 1: Rcpp::rnorm - Fast but NOT thread-safe
// Use for single-threaded operations and benchmarking
// =============================================================================

//' Fast AR Generation using Rcpp::rnorm (NOT thread-safe)
//'
//' Uses R's vectorized rnorm for maximum single-threaded performance.
//' NOT safe for parallel execution - uses R's global RNG state.
//'
//' @param n Integer, series length to generate.
//' @param phi Numeric vector, AR coefficients.
//' @param vara Double, innovation variance.
//' @return Numeric vector of length n.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::vec gen_ar_rcpp(int n, const arma::vec& phi, double vara = 1.0) {
  const int p = phi.n_elem;

  // Calculate adaptive burn-in
  int burn = (p == 0) ? 0 : calc_ar_burnin_cpp(phi, n);
  int n_total = n + burn;

  // Use R's vectorized rnorm - single call, very fast
  double sd = std::sqrt(vara);
  Rcpp::NumericVector innovations = Rcpp::rnorm(n_total, 0.0, sd);

  // Convert to arma::vec for AR recursion
  arma::vec a(innovations.begin(), n_total, false, true);

  // Handle AR(0) case - just return the innovations
  if (p == 0) {
    return a;
  }

  // AR recursion (in-place)
  for (int t = p; t < n_total; ++t) {
    for (int j = 0; j < p; ++j) {
      a[t] += phi[j] * a[t - j - 1];
    }
  }

  // Return after burn-in
  return a.tail(n);
}


// =============================================================================
// VARIANT 2: dqrng scalar - Thread-safe but slower
// Original implementation with per-element generation
// =============================================================================

//' AR Generation with dqrng scalar loop (thread-safe, slower)
//'
//' Uses dqrng's xoshiro256+ with scalar generation.
//' Thread-safe but has per-element call overhead.
//'
//' @param n Integer, series length to generate.
//' @param phi Numeric vector, AR coefficients.
//' @param vara Double, innovation variance.
//' @param rng_seed Integer, seed for this specific generation.
//' @return Numeric vector of length n.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::vec gen_ar_dqrng_scalar(int n, const arma::vec& phi, double vara,
                               uint64_t rng_seed) {
  const int p = phi.n_elem;

  // Create thread-local RNG with specific seed
  dqrng::xoshiro256plus rng(rng_seed);
  boost::random::normal_distribution<double> norm(0.0, std::sqrt(vara));

  // Calculate adaptive burn-in
  int burn = (p == 0) ? 0 : calc_ar_burnin_cpp(phi, n);
  int n_total = n + burn;

  // Generate innovations one at a time (the slow part)
  arma::vec a(n_total);
  for (int i = 0; i < n_total; ++i) {
    a[i] = norm(rng);
  }

  // Handle AR(0) case
  if (p == 0) {
    return a.head(n);
  }

  // AR recursion
  for (int t = p; t < n_total; ++t) {
    for (int j = 0; j < p; ++j) {
      a[t] += phi[j] * a[t - j - 1];
    }
  }

  return a.tail(n);
}


// =============================================================================
// VARIANT 3: Armadillo randn - Thread-safe with arma_rng::set_seed_random()
// Uses Armadillo's internal RNG which can be thread-local
// =============================================================================

//' AR Generation with Armadillo randn (thread-safe)
//'
//' Uses Armadillo's randn with explicit seeding.
//' Thread-safe when each thread uses different seed.
//'
//' @param n Integer, series length to generate.
//' @param phi Numeric vector, AR coefficients.
//' @param vara Double, innovation variance.
//' @param rng_seed Integer, seed for Armadillo's RNG.
//' @return Numeric vector of length n.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::vec gen_ar_arma(int n, const arma::vec& phi, double vara,
                       unsigned int rng_seed) {
  const int p = phi.n_elem;

  // Seed Armadillo's RNG (thread-local in recent versions)
  arma::arma_rng::set_seed(rng_seed);

  // Calculate adaptive burn-in
  int burn = (p == 0) ? 0 : calc_ar_burnin_cpp(phi, n);
  int n_total = n + burn;

  // Generate innovations - vectorized!
  double sd = std::sqrt(vara);
  arma::vec a = sd * arma::randn(n_total);

  // Handle AR(0) case
  if (p == 0) {
    return a.head(n);
  }

  // AR recursion
  for (int t = p; t < n_total; ++t) {
    for (int j = 0; j < p; ++j) {
      a[t] += phi[j] * a[t - j - 1];
    }
  }

  return a.tail(n);
}


// =============================================================================
// VARIANT 4: dqrng with Ziggurat batching
// Generate uniforms in batch, transform to normals
// =============================================================================

//' AR Generation with dqrng batched Ziggurat (thread-safe, optimized)
//'
//' Uses dqrng's fast uniform generation in batches, then applies
//' Box-Muller transform for normals. Reduces function call overhead.
//'
//' @param n Integer, series length to generate.
//' @param phi Numeric vector, AR coefficients.
//' @param vara Double, innovation variance.
//' @param rng_seed Integer, seed for this specific generation.
//' @return Numeric vector of length n.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::vec gen_ar_dqrng_boxmuller(int n, const arma::vec& phi, double vara,
                                  uint64_t rng_seed) {
  const int p = phi.n_elem;

  // Create thread-local RNG
  dqrng::xoshiro256plus rng(rng_seed);

  // Calculate adaptive burn-in
  int burn = (p == 0) ? 0 : calc_ar_burnin_cpp(phi, n);
  int n_total = n + burn;

  // Generate normals using Box-Muller with batched uniforms
  double sd = std::sqrt(vara);
  arma::vec a(n_total);

  // Box-Muller generates pairs of normals
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  const double two_pi = 2.0 * M_PI;

  int i = 0;
  while (i < n_total) {
    double u1 = unif(rng);
    double u2 = unif(rng);

    // Avoid log(0)
    while (u1 <= 1e-15) u1 = unif(rng);

    double r = sd * std::sqrt(-2.0 * std::log(u1));
    double theta = two_pi * u2;

    a[i] = r * std::cos(theta);
    if (i + 1 < n_total) {
      a[i + 1] = r * std::sin(theta);
    }
    i += 2;
  }

  // Handle AR(0) case
  if (p == 0) {
    return a.head(n);
  }

  // AR recursion
  for (int t = p; t < n_total; ++t) {
    for (int j = 0; j < p; ++j) {
      a[t] += phi[j] * a[t - j - 1];
    }
  }

  return a.tail(n);
}


// =============================================================================
// Keep original exports for backwards compatibility
// =============================================================================

//' Fast AR Generation (C++ implementation) - Original API
//'
//' @param n Integer, series length to generate.
//' @param phi Numeric vector, AR coefficients.
//' @param vara Double, innovation variance.
//' @param seed Integer, random seed for reproducibility.
//' @return Numeric vector of length n.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::vec gen_ar_cpp(int n, const arma::vec& phi, double vara = 1.0,
                      Rcpp::Nullable<int> seed = R_NilValue) {
  if (seed.isNotNull()) {
    // If seed provided, use Armadillo version for reproducibility
    return gen_ar_arma(n, phi, vara, static_cast<unsigned int>(Rcpp::as<int>(seed)));
  } else {
    // No seed: use fast Rcpp::rnorm version
    return gen_ar_rcpp(n, phi, vara);
  }
}

//' Fast AR Generation with External Seed (for parallel) - Original API
//'
//' Uses Box-Muller with dqrng - reproducible and thread-safe.
//' Note: Armadillo is faster but not reproducible across calls.
//'
//' @param n Integer, series length to generate.
//' @param phi Numeric vector, AR coefficients.
//' @param vara Double, innovation variance.
//' @param rng_seed Integer, seed for this specific generation.
//' @return Numeric vector of length n.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::vec gen_ar_seeded_cpp(int n, const arma::vec& phi, double vara,
                             uint64_t rng_seed) {
  // Use Box-Muller with dqrng - reproducible and thread-safe
  // (Armadillo is faster but not reproducible across calls)
  return gen_ar_dqrng_boxmuller(n, phi, vara, rng_seed);
}
