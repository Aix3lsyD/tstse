// gen_ar_fast_cpp.cpp - Fast AR generation with thread-safe RNG
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
// Thread-Safe AR Generation using dqrng (xoshiro256+)
// =============================================================================

//' AR Generation with dqrng scalar loop (thread-safe)
//'
//' Uses dqrng's xoshiro256+ with scalar generation.
//' Thread-safe: creates local RNG instance per call.
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

  // Generate innovations one at a time
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


//' AR Generation with dqrng Box-Muller (thread-safe, optimized)
//'
//' Uses dqrng's fast uniform generation with Box-Muller transform.
//' Thread-safe: creates local RNG instance per call.
//' Generates normals in pairs for efficiency.
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
// Primary API Functions
// =============================================================================

//' Fast AR Generation with External Seed (for parallel use)
//'
//' Uses Box-Muller with dqrng - reproducible and thread-safe.
//' This is the primary function for parallel bootstrap.
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
  return gen_ar_dqrng_boxmuller(n, phi, vara, rng_seed);
}


//' Fast AR Generation (C++ implementation)
//'
//' Generates AR process using thread-safe dqrng.
//' If seed is provided, uses it for reproducibility.
//' If no seed, generates a random seed from R's RNG.
//'
//' @param n Integer, series length to generate.
//' @param phi Numeric vector, AR coefficients.
//' @param vara Double, innovation variance.
//' @param seed Integer, random seed for reproducibility (optional).
//' @return Numeric vector of length n.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::vec gen_ar_cpp(int n, const arma::vec& phi, double vara = 1.0,
                      Rcpp::Nullable<int> seed = R_NilValue) {
  uint64_t rng_seed;

  if (seed.isNotNull()) {
    // Use provided seed
    rng_seed = static_cast<uint64_t>(Rcpp::as<int>(seed));
  } else {
    // Generate random seed from R's RNG (not thread-safe, but this is single-threaded)
    rng_seed = static_cast<uint64_t>(R::runif(0, 1) * std::numeric_limits<uint32_t>::max());
  }

  return gen_ar_dqrng_boxmuller(n, phi, vara, rng_seed);
}
