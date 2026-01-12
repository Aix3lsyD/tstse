// =============================================================================
// FILE: kernel_gen_ar.cpp
// CATEGORY: HOT PATH - Thread-safe AR generation with dqrng
// THREAD-SAFE: YES (each call creates local xoshiro256+ RNG)
//
// Generates AR processes for parallel bootstrap. Uses dqrng's xoshiro256+
// PRNG which is faster than R's Mersenne Twister and fully thread-safe.
//
// Key functions:
//   - gen_ar_into(): Zero-copy generation into workspace buffer
//   - gen_ar_seeded_cpp(): Seeded AR generation for reproducibility
//
// Called by: kernel_wbg_boot.cpp (WBGBootstrapWorker)
// =============================================================================
// [[Rcpp::depends(RcppArmadillo, dqrng, BH)]]

#include <RcppArmadillo.h>
#include <dqrng.h>
#include <dqrng_distribution.h>
#include <xoshiro.h>
#include <cmath>
#include <limits>
#include <algorithm>  // for std::max with initializer_list
#include <random>     // for std::uniform_real_distribution

// M_PI fallback for portability (not defined on MSVC without _USE_MATH_DEFINES)
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace Rcpp;

// Calculate Adaptive Burn-in for AR Process
// Internal helper: Calculate adaptive burn-in (not exported to R)
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

// AR Generation with dqrng Box-Muller (thread-safe, optimized)
// Uses dqrng's fast uniform generation with Box-Muller transform.
// Thread-safe: creates local RNG instance per call.
// Generates normals in pairs for efficiency.
// Internal: AR generation with Box-Muller (not exported to R)
// This is the primary internal implementation used by gen_ar_seeded_cpp
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
// Buffer-based AR Generation (for worker reuse)
// =============================================================================

// Helper: fill vector with Box-Muller normals
// Template to work with both arma::vec& and subviews (e.g., burn_buf.head(n))
template<typename VecType>
static void fill_normals_boxmuller(VecType&& out, dqrng::xoshiro256plus& rng, double sd) {
    const int n = out.n_elem;
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    const double two_pi = 2.0 * M_PI;

    int i = 0;
    while (i < n) {
        double u1 = unif(rng);
        double u2 = unif(rng);
        while (u1 <= 1e-15) u1 = unif(rng);

        double r = sd * std::sqrt(-2.0 * std::log(u1));
        double theta = two_pi * u2;

        out[i] = r * std::cos(theta);
        if (i + 1 < n) {
            out[i + 1] = r * std::sin(theta);
        }
        i += 2;
    }
}

// =============================================================================
// Internal: AR generation with precomputed sd and burn
// Used by parallel workers where sd and burn are constant across iterations
// =============================================================================
void gen_ar_into_precomputed(arma::vec& out, const arma::vec& phi,
                              double sd, int burn, uint64_t rng_seed) {
    const int n = out.n_elem;
    const int p = phi.n_elem;

    // Guard against buffer overflow - MAX_P=20 for stack arrays below
    constexpr int MAX_P = 20;
    if (p > MAX_P) {
        out.zeros();
        return;
    }

    dqrng::xoshiro256plus rng(rng_seed);

    if (p == 0) {
        // AR(0): fill directly with normals
        fill_normals_boxmuller(out, rng, sd);
        return;
    }

    // Pre-generate all burn-in innovations efficiently
    // Uses both Box-Muller normals (cos AND sin), halving log/sqrt calls
    arma::vec burn_innovations(burn);
    fill_normals_boxmuller(burn_innovations, rng, sd);

    // Stack buffer for burn-in (ring buffer, only keep last p values)
    double prev[MAX_P] = {0};

    // Run burn-in phase using pre-generated innovations
    for (int t = 0; t < burn; ++t) {
        double x = burn_innovations[t];
        for (int j = 0; j < p && j <= t; ++j) {
            int idx = ((t - 1 - j) % MAX_P + MAX_P) % MAX_P;
            x += phi[j] * prev[idx];
        }
        prev[t % MAX_P] = x;
    }

    // Extract last p values from ring buffer for initial conditions
    double burn_tail[MAX_P];
    for (int i = 0; i < p; ++i) {
        burn_tail[i] = prev[((burn - p + i) % MAX_P + MAX_P) % MAX_P];
    }

    // Fill output with innovations
    fill_normals_boxmuller(out, rng, sd);

    // Apply AR recursion using burn_tail for initial dependencies
    for (int t = 0; t < n; ++t) {
        for (int j = 0; j < p; ++j) {
            int lag = t - 1 - j;
            if (lag >= 0) {
                out[t] += phi[j] * out[lag];
            } else {
                // Use burn-in tail for negative lags
                out[t] += phi[j] * burn_tail[p + lag];
            }
        }
    }
}


// =============================================================================
// Workspace-aware AR generation (ZERO allocations in hot path)
// Uses pre-allocated burn_buf from workspace instead of allocating per call
// =============================================================================
void gen_ar_into_ws(arma::vec& out, const arma::vec& phi,
                    double sd, int burn, uint64_t rng_seed,
                    arma::vec& burn_buf) {
    const int n = out.n_elem;
    const int p = phi.n_elem;

    // Guard against buffer overflow - MAX_P=20 for stack arrays below
    constexpr int MAX_P = 20;
    if (p > MAX_P) {
        out.zeros();
        return;
    }

    dqrng::xoshiro256plus rng(rng_seed);

    if (p == 0) {
        // AR(0): fill directly with normals
        fill_normals_boxmuller(out, rng, sd);
        return;
    }

    // Use workspace buffer for burn-in innovations (ZERO allocation!)
    // Fill burn_buf directly using subview (head() returns a reference when used in-place)
    fill_normals_boxmuller(burn_buf.head(burn), rng, sd);

    // Stack buffer for burn-in (ring buffer, only keep last p values)
    double prev[MAX_P] = {0};

    // Run burn-in phase using workspace buffer
    for (int t = 0; t < burn; ++t) {
        double x = burn_buf[t];
        for (int j = 0; j < p && j <= t; ++j) {
            int idx = ((t - 1 - j) % MAX_P + MAX_P) % MAX_P;
            x += phi[j] * prev[idx];
        }
        prev[t % MAX_P] = x;
    }

    // Extract last p values from ring buffer for initial conditions
    double burn_tail[MAX_P];
    for (int i = 0; i < p; ++i) {
        burn_tail[i] = prev[((burn - p + i) % MAX_P + MAX_P) % MAX_P];
    }

    // Fill output with innovations
    fill_normals_boxmuller(out, rng, sd);

    // Apply AR recursion using burn_tail for initial dependencies
    for (int t = 0; t < n; ++t) {
        for (int j = 0; j < p; ++j) {
            int lag = t - 1 - j;
            if (lag >= 0) {
                out[t] += phi[j] * out[lag];
            } else {
                // Use burn-in tail for negative lags
                out[t] += phi[j] * burn_tail[p + lag];
            }
        }
    }
}


// Generate AR series directly into provided buffer (avoids copy on return)
// Public API: computes sd and burn internally (use gen_ar_into_precomputed for hot paths)
void gen_ar_into(arma::vec& out, const arma::vec& phi, double vara, uint64_t rng_seed) {
    const int n = out.n_elem;
    const int p = phi.n_elem;
    const double sd = std::sqrt(vara);
    const int burn = (p == 0) ? 0 : calc_ar_burnin_cpp(phi, n);

    gen_ar_into_precomputed(out, phi, sd, burn, rng_seed);
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
