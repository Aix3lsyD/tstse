// =============================================================================
// FILE: util_gen_arch.cpp
// CATEGORY: UTILITY
// THREAD-SAFE: NO (uses R's RNG)
//
// ARCH process generation for volatility modeling.
// Not used in WBG bootstrap (which is AR-based, not ARCH/GARCH).
//
// Exports:
//   - gen_arch_cpp(): Generate ARCH(q) process
// =============================================================================

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector gen_arch_cpp(int n,
                           double alpha0,
                           NumericVector alpha,
                           int spin = 1000,
                           Nullable<int> seed = R_NilValue) {

  if (seed.isNotNull()) {
    Environment base_env("package:base");
    Function set_seed = base_env["set.seed"];
    set_seed(as<int>(seed));
  }

  const int q = alpha.size();
  const int ntot = n + spin + q;

  // Generate standard normal innovations
  NumericVector eps = rnorm(ntot, 0.0, 1.0);

  // Initialize variance and squared innovations
  NumericVector sig2(ntot, 1.0);
  NumericVector asq(ntot, 0.0);

  // ARCH recursion: sig2[t] = alpha0 + sum(alpha[i] * asq[t-i]^2)
  for (int t = q; t < ntot; ++t) {
    double s2 = alpha0;

    // ARCH terms
    for (int i = 0; i < q; ++i) {
      s2 += alpha[i] * asq[t - i - 1] * asq[t - i - 1];
    }

    sig2[t] = s2;
    asq[t] = eps[t] * std::sqrt(s2);
  }

  // Extract after spin-up
  const int start_idx = spin + q;
  NumericVector result(n);
  for (int i = 0; i < n; ++i) {
    result[i] = asq[start_idx + i];
  }

  return result;
}
