// =============================================================================
// FILE: util_gen_garch.cpp
// CATEGORY: UTILITY
// THREAD-SAFE: NO (uses R's RNG)
//
// GARCH process generation for volatility modeling.
// Not used in WBG bootstrap (which is AR-based, not ARCH/GARCH).
//
// Exports:
//   - gen_garch_cpp(): Generate GARCH(p,q) process
// =============================================================================

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector gen_garch_cpp(int n,
                            double alpha0,
                            NumericVector alpha,
                            NumericVector beta,
                            int spin = 1000,
                            Nullable<int> seed = R_NilValue) {

  if (seed.isNotNull()) {
    Environment base_env("package:base");
    Function set_seed = base_env["set.seed"];
    set_seed(as<int>(seed));
  }

  const int q0 = alpha.size();
  const int p0 = beta.size();
  const int ntot = n + spin + q0 + p0;

  // Generate standard normal innovations
  NumericVector eps = rnorm(ntot, 0.0, 1.0);

  // Initialize
  NumericVector sig2(ntot, 1.0);
  NumericVector asq(ntot, 0.0);

  // GARCH recursion
  const int start = std::max(q0, p0);
  for (int t = start; t < ntot; ++t) {
    double s2 = alpha0;

    // ARCH terms
    for (int i = 0; i < q0; ++i) {
      s2 += alpha[i] * asq[t - i - 1] * asq[t - i - 1];
    }

    // GARCH terms
    for (int j = 0; j < p0; ++j) {
      s2 += beta[j] * sig2[t - j - 1];
    }

    sig2[t] = s2;
    asq[t] = eps[t] * std::sqrt(s2);
  }

  // Extract after spin-up
  const int start_idx = spin + q0 + p0;
  NumericVector result(n);
  for (int i = 0; i < n; ++i) {
    result[i] = asq[start_idx + i];
  }

  return result;
}
