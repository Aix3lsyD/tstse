// =============================================================================
// FILE: kernel_gen_arma.cpp
// CATEGORY: ARMA filter for gen_aruma_flex
//
// Applies ARMA recursion to pre-generated innovations, replacing arima.sim().
// Accepts innovations from any R generator (normal, t, laplace, mixture, etc.)
// and applies MA convolution + AR recursive filter + optional inverse differencing.
//
// Key function:
//   - arma_filter_cpp(): Exported to R, called by gen_aruma_flex()
//
// Replicates arima.sim() behavior exactly:
//   1. MA filter: y[t] = x[t] + ma[1]*x[t-1] + ... (sides=1, then zero first q)
//   2. AR filter: z[t] = phi[1]*z[t-1] + ... + y[t]  (recursive)
//   3. Trim burn-in
//   4. Inverse differencing (diffinv): prepend 0 + cumsum, d times
// =============================================================================

#include <Rcpp.h>
#include <algorithm>

using namespace Rcpp;

//' C++ ARMA Filter (replaces arima.sim in gen_aruma_flex)
//'
//' Applies MA convolution, AR recursion, burn-in trimming, and inverse
//' differencing to a vector of pre-generated innovations. Produces output
//' identical to arima.sim() for the same inputs.
//'
//' @param innovations Numeric vector of length n + n_start (burn-in + main).
//' @param phi Numeric vector of AR coefficients (package convention).
//' @param theta Numeric vector of MA coefficients (package/textbook convention;
//'   negated internally to match arima.sim's "+" convention).
//' @param n Integer, desired output length (before differencing).
//' @param n_start Integer, burn-in length.
//' @param d Integer, differencing order (0 = none).
//' @return Numeric vector of length n (if d=0) or n+d (if d>0).
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
NumericVector arma_filter_cpp(NumericVector innovations,
                              NumericVector phi,
                              NumericVector theta,
                              int n,
                              int n_start,
                              int d) {

  const int p = phi.size();
  const int q = theta.size();
  const int n_total = n + n_start;

  // Guard: MAX_P and MAX_Q
  constexpr int MAX_P = 20;
  constexpr int MAX_Q = 20;
  if (p > MAX_P || q > MAX_Q) {
    stop("AR order (%d) or MA order (%d) exceeds maximum (%d)", p, q,
         std::max(MAX_P, MAX_Q));
  }

  if (innovations.size() < n_total) {
    stop("innovations length (%d) < n + n_start (%d)", innovations.size(), n_total);
  }

  // --- Step 1: MA filter ---
  // R's filter(x, c(1, model$ma), sides=1) where model$ma = -theta
  // y[t] = 1*x[t] + (-theta[1])*x[t-1] + ... + (-theta[q])*x[t-q]
  // Then zero out y[0:q-1]
  //
  // IMPORTANT: MA filter reads from ORIGINAL innovations, writes to separate buffer.
  // R's filter(sides=1) creates a new vector while reading from input.

  // Working buffer: start with a copy of innovations
  NumericVector buf(n_total);

  if (q > 0) {
    // MA pass: read from innovations, write to buf
    // Zero out first q elements (matching arima.sim's x[seq_along(model$ma)] <- 0)
    for (int t = 0; t < q; ++t) {
      buf[t] = 0.0;
    }
    // Apply MA filter for t >= q
    for (int t = q; t < n_total; ++t) {
      double val = innovations[t];  // coefficient 1 * x[t]
      for (int k = 0; k < q; ++k) {
        // model$ma = -theta, so: val += (-theta[k]) * x[t-k-1]
        val -= theta[k] * innovations[t - k - 1];
      }
      buf[t] = val;
    }
  } else {
    // No MA: just copy innovations
    for (int t = 0; t < n_total; ++t) {
      buf[t] = innovations[t];
    }
  }

  // --- Step 2: AR filter ---
  // R's filter(x, phi, method="recursive")
  // z[t] = phi[1]*z[t-1] + ... + phi[p]*z[t-p] + x[t]
  // In-place: reads previously-written values from buf

  if (p > 0) {
    for (int t = 0; t < n_total; ++t) {
      for (int j = 0; j < p; ++j) {
        int lag = t - j - 1;
        if (lag >= 0) {
          buf[t] += phi[j] * buf[lag];
        }
        // lag < 0: initial conditions are 0 (matching R's filter behavior)
      }
    }
  }

  // --- Step 3: Trim burn-in ---
  // Extract buf[n_start .. n_start+n-1]
  NumericVector result(n);
  for (int i = 0; i < n; ++i) {
    result[i] = buf[n_start + i];
  }

  // --- Step 4: Inverse differencing (if d > 0) ---
  // Replicates diffinv(x, differences = d):
  //   For each level: prepend 0, then cumsum
  //   Output length grows by 1 per level
  if (d > 0) {
    for (int level = 0; level < d; ++level) {
      int cur_len = result.size();
      NumericVector expanded(cur_len + 1);
      expanded[0] = 0.0;
      double cumsum = 0.0;
      for (int i = 0; i < cur_len; ++i) {
        cumsum += result[i];
        expanded[i + 1] = cumsum;
      }
      result = expanded;
    }
  }

  return result;
}
