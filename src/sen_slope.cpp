#include <Rcpp.h>
using namespace Rcpp;

//' Compute Sen's Slope Estimator (C++ implementation)
//'
//' Computes the median of all pairwise slopes for robust trend estimation.
//' This is an O(n^2/2) operation optimized in C++ for speed.
//'
//' @param x Numeric vector of observations
//' @return The Sen's slope estimate (median of pairwise slopes)
//'
//' @noRd
// [[Rcpp::export]]
double sen_slope_cpp(NumericVector x) {
  int n = x.size();

  if (n < 2) {
    return NA_REAL;
  }

  // Number of pairwise slopes: n*(n-1)/2
int n_slopes = n * (n - 1) / 2;
  NumericVector slopes(n_slopes);

  int k = 0;
  for (int i = 0; i < n - 1; i++) {
    for (int j = i + 1; j < n; j++) {
      // slope = (x[j] - x[i]) / (j - i)
      // Note: j > i always, so (j - i) > 0
      slopes[k++] = (x[j] - x[i]) / (double)(j - i);
    }
  }

  // Use Rcpp's median function (calls R's median internally)
  return median(slopes);
}
