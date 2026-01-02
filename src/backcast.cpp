#include <Rcpp.h>
using namespace Rcpp;

//' Backcast residuals (C++ implementation)
//'
//' @param x Numeric vector, the time series.
//' @param phi Numeric vector, AR coefficients.
//' @param theta Numeric vector, MA coefficients.
//' @param n_back Integer, how far back to backcast.
//' @return Numeric vector of residuals.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
NumericVector backcast_cpp(NumericVector x,
                            NumericVector phi,
                            NumericVector theta,
                            int n_back = 50) {
   const int n = x.size();

   // Effective AR order (0 if phi is all zeros)
   int p = phi.size();
   double sum_phi_sq = sum(phi * phi);
   if (sum_phi_sq == 0.0) p = 0;

   // Effective MA order (0 if theta is all zeros)
   int q = theta.size();
   double sum_theta_sq = sum(theta * theta);
   if (sum_theta_sq == 0.0) q = 0;

   // Mean and constant
   const double xbar = mean(x);
   double cst = 1.0;
   for (int j = 0; j < p; ++j) cst -= phi[j];
   const double maconst = cst * xbar;

   // Backcast residuals
   NumericVector deltar(n, 0.0);
   const int mpq = std::max(p, q);
   const int np = n - mpq;

   for (int i = np - 1; i >= 0; --i) {
     double d = x[i];

     for (int jp = 0; jp < p; ++jp) {
       d -= phi[jp] * x[i + jp + 1];
     }
     for (int jq = 0; jq < q; ++jq) {
       d += theta[jq] * deltar[i + jq + 1];
     }

     deltar[i] = d - maconst;
   }

   // Backcast X(0), X(-1), ..., X(-n_back)
   const int n_full = n + n_back + 1;
   NumericVector xhat(n_full, 0.0);
   const int mm = n_back + 1;

   // Copy x into xhat
   for (int i = 0; i < n; ++i) {
     xhat[mm + i] = x[i];
   }

   // Backcast loop
   for (int h = 1; h <= mm; ++h) {
     const int idx = mm - h;
     double val = xhat[idx];

     for (int jp = 0; jp < p; ++jp) {
       val += phi[jp] * xhat[idx + jp + 1];
     }

     if (h <= q) {
       const double th = theta[h - 1];
       for (int jq = h; jq <= q; ++jq) {
         val -= th * deltar[jq - 1];
       }
     }

     xhat[idx] = val + maconst;
   }

   // Calculate residuals after backcasting
   NumericVector resid(n_full, 0.0);
   const int p1q1 = std::max(p + 1, q + 1);

   for (int i = p1q1 - 1; i < n; ++i) {
     xhat[mm + i] = x[i];
   }

   for (int i = p1q1 - 1; i < n_full; ++i) {
     double r = xhat[i];

     for (int jp = 0; jp < p; ++jp) {
       r -= phi[jp] * xhat[i - jp - 1];
     }
     for (int jq = 0; jq < q; ++jq) {
       r += theta[jq] * resid[i - jq - 1];
     }

     resid[i] = r - maconst;
   }

   // Extract final residuals
   NumericVector residb(n);
   for (int i = 0; i < n; ++i) {
     residb[i] = resid[mm + i];
   }

   return residb;
 }
