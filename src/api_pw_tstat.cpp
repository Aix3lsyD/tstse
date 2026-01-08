// =============================================================================
// FILE: api_pw_tstat.cpp
// CATEGORY: INTERFACE (R-facing)
// THREAD-SAFE: NO (uses Rcpp types, not for parallel use)
//
// Public Rcpp exports for Prais-Winsten t-statistic.
// For parallel bootstrap, use kernel_pw_tstat.cpp functions instead.
//
// Exports:
//   - pw_tstat_cpp(): Primary PW t-statistic API
//   - pw_full_cpp(): Full results including rho
// =============================================================================
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "kernel_types.h"
using namespace Rcpp;

// Forward declarations of pure C++ functions (defined in kernel_pw_tstat.cpp)
double pw_tstat_pure(const arma::vec& x);
PWResult pw_full_pure(const arma::vec& x);
double pw_tstat_iterative_pure(const arma::vec& x, int max_iter, double tol);
PWResult pw_full_iterative_pure(const arma::vec& x, int max_iter, double tol);

//' Prais-Winsten t-statistic (C++ implementation)
//'
//' Computes the Prais-Winsten t-statistic for testing H0: slope = 0.
//' This is a single-pass implementation using AR(1) quasi-differencing.
//'
//' @param x Numeric vector, the time series.
//' @return Double, Prais-Winsten t-statistic.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double pw_tstat_cpp(const arma::vec& x) {
    return pw_tstat_pure(x);
}

//' Prais-Winsten Full Results (C++ implementation)
//'
//' Returns full PW results including AR(1) coefficient rho.
//' For debugging and validation purposes.
//'
//' @param x Numeric vector, the time series.
//' @return List with tpw, rho, and vara.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Rcpp::List pw_full_cpp(const arma::vec& x) {
    PWResult result = pw_full_pure(x);

    return Rcpp::List::create(
        Rcpp::Named("tpw") = result.tpw,
        Rcpp::Named("rho") = result.rho,
        Rcpp::Named("vara") = result.vara
    );
}

//' Iterative Prais-Winsten t-statistic (C++ implementation)
//'
//' Computes the Prais-Winsten t-statistic using iterative estimation.
//' Iterates until rho converges or max_iter is reached.
//'
//' @param x Numeric vector, the time series.
//' @param max_iter Integer, maximum iterations (default 50).
//' @param tol Double, convergence tolerance for rho (default 1e-6).
//' @return Double, Prais-Winsten t-statistic.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double pw_tstat_iterative_cpp(const arma::vec& x, int max_iter = 50, double tol = 1e-6) {
    return pw_tstat_iterative_pure(x, max_iter, tol);
}

//' Iterative Prais-Winsten Full Results (C++ implementation)
//'
//' Returns full iterative PW results including AR(1) coefficient rho.
//'
//' @param x Numeric vector, the time series.
//' @param max_iter Integer, maximum iterations (default 50).
//' @param tol Double, convergence tolerance for rho (default 1e-6).
//' @return List with tpw, rho, and vara.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Rcpp::List pw_full_iterative_cpp(const arma::vec& x, int max_iter = 50, double tol = 1e-6) {
    PWResult result = pw_full_iterative_pure(x, max_iter, tol);

    return Rcpp::List::create(
        Rcpp::Named("tpw") = result.tpw,
        Rcpp::Named("rho") = result.rho,
        Rcpp::Named("vara") = result.vara
    );
}
