// burg_aic_pure.cpp - Pure C++ Burg algorithm with AIC selection
// Thread-safe version using C++ structs instead of Rcpp::List
// Matches R's ar.burg() with var.method=1 (recursive variance)
// [[Rcpp::depends(RcppArmadillo)]]

#include "types_pure.h"
#include <cmath>


// Pure C++ Burg with AIC/BIC selection
// Matches R's ar.burg() implementation exactly
// Thread-safe: no Rcpp types, only arma and C++ structs
BurgResult burg_aic_select_pure(const arma::vec& x, int maxp,
                                 const std::string& criterion) {
    const int n = x.n_elem;

    if (maxp >= n - 1) {
        maxp = n - 2;
    }
    if (maxp < 1) {
        maxp = 1;
    }

    // Center the series
    const double x_mean = arma::mean(x);
    arma::vec xc = x - x_mean;

    // Initialize forward and backward prediction errors
    // Note: R reverses data, but due to symmetry this doesn't affect coefficients
    arma::vec ef = xc;
    arma::vec eb = xc;

    // AR(0) variance
    double vara0 = arma::dot(xc, xc) / n;

    // Track best model
    int best_p = 0;
    double best_ic = std::numeric_limits<double>::infinity();
    arma::vec best_phi;
    double best_vara = vara0;

    // IC for AR(0): k = 1 parameter (just variance, or mean if demean=TRUE)
    // R uses: n * log(var) + 2 * p + 2 * demean = n * log(var) + 2 * (p + 1) when demean=TRUE
    // For AR(0), p=0, so IC = n * log(var0) + 2 * 1
    double ic0;
    if (criterion == "aic") {
        ic0 = n * std::log(vara0) + 2.0 * 1;
    } else if (criterion == "aicc") {
        ic0 = n * std::log(vara0) + 2.0 * 1 * n / (n - 2);
    } else { // bic
        ic0 = n * std::log(vara0) + std::log(static_cast<double>(n)) * 1;
    }
    best_ic = ic0;

    // Current variance (recursive formula like R's var1)
    double var_recursive = vara0;

    // Temporary storage for Levinson recursion
    arma::vec a_prev;

    // Single pass through all orders (like R's implementation)
    for (int p = 1; p <= maxp; ++p) {
        // Compute reflection coefficient
        double num = 0.0;
        double den = 0.0;

        for (int t = p; t < n; ++t) {
            num += ef[t] * eb[t - 1];
            den += ef[t] * ef[t] + eb[t - 1] * eb[t - 1];
        }

        double phii = (den > 1e-15) ? 2.0 * num / den : 0.0;

        // Update variance using recursive formula (R's var1)
        var_recursive = var_recursive * (1.0 - phii * phii);

        // Levinson recursion for AR coefficients
        arma::vec a_curr(p);
        a_curr(p - 1) = phii;  // Last coefficient is reflection coefficient

        if (p > 1) {
            for (int j = 0; j < p - 1; ++j) {
                a_curr(j) = a_prev(j) - phii * a_prev(p - 2 - j);
            }
        }

        a_prev = a_curr;

        // Update prediction errors in-place using reverse loop
        // Going backwards preserves correctness: at step t, we need eb(t-1)
        // which hasn't been modified yet since we process higher indices first
        for (int t = n - 1; t >= p; --t) {
            double old_ef = ef(t);  // Save before overwriting
            ef(t) = old_ef - phii * eb(t - 1);
            eb(t) = eb(t - 1) - phii * old_ef;
        }

        // Compute information criterion
        // k = p + 1 (p AR coefficients + 1 for mean/variance)
        int k = p + 1;
        double ic;

        if (var_recursive <= 0) {
            ic = std::numeric_limits<double>::infinity();
        } else if (criterion == "aic") {
            ic = n * std::log(var_recursive) + 2.0 * k;
        } else if (criterion == "aicc") {
            if (n - k - 1 > 0) {
                ic = n * std::log(var_recursive) + 2.0 * k * n / (n - k - 1);
            } else {
                ic = std::numeric_limits<double>::infinity();
            }
        } else { // bic
            ic = n * std::log(var_recursive) + std::log(static_cast<double>(n)) * k;
        }

        if (ic < best_ic) {
            best_ic = ic;
            best_p = p;
            best_phi = a_curr;
            best_vara = var_recursive;
        }
    }

    return BurgResult(best_phi, best_vara, best_p, best_ic);
}


// Standalone Burg fit for a specific order (for burg_fit_cpp compatibility)
// This version uses the direct variance estimate from residuals
void burg_fit_pure(const arma::vec& x, int p, arma::vec& phi_out, double& vara_out) {
    const int n = x.n_elem;

    if (p <= 0 || p >= n) {
        phi_out.reset();
        vara_out = arma::dot(x, x) / n;
        return;
    }

    // Center the series
    const double x_mean = arma::mean(x);
    arma::vec xc = x - x_mean;

    // Initialize forward and backward prediction errors
    arma::vec ef = xc;
    arma::vec eb = xc;

    // AR coefficients
    arma::vec a_prev;

    for (int k = 1; k <= p; ++k) {
        // Compute reflection coefficient
        double num = 0.0;
        double den = 0.0;

        for (int t = k; t < n; ++t) {
            num += ef[t] * eb[t - 1];
            den += ef[t] * ef[t] + eb[t - 1] * eb[t - 1];
        }

        double phii = (den > 1e-15) ? 2.0 * num / den : 0.0;

        // Levinson recursion
        arma::vec a_curr(k);
        a_curr(k - 1) = phii;

        if (k > 1) {
            for (int j = 0; j < k - 1; ++j) {
                a_curr(j) = a_prev(j) - phii * a_prev(k - 2 - j);
            }
        }
        a_prev = a_curr;

        // Update prediction errors in-place using reverse loop
        for (int t = n - 1; t >= k; --t) {
            double old_ef = ef(t);
            ef(t) = old_ef - phii * eb(t - 1);
            eb(t) = eb(t - 1) - phii * old_ef;
        }
    }

    phi_out = a_prev;

    // Compute residual variance from forward errors (direct estimate)
    double ss = 0.0;
    for (int t = p; t < n; ++t) {
        ss += ef[t] * ef[t];
    }
    vara_out = ss / (n - p);
}


// Rcpp export wrapper (for R interface)
// [[Rcpp::export]]
Rcpp::List burg_aic_select_pure_export(const arma::vec& x, int maxp = 5,
                                        std::string criterion = "aic") {
    BurgResult result = burg_aic_select_pure(x, maxp, criterion);

    return Rcpp::List::create(
        Rcpp::Named("phi") = result.phi,
        Rcpp::Named("vara") = result.vara,
        Rcpp::Named("p") = result.p,
        Rcpp::Named("ic") = result.ic
    );
}
