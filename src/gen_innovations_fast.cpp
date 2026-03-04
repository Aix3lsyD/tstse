// =============================================================================
// FILE: gen_innovations_fast.cpp
// CATEGORY: UTILITY - Fast innovation generators using dqrng
//
// Provides C++ implementations of innovation generators using dqrng's
// xoshiro256+ PRNG (~2-3x faster than R's Mersenne Twister).
// Each function takes an explicit uint64_t seed for reproducibility.
//
// Exports:
//   gen_norm_fast_cpp, gen_unif_fast_cpp, gen_laplace_fast_cpp,
//   gen_t_fast_cpp, gen_mixnorm_fast_cpp, gen_hetero_fast_cpp,
//   gen_garch_fast_cpp, gen_skt_fast_cpp, gen_ged_fast_cpp
// =============================================================================
// [[Rcpp::depends(RcppArmadillo, dqrng, BH)]]

#include <RcppArmadillo.h>
#include <dqrng.h>
#include <dqrng_distribution.h>
#include <xoshiro.h>
#include <cmath>

using namespace Rcpp;

// =============================================================================
// Internal helpers
// =============================================================================

// Fill arma::vec with Box-Muller normals (same pattern as kernel_gen_ar.cpp)
static void fill_normals(arma::vec& out, dqrng::xoshiro256plus& rng, double sd) {
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

// Fill arma::vec with uniform values in [lo, hi)
static void fill_uniforms(arma::vec& out, dqrng::xoshiro256plus& rng,
                           double lo, double hi) {
    std::uniform_real_distribution<double> unif(lo, hi);
    for (arma::uword i = 0; i < out.n_elem; ++i) {
        out[i] = unif(rng);
    }
}

// Marsaglia-Tsang (2000) gamma sampling: Gamma(shape, 1)
// Requires shape >= 1. For shape < 1, caller uses Gamma(shape+1) * U^(1/shape)
static double gamma_mt(double shape, dqrng::xoshiro256plus& rng) {
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    const double d = shape - 1.0 / 3.0;
    const double c = 1.0 / std::sqrt(9.0 * d);

    // Standard normal via Box-Muller (single pair)
    auto gen_normal = [&]() -> double {
        double u1, u2;
        u1 = unif(rng);
        while (u1 <= 1e-15) u1 = unif(rng);
        u2 = unif(rng);
        return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
    };

    while (true) {
        double x, v;
        do {
            x = gen_normal();
            v = 1.0 + c * x;
        } while (v <= 0.0);

        v = v * v * v;
        double u = unif(rng);

        // Squeeze test
        if (u < 1.0 - 0.0331 * (x * x) * (x * x)) return d * v;
        if (std::log(u) < 0.5 * x * x + d * (1.0 - v + std::log(v))) return d * v;
    }
}

// Gamma(shape, scale) for arbitrary shape > 0
static double gamma_sample(double shape, double scale, dqrng::xoshiro256plus& rng) {
    if (shape >= 1.0) {
        return gamma_mt(shape, rng) * scale;
    }
    // shape < 1: Gamma(shape) = Gamma(shape+1) * U^(1/shape)
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    double u = unif(rng);
    while (u <= 1e-15) u = unif(rng);
    return gamma_mt(shape + 1.0, rng) * std::pow(u, 1.0 / shape) * scale;
}


// =============================================================================
// 1. Normal generator
// =============================================================================

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::vec gen_norm_fast_cpp(int n, double sd, uint64_t seed) {
    dqrng::xoshiro256plus rng(seed);
    arma::vec out(n);
    fill_normals(out, rng, sd);
    return out;
}


// =============================================================================
// 2. Uniform generator
// =============================================================================

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::vec gen_unif_fast_cpp(int n, double half_width, uint64_t seed) {
    dqrng::xoshiro256plus rng(seed);
    arma::vec out(n);
    fill_uniforms(out, rng, -half_width, half_width);
    return out;
}


// =============================================================================
// 3. Laplace generator (inverse CDF)
// =============================================================================

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::vec gen_laplace_fast_cpp(int n, double scale, uint64_t seed) {
    dqrng::xoshiro256plus rng(seed);
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    arma::vec out(n);
    const double max_abs_u = 0.5 - 1e-15;

    for (int i = 0; i < n; ++i) {
        double u = unif(rng) - 0.5;
        // sign(u) * scale * log(1 - 2*|u|)
        // Matches R: scale * sign(u) * log(1 - 2*abs(u))
        double abs_u = std::fabs(u);
        if (abs_u >= max_abs_u) abs_u = max_abs_u;
        double sign_u = (u >= 0.0) ? 1.0 : -1.0;
        out[i] = sign_u * scale * std::log(1.0 - 2.0 * abs_u);
    }
    return out;
}


// =============================================================================
// 4. Student's t generator via Normal/Gamma
//    t(df) = Z / sqrt(V/df) where Z~N(0,1), V~Gamma(df/2, 1) scaled
// =============================================================================

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::vec gen_t_fast_cpp(int n, double df, bool scale_to_unit, uint64_t seed) {
    dqrng::xoshiro256plus rng(seed);
    arma::vec out(n);

    // Generate normals first
    fill_normals(out, rng, 1.0);

    // Scale by chi2: t = Z * sqrt(df / V) where V ~ Gamma(df/2, 2)
    const double half_df = df / 2.0;
    for (int i = 0; i < n; ++i) {
        double v = gamma_sample(half_df, 1.0, rng);  // Gamma(df/2, 1)
        out[i] *= std::sqrt(half_df / v);  // Z * sqrt(df/2 / V) = Z * sqrt(df / (2V))
        // Actually: t = Z / sqrt(V/df) = Z * sqrt(df/V)
        // V ~ Gamma(df/2, 2) so V = 2*Gamma(df/2,1)
        // t = Z * sqrt(df / (2*Gamma(df/2,1))) = Z * sqrt(half_df / Gamma(df/2,1))
    }

    if (scale_to_unit && df > 2.0) {
        // Scale to unit variance: Var(t) = df/(df-2)
        double scale_factor = std::sqrt((df - 2.0) / df);
        out *= scale_factor;
    }

    return out;
}


// =============================================================================
// 5. Mixture of normals generator
// =============================================================================

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::vec gen_mixnorm_fast_cpp(int n, double sd1, double sd2, double prob1,
                                uint64_t seed) {
    dqrng::xoshiro256plus rng(seed);
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    arma::vec out(n);

    // Generate normals with sd=1, then scale by component
    fill_normals(out, rng, 1.0);

    for (int i = 0; i < n; ++i) {
        double u = unif(rng);
        out[i] *= (u < prob1) ? sd1 : sd2;
    }

    return out;
}


// =============================================================================
// 6. Heteroscedastic generator (precomputed weights * normals)
// =============================================================================

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::vec gen_hetero_fast_cpp(int n, const arma::vec& weights, double sd,
                               uint64_t seed) {
    dqrng::xoshiro256plus rng(seed);
    arma::vec out(n);
    fill_normals(out, rng, sd);

    // Element-wise multiply by weights
    for (int i = 0; i < n; ++i) {
        out[i] *= weights[i];
    }

    return out;
}


// =============================================================================
// 7. GARCH generator (direct recursion, replaces rugarch dependency)
//    Supports sGARCH(q,p): sigma2_t = omega + sum(alpha_i * a_{t-i}^2) + sum(beta_j * sigma2_{t-j})
// =============================================================================

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::vec gen_garch_fast_cpp(int n, double omega,
                              const arma::vec& alpha, const arma::vec& beta,
                              int spin, uint64_t seed) {
    dqrng::xoshiro256plus rng(seed);

    const int q0 = alpha.n_elem;  // ARCH order
    const int p0 = beta.n_elem;   // GARCH order
    const int start = std::max(q0, p0);
    const int ntot = n + spin + start;

    // Generate all normals at once
    arma::vec eps(ntot);
    fill_normals(eps, rng, 1.0);

    // GARCH recursion
    arma::vec sig2(ntot, arma::fill::ones);   // sigma^2
    arma::vec asq(ntot, arma::fill::zeros);   // a_t = eps_t * sigma_t

    for (int t = start; t < ntot; ++t) {
        double s2 = omega;

        // ARCH terms: alpha_i * a_{t-i}^2
        for (int i = 0; i < q0; ++i) {
            double a = asq[t - i - 1];
            s2 += alpha[i] * a * a;
        }

        // GARCH terms: beta_j * sigma2_{t-j}
        for (int j = 0; j < p0; ++j) {
            s2 += beta[j] * sig2[t - j - 1];
        }

        sig2[t] = s2;
        asq[t] = eps[t] * std::sqrt(s2);
    }

    // Extract after spin-up
    const int start_idx = spin + start;
    return asq.subvec(start_idx, start_idx + n - 1);
}


// =============================================================================
// 8. Skew-t generator (Azzalini representation)
//    X = delta*|Z0| + sqrt(1-delta^2)*Z1 where Z0,Z1~N(0,1) independent
//    then X ~ skew-normal, and X * sqrt(df/V) ~ skew-t
//    delta = alpha / sqrt(1 + alpha^2)
// =============================================================================

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::vec gen_skt_fast_cpp(int n, double df, double skew_alpha,
                            bool scale_to_unit, uint64_t seed) {
    dqrng::xoshiro256plus rng(seed);

    const double delta = skew_alpha / std::sqrt(1.0 + skew_alpha * skew_alpha);
    const double sqrt_1md2 = std::sqrt(1.0 - delta * delta);
    const double half_df = df / 2.0;

    // Generate 2n normals for Z0, Z1
    arma::vec z0(n), z1(n);
    fill_normals(z0, rng, 1.0);
    fill_normals(z1, rng, 1.0);

    arma::vec out(n);
    for (int i = 0; i < n; ++i) {
        // Skew-normal: X = delta*|Z0| + sqrt(1-delta^2)*Z1
        double sn = delta * std::fabs(z0[i]) + sqrt_1md2 * z1[i];

        // Skew-t: multiply by sqrt(df/V) where V ~ Gamma(df/2, 2)
        double v = gamma_sample(half_df, 1.0, rng);
        out[i] = sn * std::sqrt(half_df / v);
    }

    if (scale_to_unit) {
        // Empirical centering and scaling (matches sn::rst behavior)
        // Theoretical mean: mu = delta * sqrt(df/pi) * Gamma((df-1)/2) / Gamma(df/2)
        //   when df > 1
        // Theoretical variance: df/(df-2) - mu^2 when df > 2
        // For robustness, use sample moments:
        double mu = arma::mean(out);
        double sigma = arma::stddev(out, 1);  // population stddev
        if (sigma > 1e-15) {
            out = (out - mu) / sigma;
        }
    }

    return out;
}


// =============================================================================
// 9. GED (Generalized Error Distribution) generator
//    |X| ~ Gamma(1/nu, 1)^(1/nu) scaled, sign from uniform
//    Specifically: X = sign * lambda * |Gamma(1/nu, 1)|^(1/nu)
//    where lambda = sqrt(Gamma(1/nu) / Gamma(3/nu)) / sd
// =============================================================================

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::vec gen_ged_fast_cpp(int n, double nu, double sd, uint64_t seed) {
    dqrng::xoshiro256plus rng(seed);
    std::uniform_real_distribution<double> unif(0.0, 1.0);

    // GED via gamma: |X| = Gamma(1/nu, 1)^(1/nu) with random sign
    // Then scale so that Var(X) = sd^2
    // Var of unscaled = Gamma(3/nu) / Gamma(1/nu)
    // So scale = sd * sqrt(Gamma(1/nu) / Gamma(3/nu))
    const double inv_nu = 1.0 / nu;
    double lgam1 = std::lgamma(inv_nu);
    double lgam3 = std::lgamma(3.0 * inv_nu);
    double scale = sd * std::exp(0.5 * (lgam1 - lgam3));

    arma::vec out(n);
    for (int i = 0; i < n; ++i) {
        // Generate |X| ~ Gamma(1/nu, 1)^(1/nu)
        double g = gamma_sample(inv_nu, 1.0, rng);
        double abs_x = std::pow(g, inv_nu);

        // Random sign
        double u = unif(rng);
        out[i] = (u < 0.5 ? -1.0 : 1.0) * scale * abs_x;
    }

    return out;
}
