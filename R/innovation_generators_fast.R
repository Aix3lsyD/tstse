#' @title Fast Innovation Generators (dqrng-accelerated)
#' @description Factory functions for creating fast innovation generators
#'   using dqrng's xoshiro256+ PRNG. Drop-in replacements for the standard
#'   generators with \code{_fast} suffix. Approximately 2-3x faster for
#'   simple distributions; 5-10x faster for GARCH (bypasses rugarch).
#' @name innovation_generators_fast
#' @keywords internal
NULL


# ==============================================================================
# Normal Innovation Generator (Fast)
# ==============================================================================

#' Create a Fast Normal Innovation Generator
#'
#' Uses dqrng Box-Muller C++ implementation.
#'
#' @param sd Standard deviation (default 1).
#' @return A function with signature \code{function(n)}.
#' @export
make_gen_norm_fast <- function(sd = 1) {
  if (sd <= 0) stop("sd must be positive")
  force(sd)

  function(n) {
    seed <- sample.int(.Machine$integer.max, 1L)
    gen_norm_fast_cpp(n, sd, seed)
  }
}


# ==============================================================================
# Student's t Innovation Generator (Fast)
# ==============================================================================

#' Create a Fast Student's t Innovation Generator
#'
#' Uses C++ Normal/Gamma representation with dqrng.
#'
#' @param df Degrees of freedom. Must be > 2 when \code{scale = TRUE}.
#' @param scale Logical. If \code{TRUE}, scale to unit variance.
#' @return A function with signature \code{function(n)}.
#' @export
make_gen_t_fast <- function(df, scale = FALSE) {
  if (df <= 0) stop("df must be positive")

  if (df <= 2 && scale) {
    warning("df <= 2 has infinite variance; scale set to FALSE")
    scale <- FALSE
  }

  force(df)
  force(scale)

  function(n) {
    seed <- sample.int(.Machine$integer.max, 1L)
    gen_t_fast_cpp(n, df, scale, seed)
  }
}


# ==============================================================================
# Skew-t Innovation Generator (Fast)
# ==============================================================================

#' Create a Fast Skew-t Innovation Generator
#'
#' Uses C++ Azzalini representation with dqrng. No sn dependency.
#'
#' @param df Degrees of freedom. Must be > 2 when \code{scale = TRUE}.
#' @param alpha Slant/skewness parameter. Default 0 (symmetric).
#' @param scale Logical. If \code{TRUE}, standardize using theoretical
#'   moments (to match \code{make_gen_skt}). Default \code{FALSE}.
#' @return A function with signature \code{function(n)}.
#' @export
make_gen_skt_fast <- function(df, alpha = 0, scale = FALSE) {
  if (df <= 0) stop("df must be positive")

  if (df <= 2 && scale) {
    warning("df <= 2 has infinite variance; scale set to FALSE")
    scale <- FALSE
  }

  if (df <= 1 && alpha != 0) {
    warning("df <= 1 with alpha != 0: mean is undefined; innovations will not be centered")
  }

  theo_mean <- 0
  theo_sd <- 1

  if (df > 1 && alpha != 0) {
    delta <- alpha / sqrt(1 + alpha^2)
    b_nu <- sqrt(df / pi) * exp(lgamma((df - 1) / 2) - lgamma(df / 2))
    theo_mean <- b_nu * delta
  }

  if (scale && df > 2) {
    theo_var <- df / (df - 2) - theo_mean^2
    if (!is.finite(theo_var) || theo_var <= 0) {
      stop("Could not compute finite theoretical variance for given alpha/df")
    }
    theo_sd <- sqrt(theo_var)
  }

  force(df)
  force(alpha)
  force(scale)
  force(theo_mean)
  force(theo_sd)

  function(n) {
    seed <- sample.int(.Machine$integer.max, 1L)
    x <- gen_skt_fast_cpp(n, df, alpha, FALSE, seed)
    (x - theo_mean) / theo_sd
  }
}


# ==============================================================================
# GED Innovation Generator (Fast)
# ==============================================================================

#' Create a Fast GED Innovation Generator
#'
#' Uses C++ Gamma-based generation with dqrng. No fGarch dependency.
#'
#' @param nu Shape parameter. \code{nu = 2} gives normal. Default 2.
#' @param sd Standard deviation. Default 1.
#' @return A function with signature \code{function(n)}.
#' @export
make_gen_ged_fast <- function(nu = 2, sd = 1) {
  if (nu <= 0) stop("nu must be positive")
  if (sd <= 0) stop("sd must be positive")

  force(nu)
  force(sd)

  function(n) {
    seed <- sample.int(.Machine$integer.max, 1L)
    gen_ged_fast_cpp(n, nu, sd, seed)
  }
}


# ==============================================================================
# Laplace Innovation Generator (Fast)
# ==============================================================================

#' Create a Fast Laplace Innovation Generator
#'
#' Uses C++ inverse CDF with dqrng uniforms.
#'
#' @param scale Scale parameter. Default \code{1/sqrt(2)} for unit variance.
#' @return A function with signature \code{function(n)}.
#' @export
make_gen_laplace_fast <- function(scale = 1 / sqrt(2)) {
  if (scale <= 0) stop("scale must be positive")
  force(scale)

  function(n) {
    seed <- sample.int(.Machine$integer.max, 1L)
    gen_laplace_fast_cpp(n, scale, seed)
  }
}


# ==============================================================================
# Uniform Innovation Generator (Fast)
# ==============================================================================

#' Create a Fast Uniform Innovation Generator
#'
#' Uses dqrng xoshiro256+ uniform generation.
#'
#' @param half_width Half the width of the uniform range.
#'   Default \code{sqrt(3)} for unit variance.
#' @return A function with signature \code{function(n)}.
#' @export
make_gen_unif_fast <- function(half_width = sqrt(3)) {
  if (half_width <= 0) stop("half_width must be positive")
  force(half_width)

  function(n) {
    seed <- sample.int(.Machine$integer.max, 1L)
    gen_unif_fast_cpp(n, half_width, seed)
  }
}


# ==============================================================================
# Mixture of Normals Innovation Generator (Fast)
# ==============================================================================

#' Create a Fast Mixture of Normals Innovation Generator
#'
#' Uses C++ dqrng component selection and normals in one pass.
#'
#' @param sd1 Standard deviation of component 1 (default 1).
#' @param sd2 Standard deviation of component 2 (default 3).
#' @param prob1 Probability of component 1 (default 0.9).
#' @return A function with signature \code{function(n)}.
#' @export
make_gen_mixnorm_fast <- function(sd1 = 1, sd2 = 3, prob1 = 0.9) {
  if (sd1 <= 0 || sd2 <= 0) stop("sd1 and sd2 must be positive")
  if (prob1 <= 0 || prob1 >= 1) stop("prob1 must be between 0 and 1 (exclusive)")

  force(sd1)
  force(sd2)
  force(prob1)

  function(n) {
    seed <- sample.int(.Machine$integer.max, 1L)
    gen_mixnorm_fast_cpp(n, sd1, sd2, prob1, seed)
  }
}


# ==============================================================================
# GARCH Innovation Generator (Fast)
# ==============================================================================

#' Create a Fast GARCH Innovation Generator
#'
#' Uses direct C++ GARCH recursion with dqrng. Bypasses rugarch entirely.
#' Only supports sGARCH with normal innovations.
#'
#' @param omega Constant term in variance equation. Must be positive.
#' @param alpha Numeric vector of ARCH coefficients.
#' @param beta Numeric vector of GARCH coefficients (default empty).
#' @param spin Burn-in length (default 1000).
#' @return A function with signature \code{function(n)}.
#' @export
make_gen_garch_fast <- function(omega, alpha, beta = numeric(0), spin = 1000L) {
  if (omega <= 0) stop("omega must be positive")
  if (any(alpha < 0)) stop("alpha coefficients must be non-negative")
  if (length(beta) > 0 && any(beta < 0)) stop("beta coefficients must be non-negative")
  if (sum(alpha) + sum(beta) >= 1) {
    warning("sum(alpha) + sum(beta) >= 1: process may be non-stationary")
  }

  alpha <- as.numeric(alpha)
  beta <- as.numeric(beta)
  force(omega)
  force(alpha)
  force(beta)
  force(spin)

  function(n) {
    seed <- sample.int(.Machine$integer.max, 1L)
    gen_garch_fast_cpp(n, omega, alpha, beta, as.integer(spin), seed)
  }
}


# ==============================================================================
# Heteroscedastic Innovation Generator (Fast)
# ==============================================================================

#' Create a Fast Heteroscedastic Innovation Generator
#'
#' Uses C++ dqrng normals with precomputed weight vector.
#' Supports the same shape specifications as \code{make_gen_hetero}.
#'
#' @param w Weight vector, function, or NULL. If NULL and shape is NULL,
#'   defaults to linear weights from 1 to 10.
#' @param sd Standard deviation of base normal. Default 1.
#' @param shape Named shape: "linear", "sqrt", "log", "power", "exp",
#'   "step", "periodic".
#' @param from Start value for shape (default 1).
#' @param to End value for shape (default 10).
#' @param power Power for "power" shape (default 2).
#' @param breaks Break points for "step" shape.
#' @param levels Levels for "step" shape.
#' @param base_w Base weight for "periodic" shape (default 1).
#' @param amplitude Amplitude for "periodic" shape (default 0.5).
#' @param period Period for "periodic" shape (default 12).
#' @param base Optional base innovation generator function. When NULL,
#'   uses fast normal base with \code{sd}. Mutually exclusive with non-default
#'   \code{sd}.
#' @return A function with signature \code{function(n)}.
#' @export
make_gen_hetero_fast <- function(w = NULL, sd = 1, shape = NULL,
                                  from = 1, to = 10, power = 2,
                                  breaks = NULL, levels = NULL,
                                  base_w = 1, amplitude = 0.5, period = 12,
                                  base = NULL) {
  if (!is.null(shape) && !is.null(w)) {
    stop("Cannot specify both 'shape' and 'w'.")
  }
  if (!is.null(base) && !identical(sd, 1)) {
    stop("Cannot specify both 'sd' (non-default) and 'base'.")
  }
  if (is.numeric(sd) && sd <= 0) stop("sd must be positive")
  if (!is.null(base) && !is.function(base)) {
    stop("'base' must be a function")
  }
  if (!is.null(w) && !is.function(w) && !is.numeric(w)) {
    stop("w must be NULL, a function, or a numeric vector")
  }

  if (!is.null(base)) {
    base_gen <- base
  } else {
    base_gen <- NULL
  }

  # Build weight function (same as make_gen_hetero)
  if (!is.null(shape)) {
    w_fun <- .hetero_shape_weights(shape, from, to, power, breaks, levels,
                                    base_w, amplitude, period)
  } else {
    w_fun <- NULL
  }

  force(w)
  force(sd)
  force(base_gen)
  force(w_fun)

  get_weights <- function(n) {
    if (!is.null(w_fun)) {
      weights <- w_fun(n)
    } else if (is.null(w)) {
      weights <- seq(1, 10, length.out = n)
    } else if (is.function(w)) {
      weights <- w(n)
      if (length(weights) != n) stop("Weight function must return exactly n values")
    } else {
      if (length(w) < n) {
        weights <- c(rep(w[1], n - length(w)), w)
      } else {
        weights <- w[seq_len(n)]
      }
    }
    weights
  }

  gen <- function(n) {
    weights <- get_weights(n)
    if (!is.null(base_gen)) {
      weights * base_gen(n)
    } else {
      seed <- sample.int(.Machine$integer.max, 1L)
      gen_hetero_fast_cpp(n, weights, sd, seed)
    }
  }

  attr(gen, "tstse_innov_kind") <- "hetero"
  attr(gen, "tstse_hetero_base_gen") <- if (!is.null(base_gen)) base_gen else function(n) {
    seed <- sample.int(.Machine$integer.max, 1L)
    gen_norm_fast_cpp(n, sd, seed)
  }
  attr(gen, "tstse_hetero_weight_builder") <- get_weights
  gen
}
