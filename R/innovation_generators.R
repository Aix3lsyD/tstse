#' @title Innovation Generators
#' @description Factory functions for creating innovation generators for time series simulation.
#' @name innovation_generators
#' @importFrom stats rt runif rbinom
#' @keywords internal
NULL


# ==============================================================================
# Normal Innovation Generator
# ==============================================================================

#' Create a Normal Innovation Generator
#'
#' Factory function that returns a normal (Gaussian) innovation generator
#' with mean zero.
#'
#' @param sd Standard deviation (default 1).
#'
#' @return A function with signature \code{function(n)} that generates
#'   \code{n} normal innovations with mean zero.
#'
#' @export
#'
#' @examples
#' # Default standard normal
#' norm_gen <- make_gen_norm()
#' innovations <- norm_gen(100)
#'
#' # Higher variance innovations
#' norm_gen_hv <- make_gen_norm(sd = 2)
#' innovations_hv <- norm_gen_hv(100)
make_gen_norm <- function(sd = 1) {
  if (sd <= 0) {
    stop("sd must be positive")
  }

  force(sd)

  function(n) {
    rnorm(n, mean = 0, sd = sd)
  }
}


# ==============================================================================
# Student's t Innovation Generator
# ==============================================================================

#' Create a Student's t Innovation Generator
#'
#' Factory function that returns a t-distributed innovation generator.
#' Optionally scales output to have unit variance for comparability
#' with normal innovations.
#'
#' @param df Degrees of freedom. Must be > 2 for finite variance when
#'   \code{scale = TRUE}.
#' @param scale Logical. If \code{TRUE}, scale to unit variance.
#'   Requires \code{df > 2}. Default is \code{FALSE}.
#'
#' @return A function with signature \code{function(n)} that generates
#'   \code{n} t-distributed innovations with mean zero.
#'
#' @details
#' The variance of a t-distribution with \code{df} degrees of freedom is
#' \code{df / (df - 2)} for \code{df > 2}. When \code{scale = TRUE}, the
#' output is divided by \code{sqrt(df / (df - 2))} to achieve unit variance.
#'
#' For \code{df <= 2}, the variance is infinite, so scaling is not possible
#' and will be disabled with a warning.
#'
#' Heavy-tailed innovations are useful for simulating financial returns
#' and other processes where extreme values occur more frequently than
#' a normal distribution predicts.
#'
#' @export
#'
#' @examples
#' # Heavy-tailed innovations (df = 5)
#' t_gen <- make_gen_t(df = 5)
#' innovations <- t_gen(1000)
#' var(innovations)
#'
#' # Scaled to unit variance
#' t_gen_scaled <- make_gen_t(df = 5, scale = TRUE)
#' var(t_gen_scaled(1000))
#'
#' # Very heavy tails (df = 3)
#' t_gen_heavy <- make_gen_t(df = 3)
make_gen_t <- function(df, scale = FALSE) {
  if (df <= 0) {
    stop("df must be positive")
  }

  if (df <= 2 && scale) {
    warning("df <= 2 has infinite variance; scale set to FALSE")
    scale <- FALSE
  }

  force(df)
  force(scale)

  function(n) {
    x <- rt(n, df = df)
    if (scale && df > 2) {
      x <- x * sqrt((df - 2) / df)
    }
    x
  }
}


# ==============================================================================
# Skew-t Innovation Generator
# ==============================================================================

#' Create a Skew-t Innovation Generator
#'
#' Factory function that returns a skew-t distributed innovation generator
#' using the \pkg{sn} package. Combines heavy tails with asymmetry.
#'
#' @param df Degrees of freedom (tail heaviness). Must be > 2 for finite
#'   variance when \code{scale = TRUE}. Called \code{nu} in some
#'   parameterizations.
#' @param alpha Slant/skewness parameter. \code{alpha = 0} gives symmetric t.
#'   Positive values give right skew, negative gives left skew. Default is 0.
#' @param scale Logical. If \code{TRUE}, standardize to zero mean and unit
#'   variance using theoretical moments. Default is \code{FALSE}.
#'
#' @return A function with signature \code{function(n)} that generates
#'   \code{n} skew-t distributed innovations with mean zero.
#'
#' @details
#' The skew-t distribution is widely used in financial econometrics to
#' model asset returns that exhibit both heavy tails and asymmetry
#' (e.g., markets tend to fall faster than they rise).
#'
#' When \code{alpha = 0}, the distribution reduces to a symmetric
#' Student's t distribution.
#'
#' \strong{Standardization approach}: When \code{scale = TRUE}, this
#' implementation standardizes using \emph{theoretical} mean and variance
#' (via \code{sn::st.cumulants}), rather than sample statistics. This
#' preserves the iid property of the innovations, which is essential for
#' valid time series simulation.
#'
#' When \code{scale = FALSE} and \code{alpha != 0}, the output is still
#' centered to zero mean using the theoretical mean, but variance is not
#' standardized.
#'
#' @seealso \code{\link[sn]{rst}} for the underlying random generator.
#' @seealso \code{\link[sn]{st.cumulants}} for theoretical cumulants.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Symmetric heavy-tailed (equivalent to t)
#' skt_sym <- make_gen_skt(df = 5, alpha = 0)
#'
#' # Left-skewed (negative shocks more extreme)
#' skt_left <- make_gen_skt(df = 5, alpha = -2)
#' innovations <- skt_left(1000)
#' hist(innovations, breaks = 50)
#'
#' # Right-skewed, standardized to unit variance
#' skt_right <- make_gen_skt(df = 5, alpha = 2, scale = TRUE)
#' innovations_scaled <- skt_right(1000)
#' var(innovations_scaled)
#' }
make_gen_skt <- function(df, alpha = 0, scale = FALSE) {

  if (!requireNamespace("sn", quietly = TRUE)) {
    stop("Package 'sn' is required. Install with: install.packages('sn')")
  }

  # Input validation
  if (!is.numeric(df) || length(df) != 1L || is.na(df) || df <= 0) {
    stop("df must be a single positive numeric value")
  }

  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha)) {
    stop("alpha must be a single numeric value")
  }
  if (!is.logical(scale) || length(scale) != 1L || is.na(scale)) {
    stop("scale must be TRUE or FALSE")
  }

  # Variance requires df > 2
  if (df <= 2 && scale) {
    warning("df <= 2 has infinite variance; scale set to FALSE")
    scale <- FALSE
  }

  # Mean requires df > 1 for centering when alpha != 0
  if (df <= 1 && alpha != 0) {
    warning("df <= 1 with alpha != 0: mean is undefined; innovations will not be centered")
  }

  # Precompute theoretical moments (avoids repeated computation in generator)
  # This is done at factory time, not at each call to the generator
  theo_mean <- 0
  theo_sd <- 1

  if (df > 1 && alpha != 0) {
    # Theoretical mean for centering
    # Formula: E[X] = omega * b_nu * delta, where omega = 1, xi = 0
    delta <- alpha / sqrt(1 + alpha^2)
    b_nu <- sqrt(df / pi) * exp(lgamma((df - 1) / 2) - lgamma(df / 2))
    theo_mean <- b_nu * delta
  }

  if (scale && df > 2) {
    # Get theoretical variance for standardization
    # Use st.cumulants which returns (mean, variance, ...)
    cumulants <- tryCatch(
      sn::st.cumulants(xi = 0, omega = 1, alpha = alpha, nu = df, n = 2),
      error = function(e) c(NA_real_, NA_real_)
    )

    theo_var <- as.numeric(cumulants[2])

    if (!is.finite(theo_var) || theo_var <= 0) {
      stop("Could not compute finite theoretical variance for given alpha/df")
    }

    theo_sd <- sqrt(theo_var)
    # Use cumulants mean for consistency (should match our formula)
    theo_mean <- as.numeric(cumulants[1])
  }

  # Freeze all values for closure
  force(df)
  force(alpha)
  force(scale)
  force(theo_mean)
  force(theo_sd)

  # Return generator function
  function(n) {
    if (!is.numeric(n) || length(n) != 1L || n < 1) {
      stop("n must be a positive integer")
    }

    # Generate raw skew-t variates
    # sn::rst uses (xi, omega, alpha, nu) parameterization
    # xi = location (0), omega = scale (1), alpha = slant, nu = df
    x <- sn::rst(n, xi = 0, omega = 1, alpha = alpha, nu = df)

    # Apply theoretical centering and scaling
    # theo_mean = 0 and theo_sd = 1 when no adjustment needed
    (x - theo_mean) / theo_sd
  }
}


# ==============================================================================
# Generalized Error Distribution (GED) Innovation Generator
# ==============================================================================

#' Create a Generalized Error Distribution (GED) Innovation Generator
#'
#' Factory function that returns a GED innovation generator using the
#' \pkg{fGarch} package. The GED allows flexible control over tail weight.
#'
#' @param nu Shape parameter controlling tail weight. \code{nu = 2} gives
#'   exactly the normal distribution. \code{nu < 2} gives heavier tails,
#'   \code{nu > 2} gives lighter tails. Must be positive. Default is 2.
#' @param sd Standard deviation. Default is 1.
#'
#' @return A function with signature \code{function(n)} that generates
#'   \code{n} GED innovations with mean zero.
#'
#' @details
#' The Generalized Error Distribution (also called Generalized Normal
#' Distribution) is a flexible family that includes:
#' \itemize{
#'   \item \code{nu = 1}: Laplace (double exponential) distribution
#'   \item \code{nu = 2}: Normal distribution
#'   \item \code{nu -> Inf}: Uniform distribution
#' }
#'
#' The GED is commonly used in GARCH modeling to capture leptokurtosis
#' (fat tails) in financial returns.
#'
#' @seealso \code{\link[fGarch]{rged}} for the underlying random generator.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Standard normal (nu = 2)
#' ged_norm <- make_gen_ged(nu = 2)
#'
#' # Heavy-tailed (nu = 1.5)
#' ged_heavy <- make_gen_ged(nu = 1.5)
#' innovations <- ged_heavy(1000)
#' hist(innovations, breaks = 50)
#'
#' # Laplace distribution (nu = 1)
#' ged_laplace <- make_gen_ged(nu = 1)
#'
#' # Light-tailed (nu = 3)
#' ged_light <- make_gen_ged(nu = 3)
#' }
make_gen_ged <- function(nu = 2, sd = 1) {
  if (!requireNamespace("fGarch", quietly = TRUE)) {
    stop("Package 'fGarch' is required. Install with: install.packages('fGarch')")
  }

  if (nu <= 0) {
    stop("nu must be positive")
  }

  if (sd <= 0) {
    stop("sd must be positive")
  }

  force(nu)
  force(sd)

  function(n) {
    # fGarch::rged generates standardized GED (mean 0, sd 1)
    fGarch::rged(n, mean = 0, sd = sd, nu = nu)
  }
}


# ==============================================================================
# Laplace (Double Exponential) Innovation Generator
# ==============================================================================

#' Create a Laplace (Double Exponential) Innovation Generator
#'
#' Factory function that returns a Laplace-distributed innovation generator.
#' The Laplace distribution has heavier tails than normal but lighter than
#' Student's t with low df.
#'
#' @param scale Scale parameter (related to standard deviation).
#'   The variance is \code{2 * scale^2}. Default is \code{1 / sqrt(2)}
#'   which gives unit variance.
#'
#' @return A function with signature \code{function(n)} that generates
#'   \code{n} Laplace innovations with mean zero.
#'
#' @details
#' The Laplace distribution has a sharper peak and heavier tails than
#' the normal distribution. It arises naturally as the difference of
#' two independent exponential random variables.
#'
#' Key properties:
#' \itemize{
#'   \item Mean: 0
#'   \item Variance: \code{2 * scale^2}
#'   \item Excess kurtosis: 3 (vs 0 for normal)
#' }
#'
#' This implementation uses base R only (no external packages).
#'
#' @export
#'
#' @examples
#' # Unit variance Laplace
#' lap_gen <- make_gen_laplace()
#' innovations <- lap_gen(1000)
#' var(innovations)
#'
#' # Higher variance
#' lap_gen_hv <- make_gen_laplace(scale = 1)
#' var(lap_gen_hv(1000))
make_gen_laplace <- function(scale = 1 / sqrt(2)) {
  if (scale <= 0) {
    stop("scale must be positive")
  }

  force(scale)

  function(n) {
    # Laplace via inverse CDF
    # Avoid exact 0.5 which gives log(0)
    u <- runif(n, min = 1e-10, max = 1 - 1e-10) - 0.5
    scale * sign(u) * log(1 - 2 * abs(u))
  }
}


# ==============================================================================
# Uniform Innovation Generator
# ==============================================================================

#' Create a Uniform Innovation Generator
#'
#' Factory function that returns a uniform innovation generator centered
#' at zero. Useful for robustness testing with bounded innovations.
#'
#' @param half_width Half-width of the uniform distribution. Innovations
#'   are generated on \code{(-half_width, half_width)}. Default is
#'   \code{sqrt(3)} which gives unit variance.
#'
#' @return A function with signature \code{function(n)} that generates
#'   \code{n} uniform innovations with mean zero.
#'
#' @details
#' The uniform distribution on \code{(-a, a)} has:
#' \itemize{
#'   \item Mean: 0
#'   \item Variance: \code{a^2 / 3}
#' }
#'
#' Setting \code{half_width = sqrt(3)} gives unit variance, matching
#' the standard normal for easier comparison.
#'
#' Uniform innovations are bounded (no extreme outliers), making them
#' useful for testing how estimators behave without heavy tails.
#'
#' @export
#'
#' @examples
#' # Unit variance uniform
#' unif_gen <- make_gen_unif()
#' innovations <- unif_gen(1000)
#' var(innovations)
#'
#' # Custom range (-2, 2)
#' unif_gen_2 <- make_gen_unif(half_width = 2)
#' innovations_2 <- unif_gen_2(1000)
#' range(innovations_2)
make_gen_unif <- function(half_width = sqrt(3)) {
  if (half_width <= 0) {
    stop("half_width must be positive")
  }

  force(half_width)

  function(n) {
    runif(n, min = -half_width, max = half_width)
  }
}


# ==============================================================================
# Mixture of Normals Innovation Generator
# ==============================================================================

#' Create a Mixture of Normals Innovation Generator
#'
#' Factory function that returns an innovation generator from a mixture
#' of two normal distributions. Useful for modeling occasional outliers
#' or regime-switching behavior.
#'
#' @param sd1 Standard deviation of the first (primary) component.
#'   Default is 1.
#' @param sd2 Standard deviation of the second (outlier) component.
#'   Default is 3.
#' @param prob1 Probability of drawing from the first component.
#'   Default is 0.9 (90% from primary, 10% outliers).
#'
#' @return A function with signature \code{function(n)} that generates
#'   \code{n} mixture-normal innovations with mean zero.
#'
#' @details
#' A mixture of normals can capture:
#' \itemize{
#'   \item Occasional large shocks (outliers) when \code{sd2 >> sd1}
#'   \item Regime-switching volatility
#'   \item Heavy tails through the mixing
#' }
#'
#' Both components are centered at zero, so the mixture has mean zero.
#' The variance is \code{prob1 * sd1^2 + (1 - prob1) * sd2^2}.
#'
#' @export
#'
#' @examples
#' # Default: 90% N(0,1), 10% N(0,9)
#' mix_gen <- make_gen_mixnorm()
#' innovations <- mix_gen(1000)
#' hist(innovations, breaks = 50)
#'
#' # More frequent outliers
#' mix_gen_30 <- make_gen_mixnorm(prob1 = 0.7, sd2 = 5)
#'
#' # Variance of default mixture
#' # 0.9 * 1 + 0.1 * 9 = 1.8
make_gen_mixnorm <- function(sd1 = 1, sd2 = 3, prob1 = 0.9) {
  if (sd1 <= 0 || sd2 <= 0) {
    stop("sd1 and sd2 must be positive")
  }

  if (prob1 <= 0 || prob1 >= 1) {
    stop("prob1 must be between 0 and 1 (exclusive)")
  }

  force(sd1)
  force(sd2)
  force(prob1)

  function(n) {
    # Draw component indicators
    component <- rbinom(n, size = 1, prob = prob1)

    # Generate only needed values from each component
    x <- numeric(n)
    idx1 <- component == 1
    n1 <- sum(idx1)
    if (n1 > 0) x[idx1] <- rnorm(n1, mean = 0, sd = sd1)
    if (n1 < n) x[!idx1] <- rnorm(n - n1, mean = 0, sd = sd2)
    x
  }
}


# ==============================================================================
# GARCH Process Generator
# ==============================================================================

#' Create a GARCH Process Generator
#'
#' Factory function that returns a GARCH process generator using the
#' \pkg{rugarch} package. Supports various GARCH variants and innovation
#' distributions.
#'
#' @param omega The constant term in the variance equation. Must be positive.
#' @param alpha Numeric vector of ARCH coefficients. Length determines the
#'   ARCH order (q).
#' @param beta Numeric vector of GARCH coefficients (default \code{NULL} for
#'   pure ARCH). Length determines the GARCH order (p).
#' @param model Character string specifying the GARCH variant. One of:
#'   \describe{
#'     \item{\code{"sGARCH"}}{Standard GARCH (default)}
#'     \item{\code{"eGARCH"}}{Exponential GARCH (Nelson, 1991)}
#'     \item{\code{"gjrGARCH"}}{GJR-GARCH with leverage (Glosten et al., 1993)}
#'     \item{\code{"apARCH"}}{Asymmetric Power ARCH (Ding et al., 1993)}
#'     \item{\code{"iGARCH"}}{Integrated GARCH}
#'     \item{\code{"csGARCH"}}{Component sGARCH}
#'   }
#' @param distribution Character string specifying the innovation distribution.
#'   One of:
#'   \describe{
#'     \item{\code{"norm"}}{Normal distribution (default)}
#'     \item{\code{"std"}}{Student's t distribution}
#'     \item{\code{"ged"}}{Generalized Error Distribution}
#'     \item{\code{"snorm"}}{Skew normal}
#'     \item{\code{"sstd"}}{Skew Student's t}
#'     \item{\code{"sged"}}{Skew GED}
#'     \item{\code{"nig"}}{Normal Inverse Gaussian}
#'     \item{\code{"jsu"}}{Johnson's SU}
#'   }
#' @param distribution_params Named list of additional distribution parameters
#'   (e.g., \code{list(shape = 5)} for Student's t degrees of freedom).
#'
#' @return A function with signature \code{function(n)} that generates
#'   \code{n} observations from the specified GARCH process.
#'
#' @details
#' The standard GARCH(q, p) model specifies the conditional variance as:
#' \deqn{\sigma_t^2 = \omega + \sum_{i=1}^{q} \alpha_i \epsilon_{t-i}^2 +
#'   \sum_{j=1}^{p} \beta_j \sigma_{t-j}^2}
#'
#' For stationarity of standard GARCH, we require
#' \code{sum(alpha) + sum(beta) < 1}.
#'
#' \strong{Important}: This generator returns \eqn{\epsilon_t = \sigma_t z_t},
#' which are \emph{not} iid---they exhibit time-varying conditional variance
#' (volatility clustering). When passed to \code{\link{gen_aruma_flex}}, the
#' result is a proper ARMA-GARCH process where the ARMA filter is applied
#' to the heteroskedastic innovations.
#'
#' @seealso \code{\link[rugarch]{ugarchspec}}, \code{\link[rugarch]{ugarchpath}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Standard GARCH(1,1) with normal innovations
#' garch11_gen <- make_gen_garch(
#'   omega = 0.1,
#'   alpha = 0.15,
#'   beta = 0.8
#' )
#' y <- garch11_gen(500)
#' plot(y, type = "l", main = "GARCH(1,1) Simulation")
#'
#' # ARCH(2) process (no GARCH terms)
#' arch2_gen <- make_gen_garch(
#'   omega = 0.2,
#'   alpha = c(0.3, 0.2)
#' )
#'
#' # GARCH(1,1) with Student's t innovations
#' garch_t_gen <- make_gen_garch(
#'   omega = 0.1,
#'   alpha = 0.1,
#'   beta = 0.85,
#'   distribution = "std",
#'   distribution_params = list(shape = 5)
#' )
#'
#' # GJR-GARCH with leverage effect
#' gjr_gen <- make_gen_garch(
#'   omega = 0.1,
#'   alpha = 0.05,
#'   beta = 0.9,
#'   model = "gjrGARCH"
#' )
#' }
make_gen_garch <- function(omega,
                           alpha,
                           beta = NULL,
                           model = c("sGARCH", "eGARCH", "gjrGARCH",
                                     "apARCH", "iGARCH", "csGARCH"),
                           distribution = c("norm", "std", "ged", "snorm",
                                            "sstd", "sged", "nig", "jsu"),
                           distribution_params = list()) {

  if (!requireNamespace("rugarch", quietly = TRUE)) {
    stop("Package 'rugarch' is required. Install with: install.packages('rugarch')")
  }

  model <- match.arg(model)
  distribution <- match.arg(distribution)

  if (omega <= 0) {
    stop("omega must be positive")
  }
  if (any(alpha < 0)) {
    stop("alpha coefficients must be non-negative")
  }
  if (!is.null(beta) && any(beta < 0)) {
    stop("beta coefficients must be non-negative")
  }

  q <- length(alpha)
  p <- if (is.null(beta)) 0L else length(beta)

  if (model == "sGARCH" && (sum(alpha) + sum(beta)) >= 1) {
    warning("sum(alpha) + sum(beta) >= 1: process may be non-stationary")
  }

  names(alpha) <- paste0("alpha", seq_len(q))
  if (p > 0) {
    names(beta) <- paste0("beta", seq_len(p))
  }

  fixed <- c(omega = omega, alpha)
  if (p > 0) {
    fixed <- c(fixed, beta)
  }

  if (length(distribution_params) > 0) {
    fixed <- c(fixed, unlist(distribution_params))
  }

  force(model)
  force(distribution)
  force(q)
  force(p)
  force(fixed)

  function(n) {
    spec <- rugarch::ugarchspec(
      variance.model = list(model = model, garchOrder = c(q, p)),
      mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
      distribution.model = distribution,
      fixed.pars = as.list(fixed)
    )

    path <- rugarch::ugarchpath(spec, n.sim = n)
    as.numeric(path@path$seriesSim)
  }
}


# ==============================================================================
# Heteroscedastic Normal Innovation Generator
# ==============================================================================

#' Create a Heteroscedastic Normal Innovation Generator
#'
#' Factory function that returns an innovation generator with time-varying
#' variance. Innovations are scaled by weights w, allowing for patterns like
#' increasing variance, periodic variance, or custom heteroscedasticity.
#'
#' @param w Weight specification. One of:
#'   \itemize{
#'     \item \code{NULL} (default): Use linear weights \code{1, 2, ..., n}
#'       (variance increases linearly with time)
#'     \item A function: Called as \code{w(n)} to generate n weights
#'     \item A numeric vector: Used directly; must have length >= n when
#'       generator is called
#'   }
#' @param sd Base standard deviation before weighting. Default is 1.
#'
#' @return A function with signature \code{function(n)} that generates
#'   \code{n} heteroscedastic normal innovations.
#'
#' @details
#' The generator produces: \code{w[1:n] * rnorm(n, 0, sd)}
#'
#' Common weight patterns:
#' \itemize{
#'   \item Linear: \code{w = NULL} or \code{w = function(n) seq_len(n)}
#'   \item Sqrt: \code{w = function(n) sqrt(seq_len(n))}
#'   \item Exponential: \code{w = function(n) exp(seq_len(n) / n)}
#'   \item Periodic: \code{w = function(n) 1 + 0.5 * sin(2 * pi * seq_len(n) / 12)}
#' }
#'
#' @export
#'
#' @examples
#' # Default linear increasing variance
#' hetero_gen <- make_gen_hetero()
#' innovations <- hetero_gen(100)
#' plot(innovations, type = "l")
#'
#' # Custom weight function (sqrt growth)
#' hetero_sqrt <- make_gen_hetero(w = function(n) sqrt(seq_len(n)))
#'
#' # Fixed weight vector
#' weights <- c(rep(1, 50), rep(3, 50))  # Jump in variance
#' hetero_jump <- make_gen_hetero(w = weights)
make_gen_hetero <- function(w = NULL, sd = 1) {
  if (sd <= 0) {
    stop("sd must be positive")
  }

  # Validate w if provided as vector
  if (!is.null(w) && !is.function(w) && !is.numeric(w)) {
    stop("w must be NULL, a function, or a numeric vector")
  }

  force(w)
  force(sd)

  function(n) {
    # Determine weights
    if (is.null(w)) {
      weights <- seq_len(n)
    } else if (is.function(w)) {
      weights <- w(n)
      if (length(weights) != n) {
        stop("Weight function must return exactly n values")
      }
    } else {
      # w is a numeric vector
      if (length(w) < n) {
        stop("Weight vector must have length >= n (", n, ")")
      }
      weights <- w[seq_len(n)]
    }

    # Generate scaled innovations
    weights * rnorm(n, mean = 0, sd = sd)
  }
}
