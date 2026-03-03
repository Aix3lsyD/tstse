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
# Heteroscedastic Innovation Generator (internal helper + factory)
# ==============================================================================

# Build a weight-generating function from a named shape specification.
# Returns function(n) that produces the weight vector.
# @keywords internal
.hetero_shape_weights <- function(shape, from, to, power, breaks, levels,
                                  base_w, amplitude, period) {
  shape <- match.arg(shape, c("linear", "sqrt", "log", "power", "exp",
                               "step", "periodic"))

  # Validate from/to for shapes that use them
  if (!shape %in% c("step", "periodic")) {
    if (!is.numeric(from) || length(from) != 1L || !is.finite(from) || from <= 0) {
      stop("'from' must be a single positive finite number")
    }
    if (!is.numeric(to) || length(to) != 1L || !is.finite(to) || to <= 0) {
      stop("'to' must be a single positive finite number")
    }
  }

  switch(shape,
    "linear" = {
      force(from); force(to)
      function(n) {
        t <- seq(0, 1, length.out = n)
        from + (to - from) * t
      }
    },
    "sqrt" = {
      force(from); force(to)
      function(n) {
        t <- seq(0, 1, length.out = n)
        from + (to - from) * sqrt(t)
      }
    },
    "log" = {
      force(from); force(to)
      function(n) {
        t <- seq(0, 1, length.out = n)
        from + (to - from) * log1p(t * (exp(1) - 1))
      }
    },
    "power" = {
      if (!is.numeric(power) || length(power) != 1L || !is.finite(power) || power <= 0) {
        stop("'power' must be a single positive finite number")
      }
      force(from); force(to); force(power)
      function(n) {
        t <- seq(0, 1, length.out = n)
        from + (to - from) * t^power
      }
    },
    "exp" = {
      force(from); force(to)
      function(n) {
        t <- seq(0, 1, length.out = n)
        from * (to / from)^t
      }
    },
    "step" = {
      if (is.null(breaks) || !is.numeric(breaks) || length(breaks) < 1L) {
        stop("'breaks' must be a numeric vector of at least length 1 for 'step' shape")
      }
      if (any(breaks <= 0) || any(breaks >= 1)) {
        stop("'breaks' values must be in (0, 1)")
      }
      if (is.unsorted(breaks)) {
        stop("'breaks' must be sorted in increasing order")
      }
      if (is.null(levels) || !is.numeric(levels)) {
        stop("'levels' must be a numeric vector for 'step' shape")
      }
      if (length(levels) != length(breaks) + 1L) {
        stop("'levels' must have length(breaks) + 1 elements")
      }
      if (any(levels <= 0)) {
        stop("All 'levels' must be positive (they represent SD weights)")
      }
      force(breaks); force(levels)
      function(n) {
        t <- seq(0, 1, length.out = n)
        cut_idx <- findInterval(t, c(0, breaks, 1), rightmost.closed = TRUE)
        cut_idx <- pmin(pmax(cut_idx, 1L), length(levels))
        levels[cut_idx]
      }
    },
    "periodic" = {
      if (!is.numeric(base_w) || length(base_w) != 1L || !is.finite(base_w) || base_w <= 0) {
        stop("'base_w' must be a single positive finite number")
      }
      if (!is.numeric(amplitude) || length(amplitude) != 1L || !is.finite(amplitude)) {
        stop("'amplitude' must be a single finite number")
      }
      if (!is.numeric(period) || length(period) != 1L || !is.finite(period) || period <= 0) {
        stop("'period' must be a single positive finite number")
      }
      if (base_w - abs(amplitude) <= 0) {
        warning("base_w - |amplitude| <= 0: weights may go non-positive, which flips sign of innovations")
      }
      force(base_w); force(amplitude); force(period)
      function(n) {
        base_w + amplitude * sin(2 * pi * seq_len(n) / period)
      }
    }
  )
}


#' Create a Heteroscedastic Innovation Generator
#'
#' Factory function that returns an innovation generator with time-varying
#' variance. Supports named shape patterns with range control, custom weight
#' functions/vectors, and composition with any base distribution.
#'
#' @param w Weight specification for legacy interface. One of:
#'   \itemize{
#'     \item \code{NULL} (default): Use linear weights from 1 to 10
#'       (when \code{shape} is also NULL)
#'     \item A function: Called as \code{w(n)} to generate n weights
#'     \item A numeric vector: Used directly; left-padded with \code{w[1]}
#'       if shorter than n
#'   }
#'   Mutually exclusive with \code{shape}.
#' @param sd Base standard deviation when using normal base distribution.
#'   Default is 1. Mutually exclusive with \code{base} (when non-default).
#' @param shape Named weight pattern. One of \code{"linear"}, \code{"sqrt"},
#'   \code{"log"}, \code{"power"}, \code{"exp"}, \code{"step"},
#'   \code{"periodic"}. Mutually exclusive with \code{w}.
#' @param from Starting weight (SD multiplier) for shape patterns. Default 1.
#' @param to Ending weight (SD multiplier) for shape patterns. Default 10.
#' @param power Exponent for \code{shape = "power"}. Default 2.
#' @param breaks Numeric vector of breakpoints in (0, 1) for
#'   \code{shape = "step"}. Proportions of the series length.
#' @param levels Numeric vector of positive SD weights for
#'   \code{shape = "step"}. Must have \code{length(breaks) + 1} elements.
#' @param base_w Base weight for \code{shape = "periodic"}. Default 1.
#' @param amplitude Amplitude for \code{shape = "periodic"}. Default 0.5.
#' @param period Period (in observations) for \code{shape = "periodic"}.
#'   Default 12.
#' @param base A base innovation generator function (e.g.,
#'   \code{make_gen_t(df = 3)}). When NULL (default), uses
#'   \code{make_gen_norm(sd = sd)}.
#'
#' @return A function with signature \code{function(n)} that generates
#'   \code{n} heteroscedastic innovations.
#'
#' @details
#' The generator produces: \code{weights * base_gen(n)} where weights
#' are determined by either \code{shape} or \code{w}, and \code{base_gen}
#' is determined by \code{base} or \code{sd}.
#'
#' \strong{Shape formulas} (using normalized time \code{t} from 0 to 1):
#' \itemize{
#'   \item \code{"linear"}: \code{from + (to - from) * t}
#'   \item \code{"sqrt"}: \code{from + (to - from) * sqrt(t)}
#'   \item \code{"log"}: \code{from + (to - from) * log1p(t * (e - 1))}
#'   \item \code{"power"}: \code{from + (to - from) * t^power}
#'   \item \code{"exp"}: \code{from * (to / from)^t}
#'   \item \code{"step"}: Piecewise constant at breakpoints
#'   \item \code{"periodic"}: \code{base_w + amplitude * sin(2*pi*i/period)}
#' }
#'
#' Note: weights scale the \strong{standard deviation}, not the variance.
#' A weight of 5 means SD is 5x, so variance is 25x.
#'
#' @export
#'
#' @seealso \code{\link{make_gen_norm}}, \code{\link{make_gen_t}},
#'   \code{\link{make_gen_skt}} for base distribution generators.
#'
#' @examples
#' # Named shape with range
#' gen <- make_gen_hetero(shape = "sqrt", from = 1, to = 5)
#' plot(gen(200), type = "l")
#'
#' # Compose with t-distributed base
#' gen_t <- make_gen_hetero(shape = "linear", from = 1, to = 5,
#'                          base = make_gen_t(df = 3))
#'
#' # Step function (regime shift at midpoint)
#' gen_step <- make_gen_hetero(shape = "step",
#'                             breaks = c(0.5), levels = c(1, 5))
#'
#' # Legacy interface still works
#' hetero_gen <- make_gen_hetero()
#' hetero_sqrt <- make_gen_hetero(w = function(n) sqrt(seq_len(n)))
#' hetero_jump <- make_gen_hetero(w = c(rep(1, 50), rep(3, 50)))
make_gen_hetero <- function(w = NULL, sd = 1, shape = NULL,
                            from = 1, to = 10, power = 2,
                            breaks = NULL, levels = NULL,
                            base_w = 1, amplitude = 0.5, period = 12,
                            base = NULL) {

  # --- Mutual exclusivity: shape vs w ---
  if (!is.null(shape) && !is.null(w)) {
    stop("Cannot specify both 'shape' and 'w'. Use 'shape' for named patterns or 'w' for custom weights.")
  }

  # --- Mutual exclusivity: sd vs base ---
  if (!is.null(base) && !identical(sd, 1)) {
    stop("Cannot specify both 'sd' (non-default) and 'base'. Use 'base' for custom distributions, or 'sd' with default normal.")
  }

  # --- Validate sd ---
  if (is.numeric(sd) && sd <= 0) {
    stop("sd must be positive")
  }

  # --- Validate base ---
  if (!is.null(base) && !is.function(base)) {
    stop("'base' must be a function (e.g., make_gen_t(df = 3))")
  }

  # --- Validate w (legacy path) ---
  if (!is.null(w) && !is.function(w) && !is.numeric(w)) {
    stop("w must be NULL, a function, or a numeric vector")
  }

  # --- Resolve base distribution ---
  if (!is.null(base)) {
    base_gen <- base
  } else {
    base_gen <- make_gen_norm(sd = sd)
  }

  # --- Resolve weight strategy ---
  if (!is.null(shape)) {
    w_fun <- .hetero_shape_weights(shape, from, to, power, breaks, levels,
                                   base_w, amplitude, period)
  } else {
    w_fun <- NULL
  }

  force(w)
  force(base_gen)
  force(w_fun)

  function(n) {
    # Determine weights
    if (!is.null(w_fun)) {
      # Shape-based path
      weights <- w_fun(n)
    } else if (is.null(w)) {
      # Default legacy: linear weights from 1 to 10
      weights <- seq(1, 10, length.out = n)
    } else if (is.function(w)) {
      weights <- w(n)
      if (length(weights) != n) {
        stop("Weight function must return exactly n values")
      }
    } else {
      # w is a numeric vector; left-pad with first weight for burn-in positions
      if (length(w) < n) {
        weights <- c(rep(w[1], n - length(w)), w)
      } else {
        weights <- w[seq_len(n)]
      }
    }

    # Generate scaled innovations
    weights * base_gen(n)
  }
}
