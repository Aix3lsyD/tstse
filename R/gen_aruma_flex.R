#' Generate Flexible ARUMA Process
#'
#' Generate a realization from an ARUMA (ARIMA with general nonstationary factor)
#' process with flexible innovation distributions.
#'
#' @importFrom utils head
#' @importFrom stats sd
#'
#' @param n Integer, length of output series.
#' @param phi Numeric vector, AR coefficients (default 0 = none).
#' @param theta Numeric vector, MA coefficients (default 0 = none).
#' @param d Integer, order of differencing (default 0).
#' @param s Integer, seasonal period for (1-B^s) factor (default 0 = none).
#' @param lambda Numeric vector, coefficients for additional nonstationary factor
#'   (default 0 = none). For factor (1 - lambda_1 B - lambda_2 B^2 - ...),
#'   pass c(lambda_1, lambda_2, ...).
#' @param innov_gen Function, innovation generator with signature \code{function(n)}
#'   returning a numeric vector of length n with mean zero. If \code{NULL} (default),
#'   uses \code{\link{make_gen_norm}(sd = sqrt(vara))}.
#' @param vara Numeric, white noise variance (default 1). Ignored if \code{innov_gen}
#'   is provided.
#' @param plot Logical or character, controls plotting:
#'   \code{TRUE} = auto-detect (ggplot2 if available, else base R),
#'   \code{FALSE} = no plot, \code{"ggplot2"} = force ggplot2,
#'   \code{"base"} = force base R.
#' @param seed Integer, random seed for reproducibility (default NULL).
#'
#' @return An object of class \code{"aruma"} containing:
#'   \describe{
#'     \item{y}{Numeric vector of length n, the generated realization}
#'     \item{n}{Length of series}
#'     \item{p}{AR order}
#'     \item{q}{MA order}
#'     \item{d}{Differencing order}
#'     \item{s}{Seasonal period}
#'     \item{phi}{AR coefficients}
#'     \item{theta}{MA coefficients}
#'     \item{lambda}{Nonstationary factor coefficients}
#'     \item{innov_type}{Description of innovation distribution}
#'     \item{plot_obj}{ggplot object (if ggplot2 used) or NULL}
#'   }
#'
#' @export
#'
#' @details
#' Generates from the ARUMA model which extends ARIMA with a general
#' nonstationary operator. The model is:
#'
#' \deqn{\lambda(B)(1-B)^d(1-B^s)\phi(B) X_t = \theta(B) a_t}
#'
#' where:
#' \itemize{
#'   \item \eqn{\phi(B)} is the AR operator
#'   \item \eqn{\theta(B)} is the MA operator
#'   \item \eqn{(1-B)^d} is the differencing operator
#'   \item \eqn{(1-B^s)} is the seasonal differencing operator
#'   \item \eqn{\lambda(B)} is a general nonstationary factor
#'   \item \eqn{a_t} are innovations from the specified distribution
#' }
#'
#' The function first generates a stationary ARMA process using the specified
#' innovation distribution, then applies the inverse of all nonstationary factors.
#'
#' \strong{Burn-in}: Uses adaptive burn-in that scales with AR/MA orders:
#' \code{n_start = max(100, 10*p, 10*q)} for ARMA stabilization, plus additional
#' burn-in for nonstationary operators.
#'
#' \strong{Important for GARCH}: When using \code{\link{make_gen_garch}}, a single
#' call generates all innovations at once to preserve the dependence structure.
#'
#' @seealso
#' \code{\link{gen_aruma}} for the simpler interface with normal innovations,
#' \code{\link{make_gen_norm}}, \code{\link{make_gen_t}}, \code{\link{make_gen_skt}},
#' \code{\link{make_gen_ged}}, \code{\link{make_gen_laplace}}, \code{\link{make_gen_unif}},
#' \code{\link{make_gen_mixnorm}}, \code{\link{make_gen_garch}} for innovation generators.
#'
#' @examples
#' # ARMA(1,1) with normal innovations (default)
#' result <- gen_aruma_flex(n = 200, phi = 0.7, theta = 0.3, seed = 42)
#' print(result)
#'
#' # AR(1) with t-distributed innovations
#' t_gen <- make_gen_t(df = 5, scale = TRUE)
#' result <- gen_aruma_flex(n = 200, phi = 0.8, innov_gen = t_gen, seed = 42)
#'
#' # ARIMA(1,1,0) with mixture innovations (outliers)
#' mix_gen <- make_gen_mixnorm(sd1 = 1, sd2 = 5, prob1 = 0.95)
#' result <- gen_aruma_flex(n = 200, phi = 0.5, d = 1, innov_gen = mix_gen, seed = 42)
#'
#' # Force base R plotting
#' result <- gen_aruma_flex(n = 100, phi = 0.5, plot = "base", seed = 42)
#'
#' # No plotting, return object only
#' result <- gen_aruma_flex(n = 100, phi = 0.5, plot = FALSE, seed = 42)
#' plot(result)  # Plot later using S3 method
gen_aruma_flex <- function(n,
                           phi = 0,
                           theta = 0,
                           d = 0L,
                           s = 0L,
                           lambda = 0,
                           innov_gen = NULL,
                           vara = 1,
                           plot = TRUE,
                           seed = NULL) {

  # Input validation
  if (!is.numeric(n) || length(n) != 1 || n < 1) {
    stop("`n` must be a positive integer")
  }
  if (!is.numeric(phi)) stop("`phi` must be numeric")
  if (!is.numeric(theta)) stop("`theta` must be numeric")
  if (!is.numeric(d) || d < 0) stop("`d` must be a non-negative integer")
  if (!is.numeric(s) || s < 0) stop("`s` must be a non-negative integer")
  if (!is.numeric(lambda)) stop("`lambda` must be numeric")
  if (!is.numeric(vara) || vara <= 0) stop("`vara` must be positive")

  # Validate plot parameter
  valid_plot <- c(TRUE, FALSE, "ggplot2", "base")
  if (!isTRUE(plot) && !isFALSE(plot) && !(plot %in% c("ggplot2", "base"))) {
    stop("`plot` must be TRUE, FALSE, 'ggplot2', or 'base'")
  }

  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)

  # Determine effective orders
  p <- if (all(phi == 0)) 0L else length(phi)
  q <- if (all(theta == 0)) 0L else length(theta)
  dlam <- if (all(lambda == 0)) 0L else length(lambda)

  # Set up innovation generator
  if (is.null(innov_gen)) {
    innov_gen <- make_gen_norm(sd = sqrt(vara))
    innov_type <- sprintf("Normal(0, %.4g)", vara)
  } else {
    if (!is.function(innov_gen)) {
      stop("`innov_gen` must be a function with signature function(n)")
    }
    innov_type <- "Custom"
  }

  # Build model list for arima.sim
  # Note: arima.sim uses opposite sign convention for MA
  model <- list(order = c(p, as.integer(d), q))
  if (p > 0) model$ar <- phi
  if (q > 0) model$ma <- -theta

  # Calculate total nonstationary order from lambda and seasonal
  dlams <- dlam + s

  # Adaptive burn-in (scales with AR/MA orders)
  n_start <- max(100L, 10L * p, 10L * q)
  spin <- max(100L, 2L * dlams)

  # Generate with spin-up
  n_sim <- n + dlams + spin

  # Generate all innovations at once (preserves GARCH dependence)
  total_innov <- n_sim + n_start
  all_innov <- innov_gen(total_innov)

  # Use arima.sim with custom innovations
  tsdata <- arima.sim(
    n = n_sim,
    model = model,
    innov = all_innov[(n_start + 1):total_innov],
    n.start = n_start,
    start.innov = all_innov[1:n_start]
  )
  y <- as.double(tsdata)

  # Build combined nonstationary operator from lambda and seasonal
  lambdas <- .build_nonstationary_operator(lambda, s, dlam)

  # Apply inverse of nonstationary operator
  x <- .apply_inverse_operator(y, lambdas, dlams, d, n, spin)

  # Handle plotting
  plot_obj <- NULL
  if (!isFALSE(plot)) {
    use_ggplot2 <- .resolve_plot_type(plot)
    if (use_ggplot2) {
      plot_obj <- .plot_aruma_flex_ggplot(x, n)
      print(plot_obj)
    } else {
      .plot_aruma_flex_base(x, n)
    }
  }

  # Build S3 object
  result <- structure(
    list(
      y = x,
      n = as.integer(n),
      p = p,
      q = q,
      d = as.integer(d),
      s = as.integer(s),
      phi = if (p > 0) phi else numeric(0),
      theta = if (q > 0) theta else numeric(0),
      lambda = if (dlam > 0) lambda else numeric(0),
      innov_type = innov_type,
      plot_obj = plot_obj
    ),
    class = "aruma"
  )

  invisible(result)
}


#' @describeIn gen_aruma_flex Print method for aruma objects
#' @param x An object of class \code{"aruma"}
#' @param ... Additional arguments (ignored)
#' @export
print.aruma <- function(x, ...) {
  cat("ARUMA Realization\n")
  cat(sprintf("n = %d, p = %d, q = %d, d = %d, s = %d\n",
              x$n, x$p, x$q, x$d, x$s))

  if (x$p > 0) {
    cat("phi =", paste(round(x$phi, 4), collapse = ", "), "\n")
  }
  if (x$q > 0) {
    cat("theta =", paste(round(x$theta, 4), collapse = ", "), "\n")
  }
  if (length(x$lambda) > 0) {
    cat("lambda =", paste(round(x$lambda, 4), collapse = ", "), "\n")
  }
  cat("Innovation:", x$innov_type, "\n")
  cat("\nFirst 6 values:", paste(round(head(x$y), 4), collapse = ", "), "\n")

  invisible(x)
}


#' @describeIn gen_aruma_flex Plot method for aruma objects
#' @param x An object of class \code{"aruma"}
#' @param type Character, plot type: \code{"auto"} (default), \code{"ggplot2"}, or \code{"base"}
#' @param ... Additional arguments passed to plotting functions
#' @export
plot.aruma <- function(x, type = c("auto", "ggplot2", "base"), ...) {
  type <- match.arg(type)

  use_ggplot2 <- switch(type,
    "auto" = requireNamespace("ggplot2", quietly = TRUE),
    "ggplot2" = TRUE,
    "base" = FALSE
  )

  if (use_ggplot2) {
    p <- .plot_aruma_flex_ggplot(x$y, x$n)
    print(p)
    invisible(p)
  } else {
    .plot_aruma_flex_base(x$y, x$n)
    invisible(NULL)
  }
}


#' @describeIn gen_aruma_flex Summary method for aruma objects
#' @param object An object of class \code{"aruma"}
#' @param ... Additional arguments (ignored)
#' @export
summary.aruma <- function(object, ...) {
  y <- object$y

  stats <- list(
    mean = mean(y),
    sd = sd(y),
    min = min(y),
    max = max(y),
    skewness = .compute_skewness(y),
    kurtosis = .compute_kurtosis(y)
  )

  result <- structure(
    list(
      n = object$n,
      model = list(p = object$p, q = object$q, d = object$d, s = object$s),
      phi = object$phi,
      theta = object$theta,
      lambda = object$lambda,
      innov_type = object$innov_type,
      stats = stats
    ),
    class = "summary.aruma"
  )

  result
}


#' @describeIn gen_aruma_flex Print method for summary.aruma objects
#' @param x An object of class \code{"summary.aruma"}
#' @param ... Additional arguments (ignored)
#' @export
print.summary.aruma <- function(x, ...) {
  cat("Summary of ARUMA Realization\n")
  cat(sprintf("n = %d, p = %d, q = %d, d = %d, s = %d\n",
              x$n, x$model$p, x$model$q, x$model$d, x$model$s))

  if (x$model$p > 0) {
    cat("phi =", paste(round(x$phi, 4), collapse = ", "), "\n")
  }
  if (x$model$q > 0) {
    cat("theta =", paste(round(x$theta, 4), collapse = ", "), "\n")
  }
  if (length(x$lambda) > 0) {
    cat("lambda =", paste(round(x$lambda, 4), collapse = ", "), "\n")
  }
  cat("Innovation:", x$innov_type, "\n")

  cat("\nSample Statistics:\n")
  cat(sprintf("  Mean:     %10.4f\n", x$stats$mean))
  cat(sprintf("  Std Dev:  %10.4f\n", x$stats$sd))
  cat(sprintf("  Min:      %10.4f\n", x$stats$min))
  cat(sprintf("  Max:      %10.4f\n", x$stats$max))
  cat(sprintf("  Skewness: %10.4f\n", x$stats$skewness))
  cat(sprintf("  Kurtosis: %10.4f\n", x$stats$kurtosis))

  invisible(x)
}


# ==============================================================================
# Internal helper functions
# ==============================================================================

#' Resolve plot type to boolean
#' @noRd
.resolve_plot_type <- function(plot) {
  if (isTRUE(plot)) {
    # Auto-detect
    requireNamespace("ggplot2", quietly = TRUE)
  } else if (plot == "ggplot2") {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Package 'ggplot2' is required for ggplot2 plotting. ",
           "Install with: install.packages('ggplot2')")
    }
    TRUE
  } else {
    # "base"
    FALSE
  }
}


#' Plot ARUMA realization using ggplot2
#' @noRd
#' @importFrom rlang .data
.plot_aruma_flex_ggplot <- function(y, n) {
  df <- data.frame(time = seq_len(n), y = y)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$time, y = .data$y)) +
    ggplot2::geom_line(color = "steelblue", linewidth = 0.5) +
    ggplot2::labs(x = "Time", y = "Value", title = "ARUMA Realization") +
    ggplot2::theme_minimal()
}


#' Plot ARUMA realization using base R
#' @noRd
.plot_aruma_flex_base <- function(y, n) {
  op <- par(mfrow = c(1, 1), mar = c(3.8, 2.5, 1.5, 1))
  on.exit(par(op), add = TRUE)

  t <- seq_len(n)
  type <- if (n < 200) "o" else "l"
  pch <- if (n < 200) 16 else NA

  plot(t, y, type = type, pch = pch, cex = 0.5,
       xlab = "Time", ylab = "", main = "ARUMA Realization")

  invisible(NULL)
}


#' Compute sample skewness
#' @noRd
.compute_skewness <- function(x) {
  n <- length(x)
  m <- mean(x)
  s <- sd(x)
  if (s == 0) return(NA_real_)
  sum((x - m)^3) / (n * s^3)
}


#' Compute sample excess kurtosis
#' @noRd
.compute_kurtosis <- function(x) {
  n <- length(x)
  m <- mean(x)
  s <- sd(x)
  if (s == 0) return(NA_real_)
  sum((x - m)^4) / (n * s^4) - 3
}
