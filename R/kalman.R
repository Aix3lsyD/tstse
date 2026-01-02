#' Kalman Filter and Smoother
#'
#' Performs Kalman filtering and smoothing for a linear Gaussian state-space
#' model. This is a convenience wrapper around [astsa::Ksmooth()].
#'
#' @param y Numeric vector. The observed time series data.
#' @param x0 Numeric scalar or vector. Initial state estimate (prior mean of
#'   the state at time 0).
#' @param P0 Numeric scalar or matrix. Initial state covariance (prior
#'   variance/covariance of the state at time 0).
#' @param Phi Numeric scalar or matrix. State transition matrix (F in some
#'   notations). Describes how the state evolves: \eqn{x_t = \Phi x_{t-1} + v_t}.
#' @param Q Numeric scalar or matrix. State noise variance (process noise
#'   covariance). The variance of \eqn{v_t}.
#' @param A Numeric scalar or matrix. Observation matrix (G or H in some
#'   notations). Relates state to observations: \eqn{y_t = A x_t + w_t}.
#' @param R Numeric scalar or matrix. Observation noise variance (measurement
#'   noise covariance). The variance of \eqn{w_t}.
#'
#' @return A list with class `"kalman"` containing:
#'   \item{y}{The original observations.}
#'   \item{predicted}{List with `state` (predicted state estimates) and
#'     `var` (predicted state variances).}
#'   \item{filtered}{List with `state` (filtered state estimates) and
#'     `var` (filtered state variances).}
#'   \item{smoothed}{List with `state` (smoothed state estimates) and
#'     `var` (smoothed state variances).}
#'   \item{result}{The raw output from [astsa::Ksmooth()] for advanced use.}
#'
#' @details
#' The function implements the Kalman filter and smoother for the linear
#' Gaussian state-space model:
#'
#' **State equation:**
#' \deqn{x_t = \Phi x_{t-1} + v_t, \quad v_t \sim N(0, Q)}
#'
#' **Observation equation:**
#' \deqn{y_t = A x_t + w_t, \quad w_t \sim N(0, R)}
#'
#' The Kalman filter provides:
#' \itemize{
#'   \item **Prediction**: Estimate of \eqn{x_t} given \eqn{y_1, \ldots, y_{t-1}}
#'   \item **Filtering**: Estimate of \eqn{x_t} given \eqn{y_1, \ldots, y_t}
#'   \item **Smoothing**: Estimate of \eqn{x_t} given all data \eqn{y_1, \ldots, y_n}
#' }
#'
#' @note
#' This function requires the `astsa` package. Install it with
#' `install.packages("astsa")` if needed.
#'
#' @references
#' Shumway, R. H., & Stoffer, D. S. (2017). *Time Series Analysis and Its
#' Applications: With R Examples* (4th ed.). Springer.
#'
#' @seealso [astsa::Ksmooth()] for the underlying implementation,
#'   [kalman_miss()] for handling missing observations.
#'
#' @examples
#' \dontrun{
#' # Local level model (random walk + noise)
#' set.seed(123)
#' n <- 100
#'
#' # True state is a random walk
#' x_true <- cumsum(rnorm(n, sd = 1))
#' # Observations are state + noise
#' y <- x_true + rnorm(n, sd = 2)
#'
#' # Fit Kalman filter
#' result <- kalman(
#'   y = y,
#'   x0 = 0,           # Initial state guess
#'   P0 = 10,          # Initial state variance (uncertain)
#'   Phi = 1,          # Random walk: x_t = x_{t-1} + v_t
#'   Q = 1,            # State noise variance
#'   A = 1,            # Direct observation: y_t = x_t + w_t
#'   R = 4             # Observation noise variance
#' )
#'
#' # Plot results
#' plot(y, type = "p", pch = 16, cex = 0.5, col = "gray")
#' lines(x_true, col = "black", lwd = 2)
#' lines(result$filtered$state, col = "blue", lwd = 2)
#' lines(result$smoothed$state, col = "red", lwd = 2)
#' legend("topleft", c("Observed", "True", "Filtered", "Smoothed"),
#'        col = c("gray", "black", "blue", "red"),
#'        pch = c(16, NA, NA, NA), lty = c(NA, 1, 1, 1), lwd = 2)
#' }
#'
#' @export
kalman <- function(y, x0, P0, Phi, Q, A, R) {
  # Check for astsa package
  if (!requireNamespace("astsa", quietly = TRUE)) {
    stop(
      "The 'astsa' package is required for kalman().\n",
      "Install it with: install.packages(\"astsa\")",
      call. = FALSE
    )
  }

  # Input validation
  if (!is.numeric(y)) {
    stop("`y` must be a numeric vector.", call. = FALSE)
  }
  n <- length(y)
  if (n < 2L) {
    stop("`y` must have at least 2 observations.", call. = FALSE)
  }

  if (!is.numeric(x0)) {
    stop("`x0` must be numeric.", call. = FALSE)
  }
  if (!is.numeric(P0)) {
    stop("`P0` must be numeric.", call. = FALSE)
  }
  if (!is.numeric(Phi)) {
    stop("`Phi` must be numeric.", call. = FALSE)
  }
  if (!is.numeric(Q)) {
    stop("`Q` must be numeric.", call. = FALSE)
  }
  if (!is.numeric(A)) {
    stop("`A` must be numeric.", call. = FALSE)
  }
  if (!is.numeric(R)) {
    stop("`R` must be numeric.", call. = FALSE)
  }

  # Validate variance parameters are non-negative
  if (any(Q < 0)) {
    stop("`Q` (state noise variance) must be non-negative.", call. = FALSE)
  }
  if (any(R < 0)) {
    stop("`R` (observation noise variance) must be non-negative.", call. = FALSE)
  }
  if (any(P0 < 0)) {
    stop("`P0` (initial state variance) must be non-negative.", call. = FALSE)
  }

  # astsa::Ksmooth takes sqrt of Q and R
  sQ <- sqrt(Q)
  sR <- sqrt(R)

  # Call astsa::Ksmooth
  # Parameter mapping:
  #   Our Phi -> astsa's Phi (state transition)
  #   Our A -> astsa's A (observation matrix)
  #   Our x0 -> astsa's mu0 (initial state mean)
  #   Our P0 -> astsa's Sigma0 (initial state covariance)
  #   Our sqrt(Q) -> astsa's sQ (sqrt of state noise cov)
  #   Our sqrt(R) -> astsa's sR (sqrt of observation noise cov)
  kf_result <- astsa::Ksmooth(
    y = y,
    A = A,
    mu0 = x0,
    Sigma0 = P0,
    Phi = Phi,
    sQ = sQ,
    sR = sR
  )

  # Extract and structure results
  # Note: astsa returns arrays; we simplify for univariate case
  result <- list(
    y = y,
    predicted = list(
      state = as.numeric(kf_result$Xp),
      var = as.numeric(kf_result$Pp)
    ),
    filtered = list(
      state = as.numeric(kf_result$Xf),
      var = as.numeric(kf_result$Pf)
    ),
    smoothed = list(
      state = as.numeric(kf_result$Xs),
      var = as.numeric(kf_result$Ps)
    ),
    result = kf_result
  )

  class(result) <- "kalman"
  result
}


#' @export
print.kalman <- function(x, ...) {
  cat("Kalman Filter/Smoother Results\n")
  cat("------------------------------\n")
  cat("Observations:", length(x$y), "\n")
  cat("\nComponents available:\n")
  cat("  $predicted  - One-step-ahead predictions\n")
  cat("  $filtered   - Filtered estimates (using data up to time t)\n")
  cat("  $smoothed   - Smoothed estimates (using all data)\n")
  cat("\nEach component contains $state and $var\n")
  invisible(x)
}
