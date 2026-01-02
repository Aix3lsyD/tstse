#' Kalman Filter and Smoother with Missing Observations
#'
#' Performs Kalman filtering and smoothing for a linear Gaussian state-space
#' model that can handle missing observations via a time-varying observation
#' matrix. This is a convenience wrapper around [astsa::Ksmooth()].
#'
#' @param y Numeric vector. The observed time series data. Missing values
#'   are handled via the time-varying observation matrix `A`.
#' @param x0 Numeric scalar or vector. Initial state estimate (prior mean of
#'   the state at time 0).
#' @param P0 Numeric scalar or matrix. Initial state covariance (prior
#'   variance/covariance of the state at time 0).
#' @param Phi Numeric scalar or matrix. State transition matrix. Describes how
#'   the state evolves: \eqn{x_t = \Phi x_{t-1} + v_t}.
#' @param Q Numeric scalar or matrix. State noise variance (process noise
#'   covariance). The variance of \eqn{v_t}.
#' @param A Array of dimension `(p, q, n)` where `p` is the observation
#'   dimension, `q` is the state dimension, and `n` is the number of time
#'   points. This is the time-varying observation matrix. To indicate a missing
#'   observation at time `t`, set `A[, , t]` to zeros.
#' @param R Numeric scalar or matrix. Observation noise variance (measurement
#'   noise covariance). The variance of \eqn{w_t}.
#'
#' @return A list with class `"kalman"` containing:
#'   \item{y}{The original observations.}
#'   \item{predicted}{List with `state` (predicted state estimates) and
#'     `var` (predicted state covariances).}
#'   \item{filtered}{List with `state` (filtered state estimates) and
#'     `var` (filtered state covariances).}
#'   \item{smoothed}{List with `state` (smoothed state estimates) and
#'     `var` (smoothed state covariances).}
#'   \item{result}{The raw output from [astsa::Ksmooth()] for advanced use.}
#'
#' @details
#' This function extends [kalman()] to handle missing observations and
#' multivariate state-space models. The model is:
#'
#' **State equation:**
#' \deqn{x_t = \Phi x_{t-1} + v_t, \quad v_t \sim N(0, Q)}
#'
#' **Observation equation:**
#' \deqn{y_t = A_t x_t + w_t, \quad w_t \sim N(0, R)}
#'
#' Note that the observation matrix \eqn{A_t} is time-varying. Missing
#' observations are handled by setting the corresponding row of \eqn{A_t}
#' to zero, which effectively removes that observation from the likelihood.
#'
#' @note
#' This function requires the `astsa` package. Install it with
#' `install.packages("astsa")` if needed.
#'
#' For simple univariate models without missing data, use [kalman()] instead.
#'
#' @references
#' Shumway, R. H., & Stoffer, D. S. (2017). *Time Series Analysis and Its
#' Applications: With R Examples* (4th ed.). Springer.
#'
#' @seealso [kalman()] for simpler models without missing data,
#'   [astsa::Ksmooth()] for the underlying implementation.
#'
#' @examples
#' \dontrun{
#' # Local level model with missing observations
#' set.seed(123)
#' n <- 100
#'
#' # True state is a random walk
#' x_true <- cumsum(rnorm(n, sd = 1))
#' # Observations with some missing (NA)
#' y <- x_true + rnorm(n, sd = 2)
#' missing_idx <- c(20, 21, 22, 50, 51)
#' y[missing_idx] <- NA
#'
#' # Create time-varying A matrix (1x1xn array)
#' # A[,,t] = 1 for observed, 0 for missing
#' A <- array(1, dim = c(1, 1, n))
#' A[, , missing_idx] <- 0
#'
#' # Replace NA with arbitrary value (will be ignored due to A=0)
#' y_filled <- y
#' y_filled[is.na(y_filled)] <- 0
#'
#' # Fit Kalman filter
#' result <- kalman_miss(
#'   y = y_filled,
#'   x0 = 0,
#'   P0 = 10,
#'   Phi = 1,
#'   Q = 1,
#'   A = A,
#'   R = 4
#' )
#'
#' # Plot results
#' plot(y, type = "p", pch = 16, cex = 0.5, col = "gray")
#' lines(x_true, col = "black", lwd = 2)
#' lines(as.numeric(result$smoothed$state), col = "red", lwd = 2)
#' }
#'
#' @export
kalman_miss <- function(y, x0, P0, Phi, Q, A, R) {
  # Check for astsa package
  if (!requireNamespace("astsa", quietly = TRUE)) {
    stop(
      "The 'astsa' package is required for kalman_miss().\n",
      "Install it with: install.packages(\"astsa\")",
      call. = FALSE
    )
  }

  # Input validation
  if (!is.numeric(y)) {
    stop("`y` must be numeric.", call. = FALSE)
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
  if (!is.numeric(A) || !is.array(A)) {
    stop("`A` must be a numeric array of dimension (p, q, n).", call. = FALSE)
  }
  if (!is.numeric(R)) {
    stop("`R` must be numeric.", call. = FALSE)
  }

  # Validate dimensions of A
  if (length(dim(A)) != 3) {
    stop("`A` must be a 3-dimensional array (p x q x n).", call. = FALSE)
  }
  if (dim(A)[3] != n) {
    stop("Third dimension of `A` must equal number of time points in `y`.",
         call. = FALSE)
  }

  # Convert scalar/vector inputs to matrices as needed for astsa
  # Match tswge behavior: sqrt for variance, chol for covariance matrix
  sQ <- sqrt(Q)
  sR <- sqrt(R)

  # If R is a matrix, use chol instead
  if (is.matrix(R)) {
    sR <- chol(R)
  }

  # Call astsa::Ksmooth exactly as tswge does
  kf_result <- astsa::Ksmooth(
    y = y,
    A = A,
    mu0 = x0,
    Sigma0 = P0,
    Phi = Phi,
    sQ = sQ,
    sR = sR
  )

  # Structure results
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
