#' Estimate FARMA Model
#'
#' Estimates a Fractional ARMA (FARMA/ARFIMA) model using grid search over
#' the fractional differencing parameter `d` and AIC selection for AR order.
#'
#' @param x Numeric vector. The time series to model.
#' @param d_range Numeric vector of length 2. Range for the fractional
#'   differencing parameter d: `c(lower, upper)`. Typically `c(0, 0.5)` for
#'   stationary long memory.
#' @param d_step Numeric. Step size for the grid search over d. Default is 0.05.
#' @param p_max Integer. Maximum AR order to consider. Default is 5.
#' @param n_back Integer. Number of backcasted values for initialization.
#'   Default is 500.
#'
#' @return A list with class `"est_farma"` containing:
#'   \item{d}{Estimated fractional differencing parameter.}
#'   \item{phi}{AR coefficients (numeric(0) if AR order is 0).}
#'   \item{p}{Selected AR order.}
#'   \item{var_a}{Estimated white noise variance.}
#'   \item{aic}{AIC value for the selected model.}
#'   \item{d_grid}{Grid of d values searched.}
#'   \item{aic_grid}{AIC values for each d in the grid.}
#'
#' @details
#' FARMA (also called ARFIMA) is a special case of GARMA where the Gegenbauer
#' frequency parameter `u = 1`, giving long memory at frequency 0. The model is:
#' \deqn{(1 - B)^d \phi(B) X_t = a_t}
#'
#' The estimation procedure:
#' \enumerate{
#'   \item Backcast to initialize the fractional differencing filter
#'   \item For each d in the grid:
#'     \itemize{
#'       \item Apply fractional differencing to get residuals
#'       \item Fit AR(p) model to residuals, selecting p by AIC
#'       \item Compute overall AIC including penalty for d
#'     }
#'   \item Select d with minimum overall AIC
#' }
#'
#' @section Relationship to GARMA:
#' FARMA is GARMA with a single factor at u = 1:
#' \itemize{
#'   \item FARMA d corresponds to GARMA d (not d/2 as in some parameterizations)
#'   \item Spectral peak is at frequency 0
#'   \item For d > 0, process has positive long-range dependence
#' }
#'
#' @seealso [est_garma()] for general GARMA estimation,
#'   [gen_geg()] for FARMA simulation (use u = 1),
#'   [fore_farma()] for FARMA forecasting.
#'
#' @examples
#' # Generate FARMA(0.3, 1, 0) process
#' set.seed(123)
#' # Using GARMA generation with u = 1
#' x <- gen_geg(n = 200, u = 1, d = 0.3/2, n_trunc = 500)
#'
#' # Estimate FARMA model
#' fit <- est_farma(x, d_range = c(0, 0.5), d_step = 0.05, p_max = 3)
#' print(fit)
#'
#' @export
est_farma <- function(x, d_range, d_step = 0.05, p_max = 5L, n_back = 500L) {
  # Input validation
  if (!is.numeric(x)) {
    stop("`x` must be a numeric vector.", call. = FALSE)
  }
  n <- length(x)
  if (n < 10L) {
    stop("`x` must have at least 10 observations.", call. = FALSE)
  }
  if (anyNA(x)) {
    stop("`x` contains NA values.", call. = FALSE)
  }

  if (!is.numeric(d_range) || length(d_range) != 2L) {
    stop("`d_range` must be a numeric vector of length 2.", call. = FALSE)
  }
  if (d_range[1] >= d_range[2]) {
    stop("`d_range[1]` must be less than `d_range[2]`.", call. = FALSE)
  }

  if (!is.numeric(d_step) || length(d_step) != 1L || d_step <= 0) {
    stop("`d_step` must be a positive number.", call. = FALSE)
  }

  if (!is.numeric(p_max) || length(p_max) != 1L || p_max < 0L) {
    stop("`p_max` must be a non-negative integer.", call. = FALSE)
  }
  p_max <- as.integer(p_max)

  if (!is.numeric(n_back) || length(n_back) != 1L || n_back < 1L) {
    stop("`n_back` must be a positive integer.", call. = FALSE)
  }
  n_back <- as.integer(n_back)

  # Internal function to compute AIC for AR model selection
  ar_aic_select <- function(w, p_max) {
    w <- w - mean(w)
    n_w <- length(w)

    # AIC for AR(0)
    var0 <- sum(w^2) / n_w
    aic0 <- log(var0) + 2 * 1 / n_w
    best_p <- 0L
    best_aic <- aic0
    best_var <- var0

    if (p_max > 0L) {
      for (p in seq_len(p_max)) {
        fit <- ar.burg(w, aic = FALSE, order.max = p)
        aic_p <- log(fit$var.pred) + 2 * (p + 1) / n_w
        if (aic_p < best_aic) {
          best_aic <- aic_p
          best_p <- p
          best_var <- fit$var.pred
        }
      }
    }

    list(p = best_p, aic = best_aic, var = best_var)
  }

  # Center the data
  x_centered <- x - mean(x)

  # Create extended sequence with backcast values
  nc <- n + n_back
  z <- numeric(nc)
  z[(n_back + 1):nc] <- x_centered

  # Fit high-order AR for backcasting
  ar_p_back <- min(20L, n - 1L)
  ar_back <- ar.yw(x_centered, order.max = ar_p_back, aic = FALSE, demean = FALSE)
  phi_back <- ar_back$ar

  # Backcast: z[n_back], z[n_back-1], ..., z[1]
  for (i in n_back:1L) {
    z[i] <- sum(phi_back * z[(i + 1):(i + ar_p_back)])
  }

  # Grid search over d
  # Note: For FARMA, the Gegenbauer parameter is d/2 (due to (1-B)^d = (1-2*1*B+B^2)^(d/2))
  # So we search d_internal = d/2 and report d = 2*d_internal
  d_internal_range <- d_range / 2
  d_grid_internal <- seq(d_internal_range[1], d_internal_range[2], by = d_step / 2)
  n_d <- length(d_grid_internal)

  # Storage for results
  results <- matrix(0, nrow = n_d, ncol = 4)
  colnames(results) <- c("d_internal", "p", "ar_aic", "total_aic")
  w_best <- NULL

  for (i in seq_len(n_d)) {
    d_int <- d_grid_internal[i]
    results[i, 1] <- d_int

    # Compute Gegenbauer coefficients for fractional differencing
    # For FARMA: u = 1, so we use gegenb(u = 1, d = -d_internal, n_coef)
    # The negative d is because we're applying the inverse filter
    C <- gegenb(u = 1, d = -d_int, n_coef = nc)

    # Apply fractional differencing filter
    w <- numeric(n)
    for (j in seq_len(n)) {
      # w[j] = sum_{k=0}^{n_back+j-2} C[k+1] * z[n_back+j-k]
      idx_z <- (n_back + j):2
      idx_c <- seq_len(n_back + j - 1)
      w[j] <- sum(C[idx_c] * z[idx_z])
    }
    w <- w - mean(w)

    # Select AR order and compute AIC
    ar_result <- ar_aic_select(w, p_max)
    results[i, 2] <- ar_result$p
    results[i, 3] <- ar_result$aic

    # Total AIC with penalty for d (if d != 0)
    d_penalty <- if (d_int != 0) 2 * 2 / n else 0
    results[i, 4] <- ar_result$aic + d_penalty

    # Store best w for final model fitting
    if (i == 1 || results[i, 4] < min(results[1:(i-1), 4])) {
      w_best <- w
    }
  }

  # Select winner
  winner_idx <- which.min(results[, 4])
  d_internal_opt <- unname(results[winner_idx, 1])
  d_opt <- 2 * d_internal_opt
  p_opt <- as.integer(results[winner_idx, 2])
  aic_opt <- unname(results[winner_idx, 4])

  # Re-fit AR model at optimal d to get coefficients
  # Need to recompute w at optimal d
  C_opt <- gegenb(u = 1, d = -d_internal_opt, n_coef = nc)
  w_opt <- numeric(n)
  for (j in seq_len(n)) {
    idx_z <- (n_back + j):2
    idx_c <- seq_len(n_back + j - 1)
    w_opt[j] <- sum(C_opt[idx_c] * z[idx_z])
  }
  w_opt <- w_opt - mean(w_opt)

  if (p_opt == 0L) {
    phi_opt <- numeric(0)
    var_opt <- mean(w_opt^2)
  } else {
    ar_fit <- ar.burg(w_opt, aic = FALSE, order.max = p_opt)
    phi_opt <- ar_fit$ar
    var_opt <- ar_fit$var.pred
  }

  # Build result
  result <- list(
    d = d_opt,
    phi = phi_opt,
    p = p_opt,
    var_a = var_opt,
    aic = aic_opt,
    d_grid = 2 * d_grid_internal,
    aic_grid = results[, 4]
  )

  class(result) <- "est_farma"
  result
}


#' @export
print.est_farma <- function(x, ...) {
  cat("FARMA Model Estimation\n")
  cat("----------------------\n")
  cat("Fractional differencing: d =", round(x$d, 4), "\n")
  cat("AR order: p =", x$p, "\n")
  if (x$p > 0) {
    cat("AR coefficients:", round(x$phi, 4), "\n")
  }
  cat("White noise variance:", round(x$var_a, 4), "\n")
  cat("AIC:", round(x$aic, 4), "\n")
  invisible(x)
}
