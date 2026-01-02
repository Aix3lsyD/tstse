#' Estimate GARMA Model Parameters
#'
#' Estimate Gegenbauer ARMA (GARMA) model parameters via grid search
#' with AIC-based AR order selection.
#'
#' @param x Numeric vector. The time series to analyze.
#' @param u_range Numeric vector of length 3: `c(low, high, increment)` defining
#'   the grid for the Gegenbauer frequency parameter u.
#' @param lambda_range Numeric vector of length 3: `c(low, high, increment)`
#'   defining the grid for the Gegenbauer persistence parameter lambda.
#' @param p_max Integer. Maximum AR order to consider. Default is 5.
#' @param n_back Integer. Number of backcasted values to prepend. Default is 500.
#' @param cores Integer. Number of cores for parallel processing.
#'   Default NULL uses `getOption("tstse.cores", 1)`.
#'   Set to 0 to use all available cores.
#'
#' @return An object of class `"est_garma"` with components:
#'   \item{u}{Estimated Gegenbauer frequency parameter}
#'   \item{lambda}{Estimated Gegenbauer persistence parameter}
#'   \item{phi}{Estimated AR coefficients (0 if p = 0)}
#'   \item{p}{Selected AR order}
#'   \item{var_a}{Estimated residual variance}
#'   \item{aic}{GARMA AIC value at the selected parameters}
#'   \item{table}{Data frame with all grid search results}
#'
#' @details
#' The function estimates parameters for a single-factor GARMA model:
#' \deqn{(1 - 2uB + B^2)^{-\lambda} \phi(B) X_t = a_t}
#'
#' The estimation algorithm:
#' \enumerate{
#'   \item Backcast the series using AR(20) Yule-Walker
#'   \item For each (u, lambda) on the grid:
#'     \itemize{
#'       \item Apply Gegenbauer filter to remove long-memory component
#'       \item Fit AR(0), AR(1), ..., AR(p_max) models
#'       \item Select AR order by AIC
#'       \item Compute GARMA AIC with 2-parameter penalty
#'     }
#'   \item Return parameters at minimum GARMA AIC
#' }
#'
#' The GARMA AIC includes a penalty for estimating both u and lambda:
#' \deqn{AIC_{GARMA} = AIC_{AR} + 4/n}
#' when lambda is non-zero.
#'
#' @section Computational Notes:
#' Grid search can be computationally expensive. Parallel processing is
#' supported via the `cores` parameter. For fine grids, consider:
#' \itemize{
#'   \item Starting with coarse grids to narrow the search region
#'   \item Using `cores = 0` to use all available cores
#' }
#'
#' @seealso [gen_garma()] for generating GARMA realizations,
#'   [fore_garma()] for forecasting,
#'   [gegenb()] for Gegenbauer polynomial coefficients.
#'
#' @examples
#' # Generate GARMA data
#' x <- gen_garma(n = 200, u = 0.8, lambda = 0.3, phi = 0.5,
#'                plot = FALSE, seed = 123)
#'
#' # Estimate parameters (coarse grid for speed)
#' fit <- est_garma(x, u_range = c(0.6, 0.9, 0.1),
#'                  lambda_range = c(0.1, 0.4, 0.1), p_max = 3)
#' print(fit)
#'
#' \donttest{
#' # Finer grid with parallel processing
#' fit2 <- est_garma(x, u_range = c(0.7, 0.9, 0.05),
#'                   lambda_range = c(0.2, 0.4, 0.05),
#'                   p_max = 5, cores = 2)
#' }
#'
#' @export
est_garma <- function(x, u_range, lambda_range, p_max = 5L,
                      n_back = 500L, cores = NULL) {

  # Input validation
  if (!is.numeric(x) || length(x) < 20L) {
    stop("`x` must be a numeric vector with at least 20 observations.", call. = FALSE)
  }
  if (anyNA(x)) {
    stop("`x` contains NA values.", call. = FALSE)
  }

  if (!is.numeric(u_range) || length(u_range) != 3L) {
    stop("`u_range` must be c(low, high, increment).", call. = FALSE)
  }
  if (u_range[1] >= u_range[2]) {
    stop("`u_range[1]` must be less than `u_range[2]`.", call. = FALSE)
  }
  if (u_range[3] <= 0) {
    stop("`u_range[3]` (increment) must be positive.", call. = FALSE)
  }
  if (u_range[1] < -1 || u_range[2] > 1) {
    stop("`u_range` values must be in [-1, 1].", call. = FALSE)
  }

  if (!is.numeric(lambda_range) || length(lambda_range) != 3L) {
    stop("`lambda_range` must be c(low, high, increment).", call. = FALSE)
  }
  if (lambda_range[1] >= lambda_range[2]) {
    stop("`lambda_range[1]` must be less than `lambda_range[2]`.", call. = FALSE)
  }
  if (lambda_range[3] <= 0) {
    stop("`lambda_range[3]` (increment) must be positive.", call. = FALSE)
  }
  if (lambda_range[1] < 0) {
    stop("`lambda_range` values must be non-negative.", call. = FALSE)
  }

  if (!is.numeric(p_max) || length(p_max) != 1L || p_max < 0L) {
    stop("`p_max` must be a non-negative integer.", call. = FALSE)
  }
  p_max <- as.integer(p_max)

  if (!is.numeric(n_back) || length(n_back) != 1L || n_back < 1L) {
    stop("`n_back` must be a positive integer.", call. = FALSE)
  }
  n_back <- as.integer(n_back)

  cores <- get_cores(cores)
  n <- length(x)

  # Build grids
  u_grid <- seq(u_range[1], u_range[2], by = u_range[3])
  lambda_grid <- seq(lambda_range[1], lambda_range[2], by = lambda_range[3])
  grid <- expand.grid(u = u_grid, lambda = lambda_grid)

  # Center data and prepare backcasted sequence
  x_centered <- x - mean(x)
  nc <- n + n_back
  z <- c(rep(0, n_back), x_centered)

  # Fit AR(20) for backcasting (or less if series is short)
  ar_order <- min(20L, n - 1L)
  ar_fit <- ar.yw(x_centered, order.max = ar_order, aic = FALSE, demean = FALSE)
  ar_phi <- ar_fit$ar

  # Backcast: z[n_back], z[n_back-1], ..., z[1]
  for (i in n_back:1) {
    z[i] <- sum(ar_phi * z[(i + 1):(i + ar_order)])
  }

  # Internal function: fit AR models and select by AIC
  ar_aic_select <- function(w, p_max) {
    w_centered <- w - mean(w)
    n_w <- length(w_centered)

    # p = 0 case (white noise)
    var0 <- sum(w_centered^2) / n_w
    aic0 <- log(var0) + 2 / n_w

    if (p_max == 0L) {
      return(list(p = 0L, var_a = var0, aic = aic0))
    }

    # p = 1, ..., p_max
    aic_values <- numeric(p_max + 1L)
    var_values <- numeric(p_max + 1L)
    aic_values[1] <- aic0
    var_values[1] <- var0

    for (p in seq_len(p_max)) {
      fit <- ar.burg(w_centered, aic = FALSE, order.max = p)
      var_values[p + 1] <- fit$var.pred
      aic_values[p + 1] <- log(fit$var.pred) + 2 * (p + 1) / n_w
    }

    best_idx <- which.min(aic_values)
    list(p = best_idx - 1L, var_a = var_values[best_idx], aic = aic_values[best_idx])
  }

  # Worker function for grid search
  fit_garma_point <- function(idx) {
    u_val <- grid$u[idx]
    d_val <- grid$lambda[idx]

    # Compute Gegenbauer coefficients (negative d to remove factor)
    C <- gegenb(u = u_val, d = -d_val, n_coef = nc)

    # Apply filter: w[j] = sum(C[1:(n_back+j-1)] * z[(n_back+j):2])
    w <- numeric(n)
    for (j in seq_len(n)) {
      # Indices for convolution
      c_idx <- seq_len(n_back + j - 1)
      z_idx <- (n_back + j):2
      w[j] <- sum(C[c_idx] * z[z_idx])
    }

    # Select AR order by AIC
    ar_result <- ar_aic_select(w, p_max)

    # GARMA AIC: add penalty for u and lambda parameters (2*2/n when d != 0)
    garma_aic <- ar_result$aic + 4 * (d_val != 0) / n

    list(
      u = u_val,
      lambda = d_val,
      p = ar_result$p,
      ar_aic = ar_result$aic,
      garma_aic = garma_aic,
      var_a = ar_result$var_a,
      w = w  # Keep filtered series for final AR fit
    )
  }

  # Fit all grid points (parallel or sequential)
  results <- pmap(seq_len(nrow(grid)), fit_garma_point, cores = cores)

  # Build results table
  result_table <- data.frame(
    u = sapply(results, `[[`, "u"),
    lambda = sapply(results, `[[`, "lambda"),
    p = sapply(results, `[[`, "p"),
    ar_aic = sapply(results, `[[`, "ar_aic"),
    garma_aic = sapply(results, `[[`, "garma_aic")
  )

  # Find best by GARMA AIC
  best_idx <- which.min(result_table$garma_aic)
  best <- results[[best_idx]]

  # Get final AR coefficients
  if (best$p == 0L) {
    phi <- 0
    var_a <- mean((best$w - mean(best$w))^2)
  } else {
    final_fit <- ar.burg(best$w - mean(best$w), aic = FALSE, order.max = best$p)
    phi <- as.numeric(final_fit$ar)
    var_a <- final_fit$var.pred
  }

  structure(
    list(
      u = best$u,
      lambda = best$lambda,
      phi = phi,
      p = best$p,
      var_a = var_a,
      aic = best$garma_aic,
      table = result_table
    ),
    class = "est_garma"
  )
}


#' @export
print.est_garma <- function(x, ...) {
  cat("\nGARMA Parameter Estimation\n")
  cat("==========================\n\n")
  cat("Gegenbauer frequency (u):", round(x$u, 4), "\n")
  cat("Gegenbauer persistence (lambda):", round(x$lambda, 4), "\n")
  cat("AR order (p):", x$p, "\n")

  if (x$p > 0) {
    cat("\nAR Coefficients:\n")
    print(round(x$phi, 4))
  }

  cat("\nResidual Variance:", round(x$var_a, 4), "\n")
  cat("GARMA AIC:", round(x$aic, 4), "\n")
  cat("\nGrid search evaluated", nrow(x$table), "parameter combinations.\n")

  invisible(x)
}
