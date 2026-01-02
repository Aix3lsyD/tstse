#' Estimate G-Lambda Parameters
#'
#' Estimate optimal lambda and offset parameters for G-lambda transformation
#' by grid search, minimizing the ACF difference between first and second
#' halves of the transformed series.
#'
#' @param x Numeric vector. The time series to analyze.
#' @param lambda_range Numeric vector of length 2: `c(low, high)` defining
#'   the range for lambda. Default is `c(0, 1)`.
#' @param offset_range Numeric vector of length 2: `c(low, high)` defining
#'   the range for offset. Default is `c(0, 100)`.
#' @param lambda_by Numeric. Increment for lambda grid. Default is 0.1.
#'
#' @return An object of class `"est_glambda"` with components:
#'   \item{lambda}{Optimal lambda value}
#'   \item{offset}{Optimal offset value}
#'   \item{q}{Q-statistic at optimal parameters}
#'   \item{table}{Data frame with all lambda/offset/Q combinations}
#'
#' @details
#' The function finds optimal G-lambda transformation parameters by:
#' \enumerate{
#'   \item For each (lambda, offset) combination on the grid:
#'     \itemize{
#'       \item Transform data to dual scale using [trans_to_dual()]
#'       \item Split transformed series in half
#'       \item Compute ACF for each half
#'       \item Calculate Q = sum of squared ACF differences
#'     }
#'   \item Select parameters that minimize Q
#' }
#'
#' A smaller Q indicates more stationary behavior in the transformed domain,
#' suggesting the transformation has successfully removed non-stationarity.
#'
#' @section Grid Search:
#' Lambda is searched from `lambda_range[1]` to `lambda_range[2]` in steps
#' of `lambda_by`. Offset is searched over all integers from `offset_range[1]`
#' to `offset_range[2]`.
#'
#' @seealso [trans_to_dual()] for the transformation,
#'   [fore_glambda()] for forecasting with estimated parameters.
#'
#' @examples
#' # Generate non-stationary data
#' set.seed(123)
#' t <- 1:200
#' x <- sin(2 * pi * t / 50) * sqrt(t) + cumsum(rnorm(200, sd = 0.5))
#'
#' # Estimate optimal parameters (use narrow range for speed)
#' fit <- est_glambda(x, lambda_range = c(0, 0.5), offset_range = c(10, 30))
#' print(fit)
#'
#' # Use estimated parameters for transformation
#' dual <- trans_to_dual(x, lambda = fit$lambda, offset = fit$offset,
#'                       plot = FALSE)
#'
#' @export
est_glambda <- function(x, lambda_range = c(0, 1), offset_range = c(0, 100),
                        lambda_by = 0.1) {

  # Input validation
  if (!is.numeric(x) || length(x) < 20L) {
    stop("`x` must be a numeric vector with at least 20 observations.", call. = FALSE)
  }
  if (anyNA(x)) {
    stop("`x` contains NA values.", call. = FALSE)
  }

  # Handle list input (as in original)
  if (is.list(x)) {
    x <- as.numeric(unlist(x))
  }

  if (!is.numeric(lambda_range) || length(lambda_range) != 2L) {
    stop("`lambda_range` must be c(low, high).", call. = FALSE)
  }
  if (lambda_range[1] > lambda_range[2]) {
    stop("`lambda_range[1]` must be <= `lambda_range[2]`.", call. = FALSE)
  }

  if (!is.numeric(offset_range) || length(offset_range) != 2L) {
    stop("`offset_range` must be c(low, high).", call. = FALSE)
  }
  if (offset_range[1] > offset_range[2]) {
    stop("`offset_range[1]` must be <= `offset_range[2]`.", call. = FALSE)
  }
  if (offset_range[1] < 0) {
    stop("`offset_range` values must be non-negative.", call. = FALSE)
  }

  if (!is.numeric(lambda_by) || length(lambda_by) != 1L || lambda_by <= 0) {
    stop("`lambda_by` must be a positive number.", call. = FALSE)
  }

  n <- length(x)

  # Internal: compute ACF-based Q statistic
  # Splits data in halves, computes ACF difference
  compute_q <- function(data) {
    n_data <- length(data)
    lag0 <- floor(n_data / 8)
    if (lag0 < 1) lag0 <- 1

    half <- floor(n_data / 2)
    acf1 <- as.numeric(acf(data[1:half], lag.max = lag0, plot = FALSE)$acf)
    acf2 <- as.numeric(acf(data[(half + 1):n_data], lag.max = lag0, plot = FALSE)$acf)

    sum((acf1 - acf2)^2)
  }

  # Internal: transform and compute Q for given parameters
  pick_q <- function(data, offset, lambda) {
    result <- trans_to_dual(data, lambda = lambda, offset = offset, plot = FALSE)
    compute_q(result$int_y)
  }

  # Build grids
  lambda_grid <- seq(lambda_range[1], lambda_range[2], by = lambda_by)
  offset_grid <- seq(offset_range[1], offset_range[2], by = 1)

  # Grid search
  results <- vector("list", length(lambda_grid))

  for (j in seq_along(lambda_grid)) {
    lam <- lambda_grid[j]

    # For this lambda, find best offset
    q_values <- sapply(offset_grid, function(off) {
      tryCatch(
        pick_q(x, offset = off, lambda = lam),
        error = function(e) Inf
      )
    })

    best_offset_idx <- which.min(q_values)
    if (length(best_offset_idx) > 1) best_offset_idx <- best_offset_idx[1]

    results[[j]] <- data.frame(
      lambda = lam,
      offset = offset_grid[best_offset_idx],
      Q = q_values[best_offset_idx]
    )
  }

  # Combine results
  q_table <- do.call(rbind, results)

  # Find overall best
  best_idx <- which.min(q_table$Q)
  if (length(best_idx) > 1) best_idx <- best_idx[1]

  best_lambda <- q_table$lambda[best_idx]
  best_offset <- q_table$offset[best_idx]
  best_q <- q_table$Q[best_idx]

  structure(
    list(
      lambda = best_lambda,
      offset = best_offset,
      q = best_q,
      table = q_table
    ),
    class = "est_glambda"
  )
}


#' @export
print.est_glambda <- function(x, ...) {
  cat("\nG-Lambda Parameter Estimation\n")
  cat("=============================\n\n")
  cat("Optimal lambda:", round(x$lambda, 4), "\n")
  cat("Optimal offset:", x$offset, "\n")
  cat("Q-statistic:", round(x$q, 6), "\n\n")
  cat("Grid search evaluated", nrow(x$table), "lambda values.\n")
  invisible(x)
}
