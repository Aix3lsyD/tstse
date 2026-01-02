#' Internal Utility Functions
#'
#' Shared internal helper functions for time series generation.
#'
#' @name utils_internal
#' @keywords internal
NULL


#' Build combined nonstationary operator
#'
#' Combines lambda and seasonal factors into a single coefficient vector
#' for the nonstationary operator.
#'
#' @param lambda Numeric vector of lambda coefficients.
#' @param s Integer seasonal period.
#' @param dlam Integer length of lambda (0 if all zeros).
#'
#' @return Numeric vector of combined coefficients, or 0 if no factors.
#'
#' @noRd
.build_nonstationary_operator <- function(lambda, s, dlam) {
  if (dlam == 0 && s == 0) {
    return(0)
  }

  # Build seasonal factor: (1 - B^s) represented as c(0,0,...,1) with s elements
  seas <- NULL
  if (s > 0) {
    seas <- rep(0, s)
    seas[s] <- 1
  }

  # Combine lambda and seasonal using mult()
  if (dlam > 0 && s > 0) {
    result <- mult(lambda, seas)
    return(result$model_coef)
  } else if (dlam > 0) {
    return(lambda)
  } else {
    return(seas)
  }
}


#' Apply inverse of nonstationary operator (integration)
#'
#' Applies the inverse filter to convert a stationary series to nonstationary.
#' Uses stats::filter(method="recursive") for C-level performance.
#'
#' The recursion is: `x[t] = y[t] + sum(lambdas[j] * x[t-j])`
#'
#' @param y Numeric vector, input series from arima.sim.
#' @param lambdas Numeric vector of filter coefficients.
#' @param dlams Integer, length of lambdas (order of nonstationary operator).
#' @param d Integer, differencing order.
#' @param n Integer, desired output length.
#' @param spin Integer, burn-in period.
#'
#' @return Numeric vector of length n.
#'
#' @noRd
.apply_inverse_operator <- function(y, lambdas, dlams, d, n, spin) {
  d1 <- d + dlams + 1L
  nd <- n + d + dlams + 1L
  ndspin <- nd + spin - 1L

  if (dlams == 0) {
    # No additional nonstationary factors, just extract
    start_idx <- spin + d1
    return(y[start_idx:(start_idx + n - 1)])
  }

  # Use stats::filter with method="recursive" (C implementation)
  # Recursion: x[t] = y[t] + lambdas[1]*x[t-1] + ... + lambdas[k]*x[t-k]
  # Apply filter to the relevant slice starting at d1
  y_slice <- y[d1:ndspin]
  filtered <- stats::filter(y_slice, filter = lambdas, method = "recursive")
  filtered <- as.numeric(filtered)  # Remove ts attributes

  # Construct full array (with leading zeros to maintain indexing)
  xfull <- c(rep(0, d1 - 1L), filtered)

  # Extract after spin-up
  start_idx <- spin + d1
  xfull[start_idx:(start_idx + n - 1)]
}
