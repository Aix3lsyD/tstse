#' Check if AR Model is Stationary
#'
#' Checks whether an AR model with given coefficients is stationary by
#' verifying all roots of the characteristic polynomial lie outside the
#' unit circle.
#'
#' @param phi Numeric vector of AR coefficients.
#' @param tol Numeric. Roots must have modulus > tol to be considered
#'   outside the unit circle. Default is 1.001 to provide a small buffer.
#'
#' @return Logical. TRUE if stationary, FALSE otherwise.
#'
#' @details
#' An AR(p) process with coefficients \eqn{\phi_1, \ldots, \phi_p} is
#' stationary if all roots of the characteristic polynomial
#' \deqn{1 - \phi_1 z - \phi_2 z^2 - \cdots - \phi_p z^p = 0}
#' lie outside the unit circle (i.e., have modulus > 1).
#'
#' The default tolerance of 1.001 provides a small buffer to handle
#' numerical precision issues and avoid near-unit-root edge cases.
#'
#' @seealso [aic_ar_mle()] for MLE AR fitting with stationarity checking.
#'
#' @examples
#' # Stationary AR(1)
#' check_stationary(0.9)        # TRUE
#' check_stationary(0.5)        # TRUE
#'
#' # Non-stationary (unit root)
#' check_stationary(1.0)        # FALSE
#' check_stationary(1.01)       # FALSE
#'
#' # Stationary AR(2)
#' check_stationary(c(0.5, 0.3))  # TRUE
#'
#' # White noise (empty phi is trivially stationary)
#' check_stationary(numeric(0))  # TRUE
#'
#' @export
check_stationary <- function(phi, tol = 1.001) {
  # White noise is trivially stationary
  if (is.null(phi) || length(phi) == 0) {
    return(TRUE)
  }

  # Characteristic polynomial: 1 - phi[1]*z - phi[2]*z^2 - ...
  # Coefficients in order: c(1, -phi[1], -phi[2], ...)
  poly_coefs <- c(1, -phi)

  # Find roots
  roots <- polyroot(poly_coefs)

  # All roots must be strictly outside unit circle
  all(Mod(roots) > tol)
}


#' Fit AR Model via MLE with Information Criterion Selection
#'
#' Fits AR models of orders 1 through p_max using maximum likelihood
#' estimation, selects the best model according to the specified information
#' criterion. Optionally verifies stationarity and falls back through top
#' candidates if the best model is non-stationary.
#'
#' @param x Numeric vector. The time series to fit.
#' @param p_max Integer. Maximum AR order to consider. Can also be a vector
#'   of orders to try (e.g., 1:5).
#' @param criterion Character. Information criterion for model selection:
#'   `"aic"` (default), `"aicc"` (corrected AIC), or `"bic"`.
#' @param stationary Logical. If TRUE (default), check that the selected model
#'   is stationary and fall back to next best if not. If FALSE, return the
#'   best model by criterion regardless of stationarity.
#' @param n_best Integer. Number of top models (by criterion) to check for
#'   stationarity before failing. Ignored if `stationary = FALSE`.
#'   Default is 3.
#' @param tol Numeric. Tolerance for stationarity check (roots must have
#'   modulus > tol). Ignored if `stationary = FALSE`. Default is 1.001.
#' @param silent Logical. If TRUE (default), suppress convergence warnings
#'   from `arima()`.
#'
#' @return A list containing:
#'   \item{phi}{Numeric vector of AR coefficients.}
#'   \item{p}{Integer. The selected AR order.}
#'   \item{criterion_value}{The value of the selection criterion.}
#'   \item{method_used}{Character. Always `"mle"` for this function.}
#'   \item{all_criteria}{Named numeric vector of criterion values for all
#'     successfully fit orders.}
#'
#' @details
#' The function attempts to fit AR(p) models for each order using
#' `arima()` with `method = "CSS-ML"` (conditional sum of squares
#' initialization followed by MLE refinement).
#'
#' Models are ranked by the specified criterion. If `stationary = TRUE`,
#' the function checks the top `n_best` models for stationarity and
#' returns the best stationary model. If no stationary model is found among
#' the top candidates, an error is raised.
#'
#' If `stationary = FALSE`, the best model by criterion is returned
#' regardless of whether it is stationary. This can be useful for
#' simulation studies.
#'
#' @seealso [check_stationary()], [aic_burg()] for Burg estimation
#'   (always stationary), [aic_ar()] for order selection with different methods.
#'
#' @examples
#' # Simulate AR(2) data
#' set.seed(123)
#' x <- arima.sim(list(ar = c(0.7, 0.2)), n = 200)
#'
#' # Fit with AIC selection
#' fit <- aic_ar_mle(x, p_max = 5)
#' fit$p    # Selected order
#' fit$phi  # Coefficients
#'
#' # Fit with BIC selection (tends to select simpler models)
#' fit_bic <- aic_ar_mle(x, p_max = 5, criterion = "bic")
#'
#' # Allow non-stationary estimates (for simulation studies)
#' fit_free <- aic_ar_mle(x, p_max = 5, stationary = FALSE)
#'
#' @export
aic_ar_mle <- function(x, p_max, criterion = c("aic", "aicc", "bic"),
                       stationary = TRUE, n_best = 3L, tol = 1.001,
                       silent = TRUE) {

  criterion <- match.arg(criterion)
  n <- length(x)

  # Handle p_max as either single value or vector
  if (length(p_max) == 1) {
    if (p_max < 1) {
      stop("No valid AR orders specified", call. = FALSE)
    }
    orders <- seq_len(p_max)
  } else {
    orders <- p_max
  }
  orders <- orders[orders >= 1]  # ensure positive orders

  if (length(orders) == 0) {
    stop("No valid AR orders specified", call. = FALSE)
  }

  # Fit models for each order
  results <- list()
  all_criteria <- c()

  for (p in orders) {
    fit <- tryCatch({
      if (silent) {
        suppressWarnings(
          arima(x, order = c(p, 0, 0), method = "CSS-ML", include.mean = TRUE)
        )
      } else {
        arima(x, order = c(p, 0, 0), method = "CSS-ML", include.mean = TRUE)
      }
    }, error = function(e) NULL)

    if (!is.null(fit)) {
      # Extract AR coefficients
      coef_names <- names(fit$coef)
      ar_idx <- grep("^ar", coef_names)
      phi <- if (length(ar_idx) > 0) as.numeric(fit$coef[ar_idx]) else numeric(0)

      # Calculate information criteria
      # k = number of parameters: AR coefficients + intercept + variance
      k <- p + 2
      loglik <- fit$loglik

      aic_val <- -2 * loglik + 2 * k
      aicc_val <- aic_val + (2 * k * (k + 1)) / max(n - k - 1, 1)
      bic_val <- -2 * loglik + k * log(n)

      crit_val <- switch(criterion,
                         aic = aic_val,
                         aicc = aicc_val,
                         bic = bic_val)

      results[[length(results) + 1]] <- list(
        p = p,
        phi = phi,
        criterion_value = crit_val,
        aic = aic_val,
        aicc = aicc_val,
        bic = bic_val
      )

      all_criteria[as.character(p)] <- crit_val
    }
  }

  if (length(results) == 0) {
    stop("MLE estimation failed for all AR orders. ",
         "The series may be too short or have numerical issues.",
         call. = FALSE)
  }

  # Sort by criterion (ascending - lower is better)
  crit_vals <- sapply(results, `[[`, "criterion_value")
  sorted_idx <- order(crit_vals)

  # If not enforcing stationarity, return best model by criterion
  if (!stationary) {
    best <- results[[sorted_idx[1]]]
    return(list(
      phi = best$phi,
      p = best$p,
      criterion_value = best$criterion_value,
      method_used = "mle",
      all_criteria = all_criteria
    ))
  }

  # Check top N for stationarity
  n_check <- min(n_best, length(results))

  for (i in seq_len(n_check)) {
    candidate <- results[[sorted_idx[i]]]
    if (check_stationary(candidate$phi, tol = tol)) {
      return(list(
        phi = candidate$phi,
        p = candidate$p,
        criterion_value = candidate$criterion_value,
        method_used = "mle",
        all_criteria = all_criteria
      ))
    }
  }

  stop("No stationary model found among top ", n_best, " candidates by ",
       toupper(criterion), ". Consider increasing n_best or p_max, ",
       "or use aic_burg() which guarantees stationarity.",
       call. = FALSE)
}
