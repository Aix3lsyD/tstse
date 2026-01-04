#' AIC/BIC Model Selection for AR Models
#'
#' Select optimal AR order using information criteria.
#'
#' @param x Numeric vector, the time series.
#' @param p Integer vector, AR orders to compare (default 1:5).
#' @param type Character, criterion to use: "aic" (default), "aicc", or "bic".
#' @param method Character, estimation method: "mle" (default), "burg", or "yw".
#' @param cores Integer, number of cores for parallel processing.
#'   Default NULL uses sequential processing (cores = 1) to prevent nested
#'   parallelization when called from within parallel contexts.
#'   Set explicitly to use parallel processing.
#'
#' @return A list with components:
#'   \item{type}{Criterion used}
#'   \item{method}{Estimation method used}
#'   \item{value}{Best criterion value}
#'   \item{p}{Selected AR order}
#'   \item{phi}{AR coefficients at selected order}
#'   \item{xbar}{Sample mean}
#'   \item{vara}{Residual variance at selected order}
#'   \item{table}{Data frame of all orders and their criterion values}
#' @export
#'
#' @examples
#' x <- gen_arma(n = 200, phi = c(1.5, -0.75), plot = FALSE, seed = 123)
#'
#' # Sequential (default)
#' aic_ar(x, p = 1:5)
#'
#' \dontrun{
#' # Parallel (uses multiple cores) - must be explicit
#' aic_ar(x, p = 1:10, cores = 2)
#' }
aic_ar <- function(x,
                   p = 1:5,
                   type = c("aic", "aicc", "bic"),
                   method = c("mle", "burg", "yw"),
                   cores = NULL) {

  # Input validation
  if (!is.numeric(x) || length(x) == 0) {
    stop("`x` must be a non-empty numeric vector")
  }
  if (!is.numeric(p) || any(p < 0)) {
    stop("`p` must be a vector of non-negative integers")
  }

  type <- match.arg(type)
  method <- match.arg(method)

  # =========================================================================
  # DEFENSIVE DEFAULT: When cores is NULL, default to sequential (1L).
  # This prevents nested parallelization when aic_ar is called from within

  # parallel workers (e.g., inside wbg_boot's bootstrap loop).
  #
  # Users who want parallel aic_ar at the top level should explicitly set
  # cores = N or cores = 0.
  # =========================================================================
  if (is.null(cores)) {
    cores <- 1L
  } else {
    cores <- get_cores(cores)
  }

  n <- length(x)
  xbar <- mean(x)
  x_centered <- x - xbar

  # Define worker function
  fit_ar_order <- function(j) {
    if (j == 0) {
      phi <- 0
      res <- x_centered
    } else {
      fit <- est_ar(x_centered, p = j, method = method, factor = FALSE)
      phi <- fit$phi
      res <- backcast(x_centered, phi = phi, theta = 0, n_back = 50)
    }

    avar <- mean(res^2)

    list(
      p = j,
      phi = phi,
      avar = avar,
      aic = log(avar) + 2 * (j + 1) / n,
      aicc = log(avar) + (n + j + 1) / (n - j - 3),
      bic = log(avar) + (j + 1) * log(n) / n
    )
  }

  # Fit models (parallel or sequential)
  results <- pmap(p, fit_ar_order, cores = cores)

  # Build comparison table
  ic_table <- data.frame(
    p = sapply(results, `[[`, "p"),
    aic = sapply(results, `[[`, "aic"),
    aicc = sapply(results, `[[`, "aicc"),
    bic = sapply(results, `[[`, "bic")
  )

  # Find best by selected criterion
  ic_values <- sapply(results, `[[`, type)
  best_idx <- which.min(ic_values)
  best <- results[[best_idx]]

  structure(
    list(
      type = type,
      method = method,
      value = best[[type]],
      p = best$p,
      phi = best$phi,
      xbar = xbar,
      vara = best$avar,
      table = ic_table
    ),
    class = "aic_ar"
  )
}


#' @export
print.aic_ar <- function(x, ...) {
  cat("\nAR Order Selection (", toupper(x$type), ", method: ", x$method, ")\n\n", sep = "")
  cat("Selected order p =", x$p, "\n")
  cat("Criterion value:", round(x$value, 4), "\n\n")
  cat("AR Coefficients:\n")
  print(round(x$phi, 4))
  cat("\nMean:", round(x$xbar, 4), "\n")
  cat("Residual Variance:", round(x$vara, 4), "\n\n")
  cat("Comparison Table:\n")
  print(round(x$table, 4))
  invisible(x)
}
