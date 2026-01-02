#' AIC/BIC Model Selection for AR Models (Burg Method)
#'
#' Select optimal AR order using Burg estimation. This is a convenience
#' wrapper around \code{\link{aic_ar}} with \code{method = "burg"}.
#'
#' @param x Numeric vector, the time series.
#' @param p Integer vector, AR orders to compare (default 1:5).
#' @param type Character, criterion to use: "aic" (default), "aicc", or "bic".
#' @param cores Integer, number of cores for parallel processing.
#'   Default is 1L (sequential). This function is often called from within
#'   parallel contexts (e.g., bootstrap loops), so it defaults to sequential
#'   execution to avoid nested parallelization.
#'
#' @return A list with components (see \code{\link{aic_ar}} for details):
#'   \item{type}{Criterion used}
#'   \item{value}{Best criterion value}
#'   \item{p}{Selected AR order}
#'   \item{phi}{AR coefficients at selected order}
#'   \item{vara}{Residual variance at selected order}
#' @export
#'
#' @seealso \code{\link{aic_ar}} for the general function with method selection.
#'
#' @examples
#' x <- gen_arma(n = 200, phi = c(1.5, -0.75), plot = FALSE, seed = 123)
#' aic_burg(x, p = 1:5)
aic_burg <- function(x, p = 1:5, type = c("aic", "aicc", "bic"), cores = 1L) {

  type <- match.arg(type)

  result <- aic_ar(x = x, p = p, type = type, method = "burg", cores = cores)

  # Match original output structure (no method, no xbar, no table)
  structure(
    list(
      type = result$type,
      value = result$value,
      p = result$p,
      phi = result$phi,
      vara = result$vara
    ),
    class = "aic_burg"
  )
}


#' @export
print.aic_burg <- function(x, ...) {
  cat("\nAR Order Selection - Burg (", toupper(x$type), ")\n\n", sep = "")
  cat("Selected order p =", x$p, "\n")
  cat("Criterion value:", round(x$value, 4), "\n\n")
  cat("AR Coefficients:\n")
  print(round(x$phi, 4))
  cat("\nResidual Variance:", round(x$vara, 4), "\n")
  invisible(x)
}
