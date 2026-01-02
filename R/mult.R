#' Multiply Polynomial Factors
#'
#' Multiply AR or MA polynomial factors together to get combined coefficients.
#'
#' @param ... Numeric vectors, each representing factor coefficients.
#'   For example, for factor (1 - 0.8B), pass \code{c(0.8)}.
#'   For factor (1 - 1.6B + 0.9B^2), pass \code{c(1.6, -0.9)}.
#'
#' @return A list with components:
#'   \item{model_coef}{Combined model coefficients (excluding leading 1)}
#'   \item{poly_coef}{Full polynomial coefficients (including leading 1)}
#'   \item{char_poly}{Polynomial object (if PolynomF available, otherwise NULL)}
#' @export
#'
#' @details
#' Takes factors of the form (1 - phi_1 B - phi_2 B^2 - ...) and multiplies
#' them together. Input each factor's coefficients (phi_1, phi_2, ...) as
#' a numeric vector.
#'
#' If the \code{PolynomF} package is installed, it will be used and a polynomial
#' object returned in \code{char_poly}. Otherwise, falls back to base R convolution.
#'
#' @examples
#' # Two AR(1) factors: (1 - 0.8B)(1 - 0.5B)
#' mult(c(0.8), c(0.5))
#'
#' # AR(1) and AR(2) factor: (1 - 0.7B)(1 - 1.2B + 0.6B^2)
#' mult(c(0.7), c(1.2, -0.6))
mult <- function(...) {

  factors <- list(...)

  # Remove NULL and zero factors
  factors <- Filter(function(f) {
    !is.null(f) && !(length(f) == 1 && f == 0)
  }, factors)

  # If all factors removed, return identity
  if (length(factors) == 0) {
    return(list(model_coef = 0, poly_coef = 1, char_poly = NULL))
  }

  # Validate inputs
  for (i in seq_along(factors)) {
    if (!is.numeric(factors[[i]])) {
      stop("All factors must be numeric vectors")
    }
  }

  # Convert to polynomial coefficients: (1, -phi_1, -phi_2, ...)
  polys <- lapply(factors, function(f) c(1, -f))

  # Use PolynomF if available, otherwise base R
  if (requireNamespace("PolynomF", quietly = TRUE)) {
    result <- .mult_polynomf(polys)
  } else {
    result <- .mult_convolve(polys)
  }

  result
}


#' Multiply using PolynomF package
#' @noRd
.mult_polynomf <- function(polys) {

  # Convert to polynomial objects
  poly_objs <- lapply(polys, PolynomF::polynom)

  # Multiply
  char_poly <- Reduce(`*`, poly_objs)

  # Extract coefficients
  poly_coef <- coef(char_poly)
  poly_coef <- round(poly_coef, 10)

  # Model coefficients (exclude leading 1, negate)
  if (length(poly_coef) > 1) {
    model_coef <- -poly_coef[-1]
  } else {
    model_coef <- 0
  }

  list(
    model_coef = model_coef,
    poly_coef = poly_coef,
    char_poly = char_poly
  )
}


#' Multiply using base R convolution
#' @noRd
.mult_convolve <- function(polys) {

  result <- polys[[1]]
  if (length(polys) > 1) {
    for (i in 2:length(polys)) {
      result <- convolve(result, rev(polys[[i]]), type = "open")
    }
  }

  # Clean up tiny numerical errors
  poly_coef <- round(result, 10)

  # Model coefficients (exclude leading 1, negate)
  if (length(poly_coef) > 1) {
    model_coef <- -poly_coef[-1]
  } else {
    model_coef <- 0
  }

  list(
    model_coef = model_coef,
    poly_coef = poly_coef,
    char_poly = NULL
  )
}
