#' Factor Component Decomposition of a Time Series
#'
#' Decomposes a time series into components based on the factors of a fitted
#' AR model using partial fraction decomposition.
#'
#' @param x Numeric vector. The time series to decompose.
#' @param p Integer. The AR order to fit. If `aic = TRUE`, this is the maximum
#'   order to consider.
#' @param n_comp Integer. Number of components to extract. If greater than the
#'   number of factors, it will be reduced to the number of factors.
#' @param aic Logical. If `TRUE`, select AR order by AIC up to `p`. If `FALSE`
#'   (default), use exactly order `p`.
#' @param plot Logical. If `TRUE` (default), plot the original series and
#'   extracted components.
#'
#' @return A list with class `"factor_comp"` containing:
#'   \item{n_comp}{Number of components extracted.}
#'   \item{n_factors}{Total number of factors in the AR model.}
#'   \item{components}{Matrix of dimension `(n_comp, n)` with each row being
#'     a component series. Values before index `p` are zero.}
#'   \item{phi}{AR coefficients from the fitted model.}
#'   \item{roots}{Complex roots of the AR polynomial, sorted by absolute value.}
#'   \item{factors}{Factor table from [factor_ts()] showing the AR factors.}
#'
#' @details
#' The function works by:
#' \enumerate{
#'   \item Fitting an AR(p) model using Burg's method ([stats::ar.burg()])
#'   \item Finding the roots of the AR characteristic polynomial
#'   \item Using partial fraction decomposition to determine how each factor
#'     contributes to the overall model
#'   \item Computing each component's contribution to the original series
#' }
#'
#' For an AR polynomial with roots \eqn{r_1, r_2, \ldots, r_p}, the transfer
#' function \eqn{1/\phi(B)} can be written as a sum of partial fractions.
#' Each factor (real root or complex conjugate pair) corresponds to one
#' component of the decomposition.
#'
#' Complex conjugate pairs are combined into a single real-valued component,
#' so the number of factors is:
#' \deqn{n_{factors} = n_{real} + n_{complex pairs}}
#'
#' @seealso [factor_ts()] for displaying AR factors,
#'   [est_ar()] for AR model estimation.
#'
#' @examples
#' # Generate AR(4) with two factor pairs
#' set.seed(123)
#' # Two complex conjugate pairs at different frequencies
#' x <- gen_arma(n = 200, phi = c(1.2, -0.85, 0.4, -0.3), plot = FALSE)
#'
#' # Decompose into 2 components
#' result <- factor_comp(x, p = 4, n_comp = 2)
#'
#' # Access components
#' dim(result$components)
#'
#' @export
factor_comp <- function(x, p, n_comp, aic = FALSE, plot = TRUE) {
  # Input validation
  if (!is.numeric(x)) {
    stop("`x` must be a numeric vector.", call. = FALSE)
  }
  n <- length(x)
  if (n < 3L) {
    stop("`x` must have at least 3 observations.", call. = FALSE)
  }
  if (anyNA(x)) {
    stop("`x` contains NA values.", call. = FALSE)
  }

  if (!is.numeric(p) || length(p) != 1L || p < 1L) {
    stop("`p` must be a positive integer.", call. = FALSE)
  }
  p <- as.integer(p)
  if (p >= n) {
    stop("`p` must be less than the length of `x`.", call. = FALSE)
  }

  if (!is.numeric(n_comp) || length(n_comp) != 1L || n_comp < 1L) {
    stop("`n_comp` must be a positive integer.", call. = FALSE)
  }
  n_comp <- as.integer(n_comp)

  if (!is.logical(aic) || length(aic) != 1L || is.na(aic)) {
    stop("`aic` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.logical(plot) || length(plot) != 1L || is.na(plot)) {
    stop("`plot` must be TRUE or FALSE.", call. = FALSE)
  }

  # Fit AR model using Burg's method
  ar_fit <- ar.burg(x, order.max = p, aic = aic)
  phi <- ar_fit$ar
  p_actual <- length(phi)

  if (p_actual == 0L) {
    stop("AR model fit resulted in order 0. Try a larger `p`.", call. = FALSE)
  }

  # Get factor decomposition for display
  factor_result <- factor_ts(phi)

  # Find roots of AR polynomial
  # Characteristic equation: 1 - phi_1*z - phi_2*z^2 - ... - phi_p*z^p = 0
  # polyroot solves: coef[1] + coef[2]*z + ... + coef[p+1]*z^p = 0
  coef_poly <- c(1, -phi)
  roots_raw <- polyroot(coef_poly)

  # Sort roots by absolute value
  root_abs <- Mod(roots_raw)
  sort_idx <- order(root_abs)
  roots <- roots_raw[sort_idx]
  nr <- length(roots)

  # Count number of factors
  # Real roots (imaginary part near zero) count as 1
  # Complex conjugate pairs count as 1 (combined)
  tol <- 1e-5
  is_real <- abs(Im(roots)) <= tol
  n_factors <- as.integer(sum(is_real) + sum(!is_real) / 2)

  # Adjust n_comp if needed
  if (n_comp > n_factors) {
    warning(paste0("`n_comp` (", n_comp, ") exceeds number of factors (",
                   n_factors, "). Using ", n_factors, " components."),
            call. = FALSE)
    n_comp <- n_factors
  }

  # Build Vandermonde-like matrix for partial fractions
  # fac[i, j] = (1/root[j])^(nr - i)
  fac_matrix <- matrix(0, nr, nr)
  for (i in seq_len(nr)) {
    fac_matrix[i, ] <- (1 / roots)^(nr - i)
  }

  # Invert to get partial fraction coefficients
  fac_inv <- solve(fac_matrix)

  # Compute coefficients for each factor
  # Real roots: use coefficient directly
  # Complex pairs: sum the two conjugate coefficients
  coefs_fac <- matrix(0, n_factors, p_actual)
  j <- 1L
  for (i in seq_len(n_factors)) {
    if (abs(Im(roots[j])) <= tol) {
      # Real root
      coefs_fac[i, ] <- fac_inv[j, ]
      j <- j + 1L
    } else {
      # Complex conjugate pair - combine
      coefs_fac[i, ] <- fac_inv[j, ] + fac_inv[j + 1L, ]
      j <- j + 2L
    }
  }

  # Compute component time series
  # x_comp[k, t] = sum_{j=1}^{p} coefs_fac[k, j] * x[t - p + j]
  x_comp <- matrix(0, n_comp, n)

  for (k in seq_len(n_comp)) {
    for (t in p_actual:n) {
      for (j in seq_len(p_actual)) {
        x_comp[k, t] <- x_comp[k, t] + coefs_fac[k, j] * x[t - p_actual + j]
      }
    }
  }

  # Take real part (imaginary should be ~0 for real data)
  x_comp <- Re(x_comp)

  # Plot if requested
  if (plot) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))

    # Determine common y-axis limits
    max_comp <- max(x_comp)
    min_comp <- min(x_comp)
    max_x <- max(x)
    min_x <- min(x)
    ylim <- c(min(min_x, min_comp), max(max_x, max_comp))

    par(mfrow = c(n_comp + 1L, 1L), mar = c(3.5, 4, 2, 1))

    # Plot original series
    plot(seq_len(n), x, type = "l", ylim = ylim,
         main = "Original Series", xlab = "", ylab = "")

    # Plot each component
    t_range <- p_actual:n
    for (k in seq_len(n_comp)) {
      plot(t_range, x_comp[k, t_range], type = "l", ylim = ylim,
           main = paste("Component", k), xlab = "", ylab = "")
    }
  }

  # Build result
  result <- list(
    n_comp = n_comp,
    n_factors = n_factors,
    components = x_comp,
    phi = phi,
    roots = roots,
    factors = factor_result
  )

  class(result) <- "factor_comp"
  invisible(result)
}


#' @export
print.factor_comp <- function(x, ...) {
  cat("Factor Component Decomposition\n")
  cat("------------------------------\n")
  cat("AR order:", length(x$phi), "\n")
  cat("Number of factors:", x$n_factors, "\n")
  cat("Components extracted:", x$n_comp, "\n")
  cat("\nAR coefficients:", round(x$phi, 4), "\n")
  cat("\nUse $components to access the component time series.\n")
  cat("Use $factors to see the factor decomposition.\n")
  invisible(x)
}
