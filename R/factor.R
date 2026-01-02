#' Factor Table for ARMA Polynomials
#'
#' Compute and display factor tables for AR and/or MA polynomials,
#' showing roots, absolute reciprocals, and system frequencies.
#'
#' @param phi Numeric vector of AR coefficients (default 0 = none).
#' @param theta Numeric vector of MA coefficients (default 0 = none).
#' @param print Logical, whether to print the factor table (default TRUE).
#'
#' @return Invisibly returns a list with components:
#'   \item{ar}{Data frame of AR factors (NULL if no AR)}
#'   \item{ma}{Data frame of MA factors (NULL if no MA)}
#'   Each data frame contains: factor1, factor2, root, abs_recip, frequency
#' @export
#'
#' @examples
#' # AR(2) process
#' factor_ts(phi = c(1.6, -0.9))
#'
#' # ARMA(1,1)
#' factor_ts(phi = 0.7, theta = 0.4)
#'
#' # Capture results
#' f <- factor_ts(phi = c(1.6, -0.9), print = FALSE)
#' f$ar
factor_ts <- function(phi = 0, theta = 0, print = TRUE) {

  # Input validation
  if (!is.numeric(phi)) stop("`phi` must be numeric")
  if (!is.numeric(theta)) stop("`theta` must be numeric")

  result <- list(ar = NULL, ma = NULL)

  # Process AR polynomial
  if (sum(phi^2) != 0) {
    coef <- c(1, -phi)
    result$ar <- .compute_factors(coef)

    if (print) {
      cat("\nCoefficients of AR polynomial:\n")
      cat(format(round(phi, 4), nsmall = 4), "\n\n")
      cat("                        AR Factor Table\n")
      cat("Factor                 Roots                Abs Recip    System Freq\n")
      .print_factors(result$ar)
      cat("\n")
    }
  }

  # Process MA polynomial
  if (sum(theta^2) != 0) {
    coef <- c(1, -theta)
    result$ma <- .compute_factors(coef)

    if (print) {
      cat("\nCoefficients of MA polynomial:\n")
      cat(format(round(theta, 4), nsmall = 4), "\n\n")
      cat("                        MA Factor Table\n")
      cat("Factor                 Roots                Abs Recip    System Freq\n")
      .print_factors(result$ma)
      cat("\n")
    }
  }

  invisible(result)
}


#' Compute factors from polynomial coefficients
#' @param coef Polynomial coefficients (leading 1)
#' @return Data frame of factors
#' @keywords internal
#' @noRd
.compute_factors <- function(coef, tol = 1e-5) {

  roots <- polyroot(coef)
  n_roots <- length(roots)

  # Sort by imaginary part to pair conjugates
  roots <- roots[order(Im(roots))]

  factors <- list()
  i <- 1
  idx <- 1

  while (i <= n_roots) {
    root <- roots[i]

    # Check if this is part of a complex conjugate pair
    is_complex <- abs(Im(root)) > tol
    has_conjugate <- (i < n_roots) &&
      (abs(Re(roots[i+1]) - Re(root)) < tol) &&
      (abs(Im(roots[i+1]) + Im(root)) < tol)

    if (is_complex && has_conjugate) {
      # Quadratic factor from conjugate pair
      abs_recip <- 1 / Mod(root)
      freq <- abs(Arg(root) / (2 * pi))
      fac1 <- -2 * Re(1 / root)
      fac2 <- 1 / Re(root * Conj(root))

      factors[[idx]] <- data.frame(
        factor1 = fac1,
        factor2 = fac2,
        root_real = Re(root),
        root_imag = abs(Im(root)),
        abs_recip = abs_recip,
        frequency = freq,
        is_quadratic = TRUE
      )

      i <- i + 2
    } else {
      # Linear factor from real root
      abs_recip <- Re(1 / Mod(root))
      freq <- abs(Arg(root) / (2 * pi))
      fac1 <- -Re(1 / root)

      factors[[idx]] <- data.frame(
        factor1 = fac1,
        factor2 = 0,
        root_real = Re(root),
        root_imag = 0,
        abs_recip = abs_recip,
        frequency = freq,
        is_quadratic = FALSE
      )

      i <- i + 1
    }
    idx <- idx + 1
  }

  result <- do.call(rbind, factors)
  result <- result[order(-result$abs_recip), ]
  rownames(result) <- NULL
  result
}


#' Print factor table
#' @param factors Data frame from .compute_factors
#' @keywords internal
#' @noRd
.print_factors <- function(factors) {

  for (i in seq_len(nrow(factors))) {
    f <- factors[i, ]

    # Format components
    fac1 <- abs(f$factor1)
    fac1_fmt <- format(round(fac1, 4), nsmall = 4, trim = TRUE)
    fac2_fmt <- format(round(abs(f$factor2), 4), nsmall = 4, trim = TRUE)
    root_r <- format(round(f$root_real, 4), nsmall = 4, trim = TRUE)
    root_i <- format(round(f$root_imag, 4), nsmall = 4, trim = TRUE)
    abs_fmt <- format(round(f$abs_recip, 4), nsmall = 4, trim = TRUE)
    freq_fmt <- format(round(f$frequency, 4), nsmall = 4, trim = TRUE)

    # Build factor string
    if (f$is_quadratic) {
      sign1 <- if (f$factor1 < 0) "+" else "-"
      sign2 <- if (f$factor2 < 0) "-" else "+"
      factor_str <- sprintf("1%s%sB%s%sB^2", sign1, fac1_fmt, sign2, fac2_fmt)
      root_str <- sprintf("%s+-%si", root_r, root_i)
    } else {
      sign1 <- if (f$factor1 < 0) "+" else "-"
      factor_str <- sprintf("1%s%sB", sign1, fac1_fmt)
      root_str <- root_r
    }

    # Print aligned
    cat(sprintf("%-22s %-20s %s       %s\n",
                factor_str, root_str, abs_fmt, freq_fmt))
  }
}
