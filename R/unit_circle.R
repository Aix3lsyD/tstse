#' Plot Roots on the Unit Circle
#'
#' Draws the unit circle and plots polynomial roots, useful for visualizing
#' the roots of AR or MA characteristic polynomials.
#'
#' @param roots Complex vector, or two numeric vectors. The roots to plot.
#'   Can be specified as:
#'   \itemize{
#'     \item A complex vector: `roots = c(0.5+0.5i, 0.8-0.3i)`
#'     \item Separate real/imaginary: `real = c(0.5, 0.8), imaginary = c(0.5, -0.3)`
#'   }
#' @param real Numeric vector. Real parts of roots (used if `roots` is NULL).
#' @param imaginary Numeric vector. Imaginary parts of roots (used if `roots`
#'   is NULL). Default is 0 (real roots).
#' @param main Character. Plot title. Default is "Roots and the Unit Circle".
#' @param show_legend Logical. If `TRUE` (default), show legend indicating
#'   roots inside/outside the unit circle.
#'
#' @return Invisibly returns a data frame with columns:
#'   \item{real}{Real part of each root.}
#'   \item{imaginary}{Imaginary part of each root.}
#'   \item{modulus}{Absolute value (distance from origin).}
#'   \item{inside}{Logical, TRUE if root is inside the unit circle.}
#'
#' @details
#' For a stationary AR process, all roots of the AR characteristic polynomial
#' must lie **outside** the unit circle (equivalently, the reciprocal roots
#' must lie inside).
#'
#' For an invertible MA process, all roots of the MA characteristic polynomial
#' must lie **outside** the unit circle.
#'
#' Roots on or very close to the unit circle indicate near-nonstationarity
#' or near-noninvertibility.
#'
#' @seealso [factor()] for factoring AR/MA polynomials and finding roots.
#'
#' @examples
#' # Plot roots using complex numbers
#' unit_circle(roots = c(0.5 + 0.5i, 0.5 - 0.5i, 1.2 + 0i))
#'
#' # Plot roots using real/imaginary parts
#' unit_circle(real = c(0.5, 0.5, 1.2), imaginary = c(0.5, -0.5, 0))
#'
#' # Single real root
#' unit_circle(real = 0.8, imaginary = 0)
#'
#' # AR(2) example: roots of 1 - 1.5z + 0.75z^2
#' roots <- polyroot(c(1, -1.5, 0.75))
#' unit_circle(roots = roots)
#'
#' @export
unit_circle <- function(roots = NULL, real = 0, imaginary = 0,
                        main = "Roots and the Unit Circle",
                        show_legend = TRUE) {
  # Handle input: either complex roots or real/imaginary parts
  if (!is.null(roots)) {
    if (!is.complex(roots) && !is.numeric(roots)) {
      stop("`roots` must be complex or numeric.", call. = FALSE)
    }
    roots <- as.complex(roots)
    real <- Re(roots)
    imaginary <- Im(roots)
  } else {
    if (!is.numeric(real)) {
      stop("`real` must be numeric.", call. = FALSE)
    }
    if (!is.numeric(imaginary)) {
      stop("`imaginary` must be numeric.", call. = FALSE)
    }
    if (length(imaginary) == 1 && length(real) > 1) {
      imaginary <- rep(imaginary, length(real))
    }
    if (length(real) != length(imaginary)) {
      stop("`real` and `imaginary` must have the same length.", call. = FALSE)
    }
  }

  n_roots <- length(real)
  if (n_roots == 0) {
    stop("No roots provided.", call. = FALSE)
  }

  # Calculate modulus (distance from origin)
  modulus <- sqrt(real^2 + imaginary^2)
  inside <- modulus < 1

  # Determine plot limits (accommodate roots outside unit circle)
  max_extent <- max(1.2, max(abs(real)) * 1.1, max(abs(imaginary)) * 1.1)
  lim <- c(-max_extent, max_extent)

  # Set up plot
  plot(NA, NA, xlim = lim, ylim = lim, asp = 1,
       xlab = expression(italic(Re)), ylab = expression(italic(Im)),
       main = main, las = 1)

  # Draw unit circle
  theta <- seq(0, 2 * pi, length.out = 361)
  lines(cos(theta), sin(theta), col = "blue", lty = 2, lwd = 2)

  # Draw axes
  abline(h = 0, col = "gray")
  abline(v = 0, col = "gray")

  # Plot roots with lines from origin
  # Color: red for inside, darkgreen for outside
  colors <- ifelse(inside, "red", "darkgreen")

  for (i in seq_len(n_roots)) {
    lines(c(0, real[i]), c(0, imaginary[i]), lwd = 2, col = colors[i])
    points(real[i], imaginary[i], pch = 19, col = colors[i], cex = 1.5)
  }

  # Add legend
  if (show_legend && any(inside) && any(!inside)) {
    legend("topright",
           legend = c("Inside unit circle", "Outside unit circle"),
           col = c("red", "darkgreen"), pch = 19, lwd = 2,
           bty = "n", cex = 0.9)
  } else if (show_legend && all(inside)) {
    legend("topright", legend = "Inside unit circle",
           col = "red", pch = 19, lwd = 2, bty = "n", cex = 0.9)
  } else if (show_legend && all(!inside)) {
    legend("topright", legend = "Outside unit circle",
           col = "darkgreen", pch = 19, lwd = 2, bty = "n", cex = 0.9)
  }

  # Build result
  result <- data.frame(
    real = real,
    imaginary = imaginary,
    modulus = modulus,
    inside = inside
  )

  invisible(result)
}
