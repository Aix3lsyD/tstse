#' Transform from Dual to Original Time Scale
#'
#' Transform a time series from the dual (equally-spaced) time scale back
#' to the original time scale using the inverse G-lambda transformation.
#'
#' @param xd Numeric vector, the dual-scale time series data.
#' @param lambda Numeric, the shape parameter (must match the value used
#'   in [trans_to_dual()]).
#' @param offset Numeric, shift parameter (must match the value used
#'   in [trans_to_dual()]).
#' @param h Numeric, scaling parameter (must match the value returned
#'   by [trans_to_dual()]).
#' @param plot Logical, whether to plot the transformed series. Default TRUE.
#'
#' @return Numeric vector of interpolated data on the original time scale.
#'
#' @details
#' This function reverses the transformation performed by [trans_to_dual()].
#' It maps the equally-spaced dual time series back to the original
#' (potentially unequally-spaced) time scale.
#'
#' The `lambda`, `offset`, and `h` parameters must match those used
#' (or returned) by the corresponding [trans_to_dual()] call.
#'
#' @seealso [trans_to_dual()] for the forward transformation.
#'
#' @examples
#' # Generate sample data
#' set.seed(123)
#' x <- cumsum(rnorm(100))
#'
#' # Transform to dual scale
#' dual <- trans_to_dual(x, lambda = 0.5, offset = 60, plot = FALSE)
#'
#' # Transform back to original scale
#' original <- trans_to_original(dual$int_y, lambda = 0.5,
#'                                offset = 60, h = dual$h, plot = FALSE)
#'
#' @export
trans_to_original <- function(xd, lambda, offset, h, plot = TRUE) {

  # Input validation
  if (!is.numeric(xd) || length(xd) == 0) {
    stop("`xd` must be a non-empty numeric vector", call. = FALSE)
  }
  if (!is.numeric(lambda) || length(lambda) != 1) {
    stop("`lambda` must be a single numeric value", call. = FALSE)
  }
  if (!is.numeric(offset) || length(offset) != 1) {
    stop("`offset` must be a single numeric value", call. = FALSE)
  }
  if (!is.numeric(h) || length(h) != 1) {
    stop("`h` must be a single numeric value", call. = FALSE)
  }

  n <- length(xd)
  start <- offset + 1
  int_start <- offset

  # Generate dual time grid
  if (lambda == 0) {
    # Geometric spacing
    m <- log(start) / log(h)
    ht <- h^seq(m, m + n - 1)
  } else {
    # Power-law spacing
    m <- (start^lambda - 1) / (h * lambda)
    kk <- seq(m, m + n - 1)
    ht <- (kk * h * lambda + 1)^(1 / lambda)
  }

  # Original time grid
  tt <- seq(int_start, int_start + n - 1)

  # Linear interpolation back to original scale
  int_data <- stats::approx(ht, xd, xout = tt, rule = 2)$y

  # Plot if requested
  if (plot) {
    op <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(op))

    graphics::par(mar = c(3.8, 2.5, 1, 1))
    t_plot <- seq_along(int_data)
    graphics::plot(t_plot, int_data, type = "o", cex = 0.5, pch = 16,
                   cex.lab = 0.75, cex.axis = 0.75, lwd = 0.75,
                   xlab = "", ylab = "")
    graphics::mtext(side = 1, text = "Time", line = 1, cex = 0.9)
    graphics::mtext(side = 3, text = "Original Scale Realization", line = 0, cex = 0.9)
  }

  int_data
}
