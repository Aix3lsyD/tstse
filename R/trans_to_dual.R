#' Transform to Dual Time Scale
#'
#' Transform a time series from original (potentially unequally-spaced) time
#' to an equally-spaced dual time scale using a power-law (G-lambda)
#' transformation.
#'
#' @param x Numeric vector, the time series data.
#' @param lambda Numeric, the shape parameter for the transformation.
#'   - `lambda = 0`: Geometric spacing
#'   - `lambda = 1`: Linear spacing (identity-like)
#'   - Other values: Power-law spacing
#' @param offset Numeric, shift parameter to avoid singularities at t=0.
#'   Default is 60.
#' @param h Numeric, scaling parameter. If <= 1 (for lambda=0) or <= 0
#'   (for lambda != 0), it is computed automatically. Default is 0.
#' @param plot Logical, whether to plot the transformed series. Default TRUE.
#'
#' @return A list containing:
#'   \item{int_x}{Dual time points (unequally spaced in original scale)}
#'   \item{int_y}{Interpolated data values on dual scale}
#'   \item{h}{The h parameter used (computed or provided)}
#'
#' @details
#' The G-lambda transformation maps time according to a power law controlled
#' by the `lambda` parameter. This is useful for:
#' \itemize{
#'   \item Analyzing data with non-linear sampling patterns
#'   \item Working with astronomical data (e.g., sunspot observations)
#'   \item Preparing data for ARMA modeling on transformed scales
#' }
#'
#' The transformation creates an equally-spaced time series in the "dual"
#' domain by interpolating the original data onto a new time grid.
#'
#' Use [trans_to_original()] to reverse the transformation.
#'
#' @references
#' Gray, H. L., Zhang, N.-F., & Woodward, W. A. (1989). On generalized
#' fractional processes. *Journal of Time Series Analysis*, 10(3), 233-257.
#'
#' @seealso [trans_to_original()] for the inverse transformation,
#'   [est_glambda()] for estimating optimal parameters.
#'
#' @examples
#' # Generate sample data
#' set.seed(123)
#' x <- cumsum(rnorm(100))
#'
#' # Transform with geometric spacing (lambda = 0)
#' result <- trans_to_dual(x, lambda = 0, offset = 60, plot = FALSE)
#' length(result$int_y)  # Transformed series length
#'
#' # Transform with power-law spacing
#' result <- trans_to_dual(x, lambda = 0.5, offset = 60, plot = FALSE)
#'
#' @export
trans_to_dual <- function(x, lambda, offset = 60, h = 0, plot = TRUE) {

 # Input validation
  if (!is.numeric(x) || length(x) == 0) {
    stop("`x` must be a non-empty numeric vector", call. = FALSE)
  }
  if (!is.numeric(lambda) || length(lambda) != 1) {
    stop("`lambda` must be a single numeric value", call. = FALSE)
  }
  if (!is.numeric(offset) || length(offset) != 1 || offset < 0) {
    stop("`offset` must be a non-negative number", call. = FALSE)
  }

  n <- length(x)
  shift <- offset

  # Helper function to compute h value
  h_value <- function(n, shift, lambda, m = n) {
    if (lambda == 0) {
      ((shift + n) / (shift + 1))^(1 / (m - 1))
    } else {
      ((shift + n)^lambda - (shift + 1)^lambda) / ((m - 1) * lambda)
    }
  }

  # Compute h if not provided or invalid
  if ((h <= 1 && lambda == 0) || (h <= 0 && lambda != 0)) {
    h <- h_value(n, shift, lambda, m = n)
  }

  # Generate original time grid
  tt <- seq(shift + 1, shift + n)

  # Generate dual time grid based on lambda
  if (lambda == 0) {
    # Geometric spacing
    k <- log(shift + 1) / log(h)
    m <- log(shift + n) / log(h)
    kk <- seq(k, m)
    ht2 <- h^kk
  } else if (lambda == 1) {
    # Linear spacing
    ht2 <- seq(shift + 1, shift + n)
  } else {
    # Power-law spacing
    k <- ((1 + shift)^lambda - 1) / (h * lambda)
    m <- ((n + shift)^lambda - 1) / (h * lambda)
    kk <- seq(k, m)
    ht2 <- (kk * h * lambda + 1)^(1 / lambda)
  }

  # Linear interpolation
  int_data <- stats::approx(tt, x, xout = ht2, rule = 2)
  int_y <- int_data$y[round(ht2, 2) >= (shift + 1)]
  int_x <- ht2[round(ht2, 2) >= (shift + 1)]

  # Plot if requested
  if (plot) {
    op <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(op))

    graphics::par(mar = c(3.8, 2.5, 1, 1))
    t_plot <- seq_along(int_y)
    graphics::plot(t_plot, int_y, type = "o", cex = 0.5, pch = 16,
                   cex.lab = 0.75, cex.axis = 0.75, lwd = 0.75,
                   xlab = "", ylab = "")
    graphics::mtext(side = 1, text = "Time", line = 1, cex = 0.9)
    graphics::mtext(side = 3, text = "Dual Scale Realization", line = 0, cex = 0.9)
  }

  list(int_x = int_x, int_y = int_y, h = h)
}
