#' Plot Time Series
#'
#' Plot a time series realization using base R or ggplot2.
#'
#' @param x Numeric vector or ts object, the time series.
#' @param style Integer, plotting style: 0 = base R (default), 1 = ggplot2.
#' @param xlab Character, x-axis label (default "Time").
#' @param ylab Character, y-axis label (default "").
#' @param main Character, plot title (default "").
#' @param col Character, line/point color (default "black").
#' @param text_size Numeric, text size for ggplot2 (default 12).
#' @param lwd Numeric, line width (default 0.75).
#' @param cex Numeric, point size (default 0.5).
#' @param cex.lab Numeric, axis label size for base R (default 0.75).
#' @param cex.axis Numeric, axis text size for base R (default 0.75).
#' @param xlim Numeric vector of length 2, x-axis limits (default NULL).
#' @param ylim Numeric vector of length 2, y-axis limits (default NULL).
#'
#' @return For style = 0, returns invisibly. For style = 1, returns a ggplot object.
#' @export
#'
#' @details
#' For style = 1 (ggplot2), the packages \code{ggplot2} and \code{zoo} must be
#' installed. If x is a ts object, the x-axis will show dates.
#'
#' For series with n <= 200, points are shown. For longer series, only lines
#' are drawn.
#'
#' @examples
#' x <- gen_arma(n = 200, phi = 0.9, plot = FALSE, seed = 123)
#'
#' # Base R plot
#' plotts(x)
#'
#' # With title and color
#' plotts(x, main = "AR(1) Process", col = "blue")
#'
#' # ggplot2 style (if installed)
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plotts(x, style = 1)
#' }
plotts <- function(x,
                   style = 0L,
                   xlab = "Time",
                   ylab = "",
                   main = "",
                   col = "black",
                   text_size = 12,
                   lwd = 0.75,
                   cex = 0.5,
                   cex.lab = 0.75,
                   cex.axis = 0.75,
                   xlim = NULL,
                   ylim = NULL) {

  # Input validation
  if (!is.numeric(x) || length(x) == 0) {
    stop("`x` must be a non-empty numeric vector or ts object")
  }
  if (!style %in% c(0, 1)) {
    stop("`style` must be 0 (base R) or 1 (ggplot2)")
  }

  n <- length(x)
  is_ts <- inherits(x, "ts")

  if (style == 0) {
    .plotts_base(x, n, is_ts, xlab, ylab, main, col, lwd, cex,
                 cex.lab, cex.axis, xlim, ylim)
    invisible(NULL)
  } else {
    .plotts_ggplot(x, n, is_ts, xlab, ylab, main, col, lwd, text_size, xlim, ylim)
  }
}


#' Base R time series plot
#' @noRd
.plotts_base <- function(x, n, is_ts, xlab, ylab, main, col, lwd, cex,
                         cex.lab, cex.axis, xlim, ylim) {

  t <- seq_len(n)
  type <- if (n <= 200) "o" else "l"
  pch <- if (n <= 200) 16 else NA

  plot(t, x, type = type, pch = pch, cex = cex,
       cex.lab = cex.lab, cex.axis = cex.axis, lwd = lwd,
       xlab = xlab, ylab = ylab, col = col, main = main,
       xlim = xlim, ylim = ylim, las = 1)
}


#' ggplot2 time series plot
#' @noRd
.plotts_ggplot <- function(x, n, is_ts, xlab, ylab, main, col, lwd, text_size, xlim, ylim) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for style = 1. Please install it.")
  }

  # Build data frame
  if (is_ts) {
    if (!requireNamespace("zoo", quietly = TRUE)) {
      stop("Package 'zoo' is required for ts objects with style = 1. Please install it.")
    }
    df <- data.frame(
      x = zoo::as.yearmon(time(x)),
      y = as.numeric(x)
    )
  } else {
    df <- data.frame(
      x = seq_len(n),
      y = as.numeric(x)
    )
  }

  # Build plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_line(color = col, linewidth = lwd) +
    ggplot2::ggtitle(main) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) +
    ggplot2::theme(text = ggplot2::element_text(size = text_size)) +
    ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)

  # Add points for short series
  if (n <= 200) {
    p <- p + ggplot2::geom_point()
  }

  p
}
