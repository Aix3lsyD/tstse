#' Sample Partial Autocorrelation Function
#'
#' Compute and plot the sample partial autocorrelations for a time series.
#'
#' @param x Numeric vector, the time series.
#' @param lag_max Integer, maximum lag (default 5).
#' @param method Character, estimation method: "yw" (default), "burg", or "mle".
#' @param plot Logical, whether to plot (default TRUE).
#' @param limits Logical, whether to show 95% confidence limits (default FALSE).
#' @param cores Integer, number of cores for parallel processing.
#'   Default NULL uses \code{getOption("tstse.cores", 1)}.
#'
#' @return A list with component:
#'   \item{pacf}{Numeric vector of partial autocorrelations for lags 1 to lag_max}
#' @export
#'
#' @details
#' Computes the sample PACF by fitting AR(k) models for k = 1, ..., lag_max
#' and extracting the last coefficient from each fit.
#'
#' The 95% confidence limits are \eqn{\pm 2/\sqrt{n}}, valid under the null
#' hypothesis that the true PACF is zero beyond some lag.
#'
#' @examples
#' # AR(2) process - PACF should cut off after lag 2
#' x <- gen_arma(n = 200, phi = c(1.5, -0.75), plot = FALSE, seed = 123)
#' pacf_ts(x, lag_max = 10)
#'
#' # With confidence limits
#' pacf_ts(x, lag_max = 10, limits = TRUE)
#'
#' # Using Burg method
#' pacf_ts(x, lag_max = 10, method = "burg")
pacf_ts <- function(x,
                    lag_max = 5L,
                    method = c("yw", "burg", "mle"),
                    plot = TRUE,
                    limits = FALSE,
                    cores = 1L) {

  # Input validation
  if (!is.numeric(x) || length(x) == 0) {
    stop("`x` must be a non-empty numeric vector")
  }
  if (!is.numeric(lag_max) || lag_max < 1) {
    stop("`lag_max` must be a positive integer")
  }

  method <- match.arg(method)
  cores <- get_cores(cores)

  n <- length(x)
  if (lag_max > n - 1) lag_max <- n - 1

  # Compute PACF by fitting AR(k) for k = 1, ..., lag_max
  compute_pacf_k <- function(k) {
    fit <- est_ar(x, p = k, method = method, factor = FALSE)
    fit$phi[k]  # Last coefficient is PACF at lag k
  }

  pacf_values <- unlist(pmap(seq_len(lag_max), compute_pacf_k, cores = cores))

  # Plot
  if (plot) {
    .plot_pacf(pacf_values, lag_max, n, method, limits)
  }

  # Return
  result <- list(pacf = pacf_values)
  if (plot) {
    invisible(result)
  } else {
    result
  }
}


#' Plot PACF
#' @noRd
.plot_pacf <- function(pacf_values, lag_max, n, method, limits) {

  op <- par(mar = c(4, 4, 3, 1))
  on.exit(par(op))

  k <- seq_len(lag_max)

  # Method label
  method_label <- switch(method,
                         yw = "Yule-Walker",
                         burg = "Burg",
                         mle = "MLE"
  )

  plot(k, pacf_values, type = "h", lwd = 2,
       xlab = "Lag", ylab = "PACF",
       ylim = c(-1, 1),
       main = paste(method_label, "Partial Autocorrelations"))
  abline(h = 0)

  # 95% confidence limits
  if (limits) {
    ul <- 2 / sqrt(n)
    abline(h = c(-ul, ul), lty = 2, col = "blue")
  }
}
