#' Generate ARIMA/SARIMA Process
#'
#' Generate a realization from an ARIMA or seasonal ARIMA process.
#'
#' @param n Integer, length of output series.
#' @param phi Numeric vector, AR coefficients (default 0 = none).
#' @param theta Numeric vector, MA coefficients (default 0 = none).
#' @param d Integer, order of differencing (default 0).
#' @param s Integer, seasonal period (default 0 = non-seasonal).
#' @param mu Numeric, theoretical mean (default 0).
#' @param vara Numeric, white noise variance (default 1).
#' @param plot Logical, whether to plot the series (default TRUE).
#' @param seed Integer, random seed for reproducibility (default NULL).
#' @param spin Integer, burn-in period (default 2000).
#'
#' @return Numeric vector of length n.
#' @export
#'
#' @examples
#' # AR(2) process
#' x <- gen_arima(n = 200, phi = c(1.6, -0.9))
#'
#' # ARMA(1,1) with mean 10
#' x <- gen_arima(n = 200, phi = 0.7, theta = 0.4, mu = 10)
#'
#' # ARIMA(1,1,1)
#' x <- gen_arima(n = 200, phi = 0.7, theta = 0.4, d = 1)
#'
#' # Seasonal (period 12)
#' x <- gen_arima(n = 200, phi = 0.7, s = 12)
gen_arima <- function(n,
                      phi = 0,
                      theta = 0,
                      d = 0L,
                      s = 0L,
                      mu = 0,
                      vara = 1,
                      plot = TRUE,
                      seed = NULL,
                      spin = 2000L) {

  # Input validation
  if (!is.numeric(n) || length(n) != 1 || n < 1) {
    stop("`n` must be a positive integer")
  }
  if (!is.numeric(phi)) stop("`phi` must be numeric")
  if (!is.numeric(theta)) stop("`theta` must be numeric")
  if (!is.numeric(d) || d < 0) stop("`d` must be a non-negative integer")
  if (!is.numeric(s) || s < 0) stop("`s` must be a non-negative integer")
  if (!is.numeric(vara) || vara <= 0) stop("`vara` must be positive")

  if (!is.null(seed)) set.seed(seed)

  # Determine effective orders
  p <- if (all(phi == 0)) 0L else length(phi)
  q <- if (all(theta == 0)) 0L else length(theta)

  # Build model list for arima.sim
  # Note: arima.sim uses opposite sign convention for MA
  model <- list(order = c(p, d, q))
  if (p > 0) model$ar <- phi
  if (q > 0) model$ma <- -theta

  # Generate with spin-up
  ngen <- n + spin + s
  sd <- sqrt(vara)
  tsdata <- arima.sim(n = ngen, model = model, sd = sd)
  y <- as.double(tsdata)

  # Apply seasonal integration if needed
  if (s > 0) {
    x <- .seasonal_integrate(y, s, n, d, spin)
  } else {
    start_idx <- spin + d + 1
    x <- y[start_idx:(start_idx + n - 1)]
  }

  # Add mean
  x <- x + mu

  if (plot) {
    .plot_realization(x, n)
    invisible(x)
  } else {
    x
  }
}


#' Apply seasonal integration
#' @keywords internal
#' @noRd
.seasonal_integrate <- function(y, s, n, d, spin) {

  ndspin <- n + d + s + spin
  xfull <- numeric(ndspin)
  d1 <- d + s + 1

  for (i in d1:ndspin) {
    xfull[i] <- y[i] + xfull[i - s]
  }

  # Extract after spin-up
  start_idx <- spin + d1
  xfull[start_idx:(start_idx + n - 1)]
}


#' Plot realization
#' @keywords internal
#' @noRd
.plot_realization <- function(x, n) {

  t <- seq_len(n)
  type <- if (n < 200) "o" else "l"
  pch <- if (n < 200) 16 else NA

  plot(t, x, type = type, pch = pch, cex = 0.5,
       xlab = "Time", ylab = "", main = "Realization")
}
