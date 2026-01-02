#' Generate ARUMA Process
#'
#' Generate a realization from an ARUMA (ARIMA with general nonstationary factor) process.
#'
#' @param n Integer, length of output series.
#' @param phi Numeric vector, AR coefficients (default 0 = none).
#' @param theta Numeric vector, MA coefficients (default 0 = none).
#' @param d Integer, order of differencing (default 0).
#' @param s Integer, seasonal period for (1-B^s) factor (default 0 = none).
#' @param lambda Numeric vector, coefficients for additional nonstationary factor
#'   (default 0 = none). For factor (1 - lambda_1 B - lambda_2 B^2 - ...),
#'   pass c(lambda_1, lambda_2, ...).
#' @param vara Numeric, white noise variance (default 1).
#' @param plot Logical, whether to plot the series (default TRUE).
#' @param seed Integer, random seed for reproducibility (default NULL).
#'
#' @return Numeric vector of length n.
#' @export
#'
#' @details
#' Generates from the ARUMA model which extends ARIMA with a general
#' nonstationary operator. The model is:
#'
#' \deqn{\lambda(B)(1-B)^d(1-B^s)\phi(B) X_t = \theta(B) a_t}
#'
#' where:
#' \itemize{
#'   \item \eqn{\phi(B)} is the AR operator
#'   \item \eqn{\theta(B)} is the MA operator
#'   \item \eqn{(1-B)^d} is the differencing operator
#'   \item \eqn{(1-B^s)} is the seasonal differencing operator
#'   \item \eqn{\lambda(B)} is a general nonstationary factor
#' }
#'
#' The function first generates a stationary ARMA process, then applies
#' the inverse of all nonstationary factors.
#'
#' @seealso [gen_arima()] for ARIMA without the lambda factor,
#'   [gen_arma()] for stationary ARMA,
#'   [mult()] for multiplying polynomial factors.
#'
#' @examples
#' # ARIMA(1,1,0) - same as gen_arima with d=1
#' x <- gen_aruma(n = 200, phi = 0.7, d = 1)
#'
#' # Seasonal ARIMA with period 12
#' x <- gen_aruma(n = 200, phi = 0.7, s = 12)
#'
#' # ARUMA with custom nonstationary factor (1 - 0.9B)
#' x <- gen_aruma(n = 200, phi = 0.5, lambda = 0.9)
#'
#' # Combined: AR(1) with differencing and seasonal
#' x <- gen_aruma(n = 200, phi = 0.5, d = 1, s = 12, seed = 42, plot = FALSE)
gen_aruma <- function(n,
                      phi = 0,
                      theta = 0,
                      d = 0L,
                      s = 0L,
                      lambda = 0,
                      vara = 1,
                      plot = TRUE,
                      seed = NULL) {

  # Input validation
  if (!is.numeric(n) || length(n) != 1 || n < 1) {
    stop("`n` must be a positive integer")
  }
  if (!is.numeric(phi)) stop("`phi` must be numeric")
  if (!is.numeric(theta)) stop("`theta` must be numeric")
  if (!is.numeric(d) || d < 0) stop("`d` must be a non-negative integer")
  if (!is.numeric(s) || s < 0) stop("`s` must be a non-negative integer")
  if (!is.numeric(lambda)) stop("`lambda` must be numeric")
  if (!is.numeric(vara) || vara <= 0) stop("`vara` must be positive")

  if (!is.null(seed)) set.seed(seed)

  # Determine effective orders
  p <- if (all(phi == 0)) 0L else length(phi)
  q <- if (all(theta == 0)) 0L else length(theta)
  dlam <- if (all(lambda == 0)) 0L else length(lambda)

  # Build model list for arima.sim
  # Note: arima.sim uses opposite sign convention for MA
  model <- list(order = c(p, as.integer(d), q))
  if (p > 0) model$ar <- phi
  if (q > 0) model$ma <- -theta

  # Calculate total nonstationary order from lambda and seasonal

  dlams <- dlam + s
  spin <- 100L

  # Generate with spin-up
  ngen <- n + dlams + spin
  sd <- sqrt(vara)
  tsdata <- arima.sim(n = ngen, model = model, sd = sd)
  y <- as.double(tsdata)

  # Build combined nonstationary operator from lambda and seasonal
  lambdas <- .build_nonstationary_operator(lambda, s, dlam)

  # Apply inverse of nonstationary operator
  x <- .apply_inverse_operator(y, lambdas, dlams, d, n, spin)

  if (plot) {
    op <- par(mfrow = c(1, 1), mar = c(3.8, 2.5, 1, 1))
    on.exit(par(op), add = TRUE)
    .plot_aruma_realization(x, n)
    invisible(x)
  } else {
    x
  }
}


#' Build combined nonstationary operator
#' @noRd
.build_nonstationary_operator <- function(lambda, s, dlam) {
  if (dlam == 0 && s == 0) {
    return(0)
  }

  # Build seasonal factor: (1 - B^s) represented as c(0,0,...,1) with s elements
  seas <- NULL
  if (s > 0) {
    seas <- rep(0, s)
    seas[s] <- 1
  }

  # Combine lambda and seasonal using mult()
  if (dlam > 0 && s > 0) {
    result <- mult(lambda, seas)
    return(result$model_coef)
  } else if (dlam > 0) {
    return(lambda)
  } else {
    return(seas)
  }
}


#' Apply inverse of nonstationary operator (integration)
#' @noRd
.apply_inverse_operator <- function(y, lambdas, dlams, d, n, spin) {
  d1 <- d + dlams + 1L
  nd <- n + d + dlams + 1L
  ndspin <- nd + spin - 1L

  xfull <- numeric(ndspin)

  if (dlams == 0) {
    # No additional nonstationary factors, just extract
    for (i in d1:ndspin) {
      xfull[i] <- y[i]
    }
  } else {
    # Apply inverse filter: x[t] = y[t] + sum(lambdas[j] * x[t-j])
    for (i in d1:ndspin) {
      xfull[i] <- y[i]
      for (j in seq_len(dlams)) {
        xfull[i] <- xfull[i] + lambdas[j] * xfull[i - j]
      }
    }
  }

  # Extract after spin-up
  start_idx <- spin + d1
  xfull[start_idx:(start_idx + n - 1)]
}


#' Plot ARUMA realization
#' @noRd
.plot_aruma_realization <- function(x, n) {
  t <- seq_len(n)
  type <- if (n < 200) "o" else "l"
  pch <- if (n < 200) 16 else NA

  plot(t, x, type = type, pch = pch, cex = 0.5,
       xlab = "Time", ylab = "", main = "Realization")
}
