#' ARUMA Forecasting
#'
#' Generate point forecasts and prediction intervals for ARUMA models.
#' ARUMA extends ARIMA with an additional lambda factor for extra AR-type roots.
#'
#' @param x Numeric vector, the time series to forecast.
#' @param phi Numeric vector, AR coefficients (default 0 = none).
#' @param theta Numeric vector, MA coefficients (default 0 = none).
#' @param d Integer, order of differencing (default 0, max 3).
#' @param s Integer, seasonal period (default 0 = none).
#' @param lambda Numeric vector, additional AR-type factor coefficients (default 0 = none).
#' @param n_ahead Integer, forecast horizon (default 5).
#' @param lastn Logical, if TRUE forecast from position n-n_ahead (holdout mode).
#' @param plot Logical, whether to plot forecasts (default TRUE).
#' @param alpha Numeric, significance level for prediction intervals (default 0.05).
#' @param limits Logical, whether to compute/show prediction limits (default TRUE).
#'
#' @return An object of class "fore_aruma" with components:
#'   \item{f}{Point forecasts (length n_ahead)}
#'   \item{ll}{Lower prediction limits}
#'   \item{ul}{Upper prediction limits}
#'   \item{resid}{Residuals (length n)}
#'   \item{wnv}{Estimated white noise variance}
#'   \item{xbar}{Sample mean of original series}
#'   \item{se}{Standard errors for each forecast horizon}
#'   \item{psi}{Psi weights used for prediction intervals}
#' @export
#'
#' @details
#' ARUMA models extend ARIMA by adding an additional AR-type factor (lambda).
#' The total AR polynomial is: phi(B) * (1-B)^d * (1-B^s) * lambda(B).
#'
#' Parameters:
#' \itemize{
#'   \item \code{phi}: AR coefficients for the stationary component
#'   \item \code{theta}: MA coefficients
#'   \item \code{d}: Order of differencing (supports 0, 1, 2, or 3)
#'   \item \code{s}: Seasonal period (e.g., 12 for monthly data)
#'   \item \code{lambda}: Additional AR-type factor (e.g., for extra unit roots)
#' }
#'
#' @seealso \code{\link{fore_arima}}, \code{\link{fore_arma}}
#'
#' @examples
#' # Generate data
#' set.seed(123)
#' x <- gen_arma(n = 100, phi = c(1.5, -0.75), plot = FALSE)
#'
#' # ARUMA forecast with lambda factor
#' fore_aruma(x, phi = 0.8, d = 1, lambda = c(0.5), n_ahead = 5)
fore_aruma <- function(x, phi = 0, theta = 0, d = 0L, s = 0L, lambda = 0,
                       n_ahead = 5L, lastn = FALSE, plot = TRUE,
                       alpha = 0.05, limits = TRUE) {

  # Input validation
  if (!is.numeric(x) || length(x) == 0) {
    stop("`x` must be a non-empty numeric vector")
  }
  if (!is.numeric(phi)) stop("`phi` must be numeric")
  if (!is.numeric(theta)) stop("`theta` must be numeric")
  if (!is.numeric(d) || length(d) != 1 || d < 0 || d > 3) {
    stop("`d` must be 0, 1, 2, or 3")
  }
  if (!is.numeric(s) || length(s) != 1 || s < 0) {
    stop("`s` must be a non-negative integer")
  }
  if (!is.numeric(lambda)) stop("`lambda` must be numeric")
  if (!is.numeric(n_ahead) || length(n_ahead) != 1 || n_ahead < 1) {
    stop("`n_ahead` must be a positive integer")
  }
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be between 0 and 1")
  }

  x <- as.double(x)
  n <- length(x)
  xbar <- mean(x)
  d <- as.integer(d)
  s <- as.integer(s)
  n_ahead <- as.integer(n_ahead)

  # Determine effective orders
  p <- if (all(phi == 0)) 0L else length(phi)
  q <- if (all(theta == 0)) 0L else length(theta)
  dlam <- if (sum(lambda^2) > 0) length(lambda) else 0L

  # Build differencing factor: (1 - B)^d
  diff_fac <- switch(as.character(d),
    "0" = 0,
    "1" = 1,
    "2" = c(2, -1),
    "3" = c(3, -3, 1)
  )

  # Build seasonal factor: (1 - B^s)
  if (s > 0) {
    seas_fac <- rep(0, s)
    seas_fac[s] <- 1
  } else {
    seas_fac <- 0
  }

  # Build lambda factor
  lambda_fac <- if (sum(lambda^2) > 0) lambda else 0

  # Combine factors for forecasting: phi * diff * seasonal * lambda
  fac1 <- if (p > 0) phi else 0
  prod_fore <- mult(fac1, diff_fac, seas_fac, lambda_fac)
  phitot_fore <- prod_fore$model_coef
  ptot_fore <- p + d + s + dlam

  # Combine factors for residuals: diff * seasonal * lambda (no AR)
  prod_res <- mult(diff_fac, seas_fac, lambda_fac)
  phitot_res <- prod_res$model_coef
  ptot_res <- d + s + dlam

  # Transform data to stationary ARMA
  if (ptot_res > 0 && !all(phitot_res == 0)) {
    y_arma <- artrans(x, phi = phitot_res, plot = FALSE)
  } else {
    y_arma <- x
  }

  # Backcast residuals on transformed data
  res_arma <- backcast(y_arma, phi = phi, theta = theta, n_back = 50L)

  # Pad residuals to full length (first ptot_res positions are zeros)
  resid <- rep(0, n)
  n_arma <- length(res_arma)
  resid[(ptot_res + 1):(ptot_res + n_arma)] <- res_arma

  # Compute constant for mean adjustment
  const <- 1 - sum(phitot_fore)

  # Forecast origin
  m <- if (lastn) n - n_ahead else n

  # Initialize forecast array with observed data
  xhat <- c(x[1:m], rep(0, n_ahead))

  # Box-Jenkins forecast recursion
  for (h in seq_len(n_ahead)) {
    idx <- m + h

    # AR contribution from total AR polynomial
    if (ptot_fore > 0 && length(phitot_fore) > 0 && !all(phitot_fore == 0)) {
      for (j in seq_len(length(phitot_fore))) {
        xhat[idx] <- xhat[idx] + phitot_fore[j] * xhat[idx - j]
      }
    }

    # MA contribution (only for horizons within MA order)
    if (q > 0 && h <= q) {
      for (j in h:q) {
        xhat[idx] <- xhat[idx] - theta[j] * resid[m + h - j]
      }
    }

    # Add mean contribution
    xhat[idx] <- xhat[idx] + mean(x[1:m]) * const
  }

  # Extract point forecasts
  f <- xhat[(m + 1):(m + n_ahead)]

  # Compute psi weights for prediction intervals
  if (ptot_fore > 0 && length(phitot_fore) > 0 && !all(phitot_fore == 0)) {
    psi <- psi_weights(phitot_fore, theta, lag_max = n_ahead)
  } else if (q > 0) {
    psi <- psi_weights(0, theta, lag_max = n_ahead)
  } else {
    psi <- rep(0, n_ahead)
  }

  # Estimate white noise variance from residuals
  p1 <- max(ptot_res + 1, q + 1)
  if (p1 <= n) {
    wnv <- sum(resid[p1:n]^2) / (n - ptot_res)
  } else {
    wnv <- var(resid[resid != 0])
  }

  # Compute cumulative variance for each horizon
  # Var(e_h) = sigma^2 * (1 + psi_1^2 + ... + psi_{h-1}^2)
  if (n_ahead == 1) {
    xisq <- 1
  } else {
    xisq <- cumsum(c(1, psi[1:(n_ahead - 1)]^2))
  }
  se <- sqrt(wnv * xisq)

  # Prediction limits
  z <- qnorm(1 - alpha / 2)
  ll <- f - z * se
  ul <- f + z * se

  # Plot if requested
  if (plot) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))

    maxh <- m + n_ahead
    t <- seq_len(n)

    # Determine y-axis limits
    if (limits) {
      minp <- min(x, ll)
      maxp <- max(x, ul)
    } else {
      if (lastn) {
        minp <- min(x, xhat[1:n])
        maxp <- max(x, xhat[1:n])
      } else {
        minp <- min(x, xhat)
        maxp <- max(x, xhat)
      }
    }

    # Plot original series
    ptype <- if (n >= 200) "l" else "o"
    plot(t, x, type = ptype, pch = 16, cex = 0.8,
         xlim = c(1, maxh), ylim = c(minp, maxp),
         xlab = "Time", ylab = "", main = "ARUMA Forecast")

    # Plot forecasts
    tf <- m:(m + n_ahead)
    fplot <- c(x[m], f)
    points(tf, fplot, type = "o", col = 2, pch = 1, cex = 0.6, lwd = 1)

    # Plot prediction limits
    if (limits) {
      llplot <- c(x[m], ll)
      ulplot <- c(x[m], ul)
      lines(tf, ulplot, lty = 3, col = "blue3", lwd = 2)
      lines(tf, llplot, lty = 3, col = "blue3", lwd = 2)
    }
  }

  # Build result
  result <- list(
    f = f,
    ll = ll,
    ul = ul,
    resid = resid,
    wnv = wnv,
    xbar = xbar,
    se = se,
    psi = psi
  )
  class(result) <- "fore_aruma"

  invisible(result)
}


#' @rdname fore_aruma
#' @param x An object of class "fore_aruma".
#' @param ... Additional arguments (ignored).
#' @export
print.fore_aruma <- function(x, ...) {

  cat("ARUMA Forecast\n")
  cat("Horizon:", length(x$f), "\n\n")

  df <- data.frame(
    h = seq_along(x$f),
    Forecast = round(x$f, 4),
    Lower = round(x$ll, 4),
    Upper = round(x$ul, 4),
    SE = round(x$se, 4)
  )
  print(df, row.names = FALSE)

  cat("\nWhite noise variance:", round(x$wnv, 6), "\n")
  cat("Series mean:", round(x$xbar, 4), "\n")

  invisible(x)
}
