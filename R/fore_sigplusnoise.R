#' Forecast Signal Plus Noise Model
#'
#' Forecasts a time series modeled as a deterministic signal (linear trend or
#' cosine function) plus AR noise.
#'
#' @param x Numeric vector. The time series to forecast.
#' @param linear Logical. If `TRUE` (default), fit a linear trend as the signal.
#'   If `FALSE`, fit a cosine function at the specified frequency.
#' @param freq Numeric. Frequency for the cosine signal when `linear = FALSE`.
#'   Ignored when `linear = TRUE`. Default is 0.
#' @param max_p Integer. Maximum AR order to consider for the noise component.
#'   Order is selected by AIC. Default is 5.
#' @param n_ahead Integer. Number of periods to forecast ahead. Default is 10.
#' @param lastn Logical. If `TRUE`, forecast the last `n_ahead` observations
#'   (for model validation). If `FALSE` (default), forecast beyond the end of
#'   the series.
#' @param method Character. Estimation method for AR model: `"mle"` (default),
#'   `"burg"`, or `"yw"` (Yule-Walker).
#' @param alpha Numeric. Significance level for prediction intervals.
#'   Default is 0.05 (95% intervals).
#' @param limits Logical. If `TRUE` (default), compute and plot prediction
#'   intervals.
#' @param plot Logical. If `TRUE` (default), plot the series and forecasts.
#'
#' @return A list with class `"fore_sigplusnoise"` containing:
#'   \item{signal}{The fitted signal component (linear trend or cosine).}
#'   \item{noise}{The noise component (residuals from signal).}
#'   \item{phi}{AR coefficients for the noise model.}
#'   \item{p}{AR order selected.}
#'   \item{f}{Point forecasts.}
#'   \item{ll}{Lower prediction limits.}
#'   \item{ul}{Upper prediction limits.}
#'   \item{resid}{Residuals from the AR model.}
#'   \item{se}{Standard errors for forecasts.}
#'   \item{wnv}{White noise variance estimate.}
#'   \item{signal_coef}{Signal coefficients: for linear, `c(intercept, slope)`;
#'     for cosine, `c(intercept, cos_coef, sin_coef)`.}
#'
#' @details
#' The model assumes:
#' \deqn{X_t = S_t + Z_t}
#' where \eqn{S_t} is a deterministic signal and \eqn{Z_t} is AR(p) noise.
#'
#' For `linear = TRUE`:
#' \deqn{S_t = \beta_0 + \beta_1 t}
#'
#' For `linear = FALSE`:
#' \deqn{S_t = \beta_0 + \beta_1 \cos(2\pi f t) + \beta_2 \sin(2\pi f t)}
#'
#' The noise AR order is selected by AIC from 0 to `max_p`.
#'
#' @seealso [gen_sigplusnoise()] for generating signal plus noise data,
#'   [fore_arma()] for ARMA forecasting.
#'
#' @examples
#' # Linear trend plus AR noise
#' set.seed(123)
#' t <- 1:100
#' signal <- 5 + 0.1 * t
#' noise <- arima.sim(model = list(ar = 0.7), n = 100)
#' x <- signal + noise
#'
#' result <- fore_sigplusnoise(x, linear = TRUE, n_ahead = 10)
#'
#' # Cosine signal plus noise
#' signal2 <- 10 + 3 * cos(2 * pi * 0.1 * t)
#' x2 <- signal2 + noise
#'
#' result2 <- fore_sigplusnoise(x2, linear = FALSE, freq = 0.1, n_ahead = 10)
#'
#' @export
fore_sigplusnoise <- function(x, linear = TRUE, freq = 0, max_p = 5L,
                              n_ahead = 10L, lastn = FALSE, method = "mle",
                              alpha = 0.05, limits = TRUE, plot = TRUE) {
  # Input validation
  if (!is.numeric(x)) {
    stop("`x` must be a numeric vector.", call. = FALSE)
  }
  n <- length(x)
  if (n < 5L) {
    stop("`x` must have at least 5 observations.", call. = FALSE)
  }
  if (anyNA(x)) {
    stop("`x` contains NA values.", call. = FALSE)
  }

  if (!is.logical(linear) || length(linear) != 1L || is.na(linear)) {
    stop("`linear` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.numeric(freq) || length(freq) != 1L) {
    stop("`freq` must be a single numeric value.", call. = FALSE)
  }

  if (!is.numeric(max_p) || length(max_p) != 1L || max_p < 0L) {
    stop("`max_p` must be a non-negative integer.", call. = FALSE)
  }
  max_p <- as.integer(max_p)

  if (!is.numeric(n_ahead) || length(n_ahead) != 1L || n_ahead < 1L) {
    stop("`n_ahead` must be a positive integer.", call. = FALSE)
  }
  n_ahead <- as.integer(n_ahead)

  if (!is.logical(lastn) || length(lastn) != 1L || is.na(lastn)) {
    stop("`lastn` must be TRUE or FALSE.", call. = FALSE)
  }

  method <- match.arg(method, c("mle", "burg", "yw"))

  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be between 0 and 1.", call. = FALSE)
  }

  if (!is.logical(limits) || length(limits) != 1L || is.na(limits)) {
    stop("`limits` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.logical(plot) || length(plot) != 1L || is.na(plot)) {
    stop("`plot` must be TRUE or FALSE.", call. = FALSE)
  }

  # Setup
  t_obs <- seq_len(n)
  n_total <- n + n_ahead

  # Fit signal component
  if (linear) {
    # Linear trend: S_t = b0 + b1 * t
    fit <- lm(x ~ t_obs)
    b0 <- unname(coef(fit)[1])
    b1 <- unname(coef(fit)[2])
    signal_coef <- c(intercept = b0, slope = b1)

    # Signal values for observed data
    signal_obs <- b0 + b1 * t_obs

    # Signal values for forecast period
    if (!lastn) {
      t_fore <- (n + 1):n_total
      signal_fore <- b0 + b1 * t_fore
    } else {
      signal_fore <- signal_obs[(n - n_ahead + 1):n]
    }
  } else {
    # Cosine signal: S_t = b0 + b1*cos(2*pi*f*t) + b2*sin(2*pi*f*t)
    cos_t <- cos(2 * pi * freq * t_obs)
    sin_t <- sin(2 * pi * freq * t_obs)
    fit <- lm(x ~ cos_t + sin_t)
    b0 <- unname(coef(fit)[1])
    b1 <- unname(coef(fit)[2])
    b2 <- unname(coef(fit)[3])
    signal_coef <- c(intercept = b0, cos_coef = b1, sin_coef = b2)

    # Signal values for observed data
    signal_obs <- b0 + b1 * cos_t + b2 * sin_t

    # Signal values for forecast period
    if (!lastn) {
      t_fore <- (n + 1):n_total
      cos_fore <- cos(2 * pi * freq * t_fore)
      sin_fore <- sin(2 * pi * freq * t_fore)
      signal_fore <- b0 + b1 * cos_fore + b2 * sin_fore
    } else {
      signal_fore <- signal_obs[(n - n_ahead + 1):n]
    }
  }

  # Extract noise component
  noise <- x - signal_obs

  # Select AR order for noise using AIC
  aic_result <- aic(noise, p = 0:max_p, q = 0)
  p <- aic_result$p

  # Fit AR model to noise
  phi <- numeric(0)
  if (p > 0L) {
    ar_fit <- est_ar(noise, p = p, method = method)
    phi <- ar_fit$phi
  }

  # Compute constant for mean adjustment
  const <- 1
  if (p > 0L) {
    const <- 1 - sum(phi)
  }
  noise_mean <- mean(noise)
  ma_const <- const * noise_mean

  # Get residuals using backcast
  resid <- backcast(noise, phi = phi, theta = 0, n_back = 50)

  # Estimate white noise variance
  wnv <- mean(resid^2)

  # Forecast noise component
  mm <- if (lastn) n - n_ahead else n

  noise_hat <- c(noise, rep(0, n_ahead))
  for (h in seq_len(n_ahead)) {
    if (p > 0L) {
      for (jp in seq_len(p)) {
        noise_hat[mm + h] <- noise_hat[mm + h] + phi[jp] * noise_hat[mm + h - jp]
      }
    }
    noise_hat[mm + h] <- noise_hat[mm + h] + ma_const
  }

  # Compute psi weights for prediction intervals
  psi <- psi_weights(phi = phi, theta = 0, lag_max = n_ahead)

  # Compute standard errors
  psi_sq_cumsum <- cumsum(c(1, psi[-length(psi)]^2))
  se <- sqrt(wnv * psi_sq_cumsum[seq_len(n_ahead)])

  # Combine signal and noise forecasts
  noise_fore <- noise_hat[(mm + 1):(mm + n_ahead)]
  f <- signal_fore + noise_fore

  # Prediction intervals
  z_alpha <- qnorm(1 - alpha / 2)
  ll <- f - z_alpha * se
  ul <- f + z_alpha * se

  # Plot if requested
  if (plot) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))

    par(mar = c(4, 4, 2, 1))

    # Determine plot limits
    max_h <- mm + n_ahead
    if (limits) {
      ylim <- range(c(x, ll, ul))
    } else {
      ylim <- range(c(x, f))
    }

    # Plot original series
    plot(t_obs, x, type = "o", pch = 16, cex = 0.8,
         xlim = c(1, max_h), ylim = ylim,
         xlab = "Time", ylab = "",
         main = paste("Signal + Noise Forecast"))

    # Plot forecasts
    t_fore_plot <- mm:(mm + n_ahead)
    f_plot <- c(x[mm], f)

    lines(t_fore_plot, f_plot, type = "o", pch = 18, lty = 3, lwd = 2, cex = 1)

    # Plot prediction limits
    if (limits) {
      ll_plot <- c(x[mm], ll)
      ul_plot <- c(x[mm], ul)
      lines(t_fore_plot, ll_plot, type = "o", pch = 18, lty = 3, cex = 0.6)
      lines(t_fore_plot, ul_plot, type = "o", pch = 18, lty = 3, cex = 0.6)
    }
  }

  # Build result
  result <- list(
    signal = signal_obs,
    noise = noise,
    phi = phi,
    p = p,
    f = f,
    ll = ll,
    ul = ul,
    resid = resid,
    se = se,
    wnv = wnv,
    signal_coef = signal_coef,
    psi = psi
  )

  class(result) <- "fore_sigplusnoise"
  invisible(result)
}


#' @export
print.fore_sigplusnoise <- function(x, ...) {
  cat("Signal Plus Noise Forecast\n")
  cat("--------------------------\n")

  if (length(x$signal_coef) == 2) {
    cat("Signal: Linear trend\n")
    cat("  Intercept:", round(x$signal_coef[1], 4), "\n")
    cat("  Slope:", round(x$signal_coef[2], 4), "\n")
  } else {
    cat("Signal: Cosine function\n")
    cat("  Intercept:", round(x$signal_coef[1], 4), "\n")
    cat("  Cos coef:", round(x$signal_coef[2], 4), "\n")
    cat("  Sin coef:", round(x$signal_coef[3], 4), "\n")
  }

  cat("\nNoise: AR(", x$p, ")\n", sep = "")
  if (x$p > 0) {
    cat("  Phi:", round(x$phi, 4), "\n")
  }

  cat("\nForecasts:", length(x$f), "periods ahead\n")
  cat("White noise variance:", round(x$wnv, 4), "\n")

  invisible(x)
}
