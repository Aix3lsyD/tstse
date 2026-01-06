#' Rolling Window RMSE for ARIMA Forecasts
#'
#' Evaluate forecast accuracy of an ARIMA/ARMA model using rolling window
#' cross-validation. Computes RMSE across multiple overlapping windows.
#'
#' @param x Numeric vector, the time series.
#' @param horizon Integer, forecast horizon (default 1).
#' @param phi Numeric vector, AR coefficients (default 0 = none).
#' @param theta Numeric vector, MA coefficients (default 0 = none).
#' @param d Integer, order of differencing (default 0).
#' @param s Integer, seasonal period (default 0 = none).
#' @param cores Integer, number of cores for parallel processing.
#'   Default NULL uses `getOption("tstse.cores", 1)`.
#'   Set to 0 to use all available cores.
#' @param plot Logical, whether to plot histogram of RMSEs (default TRUE).
#' @param verbose Logical, whether to print progress messages (default FALSE).
#'
#' @return A list with components:
#'   \item{rmse}{Mean RMSE across all windows}
#'   \item{rmse_all}{Vector of RMSE for each window}
#'   \item{n_windows}{Number of rolling windows used}
#'   \item{training_size}{Size of each training window}
#'   \item{horizon}{Forecast horizon}
#'   \item{phi}{AR coefficients used}
#'   \item{theta}{MA coefficients used}
#'   \item{d}{Differencing order}
#'   \item{s}{Seasonal period}
#'
#' @details
#' The function performs rolling window cross-validation:
#' \enumerate{
#'   \item Determines training window size based on model order
#'   \item Slides the window across the series
#'   \item For each position, fits forecast and computes error vs actual
#'   \item Returns average RMSE across all windows
#' }
#'
#' Training size is computed as:
#' \itemize{
#'   \item For ARMA (d=0, s=0): `max(p, q) + 1`
#'   \item For ARIMA/SARIMA: `p + q + d + s + 1`
#' }
#'
#' With parallelization enabled, each window is processed independently
#' across multiple cores for significant speedup on large series.
#'
#' @seealso [fore_arima()], [fore_arma()], [roll_win_rmse_nn()]
#'
#' @examples
#' \donttest{
#' # Generate AR(2) data
#' set.seed(123)
#' x <- gen_arma(n = 150, phi = c(1.5, -0.75), plot = FALSE)
#'
#' # Evaluate AR(2) model with 5-step ahead forecasts
#' result <- roll_win_rmse(x, horizon = 5, phi = c(1.5, -0.75))
#' cat("Mean RMSE:", result$rmse, "\n")
#' cat("Number of windows:", result$n_windows, "\n")
#' }
#'
#' \dontrun{
#' # With parallel processing (if multiple cores available)
#' result <- roll_win_rmse(x, horizon = 5, phi = c(1.5, -0.75), cores = 2)
#' }
#'
#' @export
roll_win_rmse <- function(x, horizon = 1L, phi = 0, theta = 0,
                          d = 0L, s = 0L,
                          cores = 1L, plot = TRUE, verbose = FALSE) {

  # Input validation
  if (!is.numeric(x) || length(x) == 0) {
    stop("`x` must be a non-empty numeric vector", call. = FALSE)
  }
  if (!is.numeric(horizon) || length(horizon) != 1 || horizon < 1) {
    stop("`horizon` must be a positive integer", call. = FALSE)
  }
  if (!is.numeric(phi)) {
    stop("`phi` must be numeric", call. = FALSE)
  }
  if (!is.numeric(theta)) {
    stop("`theta` must be numeric", call. = FALSE)
  }
  if (!is.numeric(d) || length(d) != 1 || d < 0) {
    stop("`d` must be a non-negative integer", call. = FALSE)
  }
  if (!is.numeric(s) || length(s) != 1 || s < 0) {
    stop("`s` must be a non-negative integer", call. = FALSE)
  }

  horizon <- as.integer(horizon)
  d <- as.integer(d)
  s <- as.integer(s)
  n <- length(x)
  cores <- get_cores(cores)

  # Determine effective AR/MA orders
  p <- length(phi)
  if (sum(phi^2) == 0) p <- 0L
  q <- length(theta)
  if (sum(theta^2) == 0) q <- 0L

  # Compute training size
  if (s == 0L && d == 0L) {
    # ARMA case
    training_size <- max(p, q) + 1L
  } else {
    # ARIMA/SARIMA case
    training_size <- p + q + d + s + 1L
  }

  # Ensure minimum training size (need at least 2 for forecasting)
  training_size <- max(training_size, 2L)

  # Compute number of windows
  n_windows <- n - training_size - horizon + 1L

  if (n_windows < 1L) {
    stop("Series too short for rolling window evaluation. ",
         "Need at least ", training_size + horizon, " observations.",
         call. = FALSE)
  }

  if (verbose) {
    message("Rolling Window RMSE Evaluation")
    message("  Training size: ", training_size)
    message("  Forecast horizon: ", horizon)
    message("  Number of windows: ", n_windows)
    if (cores > 1L) message("  Using ", cores, " cores")
  }

  # Store series mean for forecast (used in original tswge)
  xbar <- mean(x)

  # Define window evaluation function
  eval_window <- function(i) {
    # Extract training window
    train_start <- i
    train_end <- i + training_size - 1L
    train_data <- x[train_start:train_end]

    # Forecast
    if (s == 0L && d == 0L) {
      # Use fore_arma for pure ARMA
      fc <- fore_arma(train_data, phi = phi, theta = theta,
                      n_ahead = horizon, plot = FALSE)
    } else {
      # Use fore_arima for differencing/seasonality
      fc <- fore_arima(train_data, phi = phi, theta = theta,
                       d = d, s = s, n_ahead = horizon, plot = FALSE)
    }

    # Get actual values
    actual_start <- train_end + 1L
    actual_end <- train_end + horizon
    actual <- x[actual_start:actual_end]

    # Compute RMSE for this window
    sqrt(mean((actual - fc$f)^2))
  }

  # Run evaluation (parallel or sequential)
  if (cores > 1L && n_windows > 1L) {
    rmse_list <- pmap(seq_len(n_windows), eval_window, cores = cores)
    rmse_all <- unlist(rmse_list)
  } else {
    rmse_all <- vapply(seq_len(n_windows), eval_window, numeric(1))
  }

  # Compute mean RMSE
  rmse <- mean(rmse_all)

  # Plot histogram if requested
  if (plot && n_windows > 1L) {
    op <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(op))

    graphics::hist(rmse_all, main = "RMSE Distribution Across Windows",
                   xlab = "RMSE", col = "lightblue", border = "white")
    graphics::abline(v = rmse, col = "red", lwd = 2, lty = 2)
  }

  if (verbose) {
    message("\nResults:")
    message("  Mean RMSE: ", round(rmse, 4))
    message("  RMSE range: [", round(min(rmse_all), 4), ", ",
            round(max(rmse_all), 4), "]")
  }

  invisible(list(
    rmse          = rmse,
    rmse_all      = rmse_all,
    n_windows     = n_windows,
    training_size = training_size,
    horizon       = horizon,
    phi           = phi,
    theta         = theta,
    d             = d,
    s             = s
  ))
}
