#' Rolling Window RMSE for Neural Network Forecasts
#'
#' Evaluate forecast accuracy of an MLP model using rolling window cross-validation.
#'
#' @param x Numeric vector or ts object, the time series.
#' @param horizon Integer, forecast horizon (default 1).
#' @param model An mlp model object from the nnfor package.
#' @param plot Logical, whether to plot histogram of RMSEs (default TRUE).
#' @param verbose Logical, whether to print progress messages (default FALSE).
#'
#' @return A list with components:
#'   \item{rmse}{Mean RMSE across all windows}
#'   \item{mse}{Mean MSE across all windows}
#'   \item{rmse_all}{Vector of RMSE for each window}
#'   \item{mse_all}{Vector of MSE for each window}
#'   \item{n_windows}{Number of rolling windows used}
#'   \item{horizon}{Forecast horizon}
#'   \item{frequency}{Series frequency}
#'   \item{diff_order}{Differencing order from model}
#'
#' @details
#' For non-seasonal models, performs true rolling window evaluation across
#' multiple windows. For seasonal models (with seasonal dummies), only
#' evaluates on the last window due to MLP recalculating seasonal components.
#'
#' The training window size is determined by the model's lag structure and
#' differencing order, with a minimum of 9 observations.
#'
#' @examples
#' \dontrun{
#' # Fit an MLP model
#' library(nnfor)
#' fit <- mlp(AirPassengers)
#'
#' # Evaluate with rolling window
#' result <- roll_win_rmse_nn(AirPassengers, horizon = 12, model = fit)
#' print(result$rmse)
#' }
#'
#' @export
roll_win_rmse_nn <- function(x, horizon = 1L, model, plot = TRUE, verbose = FALSE) {

  # Check required packages
  if (!requireNamespace("nnfor", quietly = TRUE)) {
    stop("Package 'nnfor' is required.\n",
         "Install with: install.packages('nnfor')")
  }
  if (!requireNamespace("forecast", quietly = TRUE)) {
    stop("Package 'forecast' is required.\n",
         "Install with: install.packages('forecast')")
  }

  # Input validation
  if (!is.numeric(x) || length(x) == 0) {
    stop("`x` must be a non-empty numeric vector or ts object")
  }
  if (!inherits(model, "mlp")) {
    stop("`model` must be an mlp object from the nnfor package")
  }
  if (!is.numeric(horizon) || length(horizon) != 1 || horizon < 1) {
    stop("`horizon` must be a positive integer")
  }

  horizon <- as.integer(horizon)
  n <- length(x)
  s <- stats::frequency(x)

  # Extract differencing order
  d <- if (is.null(model$difforder)) 0L else max(model$difforder)

  # Branch based on seasonal dummies
  if (isFALSE(model$sdummy)) {
    # Non-seasonal: true rolling window
    lags <- model$lags
    training_size <- max(max(lags) + sum(d) + 2L, 9L)
    n_windows <- n - (training_size + horizon) + 1L

    if (n_windows < 1L) {
      stop("Series too short for rolling window evaluation")
    }

    if (verbose) {
      message("Training size: ", training_size)
      message("Number of windows: ", n_windows)
      message("Computing rolling window RMSE...")
    }

    mse_all <- numeric(n_windows)

    for (i in seq_len(n_windows)) {
      # Fit on window
      window_data <- stats::ts(x[i:(i + training_size - 1L)], frequency = s)
      fit <- nnfor::mlp(window_data, model = model)
      fit$minmax <- model$minmax
      fit$difforder <- model$difforder

      # Forecast
      forecasts <- forecast::forecast(fit, h = horizon)$mean

      # Calculate MSE
      actual <- x[(training_size + i):(training_size + i + horizon - 1L)]
      mse_all[i] <- mean((actual - forecasts)^2)
    }

  } else {
    # Seasonal: single window only
    training_size <- 20L
    n_windows <- 1L

    if (verbose) {
      message("Seasonal model: evaluating on last window only")
      message("Training size: ", training_size)
    }

    if (n < training_size + horizon) {
      stop("Series too short for evaluation")
    }

    mse_all <- numeric(1L)

    window_data <- stats::ts(x[1:(n - horizon)], frequency = s)
    fit <- nnfor::mlp(window_data, model = model)
    fit$minmax <- model$minmax
    fit$difforder <- model$difforder

    forecasts <- forecast::forecast(fit, h = horizon)$mean
    actual <- x[(n - horizon + 1L):n]
    mse_all[1L] <- mean((actual - forecasts)^2)
  }

  # Compute RMSE
  rmse_all <- sqrt(mse_all)
  rmse <- mean(rmse_all)
  mse <- mean(mse_all)

  # Plot if requested
  if (plot && n_windows > 1L) {
    op <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(op))

    graphics::hist(rmse_all, main = "RMSE Distribution Across Windows",
                   xlab = "RMSE", col = "lightblue", border = "white")
    graphics::abline(v = rmse, col = "red", lwd = 2, lty = 2)
  }

  if (verbose) {
    message("\nRolling Window Results:")
    message("  Mean RMSE: ", round(rmse, 4))
    message("  Mean MSE:  ", round(mse, 4))
  }

  invisible(list(
    rmse = rmse,
    mse = mse,
    rmse_all = rmse_all,
    mse_all = mse_all,
    n_windows = n_windows,
    horizon = horizon,
    frequency = s,
    diff_order = d
  ))
}
