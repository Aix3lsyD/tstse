#' AIC/BIC Model Selection for ARMA Models
#'
#' Select optimal ARMA(p,q) orders using information criteria.
#'
#' @param x Numeric vector, the time series.
#' @param p Integer vector, AR orders to compare (default 0:5).
#' @param q Integer vector, MA orders to compare (default 0:2).
#' @param type Character, criterion to use: "aic" (default), "aicc", or "bic".
#' @param cores Integer, number of cores for parallel processing.
#'   Default NULL uses \code{getOption("tstse.cores", 1)}.
#'   Set to 0 to use all available cores.
#'
#' @return A list with components:
#'   \item{type}{Criterion used}
#'   \item{value}{Best criterion value}
#'   \item{p}{Selected AR order}
#'   \item{q}{Selected MA order}
#'   \item{phi}{AR coefficients at selected order}
#'   \item{theta}{MA coefficients at selected order}
#'   \item{xbar}{Sample mean}
#'   \item{vara}{Residual variance at selected order}
#'   \item{table}{Data frame of all (p,q) combinations and their criterion values}
#' @export
#'
#' @examples
#' \donttest{
#' x <- gen_arma(n = 200, phi = 0.7, theta = 0.4, plot = FALSE, seed = 123)
#'
#' # Sequential
#' aic_ts(x, p = 0:3, q = 0:2)
#' }
#'
#' \dontrun{
#' # Parallel (uses multiple cores)
#' aic_ts(x, p = 0:5, q = 0:3, cores = 2)
#' }
aic_ts <- function(x,
                p = 0:5,
                q = 0:2,
                type = c("aic", "aicc", "bic"),
                cores = 1L) {

  # Input validation
  if (!is.numeric(x) || length(x) == 0) {
    stop("`x` must be a non-empty numeric vector")
  }
  if (!is.numeric(p) || any(p < 0)) {
    stop("`p` must be a vector of non-negative integers")
  }
  if (!is.numeric(q) || any(q < 0)) {
    stop("`q` must be a vector of non-negative integers")
  }

  type <- match.arg(type)
  cores <- get_cores(cores)

  n <- length(x)
  xbar <- mean(x)
  x_centered <- x - xbar

  # Build grid of (p, q) combinations
  grid <- expand.grid(p = p, q = q)

  # Define worker function
  fit_arma_order <- function(idx) {
    j <- grid$p[idx]
    k <- grid$q[idx]

    # Try to fit ARMA model
    fit <- tryCatch(
      suppressWarnings(arima(x_centered, order = c(j, 0, k))),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      return(NULL)
    }

    # Extract coefficients
    coefs <- as.vector(coef(fit))

    if (j == 0) {
      phi <- 0
    } else {
      phi <- coefs[1:j]
    }

    if (k == 0) {
      theta <- 0
    } else {
      theta <- -coefs[(j + 1):(j + k)]
    }

    # Compute residuals and variance
    res <- backcast(x_centered, phi = phi, theta = theta, n_back = 50)
    avar <- mean(res^2)

    list(
      p = j,
      q = k,
      phi = phi,
      theta = theta,
      avar = avar,
      aic = log(avar) + 2 * (j + k + 1) / n,
      aicc = log(avar) + (n + j + k + 1) / (n - j - k - 3),
      bic = log(avar) + (j + k + 1) * log(n) / n
    )
  }

  # Fit models (parallel or sequential)
  results <- pmap(seq_len(nrow(grid)), fit_arma_order, cores = cores)

  # Remove failed fits
  results <- Filter(Negate(is.null), results)

  if (length(results) == 0) {
    stop("All ARMA models failed to converge")
  }

  # Build comparison table
  ic_table <- data.frame(
    p = sapply(results, `[[`, "p"),
    q = sapply(results, `[[`, "q"),
    aic = sapply(results, `[[`, "aic"),
    aicc = sapply(results, `[[`, "aicc"),
    bic = sapply(results, `[[`, "bic")
  )

  # Find best by selected criterion
  ic_values <- sapply(results, `[[`, type)
  best_idx <- which.min(ic_values)
  best <- results[[best_idx]]

  structure(
    list(
      type = type,
      value = best[[type]],
      p = best$p,
      q = best$q,
      phi = best$phi,
      theta = best$theta,
      xbar = xbar,
      vara = best$avar,
      table = ic_table
    ),
    class = "aic_arma"
  )
}


#' @export
print.aic_arma <- function(x, ...) {
  cat("\nARMA Order Selection (", toupper(x$type), ")\n\n", sep = "")
  cat("Selected order: p =", x$p, ", q =", x$q, "\n")
  cat("Criterion value:", round(x$value, 4), "\n\n")

  if (!all(x$phi == 0)) {
    cat("AR Coefficients:\n")
    print(round(x$phi, 4))
  }

  if (!all(x$theta == 0)) {
    cat("MA Coefficients:\n")
    print(round(x$theta, 4))
  }

  cat("\nMean:", round(x$xbar, 4), "\n")
  cat("Residual Variance:", round(x$vara, 4), "\n\n")
  cat("Comparison Table:\n")
  print(round(x$table, 4))
  invisible(x)
}
