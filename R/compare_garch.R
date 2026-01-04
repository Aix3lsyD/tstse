#' Compare Multiple GARCH Model Specifications
#'
#' Fits a grid of ARCH/GARCH models to a time series and computes information
#' criteria and diagnostic statistics for model comparison and selection.
#'
#' @importFrom stats setNames
#'
#' @param data Numeric vector of returns or residuals to model.
#' @param arch_range Integer vector of ARCH orders to try (default 0:2).
#' @param garch_range Integer vector of GARCH orders to try (default 0:2).
#' @param distribution Character string specifying the innovation distribution.
#'   One of: "norm", "std", "ged", "snorm", "sged", "sstd", "nig", "ghyp", "jsu".
#'   Default is "norm".
#' @param include_mean Logical. Include a mean term in the model? Default is FALSE.
#' @param solver Character string for optimization solver. Default is "hybrid"
#'   which tries multiple solvers.
#' @param parallel Logical. Use parallel processing for model fitting? Default is TRUE.
#' @param cores Integer. Number of cores for parallel processing. If NULL, uses
#'   \code{getOption("tstse.cores")} or \code{detectCores() - 1}.
#'
#' @return A list of class \code{"garch_comparison"} containing:
#'   \describe{
#'     \item{comparison}{Data frame with model statistics}
#'     \item{fits}{Named list of uGARCHfit objects}
#'     \item{distribution}{Distribution used}
#'     \item{n}{Sample size}
#'     \item{call}{The matched call}
#'   }
#'
#' @details
#' The function fits all combinations of ARCH and GARCH orders (excluding
#' the trivial (0,0) case) and computes:
#' \itemize{
#'   \item AIC, AICc, BIC - information criteria
#'   \item WLB1, WLB2, WLB3 - Weighted Ljung-Box test p-values at different lags
#'   \item Nyblom - stability test statistic and 5\% critical value
#'   \item SignBias - sign bias test p-value
#'   \item n_sig/n_coef - count of significant coefficients
#' }
#'
#' Parallel processing uses the \code{future} package with \code{multisession}
#' backend (socket-based workers). Set \code{parallel = FALSE} to disable.
#'
#' @seealso \code{\link{table_garch_gt}}, \code{\link{table_garch_cli}} for
#'   display functions.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate GARCH data
#' garch_gen <- make_gen_garch(omega = 0.1, alpha = 0.15, beta = 0.8)
#' y <- garch_gen(1000)
#'
#' # Compare models
#' results <- compare_garch(y)
#' print(results)
#'
#' # Display as table
#' table_garch_gt(results)
#' table_garch_cli(results)
#' }
compare_garch <- function(data,
                          arch_range = 0:2,
                          garch_range = 0:2,
                          distribution = "norm",
                          include_mean = FALSE,
                          solver = "hybrid",
                          parallel = FALSE,
                          cores = NULL) {

  # Check required packages
  if (!requireNamespace("rugarch", quietly = TRUE)) {
    stop("Package 'rugarch' is required. Install with: install.packages('rugarch')")
  }

  # WeightedPortTest is OPTIONAL - warn if missing

  has_wlb <- requireNamespace("WeightedPortTest", quietly = TRUE)
  if (!has_wlb) {
    warning("Package 'WeightedPortTest' not installed. Weighted Ljung-Box tests will be skipped.",
            "\n  Install with: install.packages('WeightedPortTest')", call. = FALSE)
  }

  # Input validation
  if (!is.numeric(data) || length(data) < 10) {
    stop("data must be a numeric vector with at least 10 observations")
  }

  if (any(is.na(data))) {
    warning("NA values detected in data; these may cause fitting issues")
  }

  # Validate distribution
  valid_dists <- c("norm", "std", "ged", "snorm", "sged", "sstd", "nig", "ghyp", "jsu")
  if (!distribution %in% valid_dists) {
    stop("distribution must be one of: ", paste(valid_dists, collapse = ", "))
  }

  # Build grid of orders, excluding (0,0)
  orders <- expand.grid(arch = arch_range, garch = garch_range)
  orders <- orders[!(orders$arch == 0 & orders$garch == 0), , drop = FALSE]

  if (nrow(orders) == 0) {
    stop("No valid model orders to fit (need at least one non-zero order)")
  }

  n <- length(data)
  call <- match.call()

  # Convert orders to list for parallel processing
  order_list <- split(orders, seq_len(nrow(orders)))

  # Function to fit a single model
  fit_single_model <- function(order_row) {
    arch_order <- order_row$arch
    garch_order <- order_row$garch
    model_label <- .garch_label(arch_order, garch_order)

    tryCatch({
      spec <- rugarch::ugarchspec(
        variance.model = list(model = "sGARCH", garchOrder = c(arch_order, garch_order)),
        mean.model = list(armaOrder = c(0, 0), include.mean = include_mean),
        distribution.model = distribution
      )

      fit <- rugarch::ugarchfit(spec, data, solver = solver)
      diag <- .extract_garch_diagnostics(fit, arch_order, garch_order, n, has_wlb)

      list(
        success = TRUE,
        label = model_label,
        fit = fit,
        result = data.frame(
          Model = model_label,
          ARCH = arch_order,
          GARCH = garch_order,
          AIC = diag$aic,
          AICC = diag$aicc,
          BIC = diag$bic,
          WLB1 = diag$wlb[1],
          WLB2 = diag$wlb[2],
          WLB3 = diag$wlb[3],
          Nyblom = diag$nyblom_stat,
          Nyblom_crit = diag$nyblom_crit,
          SignBias = diag$signbias_pval,
          n_sig = diag$n_sig,
          n_coef = diag$n_coef,
          stringsAsFactors = FALSE
        )
      )
    }, error = function(e) {
      message(paste0(model_label, " failed: ", e$message))
      list(success = FALSE, label = model_label, fit = NULL, result = NULL)
    })
  }

  # Execute fitting (parallel or sequential)
  # Uses pmap() for consistent parallel handling and nesting guard
  if (parallel && length(order_list) > 1) {
    cores <- get_cores(cores)
    results_list <- pmap(order_list, fit_single_model, cores = cores)
  } else {
    # Sequential
    results_list <- lapply(order_list, fit_single_model)
  }

  # Collect successful results
  successful <- Filter(function(x) x$success, results_list)

  if (length(successful) == 0) {
    stop("All models failed to converge. Check your data or try different orders.")
  }

  # Build fits list and comparison data frame
  fits <- setNames(
    lapply(successful, function(x) x$fit),
    sapply(successful, function(x) x$label)
  )

  comparison <- do.call(rbind, lapply(successful, function(x) x$result))
  rownames(comparison) <- NULL

  structure(
    list(
      comparison = comparison,
      fits = fits,
      distribution = distribution,
      n = n,
      call = call
    ),
    class = "garch_comparison"
  )
}


# ==============================================================================
# Helper Functions
# ==============================================================================

#' Create GARCH model label
#' @noRd
.garch_label <- function(arch_order, garch_order) {
  if (garch_order == 0) {
    paste0("ARCH(", arch_order, ")")
  } else {
    paste0("GARCH(", arch_order, ",", garch_order, ")")
  }
}


#' Extract GARCH diagnostics from a fitted model
#' @noRd
.extract_garch_diagnostics <- function(fit, arch_order, garch_order, n, has_wlb = TRUE) {

  std_resid <- as.numeric(rugarch::residuals(fit, standardize = TRUE))

  # Count parameters: omega + arch coefficients + garch coefficients
  # Plus distribution params if non-normal (shape, skew, etc.)
  k <- arch_order + garch_order + 1

  # Add distribution parameters to k
 dist <- fit@model$modeldesc$distribution
  if (dist %in% c("std", "ged")) {
    k <- k + 1  # shape parameter
  } else if (dist %in% c("snorm", "sged")) {
    k <- k + 1  # skew parameter
  } else if (dist == "sstd") {
    k <- k + 2
  } else if (dist %in% c("nig", "ghyp", "jsu")) {
    k <- k + 2  # shape + skew
  }

  # Information criteria computed from likelihood directly
  # This avoids confusion about rugarch's per-observation scaling
  ll <- rugarch::likelihood(fit)
  aic <- -2 * ll + 2 * k
  bic <- -2 * ll + k * log(n)
  aicc <- aic + (2 * k^2 + 2 * k) / (n - k - 1)

  # Weighted Ljung-Box on squared standardized residuals
  # Only compute if WeightedPortTest is available
  if (has_wlb) {
    df_garch <- max(1, arch_order + garch_order)
    wlb <- tryCatch({
      box1 <- WeightedPortTest::Weighted.Box.test(
        std_resid, lag = 1,
        type = "Ljung-Box", fitdf = 0, sqrd.res = TRUE
      )
      box2 <- WeightedPortTest::Weighted.Box.test(
        std_resid,
        lag = max(2, 3 * df_garch - 1),
        type = "Ljung-Box", fitdf = df_garch, sqrd.res = TRUE
      )
      box3 <- WeightedPortTest::Weighted.Box.test(
        std_resid,
        lag = max(5, 5 * df_garch - 1),
        type = "Ljung-Box", fitdf = df_garch, sqrd.res = TRUE
      )
      c(box1$p.value, box2$p.value, box3$p.value)
    }, error = function(e) {
      c(NA_real_, NA_real_, NA_real_)
    })
  } else {
    wlb <- c(NA_real_, NA_real_, NA_real_)
  }

  # Nyblom stability test
  nyb <- tryCatch({
    rugarch::nyblom(fit)
  }, error = function(e) {
    list(JointStat = NA_real_, JointCritical = c("5%" = NA_real_))
  })
  nyblom_stat <- nyb$JointStat
  nyblom_crit <- nyb$JointCritical["5%"]

  # Sign bias test
  sb <- tryCatch({
    rugarch::signbias(fit)
  }, error = function(e) {
    matrix(NA_real_, nrow = 4, ncol = 2)
  })
  signbias_pval <- sb[4, 2]

  # Coefficient significance
  coef_mat <- fit@fit$matcoef
  n_coef <- nrow(coef_mat)
  n_sig <- sum(coef_mat[, 4] < 0.05, na.rm = TRUE)

  list(
    aic = aic,
    aicc = aicc,
    bic = bic,
    wlb = wlb,
    nyblom_stat = nyblom_stat,
    nyblom_crit = nyblom_crit,
    signbias_pval = signbias_pval,
    n_sig = n_sig,
    n_coef = n_coef
  )
}


# ==============================================================================
# S3 Methods
# ==============================================================================

#' @describeIn compare_garch Print method for garch_comparison objects
#' @param x A garch_comparison object
#' @param ... Additional arguments (ignored)
#' @export
print.garch_comparison <- function(x, ...) {
  cat("GARCH Model Comparison\n")
  cat("----------------------\n")
  cat("Distribution:", x$distribution, "\n")
  cat("Sample size:", x$n, "\n")
  cat("Models fitted:", length(x$fits), "\n\n")

  # Find best by each criterion
  comp <- x$comparison
  best_aic <- comp$Model[which.min(comp$AIC)]
  best_bic <- comp$Model[which.min(comp$BIC)]
  best_aicc <- comp$Model[which.min(comp$AICC)]

  cat("Best by AIC: ", best_aic, "\n")
  cat("Best by BIC: ", best_bic, "\n")
  cat("Best by AICc:", best_aicc, "\n")
  cat("\nUse table_garch_gt() or table_garch_cli() for full comparison.\n")

  invisible(x)
}


#' @describeIn compare_garch Summary method for garch_comparison objects
#' @param object A garch_comparison object
#' @param ... Additional arguments (ignored)
#' @export
summary.garch_comparison <- function(object, ...) {
  print(object$comparison)
  invisible(object)
}
