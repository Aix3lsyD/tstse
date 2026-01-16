#' WBG Bootstrap Test for Trend with GARCH Errors
#'
#' Performs a bootstrap hypothesis test for trend in a time series with
#' heteroscedastic (GARCH) errors. The null model is AR(p) + GARCH(1,1),
#' which accounts for volatility clustering while testing for trend.
#'
#' @param x Numeric vector. The observed time series.
#' @param stat_fn Function. A function that takes a numeric vector and returns
#'   a single numeric test statistic. Use factory functions like
#'   [make_stat_co()], [make_stat_ols_t()], etc., or create your own with
#'   signature `function(x) -> numeric(1)`.
#' @param nb Integer. Number of bootstrap replicates. Default is 399.
#' @param p_max Integer. Maximum AR order to consider when fitting the null
#'   model. Order is selected by information criterion. Default is 5.
#' @param criterion Character. Information criterion for AR order selection:
#'   `"aic"` (default), `"aicc"`, or `"bic"`.
#' @param distribution Character. Innovation distribution for GARCH model.
#'   Options include `"norm"` (default), `"std"` (Student's t), `"ged"`,
#'   `"snorm"`, `"sstd"`, `"sged"`, `"nig"`, `"jsu"`.
#' @param stationary_check Character. How to handle non-stationary AR fits:
#'   `"warn"` (default) issues a warning, `"error"` stops execution,
#'   `"none"` skips the check.
#' @param bootadj Logical. If TRUE, performs the COBA (Cochrane-Orcutt
#'   Bootstrap Adjustment) variance adjustment. Default is FALSE because
#'   fitting GARCH to each bootstrap sample is computationally expensive.
#' @param cores Integer or NULL. Number of CPU cores for parallel processing.
#'   Default 1 (sequential). Set to 0 to use all available cores.
#' @param seed Integer or NULL. Random seed for reproducibility.
#' @param verbose Logical. If TRUE (default), display progress messages.
#'
#' @return An object of class `"wbg_boot_garch"` containing:
#'   \item{obs_stat}{Observed test statistic.}
#'   \item{boot_dist}{Vector of bootstrap statistics.}
#'   \item{pvalue}{Two-sided p-value.}
#'   \item{pvalue_upper}{Upper-tail p-value.}
#'   \item{pvalue_lower}{Lower-tail p-value.}
#'   \item{ar_order}{Selected AR order under null.}
#'   \item{garch_order}{GARCH order (fixed at c(1,1)).}
#'   \item{ar_coef}{AR coefficients under null.}
#'   \item{garch_coef}{GARCH coefficients (omega, alpha1, beta1).}
#'   \item{distribution}{Innovation distribution used.}
#'   \item{n}{Sample size.}
#'   \item{nb}{Number of bootstrap replicates.}
#'   \item{ic}{Information criterion value of selected model.}
#'   \item{boot_seeds}{Vector of RNG seeds used for each bootstrap replicate.}
#'   \item{master_seed}{The seed parameter used to generate boot_seeds (for reproducibility).}
#'
#'   If `bootadj = TRUE`, also includes:
#'   \item{obs_stat_adj}{COBA-adjusted observed statistic.}
#'   \item{pvalue_adj}{COBA-adjusted two-sided p-value.}
#'   \item{adj_factor}{Variance adjustment factor.}
#'   \item{median_coef}{Median AR-GARCH coefficients from first stage.}
#'
#' @details
#' This function extends the Woodward-Bottone-Gray (WBG) bootstrap test for
#' trend to handle heteroscedastic errors. The null hypothesis is that the
#' series is AR(p) + GARCH(1,1) with no trend.
#'
#' **Procedure:**
#' 1. Compute the observed test statistic `stat_fn(x)`
#' 2. Fit AR(p) + GARCH(1,1) models for p = 0, ..., p_max via rugarch
#' 3. Select best p by information criterion (joint model selection)
#' 4. Generate `nb` bootstrap series from the fitted model using
#'    `rugarch::ugarchpath()`
#' 5. Compute test statistic for each bootstrap series
#' 6. Calculate p-value as proportion of |bootstrap stat| >= |observed stat|
#'
#' **AR Order Selection:**
#' Unlike [wbg_boot_flex()], which selects AR order from a pure AR model,
#' this function fits the full AR(p) + GARCH(1,1) model for each candidate
#' order and selects based on the joint model's information criterion.
#'
#' **Stationarity:**
#' The AR model must be stationary for bootstrap samples to be well-behaved.
#' rugarch does not enforce AR stationarity by default. Use `stationary_check`
#' to control behavior when non-stationary fits are detected.
#'
#' @references
#' Woodward, W. A., Bottone, S., and Gray, H. L. (1997). "Improved Tests for
#' Trend in Time Series Data." *Journal of Agricultural, Biological, and
#' Environmental Statistics*, 2(4), 403-416.
#'
#' @seealso [wbg_boot_flex()] for the pure AR version,
#'   [make_stat_co()] and other statistic factory functions,
#'   [compare_garch()] for GARCH model comparison.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("rugarch", quietly = TRUE)) {
#'   # Test for trend in white noise (should not reject)
#'   set.seed(123)
#'   x <- rnorm(200)
#'   stat_fn <- make_stat_co()
#'   result <- wbg_boot_garch(x, stat_fn, nb = 99, seed = 456, verbose = FALSE)
#'   print(result)
#'
#'   # Test for trend in series with actual trend (should reject)
#'   set.seed(123)
#'   x_trend <- 0.1 * (1:200) + arima.sim(list(ar = 0.7), n = 200)
#'   result <- wbg_boot_garch(x_trend, stat_fn, nb = 99, seed = 456, verbose = FALSE)
#'   result$pvalue  # Should be small
#' }
#' }
#'
#' @export
wbg_boot_garch <- function(x, stat_fn = make_stat_co(), nb = 399L, p_max = 5L,
                            criterion = c("aic", "aicc", "bic"),
                            distribution = "norm",
                            stationary_check = c("warn", "error", "none"),
                            bootadj = FALSE,
                            cores = 1L, seed = NULL,
                            verbose = TRUE) {

  criterion <- match.arg(criterion)
  stationary_check <- match.arg(stationary_check)


  # Check rugarch is installed

  if (!requireNamespace("rugarch", quietly = TRUE)) {
    stop("Package 'rugarch' is required. Install with: install.packages('rugarch')",
         call. = FALSE)
  }

  # Input validation
  if (!is.numeric(x) || length(x) == 0) {
    stop("`x` must be a non-empty numeric vector", call. = FALSE)
  }
  if (!is.function(stat_fn)) {
    stop("`stat_fn` must be a function", call. = FALSE)
  }
  if (!is.numeric(nb) || length(nb) != 1 || nb < 1) {
    stop("`nb` must be a positive integer", call. = FALSE)
  }
  if (!is.numeric(p_max) || length(p_max) != 1 || p_max < 0) {
    stop("`p_max` must be a non-negative integer", call. = FALSE)
  }

  nb <- as.integer(nb)
  p_max <- as.integer(p_max)
  n <- length(x)
  cores <- get_cores(cores)

  # Set seed and generate boot_seeds for reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
  }
  boot_seeds <- sample.int(.Machine$integer.max, nb)
  if (bootadj) {
    boot_seeds_adj <- sample.int(.Machine$integer.max, nb)
  }

  # Compute observed statistic
  obs_stat <- stat_fn(x)

  # Joint AR order selection via rugarch
  if (verbose) message("Selecting AR order and fitting AR-GARCH model...")

  model_select <- .select_ar_garch(x, p_max = p_max, criterion = criterion,
                                    distribution = distribution, verbose = verbose)

  if (is.null(model_select$fit)) {
    stop("Failed to fit any AR-GARCH model. Check data for issues.", call. = FALSE)
  }

  fit <- model_select$fit
  ar_p <- model_select$p
  best_ic <- model_select$ic

  # Extract coefficients
  fitted_coef <- rugarch::coef(fit)

  # Extract AR and GARCH coefficients separately
  ar_names <- grep("^ar", names(fitted_coef), value = TRUE)
  ar_coef <- if (length(ar_names) > 0) fitted_coef[ar_names] else numeric(0)

  garch_names <- c("omega", "alpha1", "beta1")
  garch_coef <- fitted_coef[intersect(garch_names, names(fitted_coef))]

  mu_coef <- if ("mu" %in% names(fitted_coef)) fitted_coef["mu"] else NULL

  if (verbose) {
    message("  Selected AR(", ar_p, ") + GARCH(1,1)")
    if (length(ar_coef) > 0) {
      message("  AR coefficients: [", paste(round(ar_coef, 4), collapse = ", "), "]")
    }
    message("  GARCH coefficients: omega=", round(garch_coef["omega"], 6),
            ", alpha=", round(garch_coef["alpha1"], 4),
            ", beta=", round(garch_coef["beta1"], 4))
  }

  # Post-hoc stationarity check
  if (stationary_check != "none" && length(ar_coef) > 0) {
    roots <- polyroot(c(1, -ar_coef))
    is_stationary <- all(Mod(roots) > 1)

    if (!is_stationary) {
      msg <- "Fitted AR model is non-stationary; bootstrap results may be unreliable"
      if (stationary_check == "error") {
        stop(msg, call. = FALSE)
      } else {
        warning(msg, call. = FALSE)
      }
    }
  }

  # Create fixed spec for bootstrap generation
  fixed_spec <- rugarch::ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
    mean.model = list(armaOrder = c(ar_p, 0), include.mean = TRUE),
    distribution.model = distribution,
    fixed.pars = as.list(fitted_coef)
  )

  # --- First Bootstrap ---
  if (verbose) message("Running bootstrap (", nb, " replicates)...")

  boot_fn <- function(i) {
    set.seed(boot_seeds[i])

    # Generate path from fitted model
    path <- tryCatch(
      rugarch::ugarchpath(fixed_spec, n.sim = n, m.sim = 1),
      error = function(e) NULL
    )

    if (is.null(path)) {
      return(list(stat = NA_real_, coef = NULL, phi1 = NA_real_))
    }

    xb <- as.numeric(path@path$seriesSim)
    stat <- stat_fn(xb)

    if (bootadj) {
      # Fit AR-GARCH to bootstrap sample for COBA
      ar_boot <- tryCatch({
        spec_boot <- rugarch::ugarchspec(
          variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
          mean.model = list(armaOrder = c(ar_p, 0), include.mean = TRUE),
          distribution.model = distribution
        )
        fit_boot <- rugarch::ugarchfit(spec_boot, xb, solver = "hybrid")
        coef_boot <- rugarch::coef(fit_boot)
        ar_names_boot <- grep("^ar", names(coef_boot), value = TRUE)
        ar_coefs_boot <- if (length(ar_names_boot) > 0) coef_boot[ar_names_boot] else numeric(0)
        phi1 <- 1 - sum(ar_coefs_boot)
        list(coef = coef_boot, phi1 = phi1)
      }, error = function(e) {
        list(coef = NULL, phi1 = NA_real_)
      })
      list(stat = stat, coef = ar_boot$coef, phi1 = ar_boot$phi1)
    } else {
      list(stat = stat)
    }
  }

  # Run bootstrap (parallel or sequential)
  if (cores > 1L) {
    boot_results <- pmap(seq_len(nb), boot_fn, cores = cores)
  } else {
    boot_results <- lapply(seq_len(nb), boot_fn)
  }

  boot_stats <- sapply(boot_results, `[[`, "stat")

  # Check for failed bootstrap samples
  n_failed <- sum(is.na(boot_stats))
  if (n_failed > 0) {
    if (n_failed > nb / 2) {
      stop("More than half of bootstrap samples failed to generate", call. = FALSE)
    }
    warning(n_failed, " bootstrap samples failed; using ", nb - n_failed, " samples",
            call. = FALSE)
    boot_stats <- boot_stats[!is.na(boot_stats)]
  }

  # Compute p-values (with +1 correction for proper bootstrap p-values)
  nb_actual <- length(boot_stats)
  pvalue_two <- (sum(abs(boot_stats) >= abs(obs_stat)) + 1) / (nb_actual + 1)
  pvalue_upper <- (sum(boot_stats >= obs_stat) + 1) / (nb_actual + 1)
  pvalue_lower <- (sum(boot_stats <= obs_stat) + 1) / (nb_actual + 1)

  result <- list(
    obs_stat = obs_stat,
    boot_dist = boot_stats,
    pvalue = pvalue_two,
    pvalue_upper = pvalue_upper,
    pvalue_lower = pvalue_lower,
    ar_order = ar_p,
    garch_order = c(1, 1),
    ar_coef = ar_coef,
    garch_coef = garch_coef,
    mu = mu_coef,
    distribution = distribution,
    n = n,
    nb = nb_actual,
    ic = best_ic,
    criterion = criterion,
    boot_seeds = boot_seeds,
    master_seed = seed
  )

  # --- COBA Adjustment (Second Bootstrap) ---
  if (bootadj) {
    if (verbose) message("Running COBA adjustment...")

    boot_phi1 <- sapply(boot_results, `[[`, "phi1")
    boot_coef_list <- lapply(boot_results, `[[`, "coef")

    # Filter out failed fits
    valid_idx <- !is.na(boot_phi1) & !sapply(boot_coef_list, is.null)

    if (sum(valid_idx) < nb / 2) {
      warning("More than half of bootstrap GARCH fits failed; COBA adjustment unreliable",
              call. = FALSE)
    }

    if (sum(valid_idx) > 0) {
      # Find median AR model
      valid_phi1 <- boot_phi1[valid_idx]
      valid_coef_list <- boot_coef_list[valid_idx]
      median_idx <- which.min(abs(valid_phi1 - median(valid_phi1)))
      median_coef <- valid_coef_list[[median_idx]]

      # Second stage bootstrap with median model
      median_spec <- rugarch::ugarchspec(
        variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
        mean.model = list(armaOrder = c(ar_p, 0), include.mean = TRUE),
        distribution.model = distribution,
        fixed.pars = as.list(median_coef)
      )

      boot_fn_adj <- function(i) {
        set.seed(boot_seeds_adj[i])
        path <- tryCatch(
          rugarch::ugarchpath(median_spec, n.sim = n, m.sim = 1),
          error = function(e) NULL
        )
        if (is.null(path)) return(NA_real_)
        xb <- as.numeric(path@path$seriesSim)
        stat_fn(xb)
      }

      if (cores > 1L) {
        boot_stats_adj <- unlist(pmap(seq_len(nb), boot_fn_adj, cores = cores))
      } else {
        boot_stats_adj <- vapply(seq_len(nb), boot_fn_adj, numeric(1))
      }

      boot_stats_adj <- boot_stats_adj[!is.na(boot_stats_adj)]

      if (length(boot_stats_adj) > 1 && sd(boot_stats) > 0) {
        adj_factor <- sd(boot_stats_adj) / sd(boot_stats)
        obs_stat_adj <- adj_factor * obs_stat
        pvalue_adj <- (sum(abs(boot_stats) >= abs(obs_stat_adj)) + 1) / (nb_actual + 1)

        result$obs_stat_adj <- obs_stat_adj
        result$pvalue_adj <- pvalue_adj
        result$adj_factor <- adj_factor
        result$median_coef <- median_coef
        result$boot_seeds_adj <- boot_seeds_adj
      }
    }
  }

  if (verbose) message("Done.")

  class(result) <- "wbg_boot_garch"
  result
}


#' Select AR order for AR-GARCH model via information criterion
#'
#' @param x Numeric vector. Time series data.
#' @param p_max Integer. Maximum AR order to try.
#' @param criterion Character. Information criterion ("aic", "aicc", "bic").
#' @param distribution Character. GARCH innovation distribution.
#' @param verbose Logical. Show progress.
#' @return List with p (selected order), fit (rugarch fit), ic (IC value).
#' @noRd
.select_ar_garch <- function(x, p_max, criterion, distribution, verbose = FALSE) {
  best_ic <- Inf
  best_p <- 0

  best_fit <- NULL
  n <- length(x)

  for (p in 0:p_max) {
    if (verbose) message("  Trying AR(", p, ") + GARCH(1,1)...")

    spec <- rugarch::ugarchspec(
      variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
      mean.model = list(armaOrder = c(p, 0), include.mean = TRUE),
      distribution.model = distribution
    )

    fit <- tryCatch(
      rugarch::ugarchfit(spec, x, solver = "hybrid"),
      error = function(e) NULL
    )

    if (!is.null(fit)) {
      ic <- rugarch::infocriteria(fit)

      # rugarch returns a matrix; extract the appropriate criterion
      ic_value <- switch(criterion,
        "aic" = ic["Akaike", 1],
        "bic" = ic["Bayes", 1],
        "aicc" = {
          # AICc correction: AIC + (2k^2 + 2k)/(n - k - 1)
          aic_val <- ic["Akaike", 1]
          k <- length(rugarch::coef(fit))
          if (n > k + 1) {
            aic_val + (2 * k^2 + 2 * k) / (n - k - 1)
          } else {
            Inf  # Not enough data for AICc
          }
        }
      )

      if (verbose) message("    ", toupper(criterion), " = ", round(ic_value, 4))

      if (ic_value < best_ic) {
        best_ic <- ic_value
        best_p <- p
        best_fit <- fit
      }
    } else {
      if (verbose) message("    Failed to converge")
    }
  }

  list(p = best_p, fit = best_fit, ic = best_ic)
}


#' @describeIn wbg_boot_garch Print method for wbg_boot_garch objects
#' @param x A wbg_boot_garch object
#' @param ... Additional arguments (ignored)
#' @export
print.wbg_boot_garch <- function(x, ...) {
  cat("WBG Bootstrap Test (AR-GARCH)\n")
  cat("==============================\n")
  cat("Sample size:", x$n, "\n")
  cat("Bootstrap replicates:", x$nb, "\n")
  cat("Model: AR(", x$ar_order, ") + GARCH(1,1)\n", sep = "")
  cat("Distribution:", x$distribution, "\n")
  cat("Criterion:", toupper(x$criterion), "=", round(x$ic, 4), "\n")
  cat("\n")
  if (length(x$ar_coef) > 0) {
    cat("AR coefficients:", paste(round(x$ar_coef, 4), collapse = ", "), "\n")
  }
  if (!is.null(x$mu)) {
    cat("Mean:", round(x$mu, 4), "\n")
  }
  cat("GARCH: omega=", round(x$garch_coef["omega"], 6),
      ", alpha=", round(x$garch_coef["alpha1"], 4),
      ", beta=", round(x$garch_coef["beta1"], 4), "\n", sep = "")
  cat("\n")
  cat("Observed statistic:", round(x$obs_stat, 4), "\n")
  cat("\n")
  cat("P-values:\n")
  cat("  Two-sided:", format.pval(x$pvalue, digits = 4), "\n")
  cat("  Upper-tail:", format.pval(x$pvalue_upper, digits = 4), "\n")
  cat("  Lower-tail:", format.pval(x$pvalue_lower, digits = 4), "\n")

  if (!is.null(x$pvalue_adj)) {
    cat("\nCOBA Adjustment:\n")
    cat("  Adjustment factor:", round(x$adj_factor, 4), "\n")
    cat("  Adjusted statistic:", round(x$obs_stat_adj, 4), "\n")
    cat("  Adjusted p-value:", format.pval(x$pvalue_adj, digits = 4), "\n")
  }

  invisible(x)
}


#' @describeIn wbg_boot_garch Summary method for wbg_boot_garch objects
#' @param object A wbg_boot_garch object
#' @param ... Additional arguments (ignored)
#' @export
summary.wbg_boot_garch <- function(object, ...) {
  cat("WBG Bootstrap Test Summary (AR-GARCH)\n")
  cat("-------------------------------------\n")
  cat("n =", object$n, ", nb =", object$nb, "\n")
  cat("AR(", object$ar_order, ") + GARCH(1,1) null model\n", sep = "")
  cat("Distribution:", object$distribution, "\n")
  cat("\n")
  cat("Observed:", round(object$obs_stat, 4), "\n")
  cat("Bootstrap mean:", round(mean(object$boot_dist), 4), "\n")
  cat("Bootstrap SD:", round(sd(object$boot_dist), 4), "\n")
  cat("Bootstrap range: [", round(min(object$boot_dist), 4), ", ",
      round(max(object$boot_dist), 4), "]\n", sep = "")
  cat("\n")

  # Significance interpretation
  if (object$pvalue < 0.01) {
    sig <- "*** (p < 0.01)"
  } else if (object$pvalue < 0.05) {
    sig <- "** (p < 0.05)"
  } else if (object$pvalue < 0.10) {
    sig <- "* (p < 0.10)"
  } else {
    sig <- "(not significant)"
  }
  cat("p-value:", format.pval(object$pvalue, digits = 4), sig, "\n")

  if (!is.null(object$pvalue_adj)) {
    if (object$pvalue_adj < 0.01) {
      sig_adj <- "*** (p < 0.01)"
    } else if (object$pvalue_adj < 0.05) {
      sig_adj <- "** (p < 0.05)"
    } else if (object$pvalue_adj < 0.10) {
      sig_adj <- "* (p < 0.10)"
    } else {
      sig_adj <- "(not significant)"
    }
    cat("COBA-adjusted p-value:", format.pval(object$pvalue_adj, digits = 4),
        sig_adj, "\n")
  }

  invisible(object)
}
