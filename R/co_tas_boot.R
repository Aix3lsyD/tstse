#' Bootstrap Cochrane-Orcutt Trend Test with Turner Adjusted Sample Size
#'
#' Performs a bootstrap version of the Cochrane-Orcutt trend test with
#' Turner adjusted sample size. The bootstrap p-value is computed by
#' comparing the observed test statistic to its distribution under the
#' null hypothesis (no trend).
#'
#' @param x A numeric vector containing the time series data.
#' @param nb Number of bootstrap replicates. Default is 399.
#' @param maxp Maximum AR order for model selection. Default is 5.
#' @param type Character. Information criterion for order selection:
#'   `"aic"`, `"aicc"`, or `"bic"`. Default is `"aic"`.
#' @param btest Logical. If `FALSE` (default), bootstrap using p-values.
#'   If `TRUE`, bootstrap using t-statistics (comparing absolute values).
#' @param use_fast Logical. If `TRUE` (default), use [co_tas_fast()] for
#'   improved performance. If `FALSE`, use [co_tas()]. This parameter may
#'   be deprecated in a future version when only the fast implementation
#'   is supported.
#' @param cores Integer, number of cores for parallel processing.
#'   Default is 1 (sequential). Set to 0 to use all available cores.
#' @param seed Optional integer seed for reproducibility. Default is NULL.
#'
#' @return A list with class `"co_tas_boot"` containing:
#'   \item{pvalue}{Bootstrap p-value for test of H0: no trend.}
#'   \item{pvalue_upper}{Upper-tail p-value (only if `btest=TRUE`, NA otherwise).}
#'   \item{pvalue_lower}{Lower-tail p-value (only if `btest=TRUE`, NA otherwise).}
#'   \item{tco}{t-statistic from the original data.}
#'   \item{pvalue_asymptotic}{Asymptotic p-value from [co_tas()].}
#'   \item{phi}{AR coefficients from the CO-TAS fit (differenced series model).}
#'   \item{phi_null}{AR coefficients from null model (fit to original data).}
#'   \item{vara}{Innovation variance from null model.}
#'   \item{n}{Length of input series.}
#'   \item{nb}{Number of bootstrap replicates used.}
#'   \item{btest}{Logical indicating which bootstrap method was used.}
#'   \item{boot_pvals}{Vector of bootstrap p-values (if `btest=FALSE`).}
#'   \item{boot_tco}{Vector of bootstrap t-statistics (if `btest=TRUE`).}
#'   \item{boot_seeds}{Vector of RNG seeds used for each bootstrap replicate.}
#'   \item{master_seed}{The seed parameter used to generate boot_seeds (for reproducibility).}
#'
#' @details
#' The bootstrap procedure:
#' \enumerate{
#'   \item Fit [co_tas()] to the original data to obtain the test statistic
#'   \item Fit AR(p) model directly to the original data for the null model
#'   \item Generate `nb` bootstrap series from AR(phi_null) with no trend
#'   \item For each bootstrap series:
#'     \itemize{
#'       \item If `btest=FALSE`: compute p-value using [co_tas()]
#'       \item If `btest=TRUE`: compute t-statistic using [co_tas()]
#'     }
#'   \item Compute bootstrap p-value:
#'     \itemize{
#'       \item If `btest=FALSE`: proportion of bootstrap p-values <= observed
#'       \item If `btest=TRUE`: proportion with |t| >= |observed t|
#'     }
#' }
#'
#' The t-statistic bootstrap (`btest=TRUE`) can be more powerful when the
#' null distribution of p-values is not uniform.
#'
#' By default, this function uses [co_tas_fast()] internally for improved
#' performance (~10x faster). Set `use_fast = FALSE` to use the original
#' [co_tas()] implementation. For maximum speed with parallel execution,
#' use [co_tas_boot_fast()].
#'
#' @references
#' Woodward, W. A., Gray, H. L., and Elliott, A. C. (2017).
#' *Applied Time Series Analysis with R*. CRC Press.
#'
#' @seealso [co_tas()], [co_tas_fast()] for the non-bootstrap version,
#'   [co_tas_boot_fast()] for fully parallelized C++ version,
#'   [wbg_boot()] for WBG bootstrap trend test
#'
#' @examples
#' # Generate series with trend and AR(1) errors
#' set.seed(123)
#' n <- 100
#' t <- 1:n
#' z <- arima.sim(list(ar = 0.7), n = n)
#' x <- 5 + 0.1 * t + z
#'
#' # Bootstrap trend test using p-values (default)
#' result <- co_tas_boot(x, nb = 99, seed = 42)
#' print(result)
#'
#' # Bootstrap using t-statistics
#' result_t <- co_tas_boot(x, nb = 99, btest = TRUE, seed = 42)
#' print(result_t)
#'
#' @export
co_tas_boot <- function(x, nb = 399L, maxp = 5L,
                         type = c("aic", "aicc", "bic"),
                         btest = FALSE,
                         use_fast = TRUE,
                         cores = 1L,
                         seed = NULL) {

  type <- match.arg(type)

  # Input validation
  if (!is.numeric(x) || length(x) == 0) {
    stop("`x` must be a non-empty numeric vector", call. = FALSE)
  }
  if (!is.numeric(nb) || length(nb) != 1 || nb < 1) {
    stop("`nb` must be a positive integer", call. = FALSE)
  }
  if (!is.numeric(maxp) || length(maxp) != 1 || maxp < 1) {
    stop("`maxp` must be a positive integer", call. = FALSE)
  }
  if (!is.logical(btest) || length(btest) != 1 || is.na(btest)) {
    stop("`btest` must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(use_fast) || length(use_fast) != 1 || is.na(use_fast)) {
    stop("`use_fast` must be TRUE or FALSE", call. = FALSE)
  }

  n <- length(x)
  cores <- get_cores(cores)
  if (use_fast && n < 10) {
    stop("Series must have at least 10 observations when use_fast = TRUE", call. = FALSE)
  }
  nb <- as.integer(nb)
  maxp <- as.integer(maxp)

  # Select co_tas implementation
  co_tas_fn <- if (use_fast) co_tas_fast else co_tas

  # Generate seeds upfront for reproducibility with parallel processing
  if (!is.null(seed)) {
    set.seed(seed)
  }
  boot_seeds <- sample.int(.Machine$integer.max, nb)

  # Step 1: Fit co_tas to original data for test statistic

  result <- co_tas_fn(x, maxp = maxp, type = type)
  pval_obs <- result$pvalue
  tco_obs <- result$tco
  phi <- result$phi

  # Step 2: Fit AR to ORIGINAL data for null hypothesis
  # (Not the differenced series used in co_tas)
  null_fit <- aic_burg(x, p = seq_len(maxp), type = type)
  phi_null <- null_fit$phi
  vara <- null_fit$vara

  # Step 3: Bootstrap under null hypothesis (no trend)
  # Define closure function for bootstrap replicate
  boot_fn <- function(i) {
    set.seed(boot_seeds[i])
    boot_dat <- gen_arma(n = n, phi = phi_null, vara = vara, plot = FALSE)
    boot_result <- co_tas_fn(boot_dat, maxp = maxp, type = type)
    if (btest) boot_result$tco else boot_result$pvalue
  }

  # Execute parallel or sequential
  if (cores > 1L) {
    boot_results <- pmap(seq_len(nb), boot_fn, cores = cores)
    boot_values <- unlist(boot_results)
  } else {
    boot_values <- vapply(seq_len(nb), boot_fn, numeric(1))
  }

  # Compute bootstrap p-value
  if (btest) {
    boot_tco <- boot_values
    boot_pvals <- NULL
    pvalue_boot <- sum(abs(boot_tco) >= abs(tco_obs)) / nb
    pvalue_upper <- sum(boot_tco >= tco_obs) / nb
    pvalue_lower <- sum(boot_tco <= tco_obs) / nb
  } else {
    boot_pvals <- boot_values
    boot_tco <- NULL
    pvalue_boot <- sum(boot_pvals <= pval_obs) / nb
    pvalue_upper <- NA_real_  # Not applicable for p-value bootstrap
    pvalue_lower <- NA_real_
  }

  # Return result
  result_boot <- list(
    pvalue = pvalue_boot,
    pvalue_upper = pvalue_upper,
    pvalue_lower = pvalue_lower,
    tco = tco_obs,
    pvalue_asymptotic = pval_obs,
    phi = phi,
    phi_null = phi_null,
    vara = vara,
    n = n,
    nb = nb,
    btest = btest,
    boot_pvals = boot_pvals,
    boot_tco = boot_tco,
    boot_seeds = boot_seeds,
    master_seed = seed
  )

  class(result_boot) <- "co_tas_boot"
  result_boot
}

#' @export
print.co_tas_boot <- function(x, digits = 4, ...) {
  method <- if (isTRUE(x$btest)) "t-statistic" else "p-value"
  cat("Bootstrap Cochrane-Orcutt Trend Test (Turner Adjusted Sample Size)\n")
  cat("===================================================================\n\n")
  cat("Bootstrap method:", method, "\n")
  cat("t-statistic:", round(x$tco, digits), "\n")
  cat("Bootstrap p-value:", format.pval(x$pvalue, digits = digits), "\n")
  cat("Asymptotic p-value:", format.pval(x$pvalue_asymptotic, digits = digits), "\n")
  cat("Bootstrap replicates:", x$nb, "\n")
  if (length(x$phi) > 0) {
    cat("AR coefficients (test):", paste(round(x$phi, digits), collapse = ", "), "\n")
  }
  if (length(x$phi_null) > 0) {
    cat("AR coefficients (null):", paste(round(x$phi_null, digits), collapse = ", "), "\n")
  }
  invisible(x)
}
