#' Fast Bootstrap Cochrane-Orcutt Trend Test with Turner Adjusted Sample Size
#'
#' High-performance bootstrap version of the Cochrane-Orcutt trend test with
#' Turner adjusted sample size. Uses C++ parallelization via TBB for
#' approximately 30-50x speedup over [co_tas_boot()].
#'
#' @param x A numeric vector containing the time series data.
#' @param nb Number of bootstrap replicates. Default is 399.
#' @param maxp Maximum AR order for model selection. Default is 5.
#' @param type Character. Information criterion for order selection:
#'   `"aic"`, `"aicc"`, or `"bic"`. Default is `"aic"`.
#' @param btest Logical. If `FALSE` (default), bootstrap using p-values.
#'   If `TRUE`, bootstrap using t-statistics (comparing absolute values).
#' @param seed Optional integer seed for reproducibility. Default is NULL.
#'
#' @return A list with class `"co_tas_boot"` containing:
#'   \item{pvalue}{Bootstrap p-value for test of H0: no trend.}
#'   \item{tco}{t-statistic from the original data.}
#'   \item{pvalue_asymptotic}{Asymptotic p-value from [co_tas_fast()].}
#'   \item{phi}{AR coefficients from the CO-TAS fit (differenced series model).}
#'   \item{phi_null}{AR coefficients from null model (fit to original data).}
#'   \item{vara}{Innovation variance from null model.}
#'   \item{nb}{Number of bootstrap replicates used.}
#'   \item{btest}{Logical indicating which bootstrap method was used.}
#'   \item{boot_pvals}{Vector of bootstrap p-values (if `btest=FALSE`).}
#'   \item{boot_tco}{Vector of bootstrap t-statistics (if `btest=TRUE`).}
#'   \item{boot_seeds}{Vector of RNG seeds used for each bootstrap replicate.}
#'
#' @details
#' This is a high-performance C++ implementation of [co_tas_boot()],
#' optimized for large-scale applications. Key differences:
#'
#' \itemize{
#'   \item Parallelization via TBB (works on all platforms)
#'   \item All bootstrap iterations run in C++ (no R overhead)
#'   \item Thread-safe RNG (dqrng with xoshiro256+)
#' }
#'
#' Expected speedup is 30-50x compared to [co_tas_boot()] depending on
#' series length and number of available CPU cores.
#'
#' The bootstrap procedure:
#' \enumerate{
#'   \item Fit [co_tas_fast()] to the original data to obtain the test statistic
#'   \item Fit AR(p) model directly to the original data for the null model
#'   \item Generate `nb` bootstrap series from AR(phi_null) with no trend in parallel
#'   \item For each bootstrap series:
#'     \itemize{
#'       \item If `btest=FALSE`: compute p-value using C++ CO-TAS
#'       \item If `btest=TRUE`: compute t-statistic using C++ CO-TAS
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
#' **Performance Tuning:**
#' For optimal parallel performance, you may want to disable BLAS multi-threading
#' to avoid nested parallelism. See [wbg_boot_fast()] for details.
#'
#' **Reproducibility:**
#' Results are exactly reproducible given the same `seed`, platform, compiler,
#' and math library. Minor platform differences may occur (see [wbg_boot_fast()]).
#'
#' @references
#' Turner, J. R. (1988). Effective sample size: A frequently neglected concept.
#' *Journal of Educational Statistics*, 13(2), 175-181.
#'
#' Woodward, W. A., Gray, H. L., and Elliott, A. C. (2017).
#' *Applied Time Series Analysis with R*. CRC Press.
#'
#' @seealso [co_tas_boot()] for the R implementation,
#'   [co_tas_fast()] for single-run fast version,
#'   [wbg_boot_fast()] for WBG bootstrap trend test
#'
#' @examples
#' \donttest{
#' # Generate series with trend and AR(1) errors
#' set.seed(123)
#' n <- 100
#' t <- 1:n
#' z <- arima.sim(list(ar = 0.7), n = n)
#' x <- 5 + 0.1 * t + z
#'
#' # Fast bootstrap trend test using p-values (default)
#' result <- co_tas_boot_fast(x, nb = 199, seed = 42)
#' print(result)
#'
#' # Fast bootstrap using t-statistics
#' result_t <- co_tas_boot_fast(x, nb = 199, btest = TRUE, seed = 42)
#' print(result_t)
#' }
#'
#' \dontrun{
#' # Benchmark comparison
#' library(microbenchmark)
#' x <- arima.sim(list(ar = 0.7), n = 100)
#' microbenchmark(
#'   fast = co_tas_boot_fast(x, nb = 199, seed = 42),
#'   original = co_tas_boot(x, nb = 199, seed = 42),
#'   times = 5
#' )
#' }
#'
#' @export
co_tas_boot_fast <- function(x, nb = 399L, maxp = 5L,
                              type = c("aic", "aicc", "bic"),
                              btest = FALSE,
                              seed = NULL) {

  type <- match.arg(type)

  # Input validation
  if (!is.numeric(x) || length(x) == 0) {
    stop("`x` must be a non-empty numeric vector", call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop("`x` must contain only finite values (no NA, NaN, or Inf)", call. = FALSE)
  }
  if (!is.numeric(nb) || length(nb) != 1 || nb < 1) {
    stop("`nb` must be a positive integer", call. = FALSE)
  }
  if (!is.numeric(maxp) || length(maxp) != 1 || maxp < 1) {
    stop("`maxp` must be a positive integer", call. = FALSE)
  }
  if (maxp > 20L) {
    stop("`maxp` must be <= 20 (C++ buffer limit)", call. = FALSE)
  }
  if (!is.logical(btest) || length(btest) != 1 || is.na(btest)) {
    stop("`btest` must be TRUE or FALSE", call. = FALSE)
  }

  n <- length(x)
  nb <- as.integer(nb)
  maxp <- as.integer(maxp)

  if (n < 10) {
    stop("Series must have at least 10 observations", call. = FALSE)
  }
  if (n <= maxp + 3) {
    stop("Series too short: need n > maxp + 3", call. = FALSE)
  }

  # Set seed for reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Generate bootstrap seeds upfront (R's RNG, single-threaded)
  boot_seeds <- as.numeric(sample.int(.Machine$integer.max, nb))

  # Step 1: Fit co_tas_fast to original data for test statistic
  result <- co_tas_fast(x, maxp = maxp, type = type)
  pval_obs <- result$pvalue
  tco_obs <- result$tco
  phi <- result$phi

  # Step 2: Fit AR to ORIGINAL data for null hypothesis
  # (Not the differenced series used in co_tas)
  null_fit <- burg_aic_select_cpp(x, maxp, type)
  phi_null <- as.numeric(null_fit$phi)
  vara <- null_fit$vara

  # Handle AR(0) case
  if (length(phi_null) == 0) {
    phi_null <- numeric(0)
  }

  # Step 3: Bootstrap under null hypothesis (no trend)
  num_threads <- RcppParallel::defaultNumThreads()
  grain_size <- max(16L, as.integer(nb / (4L * num_threads)))

  if (!btest) {
    # Bootstrap using p-values
    boot_pvals <- co_tas_boot_kernel_cpp(
      n = n,
      phi = phi_null,
      vara = vara,
      seeds = boot_seeds,
      maxp = maxp,
      type = type,
      grain_size = grain_size
    )

    # Bootstrap p-value = proportion of boot p-values <= observed
    pvalue_boot <- sum(boot_pvals <= pval_obs) / nb
    boot_tco <- NULL

  } else {
    # Bootstrap using t-statistics
    boot_tco <- co_tas_boot_tstat_kernel_cpp(
      n = n,
      phi = phi_null,
      vara = vara,
      seeds = boot_seeds,
      maxp = maxp,
      type = type,
      grain_size = grain_size
    )

    # Bootstrap p-value = proportion with |t| >= |observed t|
    pvalue_boot <- sum(abs(boot_tco) >= abs(tco_obs)) / nb
    boot_pvals <- NULL
  }

  # Return result with same structure as co_tas_boot for compatibility
  result_boot <- list(
    pvalue = pvalue_boot,
    tco = tco_obs,
    pvalue_asymptotic = pval_obs,
    phi = phi,
    phi_null = phi_null,
    vara = vara,
    nb = nb,
    btest = btest,
    boot_pvals = boot_pvals,
    boot_tco = boot_tco,
    boot_seeds = boot_seeds
  )

  class(result_boot) <- "co_tas_boot"
  result_boot
}
