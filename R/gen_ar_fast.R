#' Fast AR Series Generator
#'
#' Generates AR(p) realizations optimized for bootstrap applications.
#' Uses adaptive burn-in and the C-implemented `stats::filter()` for
#' efficient recursive filtering. Provides ~30-50x speedup over general
#' ARMA generators for pure AR series.
#'
#' @param n Integer. Length of output series.
#' @param phi Numeric vector. AR coefficients (same sign convention as
#'   `gen_arma()` and `arima.sim()`).
#' @param vara Numeric. White noise variance. Default is 1.
#' @param seed Integer or NULL. Random seed for reproducibility. If NULL
#'   (default), no seed is set.
#'
#' @return Numeric vector of length `n`.
#'
#' @details
#' This function is optimized for generating many AR series quickly,
#' as needed in bootstrap applications. Key optimizations include:
#'
#' \itemize{
#'   \item **Adaptive burn-in**: Adjusts based on AR persistence rather
#'     than using a fixed value. Near-unit-root models get longer burn-in.
#'   \item **C-level filtering**: Uses `stats::filter()` which is
#'     implemented in C for recursive filtering.
#'   \item **No plotting overhead**: Unlike `gen_arma()`, no plot is
#'     generated.
#'   \item **Pure AR only**: No MA/ARMA overhead.
#' }
#'
#' The adaptive burn-in strategy:
#' \itemize{
#'   \item If `sum(|phi|) >= 0.999` (near unit root): burn-in is
#'     `max(50, n/2, 500)`, capped at 2000
#'   \item If `sum(|phi|) < 0.01` (nearly white noise): burn-in is 50
#'   \item Otherwise: burn-in is based on decay rate, `3 * |log(0.001)/log(persistence)|`
#' }
#'
#' @seealso [gen_arma()] for general ARMA generation with plotting,
#'   [wbg_boot_flex()] which uses this function for bootstrap resampling.
#'
#' @examples
#' # Generate AR(1) with phi = 0.9
#' set.seed(123)
#' x <- gen_ar_fast(200, phi = 0.9)
#' plot(x, type = "l")
#'
#' # Generate AR(2)
#' x <- gen_ar_fast(200, phi = c(0.5, 0.3), seed = 456)
#'
#' # White noise (phi = 0 or empty)
#' wn <- gen_ar_fast(100, phi = numeric(0))
#'
#' # Near unit root - gets longer burn-in automatically
#' x_persistent <- gen_ar_fast(200, phi = 0.99, seed = 789)
#'
#' @export
gen_ar_fast <- function(n, phi = numeric(0), vara = 1, seed = NULL) {

 # Input validation
  if (!is.numeric(n) || length(n) != 1 || n < 1) {
    stop("`n` must be a positive integer", call. = FALSE)
  }
  if (!is.numeric(phi)) {
    stop("`phi` must be a numeric vector", call. = FALSE)
  }
  if (!is.numeric(vara) || length(vara) != 1 || vara <= 0) {
    stop("`vara` must be a positive number", call. = FALSE)
  }

  n <- as.integer(n)

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  p <- length(phi)
  sd_a <- sqrt(vara)

  # White noise case - no filtering needed
  if (p == 0 || all(phi == 0)) {
    return(rnorm(n, mean = 0, sd = sd_a))
  }

  # Calculate adaptive burn-in based on AR persistence
  # Persistence = sum of absolute AR coefficients
  persistence <- sum(abs(phi))

  if (persistence >= 0.999) {
    # Near unit root: need long burn-in
    burn <- max(50L, round(n * 0.5), 500L)
  } else if (persistence < 0.01) {
    # Nearly white noise: minimal burn-in
    burn <- 50L
  } else {
    # Normal case: calculate based on decay rate
    # We want persistence^k < 0.001 (effect decayed to < 0.1%)
    # k > log(0.001) / log(persistence)
    burn <- max(50L, round(3 * abs(log(0.001) / log(persistence))))
  }

  # Cap burn-in at reasonable maximum
  burn <- min(burn, 2000L)

  n_total <- n + burn

  # Generate white noise innovations
  a <- rnorm(n_total, mean = 0, sd = sd_a)

  # Use stats::filter for AR recursion (C implementation)
  # x[t] = a[t] + phi[1]*x[t-1] + phi[2]*x[t-2] + ...
  x <- stats::filter(a, filter = phi, method = "recursive")

  # Remove burn-in and return as numeric vector
  as.numeric(x[(burn + 1L):n_total])
}
