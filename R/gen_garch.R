#' Generate GARCH Process
#'
#' Generate a realization from a GARCH(p0, q0) process.
#'
#' @param n Integer, length of output series.
#' @param alpha0 Numeric, constant term (must be positive).
#' @param alpha Numeric vector, ARCH coefficients (q0 terms).
#' @param beta Numeric vector, GARCH coefficients (p0 terms).
#' @param spin Integer, burn-in period (default 1000).
#' @param plot Logical, whether to plot the series (default TRUE).
#' @param seed Integer, random seed for reproducibility (default NULL).
#'
#' @return Numeric vector of length n.
#' @export
#'
#' @details
#' Generates from the model:
#' \deqn{\sigma^2_t = \alpha_0 + \sum_{i=1}^{q_0} \alpha_i a_{t-i}^2 + \sum_{j=1}^{p_0} \beta_j \sigma^2_{t-j}}
#' \deqn{a_t = \epsilon_t \sigma_t}
#' where \eqn{\epsilon_t \sim N(0,1)}.
#'
#' @examples
#' # GARCH(1,1)
#' x <- gen_garch(n = 500, alpha0 = 0.1, alpha = 0.3, beta = 0.5)
#'
#' # GARCH(2,1)
#' x <- gen_garch(n = 500, alpha0 = 0.1, alpha = c(0.2, 0.1), beta = 0.6)
gen_garch <- function(n,
                      alpha0,
                      alpha,
                      beta,
                      spin = 1000L,
                      plot = TRUE,
                      seed = NULL) {

  # Input validation
  if (!is.numeric(n) || length(n) != 1 || n < 1) {
    stop("`n` must be a positive integer")
  }
  if (!is.numeric(alpha0) || length(alpha0) != 1 || alpha0 <= 0) {
    stop("`alpha0` must be a positive number")
  }
  if (!is.numeric(alpha) || length(alpha) == 0) {
    stop("`alpha` must be a non-empty numeric vector")
  }
  if (!is.numeric(beta) || length(beta) == 0) {
    stop("`beta` must be a non-empty numeric vector")
  }
  if (!is.numeric(spin) || spin < 0) {
    stop("`spin` must be a non-negative integer")
  }

  # Call Rcpp implementation
  result <- gen_garch_cpp(
    n = as.integer(n),
    alpha0 = as.double(alpha0),
    alpha = as.double(alpha),
    beta = as.double(beta),
    spin = as.integer(spin),
    seed = seed
  )

  if (plot) {
    plot(result, type = "l", xlab = "Time", ylab = "",
         main = "GARCH Process")
    invisible(result)
  } else {
    result
  }
}
