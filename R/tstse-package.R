#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel RcppParallelLibs
#' @useDynLib tstse, .registration = TRUE
## usethis namespace: end
NULL

# Base R imports ---------------------------------------------------------------

#' @importFrom grDevices gray hcl
#' @importFrom graphics abline axis image layout legend lines mtext par points segments plot
#' @importFrom stats ARMAacf ARMAtoMA HoltWinters acf approx ar.burg ar.yw arima arima.sim coef convolve fft is.ts lm na.pass pchisq predict qnorm residuals rnorm spec.pgram time ts var vcov
#' @importFrom utils capture.output
NULL

# ggplot2 .data pronoun --------------------------------------------------------

#' @importFrom rlang .data
NULL
