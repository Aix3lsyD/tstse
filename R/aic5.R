#' Display Top 5 ARMA Orders by Information Criterion
#'
#' Fit ARMA models across a grid of orders and display the 5 best by
#' the chosen information criterion.
#'
#' @param x Numeric vector, the time series.
#' @param p Integer vector, AR orders to compare (default 0:5).
#' @param q Integer vector, MA orders to compare (default 0:2).
#' @param type Character, criterion to use: "aic" (default), "aicc", or "bic".
#' @param cores Integer, number of cores for parallel processing.
#'   Default NULL uses \code{getOption("tstse.cores", 1)}.
#'   Set to 0 to use all available cores.
#'
#' @return Invisibly returns a data frame of all (p,q) combinations sorted by criterion.
#' @export
#'
#' @details
#' This is a convenience function that calls \code{\link{aic_ts}} and
#' displays the top 5 models. For programmatic use, \code{\link{aic_ts}}
#' is preferred as it returns the full results.
#'
#' @seealso \code{\link{aic_ts}} for the underlying function.
#'
#' @examples
#' \donttest{
#' x <- gen_arma(n = 200, phi = 0.7, theta = 0.4, plot = FALSE, seed = 123)
#' aic5(x, p = 0:5, q = 0:2)
#' }
#'
#' \dontrun{
#' # Parallel (uses multiple cores)
#' aic5(x, p = 0:5, q = 0:3, cores = 2)
#' }
aic5 <- function(x,
                 p = 0:5,
                 q = 0:2,
                 type = c("aic", "aicc", "bic"),
                 cores = NULL) {

  type <- match.arg(type)

  cat("---------WORKING... PLEASE WAIT...\n\n\n")

  # Use aic_ts with parallel support
  result <- tryCatch(
    aic_ts(x, p = p, q = q, type = type, cores = cores),
    error = function(e) NULL
  )

  if (is.null(result)) {
    cat("Error: All models failed to converge\n")
    return(invisible(NULL))
  }

  # Sort table by selected criterion
  tbl <- result$table
  tbl <- tbl[order(tbl[[type]]), ]

  # Prepare display table
  display_tbl <- data.frame(
    p = tbl$p,
    q = tbl$q,
    value = tbl[[type]]
  )
  colnames(display_tbl) <- c("   p", "   q", paste0("       ", type))

  cat("Five Smallest Values of ", type, "\n", sep = "")

  n_show <- min(5, nrow(display_tbl))
  print(display_tbl[1:n_show, ], row.names = FALSE)

  invisible(tbl)
}
