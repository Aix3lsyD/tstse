#' Display Top 5 AR Orders by Information Criterion
#'
#' @inheritParams aic_ar
#' @return Invisibly returns a data frame of all orders sorted by criterion.
#' @export
#'
#' @examples
#' x <- gen_arma(n = 200, phi = c(1.5, -0.75), plot = FALSE, seed = 123)
#' aic5_ar(x, p = 0:5)
#'
#' \dontrun{
#' # Parallel (uses multiple cores)
#' aic5_ar(x, p = 0:8, cores = 2)
#' }
aic5_ar <- function(x,
                    p = 0:5,
                    type = c("aic", "aicc", "bic"),
                    method = c("mle", "burg", "yw"),
                    cores = 1L) {

  type <- match.arg(type)
  method <- match.arg(method)

  cat("---------WORKING... PLEASE WAIT...\n\n\n")

  # Use aic_ar with parallel support
  capture.output({
    result <- aic_ar(x, p = p, type = type, method = method, cores = cores)
  })

  # Sort table by selected criterion
  tbl <- result$table
  tbl <- tbl[order(tbl[[type]]), ]

  # Prepare display table
  display_tbl <- data.frame(
    p = tbl$p,
    value = tbl[[type]]
  )
  colnames(display_tbl) <- c("   p", paste0("       ", type))

  cat("Five Smallest Values of ", type, "\n", sep = "")
  cat("Method=", method, "\n", sep = "")

  n_show <- min(5, nrow(display_tbl))
  print(display_tbl[1:n_show, ], row.names = FALSE)

  invisible(tbl)
}
