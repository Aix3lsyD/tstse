# R/utils-parallel.R

#' Get number of cores to use
#'
#' @param cores User-specified cores (NULL = use option, 0 = all, N = use N)
#' @return Integer number of cores
#' @noRd
get_cores <- function(cores = NULL) {

  # If not specified, check option
  if (is.null(cores)) {
    cores <- getOption("tstse.cores", default = 1L)
  }

  # 0 means all available
  if (cores == 0L) {
    cores <- parallel::detectCores(logical = FALSE)
    if (is.na(cores)) cores <- 1L
  }

  # Cap at available cores
  max_cores <- parallel::detectCores(logical = FALSE)
  if (is.na(max_cores)) max_cores <- 1L
  cores <- min(cores, max_cores)

  as.integer(max(cores, 1L))
}


#' Parallel-aware map function
#'
#' Uses mclapply on Unix, parLapply on Windows, or lapply if cores = 1.
#'
#' @param X Vector or list to iterate over
#' @param FUN Function to apply
#' @param cores Number of cores to use
#' @return List of results
#' @noRd
pmap <- function(X, FUN, cores = 1L) {

  # Prevent nested parallelization (e.g., bootstrap calling aic_ar)
  if (isTRUE(getOption("tstse.parallel_active"))) {
    cores <- 1L
  }

  if (cores <= 1L) {
    return(lapply(X, FUN))
  }

  # Set flag before spawning workers to prevent nesting
  old_opt <- getOption("tstse.parallel_active")
  options(tstse.parallel_active = TRUE)
  on.exit(options(tstse.parallel_active = old_opt), add = TRUE)

  if (.Platform$OS.type == "unix") {
    # Unix: use mclapply (fork-based, simple)
    parallel::mclapply(X, FUN, mc.cores = cores)
  } else {
    # Windows: use parLapply (socket-based, requires cluster setup)
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    # Export necessary functions to workers
    parallel::clusterExport(cl, varlist = c("backcast", "est_ar", "est_arma"),
                            envir = asNamespace("tstse"))

    parallel::parLapply(cl, X, FUN)
  }
}


#' Parallel-aware map for grid (two-variable iteration)
#'
#' @param grid Data frame with columns to iterate over
#' @param FUN Function taking a single row index
#' @param cores Number of cores
#' @return List of results
#' @noRd
pmap_grid <- function(grid, FUN, cores = 1L) {
  pmap(seq_len(nrow(grid)), FUN, cores = cores)
}
