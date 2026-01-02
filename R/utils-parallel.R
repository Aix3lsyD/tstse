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
#' Includes safeguards against nested parallelization and BLAS thread explosion.
#'
#' @param X Vector or list to iterate over
#' @param FUN Function to apply
#' @param cores Number of cores to use
#' @return List of results
#' @noRd
pmap <- function(X, FUN, cores = 1L) {

  verbose <- isTRUE(getOption("tstse.parallel_verbose", FALSE))

  # Prevent nested parallelization (e.g., bootstrap calling aic_ar)
  if (isTRUE(getOption("tstse.parallel_active"))) {
    if (verbose) message("[tstse] Nesting guard: forcing sequential execution")
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
    # Control BLAS threading to prevent thread explosion in forked children
    # On macOS, vecLib (Accelerate) is multi-threaded by default
    blas_threads <- as.integer(getOption("tstse.blas_threads", 1L))

    # Prefer RhpcBLASctl (works post-initialization, unlike env vars)
    old_blas_threads <- NULL
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      old_blas_threads <- RhpcBLASctl::blas_get_num_procs()
      RhpcBLASctl::blas_set_num_threads(blas_threads)
      on.exit(RhpcBLASctl::blas_set_num_threads(old_blas_threads), add = TRUE)
      if (verbose) {
        message("[tstse] Set BLAS threads to ", blas_threads,
                " via RhpcBLASctl (was ", old_blas_threads, ")")
      }
    } else {
      # Fallback to environment variables (only affects BLAS initialization,
      # less effective if BLAS already initialized, but doesn't hurt)
      blas_vars <- c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS",
                     "VECLIB_MAXIMUM_THREADS", "MKL_NUM_THREADS")
      old_env <- Sys.getenv(blas_vars, unset = NA)

      Sys.setenv(
        OPENBLAS_NUM_THREADS = blas_threads,
        OMP_NUM_THREADS = blas_threads,
        VECLIB_MAXIMUM_THREADS = blas_threads,
        MKL_NUM_THREADS = blas_threads
      )
      on.exit({
        for (i in seq_along(blas_vars)) {
          if (is.na(old_env[i])) {
            Sys.unsetenv(blas_vars[i])
          } else {
            do.call(Sys.setenv, stats::setNames(list(old_env[i]), blas_vars[i]))
          }
        }
      }, add = TRUE)

      if (verbose) {
        message("[tstse] Set BLAS env vars to ", blas_threads,
                " (RhpcBLASctl not available - install for better thread control)")
      }
    }

    if (verbose) {
      message("[tstse] Starting mclapply with ", cores, " cores")
    }

    # Unix: use mclapply (fork-based)
    parallel::mclapply(X, FUN, mc.cores = cores)
  } else {
    # Windows: use parLapply (socket-based, requires cluster setup)
    if (verbose) {
      message("[tstse] Starting parLapply with ", cores, " workers")
    }

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
