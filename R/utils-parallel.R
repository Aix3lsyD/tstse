# R/utils-parallel.R
# Platform-aware parallelization: mclapply (Unix/Mac) / parLapply (Windows)

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


#' Parallel map with platform-aware backend
#'
#' Uses mclapply on Unix/Mac (fork-based, preferred) and parLapply on Windows
#' (socket-based with automatic closure variable export).
#'
#' @param X Vector or list to iterate over.
#' @param FUN Function to apply to each element of X.
#' @param cores Number of cores to use (1 = sequential).
#' @param ... Additional arguments passed to FUN (currently unused).
#' @return List of results.
#' @noRd
pmap <- function(X, FUN, cores = 1L, ...) {
  # Sequential execution
 if (cores <= 1L) {
    return(lapply(X, FUN))
  }

  if (.Platform$OS.type == "unix") {
    # =========================================================================
    # Unix/Mac: fork-based parallelization (preferred)
    # - Fast: no serialization, shared memory
    # - Automatic environment capture
    # =========================================================================
    parallel::mclapply(X, FUN, mc.cores = cores)

  } else {
    # =========================================================================
    # Windows: socket-based cluster with auto-export
    # - Requires explicit variable export
    # - Workers need tstse package loaded
    # =========================================================================
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    # Load tstse on workers
    parallel::clusterEvalQ(cl, {
      if (requireNamespace("tstse", quietly = TRUE)) {
        library(tstse)
      }
    })

    # Auto-export closure variables
    # FUN is typically a closure that references variables from its enclosing env
    fn_env <- environment(FUN)
    if (!identical(fn_env, globalenv()) && !identical(fn_env, baseenv())) {
      # Get variables referenced by FUN
      fn_globals <- tryCatch(
        codetools::findGlobals(FUN, merge = FALSE)$variables,
        error = function(e) character(0)
      )

      # Filter to those that exist in the closure environment
      env_vars <- ls(envir = fn_env, all.names = TRUE)
      to_export <- intersect(fn_globals, env_vars)

      if (length(to_export) > 0) {
        parallel::clusterExport(cl, to_export, envir = fn_env)
      }
    }

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
