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


#' Check if currently running inside a future worker
#'
#' @return TRUE if inside a worker, FALSE otherwise
#' @noRd
is_inside_worker <- function() {
  # Method 1: Check for future's worker flag
  if (!is.null(future::futureSessionInfo()$worker)) {
    return(TRUE)
  }


  # Method 2: Check if we're in a worker via the plan

  # Workers have a "future.nbrOfWorkers" that's set to 1
  inherits(future::plan(), "uniprocess") &&
    isTRUE(getOption("future.fork.enable", FALSE) == FALSE)

  FALSE
}


#' Parallel map using future
#'
#' Uses future.apply::future_lapply. Backend selection:
#' - If tstse is installed: uses multisession (socket-based, safe)
#' - If using devtools::load_all(): uses multicore on Unix (fork-based)
#'
#' Workers are kept alive between calls to avoid spawn/teardown overhead.
#' Call `tstse_parallel_stop()` to explicitly shut down workers.
#'
#' @param X Vector or list to iterate over
#' @param FUN Function to apply
#' @param cores Number of cores to use (1 = sequential)
#' @return List of results
#' @noRd
pmap <- function(X, FUN, cores = 1L) {

  verbose <- isTRUE(getOption("tstse.parallel_verbose", FALSE))

  # Sequential execution
  if (cores <= 1L) {
    return(lapply(X, FUN))
  }

  # =========================================================================
  # CRITICAL: Prevent nested parallelization

  # If we're already inside a future worker, force sequential execution.
  # This prevents the 14 workers x 14 workers = 196 process explosion.
  # =========================================================================
  if (!is.null(future::futureSessionInfo()$worker)) {
    if (verbose) {
      message("[tstse] Inside worker process - forcing sequential execution")
    }
    return(lapply(X, FUN))
  }

  # Choose backend:
  #   options(tstse.parallel_backend = "multisession")  # socket-based, safe (default)
  #   options(tstse.parallel_backend = "multicore")     # fork-based, for devtools::load_all()
  backend <- getOption("tstse.parallel_backend", "multisession")

  # Multicore only works on Unix
  if (backend == "multicore" && .Platform$OS.type != "unix") {
    warning("[tstse] multicore backend not available on Windows, using multisession",
            call. = FALSE)
    backend <- "multisession"
  }

  # =========================================================================
  # ROBUST PLAN REUSE: Check worker count instead of class inheritance
  # The class hierarchy check (inherits(plan, "multisession")) is unreliable
  # because the plan object's class varies. Use nbrOfWorkers() instead.
  # =========================================================================
  current_workers <- future::nbrOfWorkers()
  need_new_plan <- current_workers < cores

  if (!need_new_plan) {
    if (verbose) {
      message("[tstse] Reusing existing plan (", current_workers, " workers)")
    }
  } else {
    if (verbose) {
      message("[tstse] Setting up ", backend, " plan with ", cores, " workers")
    }

    if (backend == "multisession") {
      future::plan(future::multisession, workers = cores)
    } else {
      future::plan(future::multicore, workers = cores)
    }

    # Mark that tstse set up this plan (for cleanup tracking)
    options(tstse.parallel_active = TRUE)
  }

  # Run the parallel operation
  # Use explicit chunking: one chunk per worker for minimal IPC overhead
  chunk_size <- max(1L, ceiling(length(X) / cores))

  if (backend == "multisession") {
    future.apply::future_lapply(X, FUN,
                                future.seed = TRUE,
                                future.chunk.size = chunk_size,
                                future.packages = "tstse")
  } else {
    future.apply::future_lapply(X, FUN,
                                future.seed = TRUE,
                                future.chunk.size = chunk_size)
  }
}


#' Stop parallel workers
#'
#' Shuts down any parallel workers that were started by tstse.
#' Call this when you're done with parallel processing to free resources.
#'
#' @return Invisible NULL
#' @export
#' @examples
#' \dontrun{
#' options(tstse.cores = 4)
#' result <- wbg_boot(x, nb = 100)  # Uses parallel workers
#' tstse_parallel_stop()            # Clean up workers
#' }
tstse_parallel_stop <- function() {
  if (isTRUE(getOption("tstse.parallel_active", FALSE))) {
    future::plan(future::sequential)
    options(tstse.parallel_active = FALSE)
    message("[tstse] Parallel workers stopped")
  }
  invisible(NULL)
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
