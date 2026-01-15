#' Launch Bootstrap Simulation Viewer
#'
#' Opens an interactive Shiny application for exploring bootstrap simulation
#' results stored in a DuckDB database.
#'
#' @param db_path Character. Path to the DuckDB database file. If NULL,
#'   opens the app without a database loaded (user can browse but won't see data).
#' @param launch.browser Logical. If TRUE (default), opens the app in a browser.
#'   If FALSE, returns the app URL for manual navigation.
#' @param port Integer. Port number for the Shiny app. Default is random.
#' @param host Character. Host address. Default is "127.0.0.1" (localhost only).
#'
#' @return Invisibly returns NULL. The function launches a Shiny app and
#'   blocks until the app is closed.
#'
#' @details
#' The viewer provides:
#' \itemize{
#'   \item Study/DGP/Method/Trial navigation and filtering
#'   \item Rejection rate summaries with standard errors
#'   \item Individual run inspection with bootstrap distribution plots
#'   \item Custom SQL query interface
#'   \item CSV export of results
#' }
#'
#' Requires the `shiny`, `DT`, `duckdb`, and `DBI` packages to be installed.
#'
#' @examples
#' \dontrun{
#' # View results from a simulation study
#' boot_db_viewer("my_simulations.duckdb")
#'
#' # Or after running simulations:
#' study <- boot_study("results.duckdb", "My Study",
#'                     dgp = list(n = 100, ar_phi = 0.7),
#'                     method = list(name = "wbg_boot", nb = 399))
#' # ... run simulations ...
#' study$complete()
#' study$end()
#'
#' # Now view the results
#' boot_db_viewer("results.duckdb")
#' }
#'
#' @seealso [boot_study()], [boot_db_connect()], [boot_db_query()]
#' @export
boot_db_viewer <- function(db_path = NULL,
                           launch.browser = TRUE,
                           port = NULL,
                           host = "127.0.0.1") {

 # Check for required packages
  required_pkgs <- c("shiny", "DT", "duckdb", "DBI")
  missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing) > 0) {
    stop(
      "The following packages are required for the viewer:\n",
      paste0("  - ", missing, collapse = "\n"), "\n",
      "Install with: install.packages(c(", paste0('"', missing, '"', collapse = ", "), "))",
      call. = FALSE
    )
  }

  # Find the app directory
  app_dir <- system.file("shiny/boot_viewer", package = "tstse")
  if (app_dir == "") {
    stop(
      "Shiny app not found. This can happen if:\n",
      "  1. The package was not properly installed\n",
      "  2. You're using devtools::load_all() and inst/ wasn't copied\n",
      "Try reinstalling with: devtools::install() or remotes::install_github(...)",
      call. = FALSE
    )
  }

  # Validate database path if provided
  if (!is.null(db_path)) {
    if (!file.exists(db_path)) {
      stop("Database file not found: ", db_path, call. = FALSE)
    }
    # Normalize to absolute path
    db_path <- normalizePath(db_path, mustWork = TRUE)
  }

  # Set the database path option (read by the app)
  old_opt <- getOption("tstse.viewer_db")
  options(tstse.viewer_db = db_path)

  # Restore option on exit
 on.exit(options(tstse.viewer_db = old_opt), add = TRUE)

  # Build runApp arguments
  run_args <- list(
    appDir = app_dir,
    launch.browser = launch.browser,
    host = host
  )
  if (!is.null(port)) {
    run_args$port <- port
  }

  # Launch the app
  do.call(shiny::runApp, run_args)

  invisible(NULL)
}
