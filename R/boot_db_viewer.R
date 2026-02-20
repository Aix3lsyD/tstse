#' Launch Bootstrap Simulation Viewer
#'
#' Opens an interactive Shiny application for exploring DuckDB-stored
#' Monte Carlo bootstrap simulation results. Provides rejection rate tables
#' with color coding and flexible ggplot2 visualizations.
#'
#' @param db_path Path to a DuckDB database file. If omitted, falls back to
#'   `getOption("tstse.viewer_db")`.
#'
#' @details
#' The viewer is read-only and does not modify the database. It queries the
#' `v_rejection_rates` and `v_rejection_rates_by_batch` views from the flat
#' two-table schema (batches + simulations). Supports pooled and per-batch
#' rejection rates, power curve plots, heatmaps, and bootstrap distribution
#' histograms.
#'
#' Required packages (all in Suggests): shiny, DT, ggplot2, duckdb, DBI.
#'
#' @return Launches the Shiny app (does not return until the app is closed).
#'
#' @examples
#' \dontrun{
#' boot_db_viewer("path/to/study.duckdb")
#'
#' # Or set the option globally
#' options(tstse.viewer_db = "path/to/study.duckdb")
#' boot_db_viewer()
#' }
#'
#' @export
boot_db_viewer <- function(db_path = NULL) {
  if (is.null(db_path)) {
    db_path <- getOption("tstse.viewer_db")
  }
  if (is.null(db_path)) {
    stop(
      "No database path provided. Pass db_path or set options(tstse.viewer_db = ...)",
      call. = FALSE
    )
  }
  # Normalize to absolute path so the Shiny app can find it regardless of working directory
  db_path <- normalizePath(db_path, mustWork = FALSE)
  if (!file.exists(db_path)) {
    stop("Database file not found: ", db_path, call. = FALSE)
  }

  pkgs <- c("shiny", "DT", "ggplot2", "duckdb", "DBI")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Required packages not installed: ", paste(missing, collapse = ", "),
      "\nInstall with: install.packages(c(",
      paste0('"', missing, '"', collapse = ", "), "))",
      call. = FALSE
    )
  }

  # Pass db_path to the app via option
  old_opt <- getOption("tstse.viewer_db")
  options(tstse.viewer_db = db_path)
  on.exit(options(tstse.viewer_db = old_opt), add = TRUE)

  app_dir <- system.file("shiny", "boot_viewer", package = "tstse")
  if (app_dir == "") {
    # Fallback for devtools::load_all() or development
    app_dir <- file.path(
      system.file(package = "tstse"),
      "..", "..", "inst", "shiny", "boot_viewer"
    )
    if (!dir.exists(app_dir)) {
      # Direct path for load_all from source tree
      app_dir <- file.path(getwd(), "inst", "shiny", "boot_viewer")
    }
  }
  if (!dir.exists(app_dir)) {
    stop("Cannot find boot_viewer app directory", call. = FALSE)
  }

  shiny::runApp(app_dir)
}
