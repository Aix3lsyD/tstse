#' Launch Capstone Simulation App
#'
#' Opens the interactive Shiny app for running and exploring Monte Carlo
#' simulations in-session.
#'
#' Required packages (all in Suggests): shiny, DT, ggplot2.
#'
#' @return Launches the Shiny app (does not return until the app is closed).
#'
#' @examples
#' \dontrun{
#' capstone_app()
#' }
#'
#' @export
capstone_app <- function() {
  pkgs <- c("shiny", "DT", "ggplot2")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Required packages not installed: ", paste(missing, collapse = ", "),
      "\nInstall with: install.packages(c(",
      paste0('"', missing, '"', collapse = ", "), "))",
      call. = FALSE
    )
  }

  # Prefer local source app when running from the project directory so the
  # newest module edits are picked up immediately.
  source_app_dir <- file.path(getwd(), "inst", "shiny", "boot_viewer")
  if (dir.exists(source_app_dir)) {
    app_dir <- source_app_dir
  } else {
    app_dir <- system.file("shiny", "boot_viewer", package = "tstse")
    if (app_dir == "") {
      # Fallback for devtools::load_all() or development
      app_dir <- file.path(
        system.file(package = "tstse"),
        "..", "..", "inst", "shiny", "boot_viewer"
      )
    }
  }
  if (!dir.exists(app_dir)) {
    stop("Cannot find boot_viewer app directory", call. = FALSE)
  }

  shiny::runApp(app_dir)
}
