# Monte Carlo Simulation Viewer
# Interactive Shiny app for exploring DuckDB-stored bootstrap rejection rates
#
# Module files in R/ are auto-sourced by shiny::loadSupport():
#   utils.R, mod_capstone.R (+ mod_capstone_*.R sub-tabs),
#   mod_innovation_comp.R, mod_benchmark.R, mod_adhoc_sim.R

library(shiny)
library(bslib)
library(thematic)
library(DT)
library(duckdb)
library(DBI)
library(ggplot2)

# thematic_shiny() is called inside the server function so it can access session

# =============================================================================
# UI
# =============================================================================

ui <- page_fluid(
  theme = bs_theme(),
  tags$head(
    tags$style(HTML('
      .nav-tabs { margin-bottom: 20px; }
      .dataTables_wrapper { font-size: 0.9em; }
      .plot-controls { background-color: var(--bs-tertiary-bg); padding: 15px; border-radius: 5px; margin-bottom: 15px; }
      /* Dark-mode overrides for DT formatStyle traffic-light cells */
      [data-bs-theme="dark"] .dataTable td[style*="background-color: rgb(255, 243, 205)"] {
        background-color: #664d03 !important;
        color: #fff3cd !important;
      }
      [data-bs-theme="dark"] .dataTable td[style*="background-color: rgb(212, 237, 218)"] {
        background-color: #1a4731 !important;
        color: #d4edda !important;
      }
      [data-bs-theme="dark"] .dataTable td[style*="background-color: rgb(248, 215, 218)"] {
        background-color: #6c2022 !important;
        color: #f8d7da !important;
      }
    '))
  ),

  div(
    class = "d-flex justify-content-between align-items-center mb-3",
    h2("Monte Carlo Simulation Viewer", style = "margin: 0;"),
    div(
      class = "d-flex align-items-center gap-3",
      input_dark_mode(id = "dark_mode", mode = "dark")
    )
  ),

  tabsetPanel(
    id = "main_tabs",
    mod_capstone_ui("cap"),
    mod_innovation_comp_ui("innovcomp"),
    mod_benchmark_ui("bench"),
    mod_adhoc_sim_ui("adhoc")
  )
)

# =============================================================================
# Server
# =============================================================================

server <- function(input, output, session) {
  thematic_shiny()

  # Database connection (read-only, NULL if no db_path provided)
  db_path <- getOption("tstse.viewer_db")
  has_db <- !is.null(db_path) && nzchar(db_path) && file.exists(db_path)

  if (has_db) {
    con <- dbConnect(duckdb(), dbdir = db_path, read_only = TRUE)
    onStop(function() dbDisconnect(con, shutdown = TRUE))
  } else {
    con <- NULL
  }

  # Trigger to refresh data after DB writes (capstone/ad-hoc modules)
  db_refresh_trigger <- reactiveVal(0L)

  # Populate filter choices (shared across modules)
  init_choices <- reactive({
    db_refresh_trigger()
    if (!has_db) return(list(n = character(0), phi = character(0), innov = character(0)))
    tryCatch({
      n_vals <- dbGetQuery(con, "SELECT DISTINCT n FROM simulations ORDER BY n")
      phi_vals <- dbGetQuery(con, "SELECT DISTINCT phi FROM simulations ORDER BY phi")
      innov_vals <- dbGetQuery(con, "SELECT DISTINCT innov_dist FROM simulations ORDER BY innov_dist")
      list(
        n = as.character(n_vals$n),
        phi = as.character(phi_vals$phi),
        innov = innov_vals$innov_dist
      )
    }, error = function(e) list(n = character(0), phi = character(0), innov = character(0)))
  })

  # Reactive values shared across modules
  rv <- reactiveValues()
  observe({ rv$dark_mode <- input$dark_mode })

  # --- Module servers ---
  mod_capstone_server("cap", con, db_path, db_refresh_trigger, init_choices)
  mod_innovation_comp_server("innovcomp", con, init_choices)
  mod_benchmark_server("bench", con)
  mod_adhoc_sim_server("adhoc", con, db_path)
}

shinyApp(ui = ui, server = server)
