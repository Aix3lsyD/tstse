# Monte Carlo Simulation Viewer
# Interactive Shiny app for exploring DuckDB-stored bootstrap rejection rates
#
# Module files in R/ are auto-sourced by shiny::loadSupport():
#   utils.R, mod_overview.R, mod_rejection_rates.R, mod_plots.R,
#   mod_pvalue.R, mod_bootstrap_dist.R, mod_analysis_grids.R,
#   mod_parallel_coords.R, mod_diagnostics.R, mod_innovation_comp.R,
#   mod_adhoc_sim.R, mod_benchmark.R

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
    input_dark_mode(id = "dark_mode", mode = "dark")
  ),

  tabsetPanel(
    id = "main_tabs",
    mod_overview_ui("overview"),
    mod_rejection_rates_ui("rr"),
    mod_plots_ui("plots"),
    mod_pvalue_ui("pval"),
    mod_bootstrap_dist_ui("bootdist"),
    mod_analysis_grids_ui("grids"),
    mod_parallel_coords_ui("parcoord"),
    mod_diagnostics_ui("diag"),
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

  # Database connection
  db_path <- getOption("tstse.viewer_db")
  con <- dbConnect(duckdb(), dbdir = db_path, read_only = TRUE)
  onStop(function() dbDisconnect(con, shutdown = TRUE))

  # Populate initial filter choices (shared across modules)
  init_choices <- tryCatch({
    n_vals <- dbGetQuery(con, "SELECT DISTINCT n FROM simulations ORDER BY n")
    phi_vals <- dbGetQuery(con, "SELECT DISTINCT phi FROM simulations ORDER BY phi")
    innov_vals <- dbGetQuery(con, "SELECT DISTINCT innov_dist FROM simulations ORDER BY innov_dist")
    list(
      n = as.character(n_vals$n),
      phi = as.character(phi_vals$phi),
      innov = innov_vals$innov_dist
    )
  }, error = function(e) list(n = character(0), phi = character(0), innov = character(0)))

  # Reactive values shared across modules
  rv <- reactiveValues()
  observe({ rv$dark_mode <- input$dark_mode })

  # Module servers
  mod_overview_server("overview", con)
  mod_rejection_rates_server("rr", con, init_choices)
  mod_plots_server("plots", con, init_choices)
  mod_pvalue_server("pval", con, init_choices)
  mod_bootstrap_dist_server("bootdist", con)
  mod_analysis_grids_server("grids", con, init_choices)
  mod_parallel_coords_server("parcoord", con, init_choices)
  mod_diagnostics_server("diag", con, init_choices)
  mod_innovation_comp_server("innovcomp", con, init_choices)
  mod_benchmark_server("bench", con)
  mod_adhoc_sim_server("adhoc", con, db_path)
}

shinyApp(ui = ui, server = server)
