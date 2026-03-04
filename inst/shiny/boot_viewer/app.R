# Monte Carlo Simulation Viewer
# Interactive Shiny app for running and exploring bootstrap simulations
#
# Module files in R/ are auto-sourced by shiny::loadSupport():
#   utils.R, mod_capstone.R (+ mod_capstone_*.R sub-tabs),
#   mod_innovation_comp.R, mod_benchmark.R, mod_adhoc_sim.R

library(shiny)
library(bslib)
library(thematic)
library(DT)
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

      /* Dark-mode fixes for shinyWidgets pickerInput */
      [data-bs-theme="dark"] .bootstrap-select > .dropdown-toggle {
        background-color: var(--bs-body-bg) !important;
        color: var(--bs-body-color) !important;
        border-color: var(--bs-border-color) !important;
      }
      [data-bs-theme="dark"] .bootstrap-select .filter-option,
      [data-bs-theme="dark"] .bootstrap-select .filter-option-inner,
      [data-bs-theme="dark"] .bootstrap-select .filter-option-inner-inner {
        color: var(--bs-body-color) !important;
      }
      [data-bs-theme="dark"] .bootstrap-select .dropdown-menu {
        background-color: var(--bs-body-bg) !important;
        border-color: var(--bs-border-color) !important;
      }
      [data-bs-theme="dark"] .bootstrap-select .dropdown-item {
        color: var(--bs-body-color) !important;
      }
      [data-bs-theme="dark"] .bootstrap-select .dropdown-item:hover,
      [data-bs-theme="dark"] .bootstrap-select .dropdown-item:focus,
      [data-bs-theme="dark"] .bootstrap-select .dropdown-item.active {
        background-color: var(--bs-secondary-bg) !important;
        color: var(--bs-body-color) !important;
      }

      /* Dark-mode fixes for shinyWidgets actionBttn "default" style */
      [data-bs-theme="dark"] .btn.btn-default,
      [data-bs-theme="dark"] .btn-material.btn-default,
      [data-bs-theme="dark"] .bttn-default {
        background-color: var(--bs-secondary-bg) !important;
        color: var(--bs-body-color) !important;
        border-color: var(--bs-border-color) !important;
      }
      [data-bs-theme="dark"] .btn.btn-default:hover,
      [data-bs-theme="dark"] .btn-material.btn-default:hover,
      [data-bs-theme="dark"] .bttn-default:hover,
      [data-bs-theme="dark"] .btn.btn-default:focus,
      [data-bs-theme="dark"] .btn-material.btn-default:focus,
      [data-bs-theme="dark"] .bttn-default:focus {
        background-color: var(--bs-tertiary-bg) !important;
        color: var(--bs-body-color) !important;
        border-color: var(--bs-border-color) !important;
      }
    '))
  ),

  div(
    class = "d-flex justify-content-between align-items-center mb-3",
    h2("Monte Carlo Simulation Viewer", style = "margin: 0;"),
    div(
      class = "d-flex align-items-center gap-3",
      actionButton(
        "app_reset", "Reset App",
        icon = icon("rotate-left"),
        class = "btn-outline-danger"
      ),
      input_dark_mode(id = "dark_mode", mode = "dark")
    )
  ),

  tabsetPanel(
    id = "main_tabs",
    mod_capstone_ui("cap"),
    mod_capstone_sanity_ui("sanity"),
    mod_innovation_comp_ui("innovcomp"),
    mod_benchmark_ui("bench"),
    mod_adhoc_sim_ui("adhoc"),
    mod_console_ui("console")
  )
)

# =============================================================================
# Server
# =============================================================================

server <- function(input, output, session) {
  thematic_shiny()

  observeEvent(input$app_reset, {
    showModal(modalDialog(
      title = "Reset App State?",
      "This will reload the app and clear all in-memory simulation results and console state.",
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_app_reset", "Reset", class = "btn-danger")
      ),
      easyClose = TRUE
    ))
  })

  observeEvent(input$confirm_app_reset, {
    removeModal()
    session$reload()
  })

  # Shared filter choices are populated from in-memory simulation results.
  init_choices <- reactive({
    list(n = character(0), phi = character(0), innov = character(0))
  })

  # Reactive values shared across modules
  rv <- reactiveValues()
  observe({ rv$dark_mode <- input$dark_mode })

  # --- Module servers ---
  mod_capstone_server("cap")
  mod_capstone_sanity_server("sanity")
  mod_innovation_comp_server("innovcomp")
  mod_benchmark_server("bench")
  mod_adhoc_sim_server("adhoc")
  mod_console_server("console")
}

shinyApp(ui = ui, server = server)
