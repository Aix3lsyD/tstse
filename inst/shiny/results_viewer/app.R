library(shiny)
library(bslib)
library(thematic)
library(DBI)
library(duckdb)
library(DT)
library(ggplot2)

resolve_db_path <- function() {
  env_path <- Sys.getenv("TSTSE_VIEWER_DB", unset = "")
  candidates <- c(
    env_path,
    "simulations_lite.duckdb",
    file.path("..", "..", "..", "simulations_lite.duckdb")
  )
  candidates <- unique(candidates[nzchar(candidates)])
  existing <- candidates[file.exists(candidates)]

  if (length(existing) == 0) {
    stop(
      "Could not find a DuckDB file. Set TSTSE_VIEWER_DB or place simulations_lite.duckdb next to app.R.",
      call. = FALSE
    )
  }

  normalizePath(existing[[1]], winslash = "/", mustWork = TRUE)
}

resolve_modules_dir <- function() {
  candidates <- c(
    file.path("..", "boot_viewer", "R"),
    file.path("inst", "shiny", "boot_viewer", "R")
  )
  existing <- candidates[dir.exists(candidates)]
  if (length(existing) == 0) {
    stop("Could not locate boot_viewer module directory.", call. = FALSE)
  }
  normalizePath(existing[[1]], winslash = "/", mustWork = TRUE)
}

source_boot_modules <- function(modules_dir) {
  files <- c(
    "utils.R",
    "mod_capstone_overview.R",
    "mod_capstone_grids.R",
    "mod_capstone_pvalue.R",
    "mod_capstone_null_diag.R"
  )

  for (f in files) {
    source(file.path(modules_dir, f), local = FALSE)
  }
}

innov_label_from_dist <- function(x) {
  if (grepl("^norm$", x)) return("Normal")
  if (grepl("^t\\(", x)) return("Student's t")
  if (grepl("^arch\\(|^garch\\(", x)) return("GARCH")
  if (grepl("^hetero", x)) return("Heteroscedastic")
  if (grepl("^laplace", x)) return("Laplace")
  x
}

load_snapshot <- function(db_path) {
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = db_path, read_only = TRUE)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  rates <- DBI::dbGetQuery(con, paste(
    "SELECT n, phi, innov_dist, n_sims, n_batches,",
    "reject_05, reject_05_se, reject_asymp_05, reject_asymp_05_se,",
    "reject_adj_05, reject_adj_05_se",
    "FROM v_rejection_rates ORDER BY innov_dist, n, phi"
  ))

  sims <- DBI::dbGetQuery(con, paste(
    "SELECT sim_id, iteration, n, phi, innov_dist, obs_stat,",
    "pvalue, pvalue_asymp, pvalue_adj, null_ar_order, null_vara, created_at",
    "FROM simulations ORDER BY innov_dist, n, phi, iteration"
  ))

  rates$innov_label <- vapply(rates$innov_dist, innov_label_from_dist, character(1))

  grid <- unique(data.frame(
    phi = rates$phi,
    n = rates$n,
    innov_label = rates$innov_label,
    innov_dist_str = rates$innov_dist,
    stringsAsFactors = FALSE
  ))
  grid <- grid[order(grid$innov_label, grid$n, grid$phi), , drop = FALSE]

  if (nrow(sims) == 0) {
    raw_map <- list()
  } else {
    sims$cell_key <- paste(sims$n, sims$phi, sims$innov_dist, sep = "|")
    split_rows <- split(sims, sims$cell_key)

    make_sim <- function(row) {
      list(
        tco_obs = as.numeric(row$obs_stat),
        obs_stat = as.numeric(row$obs_stat),
        pvalue = as.numeric(row$pvalue),
        pvalue_asymp = as.numeric(row$pvalue_asymp),
        pvalue_adj = if (is.na(row$pvalue_adj)) NULL else as.numeric(row$pvalue_adj),
        p = as.integer(row$null_ar_order),
        vara = as.numeric(row$null_vara),
        phi = numeric(0)
      )
    }

    raw_map <- lapply(split_rows, function(df) {
      results <- lapply(seq_len(nrow(df)), function(i) make_sim(df[i, , drop = FALSE]))
      list(
        results = results,
        n = as.integer(df$n[1]),
        phi = as.numeric(df$phi[1]),
        innov_dist = as.character(df$innov_dist[1]),
        innov_label = innov_label_from_dist(as.character(df$innov_dist[1])),
        innov_params = list(),
        timestamp = as.character(max(df$created_at, na.rm = TRUE))
      )
    })
  }

  list(grid = grid, rates = rates, raw_map = raw_map)
}

mod_capstone_readonly_ui <- function(id) {
  ns <- NS(id)

  tabPanel(
    "Capstone Simulation",
    br(),
    div(
      id = ns("cap_content"),
      p(
        class = "text-body-secondary",
        "Read-only view loaded from DuckDB snapshot. No simulations are run in this hosted app."
      ),
      tabsetPanel(
        id = ns("cap_tabs"),
        mod_capstone_overview_ui(ns),
        mod_capstone_grids_ui(ns),
        mod_capstone_pvalue_ui(ns),
        mod_capstone_null_diag_ui(ns)
      )
    )
  )
}

mod_capstone_readonly_server <- function(id, snapshot) {
  moduleServer(id, function(input, output, session) {
    grid_rv <- reactiveVal(snapshot$grid)
    run_results_rv <- reactiveVal(snapshot$rates)
    raw_results_rv <- reactiveVal(snapshot$raw_map)

    cap_combined_results <- reactive({
      run_results_rv()
    })

    cap_sim_data <- function(cell_key) {
      if (is.null(cell_key) || !nzchar(cell_key)) return(NULL)
      raw <- raw_results_rv()
      if (cell_key %in% names(raw)) return(raw[[cell_key]]$results)
      NULL
    }

    cap_cell_choices <- reactive({
      raw <- raw_results_rv()
      keys <- names(raw)
      if (length(keys) == 0) return(character(0))
      labels <- gsub("\\|", ", ", keys)
      setNames(keys, labels)
    })

    mod_capstone_overview_server(input, output, session,
      grid_rv = grid_rv,
      raw_results_rv = raw_results_rv
    )

    mod_capstone_grids_server(input, output, session,
      cap_combined_results = cap_combined_results,
      run_results_rv = run_results_rv
    )

    mod_capstone_pvalue_server(input, output, session,
      cap_sim_data = cap_sim_data,
      cap_cell_choices = cap_cell_choices
    )

    mod_capstone_null_diag_server(input, output, session,
      cap_sim_data = cap_sim_data,
      cap_cell_choices = cap_cell_choices
    )
  })
}

modules_dir <- resolve_modules_dir()
source_boot_modules(modules_dir)

# Read-only UI overrides: remove detailed rows and CSV download controls.
mod_capstone_overview_ui <- function(ns) {
  tabPanel("Coverage / Overview",
    br(),
    uiOutput(ns("overview_kpis")),
    p(class = "text-body-secondary",
      "Shows how many simulations exist in the loaded snapshot for each grid cell."),
    plotOutput(ns("coverage_heatmap"), height = "450px")
  )
}

mod_capstone_grids_ui <- function(ns) {
  tabPanel("Analysis Grids",
    br(),
    h4("Rejection Rate Comparison Grids"),
    p(class = "text-body-secondary",
      "Each grid fixes two dimensions and varies the third. ",
      "Rates are color-coded: ",
      span(style = "background-color:#fff3cd; color:#212529; font-weight:600; padding:2px 6px; border-radius:3px;", "< 0.03"),
      " ",
      span(style = "background-color:#d4edda; color:#212529; font-weight:600; padding:2px 6px; border-radius:3px;", "0.03 - 0.07"),
      " ",
      span(style = "background-color:#f8d7da; color:#212529; font-weight:600; padding:2px 6px; border-radius:3px;", "> 0.07")
    ),
    fluidRow(
      column(3, actionButton(ns("clear_results"), "Clear Results",
                             icon = icon("eraser"),
                             class = "btn-sm btn-outline-warning"))
    ),
    hr(),
    wellPanel(
      h5("By Sample Size"),
      p(class = "text-body-secondary", style = "margin-bottom: 10px;",
        "Fix innovation distribution and phi; rows vary by n"),
      fluidRow(
        column(4, selectInput(ns("res_grid1_innov"), "Innovation Distribution:", choices = NULL)),
        column(4, selectInput(ns("res_grid1_phi"), "Phi:", choices = NULL))
      ),
      DTOutput(ns("res_grid1_table"))
    ),
    wellPanel(
      h5("By Phi"),
      p(class = "text-body-secondary", style = "margin-bottom: 10px;",
        "Fix innovation distribution and n; rows vary by phi"),
      fluidRow(
        column(4, selectInput(ns("res_grid2_innov"), "Innovation Distribution:", choices = NULL)),
        column(4, selectInput(ns("res_grid2_n"), "Sample Size (n):", choices = NULL))
      ),
      DTOutput(ns("res_grid2_table"))
    ),
    wellPanel(
      h5("By Innovation Distribution"),
      p(class = "text-body-secondary", style = "margin-bottom: 10px;",
        "Fix phi and n; rows vary by innovation distribution"),
      fluidRow(
        column(4, selectInput(ns("res_grid3_phi"), "Phi:", choices = NULL)),
        column(4, selectInput(ns("res_grid3_n"), "Sample Size (n):", choices = NULL))
      ),
      plotOutput(ns("res_grid3_barchart"), height = "300px"),
      hr(),
      DTOutput(ns("res_grid3_table"))
    )
  )
}

db_path <- resolve_db_path()
snapshot <- load_snapshot(db_path)

ui <- page_fluid(
  theme = bs_theme(),
  tags$head(
    tags$style(HTML('
      .nav-tabs { margin-bottom: 20px; }
      .dataTables_wrapper { font-size: 0.9em; }
      .plot-controls { background-color: var(--bs-tertiary-bg); padding: 15px; border-radius: 5px; margin-bottom: 15px; }
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
      [id$="clear_results"] { display: none !important; }
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
  p(class = "text-body-secondary", sprintf("Data source: %s", db_path)),
  tabsetPanel(
    id = "main_tabs",
    mod_capstone_readonly_ui("cap")
  )
)

server <- function(input, output, session) {
  thematic_shiny()

  observeEvent(input$app_reset, {
    showModal(modalDialog(
      title = "Reset App State?",
      "This will reload the app and reset local filters/selections.",
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

  mod_capstone_readonly_server("cap", snapshot = snapshot)
}

shinyApp(ui = ui, server = server)
