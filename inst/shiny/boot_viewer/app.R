# Bootstrap Simulation Viewer
# Shiny app for exploring DuckDB-stored simulation results

library(shiny)
library(DT)
library(duckdb)
library(DBI)

# =============================================================================
# UI
# =============================================================================

ui <- fluidPage(

  tags$head(
    tags$style(HTML("
      .sidebar { background-color: #f8f9fa; padding: 15px; border-radius: 5px; }
      .info-box { background-color: #e9ecef; padding: 10px; border-radius: 5px; margin-bottom: 10px; }
      .stat-label { font-weight: bold; color: #495057; }
      .nav-tabs { margin-bottom: 20px; }
      .dataTables_wrapper { font-size: 0.9em; }
      h4 { color: #343a40; margin-top: 0; }
    "))
  ),

  titlePanel("Bootstrap Simulation Viewer"),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      class = "sidebar",

      # Database connection
      div(class = "info-box",
        h4("Database"),
        textOutput("db_path_display"),
        actionButton("refresh_db", "Refresh", icon = icon("sync"),
                     class = "btn-sm btn-outline-secondary")
      ),

      hr(),

      # Filters
      selectizeInput("study_select", "Study (multi-select)", choices = NULL,
                     multiple = TRUE, options = list(placeholder = "All studies...")),
      selectInput("error_type_select", "Error Type",
                  choices = c("All" = "", "iid" = "iid", "GARCH" = "garch",
                              "Heteroscedastic" = "hetero", "Student-t" = "t_dist")),
      selectInput("dgp_select", "DGP Configuration", choices = NULL),
      selectInput("method_select", "Method", choices = NULL),
      selectInput("trial_select", "Trial", choices = NULL),

      br(),
      actionButton("clear_filters", "Clear All Filters",
                   icon = icon("times-circle"), class = "btn-sm btn-outline-secondary"),

      hr(),

      # Summary stats
      div(class = "info-box",
        h4("Selection Summary"),
        uiOutput("selection_summary")
      ),

      hr(),

      # Export
      downloadButton("export_csv", "Export Results (CSV)", class = "btn-sm")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "main_tabs",

        # Tab 1: Overview / Rejection Rates
        tabPanel(
          "Rejection Rates",
          icon = icon("chart-bar"),
          br(),
          h4("Rejection Rates Summary"),
          p("Aggregated rejection rates across all completed trials for the selected study."),
          DTOutput("rejection_rates_table"),
          br(),
          h4("Rejection Rates by Trial"),
          DTOutput("rejection_rates_by_trial_table")
        ),

        # Tab 2: Visualizations
        tabPanel(
          "Visualizations",
          icon = icon("chart-line"),
          br(),
          h4("Rejection Rate Analysis"),
          fluidRow(
            column(6, plotOutput("plot_reject_by_n", height = "350px")),
            column(6, plotOutput("plot_method_compare", height = "350px"))
          ),
          hr(),
          h4("Study Comparison"),
          p("Compare rejection rates across selected studies (select multiple studies above)."),
          fluidRow(
            column(12, plotOutput("plot_study_compare", height = "350px"))
          ),
          hr(),
          h4("P-value Distribution"),
          p("Under the null hypothesis, p-values should be uniformly distributed."),
          fluidRow(
            column(6, plotOutput("plot_pvalue_hist", height = "300px")),
            column(6, plotOutput("plot_pvalue_qq", height = "300px"))
          )
        ),

        # Tab 3: Cross-Study Comparison
        tabPanel(
          "Cross-Study",
          icon = icon("balance-scale"),
          br(),
          h4("Rejection Rates by Error Type"),
          p("Compare CO, COB, COBA rejection rates across error types for fixed (phi, n)."),
          fluidRow(
            column(4, selectInput("compare_phi", "Phi", choices = c(0.80, 0.95, 0.99))),
            column(4, selectInput("compare_n", "Sample Size", choices = c(50, 100, 250, 500)))
          ),
          DTOutput("comparison_table"),
          br(),
          h4("Visual Comparison"),
          plotOutput("comparison_plot", height = "400px")
        ),

        # Tab 4: Studies
        tabPanel(
          "Studies",
          icon = icon("folder"),
          br(),
          h4("All Studies"),
          DTOutput("studies_table")
        ),

        # Tab 4: DGP Configs
        tabPanel(
          "DGP Configs",
          icon = icon("cogs"),
          br(),
          h4("Data Generating Process Configurations"),
          DTOutput("dgp_table")
        ),

        # Tab 5: Method Configs
        tabPanel(
          "Methods",
          icon = icon("wrench"),
          br(),
          h4("Bootstrap Method Configurations"),
          DTOutput("method_table")
        ),

        # Tab 6: Trials
        tabPanel(
          "Trials",
          icon = icon("tasks"),
          br(),
          h4("Trials"),
          p("Execution batches within the selected study."),
          DTOutput("trials_table")
        ),

        # Tab 7: Runs
        tabPanel(
          "Runs",
          icon = icon("list"),
          br(),
          h4("Individual Runs"),
          p("Select a run to view its bootstrap distribution."),
          DTOutput("runs_table"),
          br(),
          conditionalPanel(
            condition = "input.runs_table_rows_selected.length > 0",
            h4("Bootstrap Distribution"),
            fluidRow(
              column(6, plotOutput("boot_dist_hist", height = "300px")),
              column(6, plotOutput("boot_dist_ecdf", height = "300px"))
            ),
            verbatimTextOutput("run_details")
          )
        ),

        # Tab 8: Custom Query
        tabPanel(
          "SQL Query",
          icon = icon("database"),
          br(),
          h4("Custom SQL Query"),
          p("Run custom queries against the database. Available tables: studies, dgp_configs, method_configs, trials, runs."),
          p("Views: v_run_summary, v_rejection_rates, v_rejection_rates_by_trial"),
          textAreaInput("custom_sql", NULL,
                        value = "SELECT * FROM v_rejection_rates LIMIT 100",
                        rows = 4, width = "100%"),
          actionButton("run_query", "Run Query", icon = icon("play"), class = "btn-primary"),
          br(), br(),
          DTOutput("custom_query_result")
        )
      )
    )
  )
)

# =============================================================================
# Server
# =============================================================================

server <- function(input, output, session) {

  # ---------------------------------------------------------------------------
  # Database Connection
  # ---------------------------------------------------------------------------

  # Get database path from options (set by launcher function)
  db_path <- reactiveVal(getOption("tstse.viewer_db", NULL))

  # Database connection (reactive)
  con <- reactive({
    path <- db_path()
    if (is.null(path) || !file.exists(path)) {
      return(NULL)
    }
    tryCatch(
      dbConnect(duckdb(), dbdir = path, read_only = TRUE),
      error = function(e) NULL
    )
  })

  # Clean up connection on session end
 session$onSessionEnded(function() {
    conn <- isolate(con())
    if (!is.null(conn)) {
      tryCatch(dbDisconnect(conn), error = function(e) NULL)
    }
  })

  output$db_path_display <- renderText({
    path <- db_path()
    if (is.null(path)) {
      "No database loaded"
    } else {
      basename(path)
    }
  })

  # Refresh trigger
  refresh_trigger <- reactiveVal(0)
  observeEvent(input$refresh_db, {
    refresh_trigger(refresh_trigger() + 1)
  })

  # Clear all filters
  observeEvent(input$clear_filters, {
    updateSelectizeInput(session, "study_select", selected = character(0))
    updateSelectInput(session, "error_type_select", selected = "")
    updateSelectInput(session, "dgp_select", selected = "")
    updateSelectInput(session, "method_select", selected = "")
    updateSelectInput(session, "trial_select", selected = "All")
  })

  # ---------------------------------------------------------------------------
  # Data Queries (reactive)
  # ---------------------------------------------------------------------------

  # Studies
  studies_data <- reactive({
    refresh_trigger()
    conn <- con()
    if (is.null(conn)) return(data.frame())
    tryCatch(
      dbGetQuery(conn, "SELECT * FROM studies ORDER BY created_at DESC"),
      error = function(e) data.frame()
    )
  })

  # DGP Configs
  dgp_data <- reactive({
    refresh_trigger()
    conn <- con()
    if (is.null(conn)) return(data.frame())
    tryCatch(
      dbGetQuery(conn, "
        SELECT dgp_id, dgp_name, n, ar_phi, ma_theta, vara,
               innov_dist, has_trend, trend_slope, created_at
        FROM dgp_configs
        ORDER BY created_at DESC
      "),
      error = function(e) data.frame()
    )
  })

  # Method Configs
  method_data <- reactive({
    refresh_trigger()
    conn <- con()
    if (is.null(conn)) return(data.frame())
    tryCatch(
      dbGetQuery(conn, "
        SELECT method_id, method_name, nb, maxp, ar_method,
               criterion, bootadj, stat_fn_name, garch_dist, created_at
        FROM method_configs
        ORDER BY created_at DESC
      "),
      error = function(e) data.frame()
    )
  })

  # Trials (filtered by study, DGP, and method - supports multiple studies)
  trials_data <- reactive({
    refresh_trigger()
    conn <- con()
    study_ids <- input$study_select  # Now a vector (multi-select)
    dgp_id <- input$dgp_select
    method_id <- input$method_select

    if (is.null(conn)) return(data.frame())

    # Build query with optional filters
    sql <- "
      SELECT t.trial_id, t.trial_name, t.n_planned, t.n_completed,
             t.status, t.elapsed_seconds, t.started_at, t.completed_at,
             d.dgp_name, m.method_name, s.study_name
      FROM trials t
      JOIN dgp_configs d ON t.dgp_id = d.dgp_id
      JOIN method_configs m ON t.method_id = m.method_id
      JOIN studies s ON t.study_id = s.study_id
      WHERE 1=1"
    params <- list()

    # Filter by study (supports multiple)
    if (length(study_ids) > 0 && !all(study_ids == "")) {
      placeholders <- paste(rep("?", length(study_ids)), collapse = ", ")
      sql <- paste(sql, sprintf("AND t.study_id IN (%s)", placeholders))
      params <- c(params, as.list(study_ids))
    }

    if (!is.null(dgp_id) && dgp_id != "") {
      sql <- paste(sql, "AND t.dgp_id = ?")
      params <- c(params, dgp_id)
    }
    if (!is.null(method_id) && method_id != "") {
      sql <- paste(sql, "AND t.method_id = ?")
      params <- c(params, method_id)
    }

    sql <- paste(sql, "ORDER BY s.study_name, t.started_at DESC")

    tryCatch(
      dbGetQuery(conn, sql, params = params),
      error = function(e) data.frame()
    )
  })

  # Runs (filtered by trial, or by study+DGP+method when trial="All")
  runs_data <- reactive({
    refresh_trigger()
    conn <- con()
    trial_id <- input$trial_select
    study_ids <- input$study_select  # Now a vector (multi-select)
    dgp_id <- input$dgp_select
    method_id <- input$method_select

    if (is.null(conn)) return(data.frame())

    # If specific trial selected, use that (no other filters needed)
    if (!is.null(trial_id) && trial_id != "" && trial_id != "All") {
      tryCatch(
        dbGetQuery(conn, "
          SELECT run_id, iteration_num, obs_stat, pvalue, pvalue_upper,
                 pvalue_lower, pvalue_asymp, pvalue_adj, master_seed
          FROM runs
          WHERE trial_id = ?
          ORDER BY iteration_num
          LIMIT 1000
        ", params = list(trial_id)),
        error = function(e) data.frame()
      )
    } else {
      # Trial is "All" - apply study, DGP, and method filters
      sql <- "
        SELECT r.run_id, r.iteration_num, r.obs_stat, r.pvalue,
               r.pvalue_upper, r.pvalue_lower, r.pvalue_asymp,
               r.pvalue_adj, r.master_seed, t.trial_name, s.study_name
        FROM runs r
        JOIN trials t ON r.trial_id = t.trial_id
        JOIN studies s ON t.study_id = s.study_id
        WHERE 1=1"
      params <- list()

      # Filter by study (supports multiple)
      if (length(study_ids) > 0 && !all(study_ids == "")) {
        placeholders <- paste(rep("?", length(study_ids)), collapse = ", ")
        sql <- paste(sql, sprintf("AND t.study_id IN (%s)", placeholders))
        params <- c(params, as.list(study_ids))
      }

      if (!is.null(dgp_id) && dgp_id != "") {
        sql <- paste(sql, "AND t.dgp_id = ?")
        params <- c(params, dgp_id)
      }
      if (!is.null(method_id) && method_id != "") {
        sql <- paste(sql, "AND t.method_id = ?")
        params <- c(params, method_id)
      }

      sql <- paste(sql, "ORDER BY s.study_name, t.started_at DESC, r.iteration_num LIMIT 1000")

      tryCatch(
        dbGetQuery(conn, sql, params = params),
        error = function(e) data.frame()
      )
    }
  })

  # Rejection rates (filtered by study, error_type, DGP, and method - supports multiple studies)
  rejection_rates_data <- reactive({
    refresh_trigger()
    conn <- con()
    study_ids <- input$study_select  # Now a vector (multi-select)
    error_type <- input$error_type_select
    dgp_id <- input$dgp_select
    method_id <- input$method_select

    if (is.null(conn)) return(data.frame())

    # Build query with filters
    sql <- "SELECT * FROM v_rejection_rates WHERE 1=1"
    params <- list()

    # Filter by study (supports multiple)
    if (length(study_ids) > 0 && !all(study_ids == "")) {
      studies <- studies_data()
      study_names <- studies$study_name[studies$study_id %in% study_ids]
      if (length(study_names) > 0) {
        placeholders <- paste(rep("?", length(study_names)), collapse = ", ")
        sql <- paste(sql, sprintf("AND study_name IN (%s)", placeholders))
        params <- c(params, as.list(study_names))
      }
    }

    # Filter by error_type
    if (!is.null(error_type) && error_type != "") {
      sql <- paste(sql, "AND error_type = ?")
      params <- c(params, error_type)
    }

    # Filter by DGP (need to get dgp_name from dgp_id)
    if (!is.null(dgp_id) && dgp_id != "") {
      dgp_info <- tryCatch(
        dbGetQuery(conn, "SELECT dgp_name FROM dgp_configs WHERE dgp_id = ?",
                   params = list(dgp_id)),
        error = function(e) data.frame()
      )
      if (nrow(dgp_info) > 0) {
        sql <- paste(sql, "AND dgp_name = ?")
        params <- c(params, dgp_info$dgp_name[1])
      }
    }

    # Filter by method (need to get method_name from method_id)
    if (!is.null(method_id) && method_id != "") {
      method_info <- tryCatch(
        dbGetQuery(conn, "SELECT method_name FROM method_configs WHERE method_id = ?",
                   params = list(method_id)),
        error = function(e) data.frame()
      )
      if (nrow(method_info) > 0) {
        sql <- paste(sql, "AND method_name = ?")
        params <- c(params, method_info$method_name[1])
      }
    }

    tryCatch(
      dbGetQuery(conn, sql, params = params),
      error = function(e) data.frame()
    )
  })

  # Rejection rates by trial (filtered by study, error_type, DGP, method - supports multiple studies)
  rejection_rates_by_trial_data <- reactive({
    refresh_trigger()
    conn <- con()
    study_ids <- input$study_select  # Now a vector (multi-select)
    error_type <- input$error_type_select
    dgp_id <- input$dgp_select
    method_id <- input$method_select

    if (is.null(conn)) return(data.frame())

    # Build query with filters
    sql <- "SELECT * FROM v_rejection_rates_by_trial WHERE 1=1"
    params <- list()

    # Filter by study (supports multiple)
    if (length(study_ids) > 0 && !all(study_ids == "")) {
      studies <- studies_data()
      study_names <- studies$study_name[studies$study_id %in% study_ids]
      if (length(study_names) > 0) {
        placeholders <- paste(rep("?", length(study_names)), collapse = ", ")
        sql <- paste(sql, sprintf("AND study_name IN (%s)", placeholders))
        params <- c(params, as.list(study_names))
      }
    }

    # Filter by error_type
    if (!is.null(error_type) && error_type != "") {
      sql <- paste(sql, "AND error_type = ?")
      params <- c(params, error_type)
    }

    # Filter by DGP
    if (!is.null(dgp_id) && dgp_id != "") {
      dgp_info <- tryCatch(
        dbGetQuery(conn, "SELECT dgp_name FROM dgp_configs WHERE dgp_id = ?",
                   params = list(dgp_id)),
        error = function(e) data.frame()
      )
      if (nrow(dgp_info) > 0) {
        sql <- paste(sql, "AND dgp_name = ?")
        params <- c(params, dgp_info$dgp_name[1])
      }
    }

    # Filter by method
    if (!is.null(method_id) && method_id != "") {
      method_info <- tryCatch(
        dbGetQuery(conn, "SELECT method_name FROM method_configs WHERE method_id = ?",
                   params = list(method_id)),
        error = function(e) data.frame()
      )
      if (nrow(method_info) > 0) {
        sql <- paste(sql, "AND method_name = ?")
        params <- c(params, method_info$method_name[1])
      }
    }

    tryCatch(
      dbGetQuery(conn, sql, params = params),
      error = function(e) data.frame()
    )
  })

  # ---------------------------------------------------------------------------
  # Update Filter Dropdowns
  # ---------------------------------------------------------------------------

  # Update study dropdown (multi-select)
  observe({
    studies <- studies_data()
    if (nrow(studies) == 0) {
      choices <- c("No studies found" = "")
    } else {
      choices <- setNames(studies$study_id, studies$study_name)
    }
    updateSelectizeInput(session, "study_select", choices = choices,
                         server = FALSE)
  })

  # Update DGP dropdown (filtered by selected studies - supports multiple)
  observe({
    refresh_trigger()
    conn <- con()
    study_ids <- input$study_select  # Now a vector

    if (is.null(conn) || length(study_ids) == 0 || all(study_ids == "")) {
      # No study selected - show all DGPs
      dgps <- dgp_data()
    } else {
      # Show DGPs used in any of the selected studies (union)
      placeholders <- paste(rep("?", length(study_ids)), collapse = ", ")
      dgps <- tryCatch(
        dbGetQuery(conn, sprintf("
          SELECT DISTINCT d.dgp_id, d.dgp_name
          FROM dgp_configs d
          JOIN trials t ON d.dgp_id = t.dgp_id
          WHERE t.study_id IN (%s)
          ORDER BY d.dgp_name
        ", placeholders), params = as.list(study_ids)),
        error = function(e) data.frame()
      )
    }

    if (nrow(dgps) == 0) {
      choices <- c("No DGPs found" = "")
    } else {
      choices <- c("All" = "", setNames(dgps$dgp_id, dgps$dgp_name))
    }
    updateSelectInput(session, "dgp_select", choices = choices)
  })

  # Update method dropdown (filtered by selected studies - supports multiple)
  observe({
    refresh_trigger()
    conn <- con()
    study_ids <- input$study_select  # Now a vector

    if (is.null(conn) || length(study_ids) == 0 || all(study_ids == "")) {
      # No study selected - show all methods
      methods <- method_data()
    } else {
      # Show methods used in any of the selected studies (union)
      placeholders <- paste(rep("?", length(study_ids)), collapse = ", ")
      methods <- tryCatch(
        dbGetQuery(conn, sprintf("
          SELECT DISTINCT m.method_id, m.method_name, m.nb
          FROM method_configs m
          JOIN trials t ON m.method_id = t.method_id
          WHERE t.study_id IN (%s)
          ORDER BY m.method_name
        ", placeholders), params = as.list(study_ids)),
        error = function(e) data.frame()
      )
    }

    if (nrow(methods) == 0) {
      choices <- c("No methods found" = "")
    } else {
      labels <- paste0(methods$method_name, " (nb=", methods$nb, ")")
      choices <- c("All" = "", setNames(methods$method_id, labels))
    }
    updateSelectInput(session, "method_select", choices = choices)
  })

  # Update trial dropdown based on selected study
  observe({
    trials <- trials_data()
    if (nrow(trials) == 0) {
      choices <- c("No trials found" = "")
    } else {
      # Use first 5 chars of trial_id as prefix for differentiation
      id_prefix <- substr(trials$trial_id, 1, 5)
      base_name <- ifelse(
        is.na(trials$trial_name) | trials$trial_name == "",
        paste0("Trial ", seq_len(nrow(trials))),
        trials$trial_name
      )
      labels <- paste0("(", id_prefix, ") ", base_name, " [", trials$status, "]")
      choices <- c("All" = "All", setNames(trials$trial_id, labels))
    }
    updateSelectInput(session, "trial_select", choices = choices)
  })

  # ---------------------------------------------------------------------------
  # Selection Summary
  # ---------------------------------------------------------------------------

  output$selection_summary <- renderUI({
    runs <- runs_data()
    trials <- trials_data()

    n_runs <- nrow(runs)
    n_trials <- nrow(trials)
    n_completed <- sum(trials$status == "completed", na.rm = TRUE)

    tagList(
      div(class = "stat-label", "Trials: ",
          span(paste0(n_completed, "/", n_trials, " completed"))),
      div(class = "stat-label", "Runs shown: ", span(n_runs)),
      if (n_runs >= 1000) div(em("(limited to 1000 rows)"))
    )
  })

  # ---------------------------------------------------------------------------
  # Render Tables
  # ---------------------------------------------------------------------------

  # Studies table
  output$studies_table <- renderDT({
    data <- studies_data()
    if (nrow(data) == 0) return(datatable(data.frame(Message = "No studies found")))

    datatable(
      data[, c("study_name", "study_type", "description", "created_at")],
      selection = "single",
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE
    )
  })

  # DGP table
  output$dgp_table <- renderDT({
    data <- dgp_data()
    if (nrow(data) == 0) return(datatable(data.frame(Message = "No DGP configs found")))

    # Format ar_phi for display
    data$ar_phi_str <- sapply(data$ar_phi, function(x) {
      if (is.null(x) || length(x) == 0) return("-")
      paste(round(x, 3), collapse = ", ")
    })

    datatable(
      data[, c("dgp_name", "n", "ar_phi_str", "innov_dist", "has_trend", "trend_slope")],
      selection = "single",
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE,
      colnames = c("Name", "n", "AR(phi)", "Innovation", "Trend?", "Slope")
    )
  })

  # Method table
  output$method_table <- renderDT({
    data <- method_data()
    if (nrow(data) == 0) return(datatable(data.frame(Message = "No method configs found")))

    datatable(
      data[, c("method_name", "nb", "maxp", "ar_method", "criterion", "bootadj")],
      selection = "single",
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE,
      colnames = c("Method", "B", "maxp", "AR Method", "Criterion", "COBA?")
    )
  })

  # Trials table
  output$trials_table <- renderDT({
    data <- trials_data()
    if (nrow(data) == 0) return(datatable(data.frame(Message = "No trials found")))

    display_data <- data[, c("trial_name", "dgp_name", "method_name",
                              "n_planned", "n_completed", "status", "elapsed_seconds")]
    display_data$elapsed_seconds <- round(display_data$elapsed_seconds, 1)

    datatable(
      display_data,
      selection = "single",
      options = list(pageLength = 15, scrollX = TRUE),
      rownames = FALSE,
      colnames = c("Trial", "DGP", "Method", "Planned", "Completed", "Status", "Time (s)")
    ) |>
      formatStyle("status",
        backgroundColor = styleEqual(
          c("completed", "running", "failed"),
          c("#d4edda", "#fff3cd", "#f8d7da")
        )
      )
  })

  # Runs table
  output$runs_table <- renderDT({
    data <- runs_data()
    if (nrow(data) == 0) return(datatable(data.frame(Message = "No runs found")))

    # Select columns to display
    cols <- intersect(
      c("iteration_num", "trial_name", "obs_stat", "pvalue", "pvalue_asymp", "pvalue_adj"),
      names(data)
    )
    display_data <- data[, cols, drop = FALSE]

    # Round numeric columns
    num_cols <- sapply(display_data, is.numeric)
    display_data[num_cols] <- lapply(display_data[num_cols], function(x) round(x, 4))

    datatable(
      display_data,
      selection = "single",
      options = list(pageLength = 20, scrollX = TRUE),
      rownames = FALSE
    ) |>
      formatStyle("pvalue",
        backgroundColor = styleInterval(c(0.05, 0.10), c("#f8d7da", "#fff3cd", "#ffffff"))
      )
  })

  # Rejection rates table
  output$rejection_rates_table <- renderDT({
    data <- rejection_rates_data()
    if (nrow(data) == 0) return(datatable(data.frame(Message = "No rejection rates available")))

    # Format ar_phi for display
    data$ar_phi_str <- sapply(data$ar_phi, function(x) {
      if (is.null(x) || length(x) == 0) return("-")
      paste(round(x, 3), collapse = ", ")
    })

    display_data <- data.frame(
      Study = data$study_name,
      DGP = data$dgp_name,
      n = data$n,
      `AR(phi)` = data$ar_phi_str,
      Method = data$method_name,
      `N Runs` = data$n_runs,
      `CO (0.05)` = sprintf("%.3f (%.3f)", data$reject_co_05, data$reject_co_05_se),
      `COB (0.05)` = sprintf("%.3f (%.3f)", data$reject_cob_05, data$reject_cob_05_se),
      `COBA (0.05)` = sprintf("%.3f (%.3f)", data$reject_coba_05, data$reject_coba_05_se),
      check.names = FALSE
    )

    datatable(
      display_data,
      selection = "none",
      options = list(pageLength = 20, scrollX = TRUE),
      rownames = FALSE
    )
  })

  # Rejection rates by trial table
  output$rejection_rates_by_trial_table <- renderDT({
    data <- rejection_rates_by_trial_data()
    if (nrow(data) == 0) return(datatable(data.frame(Message = "No data available")))

    display_data <- data.frame(
      Trial = data$trial_name,
      DGP = data$dgp_name,
      Method = data$method_name,
      `N Runs` = data$n_runs,
      `CO (0.05)` = round(data$reject_co_05, 3),
      `COB (0.05)` = round(data$reject_cob_05, 3),
      `COBA (0.05)` = round(data$reject_coba_05, 3),
      check.names = FALSE
    )

    datatable(
      display_data,
      selection = "none",
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE
    )
  })

  # ---------------------------------------------------------------------------
  # Bootstrap Distribution Plots
  # ---------------------------------------------------------------------------

  selected_run_data <- reactive({
    runs <- runs_data()
    sel <- input$runs_table_rows_selected
    if (is.null(sel) || length(sel) == 0 || nrow(runs) == 0) return(NULL)

    run_id <- runs$run_id[sel]
    conn <- con()
    if (is.null(conn)) return(NULL)

    tryCatch({
      result <- dbGetQuery(conn, "
        SELECT obs_stat, boot_dist, pvalue, pvalue_asymp
        FROM runs WHERE run_id = ?
      ", params = list(run_id))

      if (nrow(result) == 0) return(NULL)

      list(
        run_id = run_id,
        obs_stat = result$obs_stat[1],
        boot_dist = result$boot_dist[[1]],
        pvalue = result$pvalue[1],
        pvalue_asymp = result$pvalue_asymp[1]
      )
    }, error = function(e) NULL)
  })

  output$boot_dist_hist <- renderPlot({
    data <- selected_run_data()
    if (is.null(data) || is.null(data$boot_dist) || length(data$boot_dist) == 0) {
      plot.new()
      text(0.5, 0.5, "No bootstrap distribution stored", cex = 1.2)
      return()
    }

    hist(data$boot_dist, breaks = 30, col = "steelblue", border = "white",
         main = "Bootstrap Distribution", xlab = "Test Statistic",
         freq = FALSE)
    abline(v = data$obs_stat, col = "red", lwd = 2, lty = 2)
    abline(v = -data$obs_stat, col = "red", lwd = 2, lty = 2)
    legend("topright", legend = c(sprintf("Observed: %.3f", data$obs_stat)),
           col = "red", lty = 2, lwd = 2, bty = "n")
  })

  output$boot_dist_ecdf <- renderPlot({
    data <- selected_run_data()
    if (is.null(data) || is.null(data$boot_dist) || length(data$boot_dist) == 0) {
      plot.new()
      text(0.5, 0.5, "No bootstrap distribution stored", cex = 1.2)
      return()
    }

    plot(ecdf(data$boot_dist), main = "Empirical CDF",
         xlab = "Test Statistic", ylab = "F(x)")
    abline(v = data$obs_stat, col = "red", lwd = 2, lty = 2)
    abline(v = -data$obs_stat, col = "red", lwd = 2, lty = 2)
  })

  output$run_details <- renderPrint({
    data <- selected_run_data()
    if (is.null(data)) return(cat("Select a run to view details"))

    cat("Run ID:", data$run_id, "\n")
    cat("Observed Statistic:", round(data$obs_stat, 4), "\n")
    cat("Bootstrap P-value:", round(data$pvalue, 4), "\n")
    cat("Asymptotic P-value:", round(data$pvalue_asymp, 4), "\n")
    if (!is.null(data$boot_dist)) {
      cat("Bootstrap Replicates:", length(data$boot_dist), "\n")
      cat("Boot Dist Summary: min=", round(min(data$boot_dist), 3),
          ", median=", round(median(data$boot_dist), 3),
          ", max=", round(max(data$boot_dist), 3), "\n")
    }
  })

  # ---------------------------------------------------------------------------
  # Visualization Plots
  # ---------------------------------------------------------------------------

  # Plot 1: Rejection Rate vs Sample Size
  output$plot_reject_by_n <- renderPlot({
    data <- rejection_rates_data()
    if (nrow(data) == 0) {
      plot.new()
      text(0.5, 0.5, "No data available", cex = 1.2)
      return()
    }

    # Extract phi values from ar_phi column
    data$phi <- sapply(data$ar_phi, function(x) {
      if (is.null(x) || length(x) == 0) return(NA)
      x[1]  # First AR coefficient
    })

    # Filter to rows with valid n and phi
    data <- data[!is.na(data$n) & !is.na(data$phi), ]
    if (nrow(data) == 0) {
      plot.new()
      text(0.5, 0.5, "No AR(1) data available", cex = 1.2)
      return()
    }

    # Get unique phi values for colors
    phi_vals <- sort(unique(data$phi))
    colors <- rainbow(length(phi_vals), s = 0.7, v = 0.8)
    names(colors) <- as.character(phi_vals)

    # Set up plot
    op <- par(mar = c(4, 4, 3, 8), xpd = TRUE)
    on.exit(par(op))

    ylim_max <- max(c(data$reject_cob_05, 0.15), na.rm = TRUE)
    plot(NULL, xlim = range(data$n), ylim = c(0, ylim_max),
         xlab = "Sample Size (n)", ylab = "Rejection Rate",
         main = "COB Rejection Rate vs Sample Size")

    # Horizontal line at nominal level
    abline(h = 0.05, col = "gray50", lty = 2, lwd = 1.5)

    # Plot lines for each phi
    for (phi in phi_vals) {
      sub <- data[data$phi == phi, ]
      sub <- sub[order(sub$n), ]
      if (nrow(sub) > 1) {
        lines(sub$n, sub$reject_cob_05, col = colors[as.character(phi)], lwd = 2)
        points(sub$n, sub$reject_cob_05, col = colors[as.character(phi)], pch = 19)
      } else if (nrow(sub) == 1) {
        points(sub$n, sub$reject_cob_05, col = colors[as.character(phi)], pch = 19, cex = 1.5)
      }
    }

    # Legend
    legend("topright", inset = c(-0.25, 0),
           legend = paste0("phi=", phi_vals),
           col = colors, lwd = 2, pch = 19, bty = "n", cex = 0.9)
  })

  # Plot 2: Method Comparison Bar Chart
  output$plot_method_compare <- renderPlot({
    data <- rejection_rates_data()
    if (nrow(data) == 0) {
      plot.new()
      text(0.5, 0.5, "No data available", cex = 1.2)
      return()
    }

    # Reshape data for grouped bar chart
    # Aggregate by DGP name (combining across sample sizes if needed)
    dgp_names <- unique(data$dgp_name)
    if (length(dgp_names) > 10) dgp_names <- dgp_names[1:10]  # Limit for readability

    # Prepare matrix for barplot
    methods <- c("CO", "COB", "COBA")
    mat <- matrix(NA, nrow = 3, ncol = length(dgp_names),
                  dimnames = list(methods, dgp_names))

    for (i in seq_along(dgp_names)) {
      sub <- data[data$dgp_name == dgp_names[i], ]
      if (nrow(sub) > 0) {
        mat["CO", i] <- mean(sub$reject_co_05, na.rm = TRUE)
        mat["COB", i] <- mean(sub$reject_cob_05, na.rm = TRUE)
        mat["COBA", i] <- mean(sub$reject_coba_05, na.rm = TRUE)
      }
    }

    # Shorten DGP names for display
    short_names <- gsub("AR\\(1\\) phi=", "phi=", dgp_names)
    short_names <- gsub(" \\+ ARCH\\(8\\)", "+ARCH", short_names)
    short_names <- gsub(" \\+ Hetero", "+Het", short_names)
    short_names <- gsub(" n=", " n", short_names)
    colnames(mat) <- short_names

    op <- par(mar = c(7, 4, 3, 1))
    on.exit(par(op))

    barplot(mat, beside = TRUE, col = c("#e74c3c", "#3498db", "#2ecc71"),
            las = 2, ylab = "Rejection Rate",
            main = "Method Comparison by DGP",
            ylim = c(0, max(mat, na.rm = TRUE) * 1.2))
    abline(h = 0.05, col = "gray50", lty = 2, lwd = 1.5)
    legend("topright", legend = methods, fill = c("#e74c3c", "#3498db", "#2ecc71"),
           bty = "n", cex = 0.9)
  })

  # Plot 3: Study Comparison (for multi-study selection)
  output$plot_study_compare <- renderPlot({
    data <- rejection_rates_data()
    study_ids <- input$study_select

    if (nrow(data) == 0) {
      plot.new()
      text(0.5, 0.5, "No data available", cex = 1.2)
      return()
    }

    if (length(study_ids) < 2) {
      plot.new()
      text(0.5, 0.5, "Select 2+ studies above to compare", cex = 1.2)
      return()
    }

    # Get unique studies and DGPs
    study_names <- unique(data$study_name)
    dgp_names <- unique(data$dgp_name)

    # Limit DGPs for readability
    if (length(dgp_names) > 6) dgp_names <- dgp_names[1:6]

    # Prepare matrix: rows = studies, cols = DGPs
    mat <- matrix(NA, nrow = length(study_names), ncol = length(dgp_names),
                  dimnames = list(study_names, dgp_names))

    for (i in seq_along(study_names)) {
      for (j in seq_along(dgp_names)) {
        sub <- data[data$study_name == study_names[i] & data$dgp_name == dgp_names[j], ]
        if (nrow(sub) > 0) {
          mat[i, j] <- mean(sub$reject_cob_05, na.rm = TRUE)
        }
      }
    }

    # Shorten names for display
    short_dgp <- gsub("AR\\(1\\) phi=", "phi=", dgp_names)
    short_dgp <- gsub(" \\+ ARCH\\(8\\)", "+ARCH", short_dgp)
    short_dgp <- gsub(" \\+ Hetero", "+Het", short_dgp)
    short_dgp <- gsub(" n=", " n", short_dgp)
    colnames(mat) <- short_dgp

    short_study <- gsub("Scenario [0-9]+: ", "", study_names)
    rownames(mat) <- short_study

    op <- par(mar = c(7, 4, 3, 8), xpd = TRUE)
    on.exit(par(op))

    # Use different colors for each study
    colors <- rainbow(length(study_names), s = 0.7, v = 0.8)

    barplot(mat, beside = TRUE, col = colors,
            las = 2, ylab = "COB Rejection Rate",
            main = "Cross-Study Comparison",
            ylim = c(0, max(mat, na.rm = TRUE) * 1.2))
    abline(h = 0.05, col = "gray50", lty = 2, lwd = 1.5)
    legend("topright", inset = c(-0.15, 0),
           legend = rownames(mat), fill = colors,
           bty = "n", cex = 0.8, title = "Study")
  })

  # Plot 4: P-value Histogram (now Plot 4 after adding study comparison)
  output$plot_pvalue_hist <- renderPlot({
    data <- runs_data()
    if (nrow(data) == 0 || !"pvalue" %in% names(data)) {
      plot.new()
      text(0.5, 0.5, "No p-values available", cex = 1.2)
      return()
    }

    pvals <- data$pvalue[!is.na(data$pvalue)]
    if (length(pvals) == 0) {
      plot.new()
      text(0.5, 0.5, "No p-values available", cex = 1.2)
      return()
    }

    hist(pvals, breaks = 20, col = "steelblue", border = "white",
         main = "Bootstrap P-value Distribution",
         xlab = "P-value", freq = FALSE, xlim = c(0, 1))
    abline(h = 1, col = "red", lty = 2, lwd = 2)  # Uniform reference
    legend("topright", legend = "Uniform(0,1)", col = "red", lty = 2, lwd = 2, bty = "n")
  })

  # Plot 4: P-value QQ Plot
  output$plot_pvalue_qq <- renderPlot({
    data <- runs_data()
    if (nrow(data) == 0 || !"pvalue" %in% names(data)) {
      plot.new()
      text(0.5, 0.5, "No p-values available", cex = 1.2)
      return()
    }

    pvals <- data$pvalue[!is.na(data$pvalue)]
    if (length(pvals) < 2) {
      plot.new()
      text(0.5, 0.5, "Need at least 2 p-values for QQ plot", cex = 1.2)
      return()
    }

    # QQ plot against Uniform(0,1)
    n <- length(pvals)
    theoretical <- (1:n - 0.5) / n
    observed <- sort(pvals)

    plot(theoretical, observed, pch = 19, col = "steelblue",
         xlab = "Theoretical Quantiles (Uniform)",
         ylab = "Observed P-values",
         main = "P-value QQ Plot",
         xlim = c(0, 1), ylim = c(0, 1))
    abline(0, 1, col = "red", lwd = 2, lty = 2)
    legend("bottomright", legend = "y = x", col = "red", lty = 2, lwd = 2, bty = "n")
  })

  # ---------------------------------------------------------------------------
  # Cross-Study Comparison
  # ---------------------------------------------------------------------------

  comparison_data <- reactive({
    conn <- con()
    phi <- input$compare_phi
    n <- input$compare_n

    if (is.null(conn)) return(data.frame())

    # Query rejection rates grouped by error_type for fixed phi and n
    tryCatch({
      dbGetQuery(conn, "
        SELECT error_type,
               ROUND(reject_co_05, 4) AS CO,
               ROUND(reject_cob_05, 4) AS COB,
               ROUND(reject_coba_05, 4) AS COBA,
               n_runs
        FROM v_rejection_rates
        WHERE phi_value = ? AND n = ?
        ORDER BY error_type
      ", params = list(as.numeric(phi), as.integer(n)))
    }, error = function(e) data.frame())
  })

  output$comparison_table <- renderDT({
    data <- comparison_data()
    if (nrow(data) == 0) {
      return(datatable(data.frame(Message = "No data available. Run migration first.")))
    }

    datatable(data,
              options = list(dom = 't', pageLength = 10),
              rownames = FALSE) |>
      formatStyle('COB',
                  backgroundColor = styleInterval(c(0.04, 0.06),
                                                  c('#d4edda', '#fff3cd', '#f8d7da')))
  })

  output$comparison_plot <- renderPlot({
    data <- comparison_data()
    if (nrow(data) == 0 || !all(c("CO", "COB", "COBA") %in% names(data))) {
      plot.new()
      text(0.5, 0.5, "No data available.\nRun boot_db_migrate_error_type() first.", cex = 1.2)
      return()
    }

    mat <- as.matrix(data[, c("CO", "COB", "COBA")])
    rownames(mat) <- data$error_type

    op <- par(mar = c(7, 4, 3, 2))
    on.exit(par(op))

    barplot(t(mat), beside = TRUE, col = c("#e74c3c", "#3498db", "#2ecc71"),
            las = 2, ylab = "Rejection Rate",
            main = sprintf("phi = %s, n = %s", input$compare_phi, input$compare_n),
            ylim = c(0, max(mat, na.rm = TRUE) * 1.3))
    abline(h = 0.05, lty = 2, col = "gray50", lwd = 2)
    legend("topright", legend = c("CO", "COB", "COBA"),
           fill = c("#e74c3c", "#3498db", "#2ecc71"), bty = "n")
  })

  # ---------------------------------------------------------------------------
  # Custom SQL Query
  # ---------------------------------------------------------------------------

  custom_query_result <- eventReactive(input$run_query, {
    conn <- con()
    sql <- input$custom_sql
    if (is.null(conn) || is.null(sql) || sql == "") return(data.frame())

    tryCatch(
      dbGetQuery(conn, sql),
      error = function(e) {
        data.frame(Error = conditionMessage(e))
      }
    )
  })

  output$custom_query_result <- renderDT({
    data <- custom_query_result()
    datatable(data, options = list(pageLength = 20, scrollX = TRUE), rownames = FALSE)
  })

  # ---------------------------------------------------------------------------
  # Export
  # ---------------------------------------------------------------------------

  output$export_csv <- downloadHandler(
    filename = function() {
      paste0("simulation_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      data <- runs_data()
      write.csv(data, file, row.names = FALSE)
    }
  )
}

# =============================================================================
# Run App
# =============================================================================

shinyApp(ui = ui, server = server)
