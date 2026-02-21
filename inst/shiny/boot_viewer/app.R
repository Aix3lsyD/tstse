# Monte Carlo Simulation Viewer
# Interactive Shiny app for exploring DuckDB-stored bootstrap rejection rates

library(shiny)
library(DT)
library(duckdb)
library(DBI)
library(ggplot2)

# =============================================================================
# UI
# =============================================================================

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      .sidebar { background-color: #f8f9fa; padding: 15px; border-radius: 5px; }
      .nav-tabs { margin-bottom: 20px; }
      .dataTables_wrapper { font-size: 0.9em; }
      .plot-controls { background-color: #f8f9fa; padding: 15px; border-radius: 5px; margin-bottom: 15px; }
    "))
  ),

  titlePanel("Monte Carlo Simulation Viewer"),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      class = "sidebar",

      h4("Filters"),

      selectizeInput("n_filter", "Sample Size (n)",
                     choices = NULL, multiple = TRUE,
                     options = list(placeholder = "All...")),

      selectizeInput("phi_filter", "Phi Value",
                     choices = NULL, multiple = TRUE,
                     options = list(placeholder = "All...")),

      selectizeInput("innov_filter", "Innovation Distribution",
                     choices = NULL, multiple = TRUE,
                     options = list(placeholder = "All...")),

      hr(),

      checkboxInput("compare_batches", "Compare Batches", value = FALSE),

      br(),
      actionButton("clear_filters", "Clear All Filters",
                   icon = icon("times-circle"), class = "btn-sm btn-outline-secondary"),

      hr(),
      div(style = "font-size: 0.85em; color: #6c757d;",
        textOutput("status_text")
      )
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "main_tabs",

        # =====================================================================
        # Tab 1: Study Overview
        # =====================================================================
        tabPanel(
          "Study Overview",
          br(),
          h4("Study Overview"),
          fluidRow(
            column(3,
              wellPanel(style = "text-align: center;",
                h6(style = "color: #6c757d; margin-bottom: 5px;", "Total Simulations"),
                h3(style = "margin: 0;", textOutput("overview_total_sims", inline = TRUE))
              )
            ),
            column(3,
              wellPanel(style = "text-align: center;",
                h6(style = "color: #6c757d; margin-bottom: 5px;", "Configurations"),
                h3(style = "margin: 0;", textOutput("overview_n_configs", inline = TRUE))
              )
            ),
            column(3,
              wellPanel(style = "text-align: center;",
                h6(style = "color: #6c757d; margin-bottom: 5px;", "Batches"),
                h3(style = "margin: 0;", textOutput("overview_n_batches", inline = TRUE))
              )
            ),
            column(3,
              wellPanel(style = "text-align: center;",
                h6(style = "color: #6c757d; margin-bottom: 5px;", "Date Range"),
                div(style = "font-size: 0.95em;",
                  textOutput("overview_date_range", inline = TRUE)
                )
              )
            )
          ),
          hr(),
          h5("Coverage Matrix: Simulations per Configuration"),
          p(style = "color: #6c757d;",
            "Each cell shows the number of simulations for that (n, phi) combination. ",
            "Faceted by innovation distribution. Empty/missing cells indicate gaps in coverage."
          ),
          plotOutput("coverage_matrix", height = "450px")
        ),

        # =====================================================================
        # Tab 2: Rejection Rates Table
        # =====================================================================
        tabPanel(
          "Rejection Rates",
          br(),
          fluidRow(
            column(9, h4("Rejection Rate Summary")),
            column(3, downloadButton("export_csv", "Export CSV",
                                     class = "btn-sm btn-outline-secondary pull-right"))
          ),
          DTOutput("rejection_table")
        ),

        # =====================================================================
        # Tab 3: Plots
        # =====================================================================
        tabPanel(
          "Plots",
          br(),
          div(class = "plot-controls",
            fluidRow(
              column(3,
                selectInput("plot_type", "Plot Type",
                            choices = c("Power Curve" = "power_curve",
                                        "Heatmap" = "heatmap",
                                        "Deviation from Nominal" = "deviation",
                                        "Rate vs Sample Size" = "rate_vs_n",
                                        "Method Forest Plot" = "forest"))
              ),
              column(3,
                conditionalPanel(
                  condition = "input.plot_type == 'power_curve' || input.plot_type == 'heatmap' || input.plot_type == 'rate_vs_n'",
                  selectInput("rate_type", "Rate Type",
                              choices = c("Bootstrap" = "reject_05",
                                          "Asymptotic" = "reject_asymp_05",
                                          "COBA" = "reject_adj_05"))
                )
              ),
              column(3,
                conditionalPanel(
                  condition = "input.plot_type == 'power_curve'",
                  selectInput("facet_by", "Facet By",
                              choices = c("Sample Size (n)" = "n",
                                          "Innovation" = "innov_dist",
                                          "None" = "none"))
                ),
                conditionalPanel(
                  condition = "input.plot_type == 'rate_vs_n'",
                  selectInput("rate_vs_n_facet", "Facet By",
                              choices = c("Innovation" = "innov_dist",
                                          "Phi" = "phi",
                                          "None" = "none"))
                )
              )
            )
          ),
          plotOutput("main_plot", height = "500px"),
          br(),
          downloadButton("export_plot", "Download Plot (PNG)",
                         class = "btn-sm btn-outline-secondary")
        ),

        # =====================================================================
        # Tab 4: P-value Diagnostics
        # =====================================================================
        tabPanel(
          "P-value Diagnostics",
          br(),
          h4("P-value Calibration Diagnostics"),
          p(style = "color: #6c757d;",
            "Under H0 (phi = 0), p-values should be uniformly distributed on [0, 1]. ",
            "Non-uniformity indicates size distortion. Filter to phi = 0 configs to check calibration."
          ),

          div(class = "plot-controls",
            fluidRow(
              column(3,
                selectizeInput("pval_n_filter", "Sample Size (n)",
                               choices = NULL, multiple = TRUE,
                               options = list(placeholder = "All..."))
              ),
              column(3,
                selectizeInput("pval_phi_filter", "Phi Value",
                               choices = NULL, multiple = TRUE,
                               options = list(placeholder = "All..."))
              ),
              column(3,
                selectizeInput("pval_innov_filter", "Innovation Distribution",
                               choices = NULL, multiple = TRUE,
                               options = list(placeholder = "All..."))
              ),
              column(3,
                div(style = "margin-top: 25px;",
                  textOutput("pval_status")
                )
              )
            )
          ),

          # P-value histogram
          wellPanel(
            h5("P-value Histograms"),
            p(style = "color: #6c757d; margin-bottom: 10px;",
              "Histograms of p-values for each method. ",
              "Under H0, these should be approximately uniform (flat)."),
            plotOutput("pval_hist", height = "300px")
          ),

          # QQ plot
          wellPanel(
            h5("P-value QQ Plot"),
            p(style = "color: #6c757d; margin-bottom: 10px;",
              "Quantile-quantile plot against Uniform(0, 1). ",
              "Points on the diagonal indicate well-calibrated p-values."),
            plotOutput("pval_qq", height = "350px")
          ),

          # Method comparison scatter
          wellPanel(
            h5("P-value Method Comparison"),
            p(style = "color: #6c757d; margin-bottom: 10px;",
              "Scatter plots comparing p-values between methods. ",
              "Points on the diagonal indicate agreement."),
            plotOutput("pval_scatter", height = "350px")
          )
        ),

        # =====================================================================
        # Tab 5: Bootstrap Distribution
        # =====================================================================
        tabPanel(
          "Bootstrap Distribution",
          br(),
          fluidRow(
            column(4,
              numericInput("sim_id_input", "Simulation ID (sim_id):",
                           value = 1, min = 1, step = 1)
            ),
            column(4,
              actionButton("load_dist", "Load Distribution",
                           class = "btn-primary")
            )
          ),
          br(),
          plotOutput("dist_plot", height = "400px"),
          verbatimTextOutput("dist_summary")
        ),

        # =====================================================================
        # Tab 6: Analysis Grids
        # =====================================================================
        tabPanel(
          "Analysis Grids",
          br(),
          h4("Rejection Rate Comparison Grids"),
          p(style = "color: #6c757d;",
            "Each grid fixes two dimensions and varies the third. ",
            "Rates are color-coded: ",
            span(style = "background-color:#fff3cd; padding:2px 6px; border-radius:3px;", "< 0.03"),
            " ",
            span(style = "background-color:#d4edda; padding:2px 6px; border-radius:3px;", "0.03 \u2013 0.07"),
            " ",
            span(style = "background-color:#f8d7da; padding:2px 6px; border-radius:3px;", "> 0.07")
          ),

          # Grid 1: By Sample Size
          wellPanel(
            h5("By Sample Size"),
            p(style = "color: #6c757d; margin-bottom: 10px;",
              "Fix innovation distribution and phi; rows vary by n"),
            fluidRow(
              column(3,
                selectInput("grid1_innov", "Innovation Distribution:",
                            choices = NULL)
              ),
              column(3,
                selectInput("grid1_phi", "Phi:",
                            choices = NULL)
              ),
              column(3,
                selectInput("grid1_batch", "Batch:",
                            choices = c("All (Pooled)" = "all"))
              )
            ),
            DTOutput("grid1_table")
          ),

          # Grid 2: By Phi
          wellPanel(
            h5("By Phi"),
            p(style = "color: #6c757d; margin-bottom: 10px;",
              "Fix innovation distribution and n; rows vary by phi"),
            fluidRow(
              column(3,
                selectInput("grid2_innov", "Innovation Distribution:",
                            choices = NULL)
              ),
              column(3,
                selectInput("grid2_n", "Sample Size (n):",
                            choices = NULL)
              ),
              column(3,
                selectInput("grid2_batch", "Batch:",
                            choices = c("All (Pooled)" = "all"))
              )
            ),
            DTOutput("grid2_table")
          ),

          # Grid 3: By Distribution
          wellPanel(
            h5("By Innovation Distribution"),
            p(style = "color: #6c757d; margin-bottom: 10px;",
              "Fix phi and n; rows vary by innovation distribution"),
            fluidRow(
              column(3,
                selectInput("grid3_phi", "Phi:",
                            choices = NULL)
              ),
              column(3,
                selectInput("grid3_n", "Sample Size (n):",
                            choices = NULL)
              ),
              column(3,
                selectInput("grid3_batch", "Batch:",
                            choices = c("All (Pooled)" = "all"))
              )
            ),
            DTOutput("grid3_table")
          )
        ),

        # =====================================================================
        # Tab 7: Parallel Coordinates
        # =====================================================================
        tabPanel(
          "Parallel Coordinates",
          br(),
          h4("Parallel Coordinates Plots"),
          p(style = "color: #6c757d;",
            "Each line represents one (n, phi, innov_dist) configuration. ",
            "Axes are scaled to [0, 1] with original values shown. ",
            "Use the filters below to reduce clutter."
          ),

          # Shared filters for this tab
          div(class = "plot-controls",
            fluidRow(
              column(3,
                selectizeInput("pc_n_filter", "Sample Size (n)",
                               choices = NULL, multiple = TRUE,
                               options = list(placeholder = "All..."))
              ),
              column(3,
                selectizeInput("pc_phi_filter", "Phi Value",
                               choices = NULL, multiple = TRUE,
                               options = list(placeholder = "All..."))
              ),
              column(3,
                selectizeInput("pc_innov_filter", "Innovation Distribution",
                               choices = NULL, multiple = TRUE,
                               options = list(placeholder = "All..."))
              )
            )
          ),

          # Plot 1: Full Profile
          wellPanel(
            h5("Full Profile"),
            p(style = "color: #6c757d; margin-bottom: 10px;",
              "Axes: n \u2192 phi \u2192 CO \u2192 COB \u2192 COBA. ",
              "Color by innovation distribution."),
            plotOutput("pc_full_plot", height = "350px")
          ),

          # Plot 2: Method Comparison
          wellPanel(
            h5("Method Comparison"),
            p(style = "color: #6c757d; margin-bottom: 10px;",
              "Axes: CO \u2192 COB \u2192 COBA. ",
              "Focused view comparing how the three methods relate."),
            fluidRow(
              column(3,
                selectInput("pc_color_by", "Color by:",
                            choices = c("Innovation" = "innov_dist",
                                        "Sample Size (n)" = "n",
                                        "Phi" = "phi"))
              )
            ),
            plotOutput("pc_method_plot", height = "350px")
          ),

          # Plot 3: Single-Method Sensitivity
          wellPanel(
            h5("Single-Method Sensitivity"),
            p(style = "color: #6c757d; margin-bottom: 10px;",
              "Axes: n \u2192 phi \u2192 Rate. ",
              "How one method responds to DGP parameters."),
            fluidRow(
              column(3,
                selectInput("pc_method", "Method:",
                            choices = c("Bootstrap (COB)" = "reject_05",
                                        "Asymptotic (CO)" = "reject_asymp_05",
                                        "COBA" = "reject_adj_05"))
              )
            ),
            plotOutput("pc_sensitivity_plot", height = "350px")
          )
        ),

        # =====================================================================
        # Tab 8: Diagnostics
        # =====================================================================
        tabPanel(
          "Diagnostics",
          br(),
          h4("Simulation Diagnostics"),
          p(style = "color: #6c757d;",
            "Validates the Monte Carlo machinery: null model fits, ",
            "convergence of rejection rates, batch-to-batch consistency, ",
            "and test statistic distributions."
          ),

          div(class = "plot-controls",
            fluidRow(
              column(3,
                selectizeInput("diag_n_filter", "Sample Size (n)",
                               choices = NULL, multiple = TRUE,
                               options = list(placeholder = "All..."))
              ),
              column(3,
                selectizeInput("diag_phi_filter", "Phi Value",
                               choices = NULL, multiple = TRUE,
                               options = list(placeholder = "All..."))
              ),
              column(3,
                selectizeInput("diag_innov_filter", "Innovation Distribution",
                               choices = NULL, multiple = TRUE,
                               options = list(placeholder = "All..."))
              ),
              column(3,
                div(style = "margin-top: 25px;",
                  textOutput("diag_status")
                )
              )
            )
          ),

          # AR order distribution
          wellPanel(
            h5("Fitted AR Order Distribution"),
            p(style = "color: #6c757d; margin-bottom: 10px;",
              "How often each AR order is selected under the null. ",
              "For AR(1) DGPs, order 1 should dominate."),
            plotOutput("diag_ar_order_plot", height = "300px")
          ),

          # Estimated phi distribution
          wellPanel(
            h5("Estimated Null AR Coefficient Distribution"),
            p(style = "color: #6c757d; margin-bottom: 10px;",
              "Distribution of the first fitted AR coefficient across simulations. ",
              "Should cluster near the true phi value."),
            plotOutput("diag_phi_plot", height = "300px")
          ),

          # MC Convergence
          wellPanel(
            h5("Monte Carlo Convergence"),
            p(style = "color: #6c757d; margin-bottom: 10px;",
              "Running rejection rate as simulations accumulate. ",
              "The rate should stabilize, confirming sufficient MC sample size."),
            fluidRow(
              column(4,
                selectInput("diag_conv_method", "Method:",
                            choices = c("Bootstrap (COB)" = "pvalue",
                                        "Asymptotic (CO)" = "pvalue_asymp",
                                        "COBA" = "pvalue_adj"))
              )
            ),
            plotOutput("diag_convergence_plot", height = "350px")
          ),

          # Batch Consistency
          wellPanel(
            h5("Batch-to-Batch Consistency"),
            p(style = "color: #6c757d; margin-bottom: 10px;",
              "Per-batch rejection rates for each configuration. ",
              "Low spread indicates stable Monte Carlo estimates."),
            fluidRow(
              column(4,
                selectInput("diag_batch_method", "Method:",
                            choices = c("Bootstrap (COB)" = "reject_05",
                                        "Asymptotic (CO)" = "reject_asymp_05",
                                        "COBA" = "reject_adj_05"))
              )
            ),
            plotOutput("diag_batch_plot", height = "350px")
          ),

          # Test statistic distribution
          wellPanel(
            h5("Test Statistic Distribution"),
            p(style = "color: #6c757d; margin-bottom: 10px;",
              "Distribution of the observed test statistic across simulations. ",
              "Helps identify outliers or unexpected distributional shapes."),
            plotOutput("diag_tstat_plot", height = "300px")
          )
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
  # Database connection
  # ---------------------------------------------------------------------------
  db_path <- getOption("tstse.viewer_db")
  con <- dbConnect(duckdb(), dbdir = db_path, read_only = TRUE)
  onStop(function() dbDisconnect(con, shutdown = TRUE))

  # ---------------------------------------------------------------------------
  # Populate initial filter choices
  # ---------------------------------------------------------------------------
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

  updateSelectizeInput(session, "n_filter", choices = init_choices$n, server = FALSE)
  updateSelectizeInput(session, "phi_filter", choices = init_choices$phi, server = FALSE)
  updateSelectizeInput(session, "innov_filter", choices = init_choices$innov, server = FALSE)

  # Populate Analysis Grid dropdowns (single-select, default to first value)
  updateSelectInput(session, "grid1_innov", choices = init_choices$innov,
                    selected = init_choices$innov[1])
  updateSelectInput(session, "grid1_phi", choices = init_choices$phi,
                    selected = init_choices$phi[1])
  updateSelectInput(session, "grid2_innov", choices = init_choices$innov,
                    selected = init_choices$innov[1])
  updateSelectInput(session, "grid2_n", choices = init_choices$n,
                    selected = init_choices$n[1])
  updateSelectInput(session, "grid3_phi", choices = init_choices$phi,
                    selected = init_choices$phi[1])
  updateSelectInput(session, "grid3_n", choices = init_choices$n,
                    selected = init_choices$n[1])

  # Populate Parallel Coordinates filters
  updateSelectizeInput(session, "pc_n_filter", choices = init_choices$n, server = FALSE)
  updateSelectizeInput(session, "pc_phi_filter", choices = init_choices$phi, server = FALSE)
  updateSelectizeInput(session, "pc_innov_filter", choices = init_choices$innov, server = FALSE)

  # Populate P-value Diagnostics filters
  updateSelectizeInput(session, "pval_n_filter", choices = init_choices$n, server = FALSE)
  updateSelectizeInput(session, "pval_phi_filter", choices = init_choices$phi, server = FALSE)
  updateSelectizeInput(session, "pval_innov_filter", choices = init_choices$innov, server = FALSE)

  # Populate Diagnostics filters
  updateSelectizeInput(session, "diag_n_filter", choices = init_choices$n, server = FALSE)
  updateSelectizeInput(session, "diag_phi_filter", choices = init_choices$phi, server = FALSE)
  updateSelectizeInput(session, "diag_innov_filter", choices = init_choices$innov, server = FALSE)

  # ---------------------------------------------------------------------------
  # Build query with filters
  # ---------------------------------------------------------------------------
  build_query <- function(view_name) {
    sql <- paste0("SELECT * FROM ", view_name, " WHERE 1=1")
    params <- list()

    if (length(input$n_filter) > 0) {
      placeholders <- paste(rep("?", length(input$n_filter)), collapse = ", ")
      sql <- paste0(sql, " AND n IN (", placeholders, ")")
      params <- c(params, as.list(as.integer(input$n_filter)))
    }
    if (length(input$phi_filter) > 0) {
      placeholders <- paste(rep("?", length(input$phi_filter)), collapse = ", ")
      sql <- paste0(sql, " AND phi IN (", placeholders, ")")
      params <- c(params, as.list(as.numeric(input$phi_filter)))
    }
    if (length(input$innov_filter) > 0) {
      placeholders <- paste(rep("?", length(input$innov_filter)), collapse = ", ")
      sql <- paste0(sql, " AND innov_dist IN (", placeholders, ")")
      params <- c(params, as.list(input$innov_filter))
    }

    sql <- paste0(sql, " ORDER BY innov_dist, n, phi")
    list(sql = sql, params = params)
  }

  # ---------------------------------------------------------------------------
  # Grid query helper (independent of sidebar filters)
  # ---------------------------------------------------------------------------
  grid_query <- function(innov_dist = NULL, phi = NULL, n = NULL,
                         order_col = "n", batch_id = NULL) {
    if (!is.null(batch_id) && batch_id != "all") {
      sql <- "SELECT * FROM v_rejection_rates_by_batch WHERE 1=1"
    } else {
      sql <- "SELECT * FROM v_rejection_rates WHERE 1=1"
    }
    params <- list()

    if (!is.null(batch_id) && batch_id != "all") {
      sql <- paste0(sql, " AND batch_id = ?")
      params <- c(params, list(as.integer(batch_id)))
    }
    if (!is.null(innov_dist)) {
      sql <- paste0(sql, " AND innov_dist = ?")
      params <- c(params, list(innov_dist))
    }
    if (!is.null(phi)) {
      sql <- paste0(sql, " AND phi = ?")
      params <- c(params, list(as.numeric(phi)))
    }
    if (!is.null(n)) {
      sql <- paste0(sql, " AND n = ?")
      params <- c(params, list(as.integer(n)))
    }

    sql <- paste0(sql, " ORDER BY ", order_col)
    tryCatch(
      dbGetQuery(con, sql, params = params),
      error = function(e) data.frame()
    )
  }

  # ---------------------------------------------------------------------------
  # Grid batch choices helper: returns batches matching a config filter
  # ---------------------------------------------------------------------------
  grid_batch_choices <- function(innov_dist = NULL, phi = NULL, n = NULL) {
    sql <- paste0(
      "SELECT DISTINCT b.batch_id, b.label ",
      "FROM simulations s JOIN batches b ON s.batch_id = b.batch_id WHERE 1=1"
    )
    params <- list()
    if (!is.null(innov_dist)) {
      sql <- paste0(sql, " AND s.innov_dist = ?")
      params <- c(params, list(innov_dist))
    }
    if (!is.null(phi)) {
      sql <- paste0(sql, " AND s.phi = ?")
      params <- c(params, list(as.numeric(phi)))
    }
    if (!is.null(n)) {
      sql <- paste0(sql, " AND s.n = ?")
      params <- c(params, list(as.integer(n)))
    }
    sql <- paste0(sql, " ORDER BY b.batch_id")
    batches <- tryCatch(dbGetQuery(con, sql, params = params),
                        error = function(e) data.frame())
    if (nrow(batches) == 0) return(c("All (Pooled)" = "all"))
    choices <- setNames(
      as.character(batches$batch_id),
      ifelse(is.na(batches$label) | batches$label == "",
             paste("Batch", batches$batch_id),
             paste0("Batch ", batches$batch_id, ": ", batches$label))
    )
    c("All (Pooled)" = "all", choices)
  }

  # Update per-grid batch dropdowns when their config filters change
  observe({
    req(input$grid1_innov, input$grid1_phi)
    choices <- grid_batch_choices(innov_dist = input$grid1_innov, phi = input$grid1_phi)
    sel <- if (input$grid1_batch %in% choices) input$grid1_batch else "all"
    updateSelectInput(session, "grid1_batch", choices = choices, selected = sel)
  })

  observe({
    req(input$grid2_innov, input$grid2_n)
    choices <- grid_batch_choices(innov_dist = input$grid2_innov, n = input$grid2_n)
    sel <- if (input$grid2_batch %in% choices) input$grid2_batch else "all"
    updateSelectInput(session, "grid2_batch", choices = choices, selected = sel)
  })

  observe({
    req(input$grid3_phi, input$grid3_n)
    choices <- grid_batch_choices(phi = input$grid3_phi, n = input$grid3_n)
    sel <- if (input$grid3_batch %in% choices) input$grid3_batch else "all"
    updateSelectInput(session, "grid3_batch", choices = choices, selected = sel)
  })

  # ---------------------------------------------------------------------------
  # Grid DT formatting helper
  # ---------------------------------------------------------------------------
  format_grid_dt <- function(df, row_label_col, row_label_name) {
    if (nrow(df) == 0) {
      return(datatable(data.frame(Message = "No data for this combination"),
                       rownames = FALSE, options = list(dom = "t")))
    }

    display_df <- data.frame(
      Label    = df[[row_label_col]],
      n_sims   = df$n_sims,
      CO       = df$reject_asymp_05,
      CO_SE    = df$reject_asymp_05_se,
      COB      = df$reject_05,
      COB_SE   = df$reject_05_se,
      COBA     = df$reject_adj_05,
      COBA_SE  = df$reject_adj_05_se
    )

    col_names <- c(row_label_name, "Sims", "CO Rate", "CO SE",
                   "COB Rate", "COB SE", "COBA Rate", "COBA SE")

    dt <- datatable(
      display_df,
      colnames = col_names,
      rownames = FALSE,
      options = list(
        dom = "t",
        paging = FALSE,
        searching = FALSE,
        ordering = FALSE,
        columnDefs = list(
          list(targets = c(3, 5, 7), className = "dt-body-right")
        )
      )
    )

    rate_cols <- c("CO", "COB", "COBA")
    se_cols <- c("CO_SE", "COB_SE", "COBA_SE")

    dt <- formatRound(dt, columns = rate_cols, digits = 4)
    dt <- formatRound(dt, columns = se_cols, digits = 4)

    for (col in rate_cols) {
      dt <- formatStyle(dt, col,
        backgroundColor = styleInterval(
          c(0.03, 0.07),
          c("#fff3cd", "#d4edda", "#f8d7da")
        )
      )
    }

    dt
  }

  # ---------------------------------------------------------------------------
  # Curated data reactive
  # ---------------------------------------------------------------------------
  curated_data <- reactive({
    view_name <- if (input$compare_batches) {
      "v_rejection_rates_by_batch"
    } else {
      "v_rejection_rates"
    }
    q <- build_query(view_name)
    tryCatch(
      dbGetQuery(con, q$sql, params = q$params),
      error = function(e) data.frame()
    )
  })

  # ---------------------------------------------------------------------------
  # Cascading filter updates
  # ---------------------------------------------------------------------------
  observe({
    df <- curated_data()
    if (nrow(df) == 0) return()

    n_choices <- sort(unique(df$n))
    phi_choices <- sort(unique(df$phi))
    innov_choices <- sort(unique(df$innov_dist))

    updateSelectizeInput(session, "n_filter", choices = as.character(n_choices),
                         selected = intersect(input$n_filter, as.character(n_choices)),
                         server = FALSE)
    updateSelectizeInput(session, "phi_filter", choices = as.character(phi_choices),
                         selected = intersect(input$phi_filter, as.character(phi_choices)),
                         server = FALSE)
    updateSelectizeInput(session, "innov_filter", choices = innov_choices,
                         selected = intersect(input$innov_filter, innov_choices),
                         server = FALSE)
  })

  # Clear all filters
  observeEvent(input$clear_filters, {
    updateSelectizeInput(session, "n_filter", selected = character(0))
    updateSelectizeInput(session, "phi_filter", selected = character(0))
    updateSelectizeInput(session, "innov_filter", selected = character(0))
    updateCheckboxInput(session, "compare_batches", value = FALSE)
  })

  # Status text
  output$status_text <- renderText({
    df <- curated_data()
    if (nrow(df) == 0) return("No data matching filters")
    total_sims <- sum(df$n_sims, na.rm = TRUE)
    paste0(nrow(df), " scenarios, ", format(total_sims, big.mark = ","), " total sims")
  })

  # ---------------------------------------------------------------------------
  # Tab 1: Study Overview
  # ---------------------------------------------------------------------------

  overview_stats <- reactive({
    tryCatch({
      dbGetQuery(con, "
        SELECT COUNT(*) as total_sims,
               COUNT(DISTINCT batch_id) as n_batches,
               MIN(created_at) as first_sim,
               MAX(created_at) as last_sim
        FROM simulations
      ")
    }, error = function(e) data.frame(total_sims = 0, n_batches = 0,
                                       first_sim = NA, last_sim = NA))
  })

  overview_coverage <- reactive({
    tryCatch({
      dbGetQuery(con, "
        SELECT n, phi, innov_dist, COUNT(*) as n_sims
        FROM simulations
        GROUP BY n, phi, innov_dist
        ORDER BY innov_dist, n, phi
      ")
    }, error = function(e) data.frame())
  })

  output$overview_total_sims <- renderText({
    stats <- overview_stats()
    format(stats$total_sims[1], big.mark = ",")
  })

  output$overview_n_configs <- renderText({
    df <- overview_coverage()
    as.character(nrow(df))
  })

  output$overview_n_batches <- renderText({
    stats <- overview_stats()
    as.character(stats$n_batches[1])
  })

  output$overview_date_range <- renderText({
    stats <- overview_stats()
    if (is.na(stats$first_sim[1])) return("No data")
    paste0(format(stats$first_sim[1], "%Y-%m-%d"), " to ",
           format(stats$last_sim[1], "%Y-%m-%d"))
  })

  output$coverage_matrix <- renderPlot({
    df <- overview_coverage()
    if (nrow(df) == 0) {
      plot.new()
      text(0.5, 0.5, "No data", cex = 1.2, col = "grey50")
      return()
    }

    df$n <- factor(df$n)
    df$phi <- factor(df$phi)

    p <- ggplot(df, aes(x = phi, y = n, fill = n_sims)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = format(n_sims, big.mark = ",")),
                size = 3.5, color = "black") +
      scale_fill_gradient(low = "#deebf7", high = "#2171b5",
                          name = "Simulations") +
      labs(x = expression(phi), y = "Sample Size (n)") +
      theme_minimal(base_size = 14) +
      theme(
        panel.grid = element_blank(),
        legend.position = "right",
        strip.text = element_text(face = "bold")
      ) +
      facet_wrap(~ innov_dist)

    print(p)
  })

  # ---------------------------------------------------------------------------
  # Tab 2: Rejection Rates Table
  # ---------------------------------------------------------------------------

  output$rejection_table <- renderDT({
    df <- curated_data()
    if (nrow(df) == 0) return(datatable(data.frame(Message = "No data")))

    # Select display columns based on mode
    if (input$compare_batches) {
      display_cols <- c("batch_id", "batch_label", "innov_dist", "n", "phi",
                        "n_sims", "reject_05", "reject_05_se",
                        "reject_asymp_05", "reject_asymp_05_se",
                        "reject_adj_05", "reject_adj_05_se")
      pretty_names <- c("Batch", "Label", "Innovation", "n", "Phi",
                        "Sims", "Boot Rate", "Boot SE",
                        "Asymp Rate", "Asymp SE",
                        "COBA Rate", "COBA SE")
    } else {
      display_cols <- c("innov_dist", "n", "phi", "n_sims", "n_batches",
                        "reject_05", "reject_05_se",
                        "reject_asymp_05", "reject_asymp_05_se",
                        "reject_adj_05", "reject_adj_05_se")
      pretty_names <- c("Innovation", "n", "Phi", "Sims", "Batches",
                        "Boot Rate", "Boot SE",
                        "Asymp Rate", "Asymp SE",
                        "COBA Rate", "COBA SE")
    }

    available <- intersect(display_cols, names(df))
    df_show <- df[, available, drop = FALSE]
    col_names <- pretty_names[display_cols %in% available]

    dt <- datatable(
      df_show,
      colnames = col_names,
      rownames = FALSE,
      filter = "top",
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        order = list()
      )
    )

    # Round rate and SE columns
    rate_cols <- intersect(c("reject_05", "reject_asymp_05", "reject_adj_05"), available)
    se_cols <- intersect(c("reject_05_se", "reject_asymp_05_se", "reject_adj_05_se"), available)

    if (length(rate_cols) > 0) {
      dt <- formatRound(dt, columns = rate_cols, digits = 4)
    }
    if (length(se_cols) > 0) {
      dt <- formatRound(dt, columns = se_cols, digits = 4)
    }

    # Color-code rejection rate columns
    for (col in rate_cols) {
      dt <- formatStyle(
        dt, col,
        backgroundColor = styleInterval(
          c(0.03, 0.07),
          c("#fff3cd", "#d4edda", "#f8d7da")  # yellow / green / red
        )
      )
    }

    dt
  })

  # CSV export
  output$export_csv <- downloadHandler(
    filename = function() {
      paste0("rejection_rates_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      write.csv(curated_data(), file, row.names = FALSE)
    }
  )

  # ---------------------------------------------------------------------------
  # Tab 3: Plots
  # ---------------------------------------------------------------------------

  current_plot <- reactive({
    df <- curated_data()
    if (nrow(df) == 0) return(NULL)

    rate_col <- input$rate_type
    if (!rate_col %in% names(df)) return(NULL)

    rate_labels <- c(
      reject_05 = "Bootstrap Rejection Rate",
      reject_asymp_05 = "Asymptotic Rejection Rate",
      reject_adj_05 = "COBA Rejection Rate"
    )

    if (input$plot_type == "power_curve") {
      # Power curve: rejection rate vs phi
      df$phi <- as.numeric(df$phi)

      p <- ggplot(df, aes(x = phi, y = .data[[rate_col]],
                          color = factor(innov_dist),
                          group = factor(innov_dist))) +
        geom_line(linewidth = 0.8) +
        geom_point(size = 2) +
        geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey40") +
        labs(
          x = expression(phi),
          y = "Rejection Rate",
          color = "Innovation",
          title = rate_labels[rate_col]
        ) +
        theme_minimal(base_size = 14) +
        theme(legend.position = "bottom")

      # SE error bars
      se_col <- paste0(rate_col, "_se")
      if (se_col %in% names(df)) {
        p <- p + geom_errorbar(
          aes(ymin = .data[[rate_col]] - .data[[se_col]],
              ymax = .data[[rate_col]] + .data[[se_col]]),
          width = 0.02, alpha = 0.5
        )
      }

      # Faceting
      if (input$facet_by == "n") {
        p <- p + facet_wrap(~ n, scales = "free_x",
                            labeller = labeller(n = function(x) paste0("n = ", x)))
      } else if (input$facet_by == "innov_dist") {
        p <- p + facet_wrap(~ innov_dist)
      }

      p

    } else if (input$plot_type == "heatmap") {
      # Heatmap: n vs phi, fill = rejection rate
      df$n <- factor(df$n)
      df$phi <- factor(df$phi)

      p <- ggplot(df, aes(x = phi, y = n, fill = .data[[rate_col]])) +
        geom_tile(color = "white", linewidth = 0.5) +
        scale_fill_gradient2(
          low = "#2166ac", mid = "#f7f7f7", high = "#b2182b",
          midpoint = 0.05, name = "Rejection\nRate"
        ) +
        geom_text(aes(label = sprintf("%.3f", .data[[rate_col]])),
                  size = 3, color = "black") +
        labs(
          x = expression(phi),
          y = "Sample Size (n)",
          title = rate_labels[rate_col]
        ) +
        theme_minimal(base_size = 14) +
        theme(
          panel.grid = element_blank(),
          legend.position = "right"
        ) +
        facet_wrap(~ innov_dist)

      p

    } else if (input$plot_type == "deviation") {
      # Deviation from nominal: dot plot of (rate - 0.05) for all 3 methods
      # Build long-format data with one row per (config, method)
      dev_rows <- list()
      for (i in seq_len(nrow(df))) {
        config_label <- sprintf("%s / n=%s / phi=%s", df$innov_dist[i], df$n[i], df$phi[i])
        if ("reject_asymp_05" %in% names(df)) {
          dev_rows <- c(dev_rows, list(data.frame(
            config = config_label, method = "CO",
            deviation = df$reject_asymp_05[i] - 0.05,
            se = if ("reject_asymp_05_se" %in% names(df)) df$reject_asymp_05_se[i] else NA_real_,
            stringsAsFactors = FALSE
          )))
        }
        if ("reject_05" %in% names(df)) {
          dev_rows <- c(dev_rows, list(data.frame(
            config = config_label, method = "COB",
            deviation = df$reject_05[i] - 0.05,
            se = if ("reject_05_se" %in% names(df)) df$reject_05_se[i] else NA_real_,
            stringsAsFactors = FALSE
          )))
        }
        if ("reject_adj_05" %in% names(df)) {
          dev_rows <- c(dev_rows, list(data.frame(
            config = config_label, method = "COBA",
            deviation = df$reject_adj_05[i] - 0.05,
            se = if ("reject_adj_05_se" %in% names(df)) df$reject_adj_05_se[i] else NA_real_,
            stringsAsFactors = FALSE
          )))
        }
      }
      dev_df <- do.call(rbind, dev_rows)
      dev_df$method <- factor(dev_df$method, levels = c("CO", "COB", "COBA"))

      p <- ggplot(dev_df, aes(x = deviation, y = config, color = method)) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
        geom_vline(xintercept = c(-0.02, 0.02), linetype = "dotted", color = "grey70") +
        geom_point(size = 2.5, position = position_dodge(width = 0.6)) +
        {
          if (!all(is.na(dev_df$se))) {
            geom_errorbarh(
              aes(xmin = deviation - 1.96 * se, xmax = deviation + 1.96 * se),
              height = 0.2, position = position_dodge(width = 0.6), alpha = 0.5
            )
          }
        } +
        scale_color_manual(values = c(CO = "#e41a1c", COB = "#377eb8", COBA = "#4daf4a")) +
        labs(
          x = "Deviation from Nominal (Rate - 0.05)",
          y = NULL,
          color = "Method",
          title = "Deviation from Nominal 5% Rate"
        ) +
        theme_minimal(base_size = 13) +
        theme(legend.position = "bottom")

      p

    } else if (input$plot_type == "rate_vs_n") {
      # Rate vs sample size: line plot of rate vs n
      df$n <- as.numeric(df$n)

      p <- ggplot(df, aes(x = n, y = .data[[rate_col]],
                          color = factor(innov_dist),
                          group = factor(innov_dist))) +
        geom_line(linewidth = 0.8) +
        geom_point(size = 2) +
        geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey40") +
        labs(
          x = "Sample Size (n)",
          y = "Rejection Rate",
          color = "Innovation",
          title = paste(rate_labels[rate_col], "vs Sample Size")
        ) +
        theme_minimal(base_size = 14) +
        theme(legend.position = "bottom")

      # SE error bars
      se_col <- paste0(rate_col, "_se")
      if (se_col %in% names(df)) {
        p <- p + geom_errorbar(
          aes(ymin = .data[[rate_col]] - .data[[se_col]],
              ymax = .data[[rate_col]] + .data[[se_col]]),
          width = 10, alpha = 0.5
        )
      }

      # Faceting
      facet_choice <- input$rate_vs_n_facet
      if (!is.null(facet_choice) && facet_choice == "innov_dist") {
        p <- p + facet_wrap(~ innov_dist)
      } else if (!is.null(facet_choice) && facet_choice == "phi") {
        p <- p + facet_wrap(~ phi,
                            labeller = labeller(phi = function(x) paste0("phi = ", x)))
      }

      p

    } else if (input$plot_type == "forest") {
      # Method comparison forest plot: 3 methods per config, horizontal
      forest_rows <- list()
      for (i in seq_len(nrow(df))) {
        config_label <- sprintf("%s / n=%s / phi=%s", df$innov_dist[i], df$n[i], df$phi[i])
        if ("reject_asymp_05" %in% names(df)) {
          forest_rows <- c(forest_rows, list(data.frame(
            config = config_label, method = "CO",
            rate = df$reject_asymp_05[i],
            se = if ("reject_asymp_05_se" %in% names(df)) df$reject_asymp_05_se[i] else NA_real_,
            stringsAsFactors = FALSE
          )))
        }
        if ("reject_05" %in% names(df)) {
          forest_rows <- c(forest_rows, list(data.frame(
            config = config_label, method = "COB",
            rate = df$reject_05[i],
            se = if ("reject_05_se" %in% names(df)) df$reject_05_se[i] else NA_real_,
            stringsAsFactors = FALSE
          )))
        }
        if ("reject_adj_05" %in% names(df)) {
          forest_rows <- c(forest_rows, list(data.frame(
            config = config_label, method = "COBA",
            rate = df$reject_adj_05[i],
            se = if ("reject_adj_05_se" %in% names(df)) df$reject_adj_05_se[i] else NA_real_,
            stringsAsFactors = FALSE
          )))
        }
      }
      forest_df <- do.call(rbind, forest_rows)
      forest_df$method <- factor(forest_df$method, levels = c("CO", "COB", "COBA"))

      p <- ggplot(forest_df, aes(x = rate, y = config, color = method)) +
        geom_vline(xintercept = 0.05, linetype = "dashed", color = "grey40") +
        geom_point(size = 2.5, position = position_dodge(width = 0.6)) +
        {
          if (!all(is.na(forest_df$se))) {
            geom_errorbarh(
              aes(xmin = rate - 1.96 * se, xmax = rate + 1.96 * se),
              height = 0.2, position = position_dodge(width = 0.6), alpha = 0.5
            )
          }
        } +
        scale_color_manual(values = c(CO = "#e41a1c", COB = "#377eb8", COBA = "#4daf4a")) +
        labs(
          x = "Rejection Rate",
          y = NULL,
          color = "Method",
          title = "Method Comparison: Rejection Rate with 95% CI"
        ) +
        theme_minimal(base_size = 13) +
        theme(legend.position = "bottom")

      p
    }
  })

  output$main_plot <- renderPlot({
    p <- current_plot()
    if (is.null(p)) {
      plot.new()
      text(0.5, 0.5, "No data matching filters", cex = 1.2, col = "grey50")
      return()
    }
    print(p)
  })

  # PNG download
  output$export_plot <- downloadHandler(
    filename = function() {
      paste0("mc_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      p <- current_plot()
      if (!is.null(p)) {
        ggsave(file, plot = p, width = 10, height = 6, dpi = 300)
      }
    }
  )

  # ---------------------------------------------------------------------------
  # Tab 5: Bootstrap Distribution
  # ---------------------------------------------------------------------------

  dist_data <- eventReactive(input$load_dist, {
    sim_id <- input$sim_id_input
    tryCatch({
      result <- dbGetQuery(con,
        "SELECT sim_id, obs_stat, boot_dist, pvalue, n, phi, innov_dist
         FROM simulations WHERE sim_id = ?",
        params = list(as.integer(sim_id)))
      if (nrow(result) == 0) return(NULL)
      result
    }, error = function(e) NULL)
  })

  output$dist_plot <- renderPlot({
    d <- dist_data()
    if (is.null(d)) {
      plot.new()
      text(0.5, 0.5, "Enter a sim_id and click Load", cex = 1.2, col = "grey50")
      return()
    }

    boot_dist <- d$boot_dist[[1]]
    obs_stat <- d$obs_stat[1]
    pval <- d$pvalue[1]

    hist_data <- data.frame(t_stat = boot_dist)

    p <- ggplot(hist_data, aes(x = t_stat)) +
      geom_histogram(bins = 40, fill = "#4292c6", color = "white", alpha = 0.8) +
      geom_vline(xintercept = obs_stat, color = "red", linewidth = 1.2, linetype = "solid") +
      annotate("text", x = obs_stat, y = Inf, vjust = 2, hjust = -0.1,
               label = sprintf("T_obs = %.3f\np = %.4f", obs_stat, pval),
               color = "red", size = 4, fontface = "bold") +
      labs(
        x = "Bootstrap t-statistic",
        y = "Count",
        title = sprintf("Bootstrap Distribution (sim_id = %d, n=%d, phi=%.2f, %s)",
                        d$sim_id[1], d$n[1], d$phi[1], d$innov_dist[1])
      ) +
      theme_minimal(base_size = 14)

    print(p)
  })

  output$dist_summary <- renderPrint({
    d <- dist_data()
    if (is.null(d)) return(cat("No simulation loaded"))

    boot_dist <- d$boot_dist[[1]]
    cat(sprintf("Simulation ID: %d\n", d$sim_id[1]))
    cat(sprintf("Config: n=%d, phi=%.2f, innov_dist=%s\n", d$n[1], d$phi[1], d$innov_dist[1]))
    cat(sprintf("Observed t-stat: %.4f\n", d$obs_stat[1]))
    cat(sprintf("Bootstrap replicates: %d\n", length(boot_dist)))
    cat(sprintf("Bootstrap mean: %.4f, sd: %.4f\n", mean(boot_dist), sd(boot_dist)))
    cat(sprintf("P-value (stored): %.4f\n", d$pvalue[1]))

    # Recalculate
    nb <- length(boot_dist)
    recalc_p <- (sum(abs(boot_dist) >= abs(d$obs_stat[1])) + 1) / (nb + 1)
    cat(sprintf("P-value (recalculated): %.4f\n", recalc_p))
  })

  # ---------------------------------------------------------------------------
  # Tab 4: P-value Diagnostics
  # ---------------------------------------------------------------------------

  # Filtered p-value data (independent of sidebar)
  pval_data <- reactive({
    sql <- "SELECT pvalue, pvalue_asymp, pvalue_adj, n, phi, innov_dist FROM simulations WHERE 1=1"
    params <- list()

    if (length(input$pval_n_filter) > 0) {
      ph <- paste(rep("?", length(input$pval_n_filter)), collapse = ", ")
      sql <- paste0(sql, " AND n IN (", ph, ")")
      params <- c(params, as.list(as.integer(input$pval_n_filter)))
    }
    if (length(input$pval_phi_filter) > 0) {
      ph <- paste(rep("?", length(input$pval_phi_filter)), collapse = ", ")
      sql <- paste0(sql, " AND phi IN (", ph, ")")
      params <- c(params, as.list(as.numeric(input$pval_phi_filter)))
    }
    if (length(input$pval_innov_filter) > 0) {
      ph <- paste(rep("?", length(input$pval_innov_filter)), collapse = ", ")
      sql <- paste0(sql, " AND innov_dist IN (", ph, ")")
      params <- c(params, as.list(input$pval_innov_filter))
    }

    tryCatch(
      dbGetQuery(con, sql, params = params),
      error = function(e) data.frame()
    )
  })

  output$pval_status <- renderText({
    df <- pval_data()
    if (nrow(df) == 0) return("No data")
    n_configs <- nrow(unique(df[, c("n", "phi", "innov_dist"), drop = FALSE]))
    paste0(format(nrow(df), big.mark = ","), " sims, ", n_configs, " configs")
  })

  # P-value histogram: side-by-side for 3 methods
  output$pval_hist <- renderPlot({
    df <- pval_data()
    if (nrow(df) == 0) {
      plot.new()
      text(0.5, 0.5, "No data matching filters", cex = 1.2, col = "grey50")
      return()
    }

    # Build long-format: one row per (sim, method)
    long_rows <- list()
    if ("pvalue" %in% names(df) && !all(is.na(df$pvalue))) {
      long_rows <- c(long_rows, list(data.frame(
        pval = df$pvalue[!is.na(df$pvalue)], method = "COB (Bootstrap)",
        stringsAsFactors = FALSE
      )))
    }
    if ("pvalue_asymp" %in% names(df) && !all(is.na(df$pvalue_asymp))) {
      long_rows <- c(long_rows, list(data.frame(
        pval = df$pvalue_asymp[!is.na(df$pvalue_asymp)], method = "CO (Asymptotic)",
        stringsAsFactors = FALSE
      )))
    }
    if ("pvalue_adj" %in% names(df) && !all(is.na(df$pvalue_adj))) {
      long_rows <- c(long_rows, list(data.frame(
        pval = df$pvalue_adj[!is.na(df$pvalue_adj)], method = "COBA (Adjusted)",
        stringsAsFactors = FALSE
      )))
    }

    if (length(long_rows) == 0) {
      plot.new()
      text(0.5, 0.5, "No p-value data available", cex = 1.2, col = "grey50")
      return()
    }

    pval_long <- do.call(rbind, long_rows)

    p <- ggplot(pval_long, aes(x = pval)) +
      geom_histogram(bins = 20, fill = "#4292c6", color = "white", alpha = 0.8,
                     boundary = 0) +
      facet_wrap(~ method, scales = "free_y") +
      labs(x = "P-value", y = "Count",
           title = "P-value Distributions (uniform = well-calibrated)") +
      theme_minimal(base_size = 13) +
      theme(strip.text = element_text(face = "bold"))

    # Add uniform reference line (expected count per bin)
    n_per_method <- table(pval_long$method)
    ref_df <- data.frame(
      method = names(n_per_method),
      yint = as.numeric(n_per_method) / 20,
      stringsAsFactors = FALSE
    )
    p <- p + geom_hline(data = ref_df, aes(yintercept = yint),
                        color = "red", linetype = "dashed", linewidth = 0.6)

    print(p)
  })

  # P-value QQ plot against U(0,1)
  output$pval_qq <- renderPlot({
    df <- pval_data()
    if (nrow(df) == 0) {
      plot.new()
      text(0.5, 0.5, "No data matching filters", cex = 1.2, col = "grey50")
      return()
    }

    # Build long-format for QQ
    qq_rows <- list()
    methods <- c(pvalue = "COB (Bootstrap)", pvalue_asymp = "CO (Asymptotic)",
                 pvalue_adj = "COBA (Adjusted)")

    for (col in names(methods)) {
      if (col %in% names(df)) {
        vals <- df[[col]][!is.na(df[[col]])]
        if (length(vals) > 0) {
          n_vals <- length(vals)
          theoretical <- (seq_len(n_vals) - 0.5) / n_vals
          qq_rows <- c(qq_rows, list(data.frame(
            theoretical = theoretical,
            observed = sort(vals),
            method = methods[col],
            stringsAsFactors = FALSE
          )))
        }
      }
    }

    if (length(qq_rows) == 0) {
      plot.new()
      text(0.5, 0.5, "No p-value data available", cex = 1.2, col = "grey50")
      return()
    }

    qq_df <- do.call(rbind, qq_rows)

    p <- ggplot(qq_df, aes(x = theoretical, y = observed)) +
      geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
      geom_point(size = 0.8, alpha = 0.5, color = "#2171b5") +
      facet_wrap(~ method) +
      labs(x = "Theoretical Quantiles (Uniform)",
           y = "Observed P-value Quantiles",
           title = "QQ Plot: P-values vs Uniform(0, 1)") +
      coord_equal() +
      theme_minimal(base_size = 13) +
      theme(strip.text = element_text(face = "bold"))

    print(p)
  })

  # P-value scatter: bootstrap vs asymptotic, bootstrap vs COBA
  output$pval_scatter <- renderPlot({
    df <- pval_data()
    if (nrow(df) == 0) {
      plot.new()
      text(0.5, 0.5, "No data matching filters", cex = 1.2, col = "grey50")
      return()
    }

    has_asymp <- "pvalue_asymp" %in% names(df) && !all(is.na(df$pvalue_asymp))
    has_adj <- "pvalue_adj" %in% names(df) && !all(is.na(df$pvalue_adj))
    has_boot <- "pvalue" %in% names(df) && !all(is.na(df$pvalue))

    if (!has_boot || (!has_asymp && !has_adj)) {
      plot.new()
      text(0.5, 0.5, "Need at least 2 methods with p-values", cex = 1.2, col = "grey50")
      return()
    }

    scatter_rows <- list()
    if (has_asymp) {
      idx <- !is.na(df$pvalue) & !is.na(df$pvalue_asymp)
      scatter_rows <- c(scatter_rows, list(data.frame(
        boot_pval = df$pvalue[idx],
        other_pval = df$pvalue_asymp[idx],
        comparison = "COB vs CO (Asymptotic)",
        stringsAsFactors = FALSE
      )))
    }
    if (has_adj) {
      idx <- !is.na(df$pvalue) & !is.na(df$pvalue_adj)
      scatter_rows <- c(scatter_rows, list(data.frame(
        boot_pval = df$pvalue[idx],
        other_pval = df$pvalue_adj[idx],
        comparison = "COB vs COBA (Adjusted)",
        stringsAsFactors = FALSE
      )))
    }

    scatter_df <- do.call(rbind, scatter_rows)

    p <- ggplot(scatter_df, aes(x = boot_pval, y = other_pval)) +
      geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
      geom_point(size = 0.6, alpha = 0.2, color = "#2171b5") +
      facet_wrap(~ comparison) +
      labs(x = "Bootstrap P-value (COB)",
           y = "Comparison Method P-value",
           title = "P-value Agreement Between Methods") +
      coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
      theme_minimal(base_size = 13) +
      theme(strip.text = element_text(face = "bold"))

    print(p)
  })

  # ---------------------------------------------------------------------------
  # Tab 6: Analysis Grids
  # ---------------------------------------------------------------------------

  # Grid 1: By Sample Size (fix innov_dist + phi, vary n)
  grid1_data <- reactive({
    req(input$grid1_innov, input$grid1_phi)
    grid_query(innov_dist = input$grid1_innov,
               phi = input$grid1_phi,
               order_col = "n",
               batch_id = input$grid1_batch)
  })

  output$grid1_table <- renderDT({
    format_grid_dt(grid1_data(), row_label_col = "n", row_label_name = "n")
  })

  # Grid 2: By Phi (fix innov_dist + n, vary phi)
  grid2_data <- reactive({
    req(input$grid2_innov, input$grid2_n)
    grid_query(innov_dist = input$grid2_innov,
               n = input$grid2_n,
               order_col = "phi",
               batch_id = input$grid2_batch)
  })

  output$grid2_table <- renderDT({
    format_grid_dt(grid2_data(), row_label_col = "phi", row_label_name = "Phi")
  })

  # Grid 3: By Distribution (fix phi + n, vary innov_dist)
  grid3_data <- reactive({
    req(input$grid3_phi, input$grid3_n)
    grid_query(phi = input$grid3_phi,
               n = input$grid3_n,
               order_col = "innov_dist",
               batch_id = input$grid3_batch)
  })

  output$grid3_table <- renderDT({
    format_grid_dt(grid3_data(), row_label_col = "innov_dist",
                   row_label_name = "Distribution")
  })

  # ---------------------------------------------------------------------------
  # Tab 7: Parallel Coordinates
  # ---------------------------------------------------------------------------

  # Filtered data for PC plots (independent of sidebar)
  pc_data <- reactive({
    sql <- "SELECT * FROM v_rejection_rates WHERE 1=1"
    params <- list()

    if (length(input$pc_n_filter) > 0) {
      ph <- paste(rep("?", length(input$pc_n_filter)), collapse = ", ")
      sql <- paste0(sql, " AND n IN (", ph, ")")
      params <- c(params, as.list(as.integer(input$pc_n_filter)))
    }
    if (length(input$pc_phi_filter) > 0) {
      ph <- paste(rep("?", length(input$pc_phi_filter)), collapse = ", ")
      sql <- paste0(sql, " AND phi IN (", ph, ")")
      params <- c(params, as.list(as.numeric(input$pc_phi_filter)))
    }
    if (length(input$pc_innov_filter) > 0) {
      ph <- paste(rep("?", length(input$pc_innov_filter)), collapse = ", ")
      sql <- paste0(sql, " AND innov_dist IN (", ph, ")")
      params <- c(params, as.list(input$pc_innov_filter))
    }

    sql <- paste0(sql, " ORDER BY innov_dist, n, phi")
    tryCatch(
      dbGetQuery(con, sql, params = params),
      error = function(e) data.frame()
    )
  })

  # Helper: build a parallel coordinates ggplot
  # cols = named list of column_name -> axis_label
  # color_col = column to use for color
  # color_label = legend title
  build_parcoord <- function(df, cols, color_col, color_label, title = "") {
    if (nrow(df) == 0) return(NULL)

    col_names <- names(cols)
    axis_labels <- unname(cols)
    n_axes <- length(col_names)

    # Assign config ID to each row
    df$config_id <- seq_len(nrow(df))

    # Scale each column to [0, 1] and store original range
    ranges <- list()
    scaled <- data.frame(config_id = df$config_id)

    for (col in col_names) {
      vals <- as.numeric(df[[col]])
      rng <- range(vals, na.rm = TRUE)
      ranges[[col]] <- rng
      if (rng[1] == rng[2]) {
        scaled[[col]] <- 0.5
      } else {
        scaled[[col]] <- (vals - rng[1]) / (rng[2] - rng[1])
      }
    }

    # Build long-format data for geom_line
    long_rows <- vector("list", nrow(df) * n_axes)
    idx <- 1
    for (i in seq_len(nrow(df))) {
      for (j in seq_len(n_axes)) {
        long_rows[[idx]] <- data.frame(
          config_id = i,
          axis_num = j,
          axis_label = axis_labels[j],
          value = scaled[[col_names[j]]][i],
          color_var = as.character(df[[color_col]][i]),
          stringsAsFactors = FALSE
        )
        idx <- idx + 1
      }
    }
    long_df <- do.call(rbind, long_rows)
    long_df$axis_num <- as.integer(long_df$axis_num)

    # Build annotation data for axis tick labels (min/max at each axis)
    tick_labels <- unlist(lapply(col_names, function(col) {
      rng <- ranges[[col]]
      c(format(rng[1], digits = 3), format(rng[2], digits = 3))
    }))
    tick_df <- data.frame(
      axis_num = rep(seq_len(n_axes), each = 2),
      y = rep(c(0, 1), times = n_axes),
      label = tick_labels,
      stringsAsFactors = FALSE
    )

    p <- ggplot(long_df, aes(x = axis_num, y = value,
                              group = config_id, color = color_var)) +
      geom_line(alpha = 0.6, linewidth = 0.7) +
      geom_point(size = 1.5, alpha = 0.7) +
      # Reference line at nominal 0.05 for rate axes
      scale_x_continuous(breaks = seq_len(n_axes), labels = axis_labels,
                         limits = c(0.5, n_axes + 0.5)) +
      scale_y_continuous(limits = c(-0.08, 1.08), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
      # Axis tick labels showing original values
      geom_text(data = tick_df,
                aes(x = axis_num, y = y, label = label),
                inherit.aes = FALSE, size = 2.8, color = "grey40",
                hjust = 1.2, fontface = "italic") +
      # Vertical axis lines
      geom_vline(xintercept = seq_len(n_axes), color = "grey80", linewidth = 0.3) +
      labs(x = NULL, y = "Scaled Value [0, 1]", color = color_label, title = title) +
      theme_minimal(base_size = 13) +
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom"
      )

    p
  }

  # Plot 1: Full Profile (n -> phi -> CO -> COB -> COBA, color by innov_dist)
  output$pc_full_plot <- renderPlot({
    df <- pc_data()
    if (nrow(df) == 0) {
      plot.new()
      text(0.5, 0.5, "No data matching filters", cex = 1.2, col = "grey50")
      return()
    }

    p <- build_parcoord(
      df,
      cols = c(n = "n", phi = "phi",
               reject_asymp_05 = "CO Rate",
               reject_05 = "COB Rate",
               reject_adj_05 = "COBA Rate"),
      color_col = "innov_dist",
      color_label = "Innovation",
      title = "Full Profile: n \u2192 phi \u2192 CO \u2192 COB \u2192 COBA"
    )

    if (!is.null(p)) print(p)
  })

  # Plot 2: Method Comparison (CO -> COB -> COBA, color by user selection)
  output$pc_method_plot <- renderPlot({
    df <- pc_data()
    if (nrow(df) == 0) {
      plot.new()
      text(0.5, 0.5, "No data matching filters", cex = 1.2, col = "grey50")
      return()
    }

    color_col <- input$pc_color_by
    color_labels <- c(innov_dist = "Innovation", n = "Sample Size", phi = "Phi")
    # Ensure color column is character for consistent legend
    df[[color_col]] <- as.character(df[[color_col]])

    p <- build_parcoord(
      df,
      cols = c(reject_asymp_05 = "CO Rate",
               reject_05 = "COB Rate",
               reject_adj_05 = "COBA Rate"),
      color_col = color_col,
      color_label = color_labels[color_col],
      title = "Method Comparison: CO \u2192 COB \u2192 COBA"
    )

    if (!is.null(p)) print(p)
  })

  # Plot 3: Single-Method Sensitivity (n -> phi -> rate, color by innov_dist)
  output$pc_sensitivity_plot <- renderPlot({
    df <- pc_data()
    if (nrow(df) == 0) {
      plot.new()
      text(0.5, 0.5, "No data matching filters", cex = 1.2, col = "grey50")
      return()
    }

    rate_col <- input$pc_method
    rate_labels <- c(reject_05 = "COB Rate",
                     reject_asymp_05 = "CO Rate",
                     reject_adj_05 = "COBA Rate")

    cols <- c("n", "phi", rate_col)
    names(cols) <- NULL
    named_cols <- setNames(cols, cols)
    named_cols <- c(n = "n", phi = "phi")
    named_cols[[rate_col]] <- rate_labels[rate_col]

    p <- build_parcoord(
      df,
      cols = named_cols,
      color_col = "innov_dist",
      color_label = "Innovation",
      title = paste0("Sensitivity: n \u2192 phi \u2192 ", rate_labels[rate_col])
    )

    if (!is.null(p)) print(p)
  })

  # ---------------------------------------------------------------------------
  # Tab 8: Diagnostics
  # ---------------------------------------------------------------------------

  # Filtered diagnostics data (independent of sidebar)
  diag_data <- reactive({
    sql <- paste0(
      "SELECT null_ar_order, null_ar_phi[1] AS null_phi1, n, phi, innov_dist ",
      "FROM simulations WHERE null_ar_order IS NOT NULL"
    )
    params <- list()

    if (length(input$diag_n_filter) > 0) {
      ph <- paste(rep("?", length(input$diag_n_filter)), collapse = ", ")
      sql <- paste0(sql, " AND n IN (", ph, ")")
      params <- c(params, as.list(as.integer(input$diag_n_filter)))
    }
    if (length(input$diag_phi_filter) > 0) {
      ph <- paste(rep("?", length(input$diag_phi_filter)), collapse = ", ")
      sql <- paste0(sql, " AND phi IN (", ph, ")")
      params <- c(params, as.list(as.numeric(input$diag_phi_filter)))
    }
    if (length(input$diag_innov_filter) > 0) {
      ph <- paste(rep("?", length(input$diag_innov_filter)), collapse = ", ")
      sql <- paste0(sql, " AND innov_dist IN (", ph, ")")
      params <- c(params, as.list(input$diag_innov_filter))
    }

    tryCatch(
      dbGetQuery(con, sql, params = params),
      error = function(e) data.frame()
    )
  })

  # Status text
  output$diag_status <- renderText({
    df <- diag_data()
    if (nrow(df) == 0) return("No data")
    n_configs <- nrow(unique(df[, c("n", "phi", "innov_dist"), drop = FALSE]))
    paste0(format(nrow(df), big.mark = ","), " sims, ", n_configs, " configs")
  })

  # AR order distribution bar chart
  output$diag_ar_order_plot <- renderPlot({
    df <- diag_data()
    if (nrow(df) == 0) {
      plot.new()
      text(0.5, 0.5, "No data matching filters", cex = 1.2, col = "grey50")
      return()
    }

    df$null_ar_order <- factor(df$null_ar_order)

    # Build config label for faceting
    df$config <- sprintf("n=%s, phi=%s, %s", df$n, df$phi, df$innov_dist)

    # If many configs, facet; otherwise single plot
    n_configs <- length(unique(df$config))

    p <- ggplot(df, aes(x = null_ar_order)) +
      geom_bar(fill = "#4292c6", color = "white", alpha = 0.8) +
      labs(
        x = "Fitted AR Order",
        y = "Count",
        title = "Distribution of Fitted AR Orders Under the Null"
      ) +
      theme_minimal(base_size = 14)

    if (n_configs > 1 && n_configs <= 12) {
      p <- p + facet_wrap(~ config, scales = "free_y")
    } else if (n_configs > 12) {
      p <- p + facet_wrap(~ config, scales = "free_y", ncol = 4)
    }

    print(p)
  })

  # Estimated phi distribution
  output$diag_phi_plot <- renderPlot({
    df <- diag_data()
    if (nrow(df) == 0 || all(is.na(df$null_phi1))) {
      plot.new()
      text(0.5, 0.5, "No data matching filters", cex = 1.2, col = "grey50")
      return()
    }

    df <- df[!is.na(df$null_phi1), ]
    df$true_phi_label <- paste0("True phi = ", df$phi)

    p <- ggplot(df, aes(x = null_phi1)) +
      geom_histogram(bins = 50, fill = "#4292c6", color = "white", alpha = 0.8) +
      geom_vline(aes(xintercept = phi), color = "red", linewidth = 1, linetype = "dashed") +
      labs(
        x = "Estimated First AR Coefficient",
        y = "Count",
        title = "Distribution of Estimated Null AR(1) Coefficient"
      ) +
      theme_minimal(base_size = 14)

    n_configs <- length(unique(df$true_phi_label))
    if (n_configs > 1) {
      p <- p + facet_wrap(~ true_phi_label, scales = "free")
    }

    print(p)
  })

  # MC Convergence: running rejection rate as sims accumulate
  output$diag_convergence_plot <- renderPlot({
    pval_col <- input$diag_conv_method
    method_labels <- c(pvalue = "Bootstrap (COB)",
                       pvalue_asymp = "Asymptotic (CO)",
                       pvalue_adj = "COBA")

    # Query raw sims ordered by sim_id, with diag filters
    sql <- paste0("SELECT sim_id, ", pval_col, " AS pval, n, phi, innov_dist ",
                  "FROM simulations WHERE ", pval_col, " IS NOT NULL")
    params <- list()

    if (length(input$diag_n_filter) > 0) {
      ph <- paste(rep("?", length(input$diag_n_filter)), collapse = ", ")
      sql <- paste0(sql, " AND n IN (", ph, ")")
      params <- c(params, as.list(as.integer(input$diag_n_filter)))
    }
    if (length(input$diag_phi_filter) > 0) {
      ph <- paste(rep("?", length(input$diag_phi_filter)), collapse = ", ")
      sql <- paste0(sql, " AND phi IN (", ph, ")")
      params <- c(params, as.list(as.numeric(input$diag_phi_filter)))
    }
    if (length(input$diag_innov_filter) > 0) {
      ph <- paste(rep("?", length(input$diag_innov_filter)), collapse = ", ")
      sql <- paste0(sql, " AND innov_dist IN (", ph, ")")
      params <- c(params, as.list(input$diag_innov_filter))
    }

    sql <- paste0(sql, " ORDER BY n, phi, innov_dist, sim_id")
    df <- tryCatch(dbGetQuery(con, sql, params = params), error = function(e) data.frame())

    if (nrow(df) == 0) {
      plot.new()
      text(0.5, 0.5, "No data matching filters", cex = 1.2, col = "grey50")
      return()
    }

    # Compute running rejection rate per config
    df$config <- sprintf("%s / n=%s / phi=%s", df$innov_dist, df$n, df$phi)
    configs <- unique(df$config)

    conv_rows <- list()
    for (cfg in configs) {
      sub <- df[df$config == cfg, ]
      rejected <- sub$pval < 0.05
      cum_rate <- cumsum(rejected) / seq_along(rejected)
      conv_rows <- c(conv_rows, list(data.frame(
        sim_num = seq_along(cum_rate),
        cum_reject_rate = cum_rate,
        config = cfg,
        stringsAsFactors = FALSE
      )))
    }
    conv_df <- do.call(rbind, conv_rows)

    p <- ggplot(conv_df, aes(x = sim_num, y = cum_reject_rate, color = config)) +
      geom_line(alpha = 0.7, linewidth = 0.6) +
      geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey40") +
      labs(
        x = "Simulation Number",
        y = "Cumulative Rejection Rate",
        color = "Configuration",
        title = paste0("MC Convergence: ", method_labels[pval_col])
      ) +
      theme_minimal(base_size = 13) +
      theme(legend.position = "bottom",
            legend.text = element_text(size = 8))

    print(p)
  })

  # Batch consistency: boxplot of per-batch rejection rates
  output$diag_batch_plot <- renderPlot({
    rate_col <- input$diag_batch_method
    rate_labels <- c(reject_05 = "Bootstrap (COB)",
                     reject_asymp_05 = "Asymptotic (CO)",
                     reject_adj_05 = "COBA")

    # Query per-batch view with diag filters
    sql <- "SELECT * FROM v_rejection_rates_by_batch WHERE 1=1"
    params <- list()

    if (length(input$diag_n_filter) > 0) {
      ph <- paste(rep("?", length(input$diag_n_filter)), collapse = ", ")
      sql <- paste0(sql, " AND n IN (", ph, ")")
      params <- c(params, as.list(as.integer(input$diag_n_filter)))
    }
    if (length(input$diag_phi_filter) > 0) {
      ph <- paste(rep("?", length(input$diag_phi_filter)), collapse = ", ")
      sql <- paste0(sql, " AND phi IN (", ph, ")")
      params <- c(params, as.list(as.numeric(input$diag_phi_filter)))
    }
    if (length(input$diag_innov_filter) > 0) {
      ph <- paste(rep("?", length(input$diag_innov_filter)), collapse = ", ")
      sql <- paste0(sql, " AND innov_dist IN (", ph, ")")
      params <- c(params, as.list(input$diag_innov_filter))
    }

    sql <- paste0(sql, " ORDER BY innov_dist, n, phi")
    df <- tryCatch(dbGetQuery(con, sql, params = params), error = function(e) data.frame())

    if (nrow(df) == 0 || !rate_col %in% names(df)) {
      plot.new()
      text(0.5, 0.5, "No per-batch data available", cex = 1.2, col = "grey50")
      return()
    }

    df$config <- sprintf("%s / n=%s / phi=%s", df$innov_dist, df$n, df$phi)

    p <- ggplot(df, aes(x = config, y = .data[[rate_col]])) +
      geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey40") +
      geom_boxplot(fill = "#d4edda", alpha = 0.6, outlier.shape = NA) +
      geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "#2171b5") +
      labs(
        x = NULL,
        y = "Rejection Rate",
        title = paste0("Batch-to-Batch Consistency: ", rate_labels[rate_col])
      ) +
      theme_minimal(base_size = 13) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

    print(p)
  })

  # Test statistic distribution: histogram of obs_stat per config
  output$diag_tstat_plot <- renderPlot({
    sql <- "SELECT obs_stat, n, phi, innov_dist FROM simulations WHERE obs_stat IS NOT NULL"
    params <- list()

    if (length(input$diag_n_filter) > 0) {
      ph <- paste(rep("?", length(input$diag_n_filter)), collapse = ", ")
      sql <- paste0(sql, " AND n IN (", ph, ")")
      params <- c(params, as.list(as.integer(input$diag_n_filter)))
    }
    if (length(input$diag_phi_filter) > 0) {
      ph <- paste(rep("?", length(input$diag_phi_filter)), collapse = ", ")
      sql <- paste0(sql, " AND phi IN (", ph, ")")
      params <- c(params, as.list(as.numeric(input$diag_phi_filter)))
    }
    if (length(input$diag_innov_filter) > 0) {
      ph <- paste(rep("?", length(input$diag_innov_filter)), collapse = ", ")
      sql <- paste0(sql, " AND innov_dist IN (", ph, ")")
      params <- c(params, as.list(input$diag_innov_filter))
    }

    df <- tryCatch(dbGetQuery(con, sql, params = params), error = function(e) data.frame())

    if (nrow(df) == 0) {
      plot.new()
      text(0.5, 0.5, "No data matching filters", cex = 1.2, col = "grey50")
      return()
    }

    df$config <- sprintf("%s / n=%s / phi=%s", df$innov_dist, df$n, df$phi)
    n_configs <- length(unique(df$config))

    p <- ggplot(df, aes(x = obs_stat)) +
      geom_histogram(bins = 50, fill = "#4292c6", color = "white", alpha = 0.8) +
      labs(
        x = "Observed Test Statistic",
        y = "Count",
        title = "Distribution of Observed Test Statistics"
      ) +
      theme_minimal(base_size = 14)

    if (n_configs > 1 && n_configs <= 12) {
      p <- p + facet_wrap(~ config, scales = "free")
    } else if (n_configs > 12) {
      p <- p + facet_wrap(~ config, scales = "free", ncol = 4)
    }

    print(p)
  })
}

# =============================================================================
# Run
# =============================================================================

shinyApp(ui = ui, server = server)
