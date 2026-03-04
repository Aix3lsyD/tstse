# Capstone Sub-Tab: Analysis Grids (Tab 4)
# Three rejection rate comparison grids + bar chart

mod_capstone_grids_ui <- function(ns) {
  tabPanel("Analysis Grids",
    br(),
    h4("Rejection Rate Comparison Grids"),
    p(class = "text-body-secondary",
      "Each grid fixes two dimensions and varies the third. ",
      "Rates are color-coded: ",
      span(style = "background-color:#fff3cd; color:#212529; font-weight:600; padding:2px 6px; border-radius:3px;", "< 0.03"),
      " ",
      span(style = "background-color:#d4edda; color:#212529; font-weight:600; padding:2px 6px; border-radius:3px;", "0.03 \u2013 0.07"),
      " ",
      span(style = "background-color:#f8d7da; color:#212529; font-weight:600; padding:2px 6px; border-radius:3px;", "> 0.07")
    ),
    fluidRow(
      column(3, actionButton(ns("clear_results"), "Clear Results",
                             icon = icon("eraser"),
                             class = "btn-sm btn-outline-warning")),
      column(3, downloadButton(ns("download_results_csv"), "Download CSV",
                                class = "btn-sm btn-outline-secondary"))
    ),
    hr(),

    # Grid 1: By Sample Size
    wellPanel(
      h5("By Sample Size"),
      p(class = "text-body-secondary", style = "margin-bottom: 10px;",
        "Fix innovation distribution and phi; rows vary by n"),
      fluidRow(
        column(4, selectInput(ns("res_grid1_innov"), "Innovation Distribution:",
                              choices = NULL)),
        column(4, selectInput(ns("res_grid1_phi"), "Phi:", choices = NULL))
      ),
      DTOutput(ns("res_grid1_table"))
    ),

    # Grid 2: By Phi
    wellPanel(
      h5("By Phi"),
      p(class = "text-body-secondary", style = "margin-bottom: 10px;",
        "Fix innovation distribution and n; rows vary by phi"),
      fluidRow(
        column(4, selectInput(ns("res_grid2_innov"), "Innovation Distribution:",
                              choices = NULL)),
        column(4, selectInput(ns("res_grid2_n"), "Sample Size (n):", choices = NULL))
      ),
      DTOutput(ns("res_grid2_table"))
    ),

    # Grid 3: By Distribution
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

mod_capstone_grids_server <- function(input, output, session,
                                      cap_combined_results, run_results_rv) {

  # Populate filter dropdowns
  observe({
    df <- cap_combined_results()
    if (nrow(df) == 0) {
      empty <- character(0)
      updateSelectInput(session, "res_grid1_innov", choices = empty)
      updateSelectInput(session, "res_grid1_phi", choices = empty)
      updateSelectInput(session, "res_grid2_innov", choices = empty)
      updateSelectInput(session, "res_grid2_n", choices = empty)
      updateSelectInput(session, "res_grid3_phi", choices = empty)
      updateSelectInput(session, "res_grid3_n", choices = empty)
      return()
    }
    innov_ch <- sort(unique(df$innov_dist))
    phi_ch <- as.character(sort(unique(df$phi)))
    n_ch <- as.character(sort(unique(df$n)))

    updateSelectInput(session, "res_grid1_innov", choices = innov_ch,
                      selected = innov_ch[1])
    updateSelectInput(session, "res_grid1_phi", choices = phi_ch,
                      selected = phi_ch[1])
    updateSelectInput(session, "res_grid2_innov", choices = innov_ch,
                      selected = innov_ch[1])
    updateSelectInput(session, "res_grid2_n", choices = n_ch,
                      selected = n_ch[1])
    updateSelectInput(session, "res_grid3_phi", choices = phi_ch,
                      selected = phi_ch[1])
    updateSelectInput(session, "res_grid3_n", choices = n_ch,
                      selected = n_ch[1])
  })

  # Helper: filter combined results
  .mem_grid_query <- function(df, innov_dist = NULL, phi = NULL, n = NULL) {
    if (nrow(df) == 0) return(df)
    if (!is.null(innov_dist) && nzchar(innov_dist))
      df <- df[df$innov_dist == innov_dist, , drop = FALSE]
    if (!is.null(phi) && nzchar(phi))
      df <- df[abs(df$phi - as.numeric(phi)) < 1e-9, , drop = FALSE]
    if (!is.null(n) && nzchar(n))
      df <- df[df$n == as.integer(n), , drop = FALSE]
    df
  }

  # Grid 1: By Sample Size
  res_grid1_data <- reactive({
    req(input$res_grid1_innov, input$res_grid1_phi)
    df <- .mem_grid_query(cap_combined_results(),
                           innov_dist = input$res_grid1_innov,
                           phi = input$res_grid1_phi)
    df[order(df$n), , drop = FALSE]
  })
  output$res_grid1_table <- renderDT({
    format_grid_dt(res_grid1_data(), row_label_col = "n", row_label_name = "n")
  })

  # Grid 2: By Phi
  res_grid2_data <- reactive({
    req(input$res_grid2_innov, input$res_grid2_n)
    df <- .mem_grid_query(cap_combined_results(),
                           innov_dist = input$res_grid2_innov,
                           n = input$res_grid2_n)
    df[order(df$phi), , drop = FALSE]
  })
  output$res_grid2_table <- renderDT({
    format_grid_dt(res_grid2_data(), row_label_col = "phi", row_label_name = "Phi")
  })

  # Grid 3: By Distribution
  res_grid3_data <- reactive({
    req(input$res_grid3_phi, input$res_grid3_n)
    df <- .mem_grid_query(cap_combined_results(),
                           phi = input$res_grid3_phi,
                           n = input$res_grid3_n)
    df[order(df$innov_dist), , drop = FALSE]
  })

  output$res_grid3_barchart <- renderPlot(bg = "transparent", {
    df <- res_grid3_data()
    if (nrow(df) == 0) {
      plot.new()
      text(0.5, 0.5, "No data for this combination", cex = 1.2,
           col = viewer_plot_fg())
      return()
    }
    p <- build_rejection_barchart(df, label_col = "innov_dist",
                                   title = "Rejection Rates by Innovation Distribution")
    if (!is.null(p)) print(p)
  })

  output$res_grid3_table <- renderDT({
    format_grid_dt(res_grid3_data(), row_label_col = "innov_dist",
                   row_label_name = "Distribution")
  })

  # Clear results
  observeEvent(input$clear_results, {
    run_results_rv(data.frame(
      n = integer(0), phi = numeric(0), innov_dist = character(0),
      innov_label = character(0), n_sims = integer(0),
      reject_05 = numeric(0), reject_05_se = numeric(0),
      reject_asymp_05 = numeric(0), reject_asymp_05_se = numeric(0),
      reject_adj_05 = numeric(0), reject_adj_05_se = numeric(0),
      stringsAsFactors = FALSE
    ))
    showNotification("In-memory results cleared.", type = "message", duration = 3)
  })

  # CSV download
  output$download_results_csv <- downloadHandler(
    filename = function() {
      paste0("capstone_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      write.csv(cap_combined_results(), file, row.names = FALSE)
    }
  )
}
