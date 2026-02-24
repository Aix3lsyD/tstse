# Module: Analysis Grids (Tab 6)
# Extracted from the monolithic app.R

mod_analysis_grids_ui <- function(id) {
  ns <- NS(id)

  tabPanel(
    "Analysis Grids",
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

    # Grid 1: By Sample Size
    wellPanel(
      h5("By Sample Size"),
      p(class = "text-body-secondary", style = "margin-bottom: 10px;",
        "Fix innovation distribution and phi; rows vary by n"),
      fluidRow(
        column(3,
          selectInput(ns("grid1_innov"), "Innovation Distribution:",
                      choices = NULL)
        ),
        column(3,
          selectInput(ns("grid1_phi"), "Phi:",
                      choices = NULL)
        ),
        column(3,
          selectInput(ns("grid1_batch"), "Batch:",
                      choices = c("All (Pooled)" = "all"))
        )
      ),
      DTOutput(ns("grid1_table"))
    ),

    # Grid 2: By Phi
    wellPanel(
      h5("By Phi"),
      p(class = "text-body-secondary", style = "margin-bottom: 10px;",
        "Fix innovation distribution and n; rows vary by phi"),
      fluidRow(
        column(3,
          selectInput(ns("grid2_innov"), "Innovation Distribution:",
                      choices = NULL)
        ),
        column(3,
          selectInput(ns("grid2_n"), "Sample Size (n):",
                      choices = NULL)
        ),
        column(3,
          selectInput(ns("grid2_batch"), "Batch:",
                      choices = c("All (Pooled)" = "all"))
        )
      ),
      DTOutput(ns("grid2_table"))
    ),

    # Grid 3: By Distribution
    wellPanel(
      h5("By Innovation Distribution"),
      p(class = "text-body-secondary", style = "margin-bottom: 10px;",
        "Fix phi and n; rows vary by innovation distribution"),
      fluidRow(
        column(3,
          selectInput(ns("grid3_phi"), "Phi:",
                      choices = NULL)
        ),
        column(3,
          selectInput(ns("grid3_n"), "Sample Size (n):",
                      choices = NULL)
        ),
        column(3,
          selectInput(ns("grid3_batch"), "Batch:",
                      choices = c("All (Pooled)" = "all"))
        )
      ),
      plotOutput(ns("grid3_barchart"), height = "300px"),
      hr(),
      DTOutput(ns("grid3_table"))
    )
  )
}

mod_analysis_grids_server <- function(id, con, init_choices) {
  moduleServer(id, function(input, output, session) {

    # Populate filter dropdowns from initial choices
    updateSelectInput(session, "grid1_innov",
                      choices = init_choices$innov,
                      selected = init_choices$innov[1])
    updateSelectInput(session, "grid1_phi",
                      choices = init_choices$phi,
                      selected = init_choices$phi[1])
    updateSelectInput(session, "grid2_innov",
                      choices = init_choices$innov,
                      selected = init_choices$innov[1])
    updateSelectInput(session, "grid2_n",
                      choices = init_choices$n,
                      selected = init_choices$n[1])
    updateSelectInput(session, "grid3_phi",
                      choices = init_choices$phi,
                      selected = init_choices$phi[1])
    updateSelectInput(session, "grid3_n",
                      choices = init_choices$n,
                      selected = init_choices$n[1])

    # Update per-grid batch dropdowns when their config filters change
    observe({
      req(input$grid1_innov, input$grid1_phi)
      choices <- grid_batch_choices(con, innov_dist = input$grid1_innov, phi = input$grid1_phi)
      sel <- if (input$grid1_batch %in% choices) input$grid1_batch else "all"
      updateSelectInput(session, "grid1_batch", choices = choices, selected = sel)
    })

    observe({
      req(input$grid2_innov, input$grid2_n)
      choices <- grid_batch_choices(con, innov_dist = input$grid2_innov, n = input$grid2_n)
      sel <- if (input$grid2_batch %in% choices) input$grid2_batch else "all"
      updateSelectInput(session, "grid2_batch", choices = choices, selected = sel)
    })

    observe({
      req(input$grid3_phi, input$grid3_n)
      choices <- grid_batch_choices(con, phi = input$grid3_phi, n = input$grid3_n)
      sel <- if (input$grid3_batch %in% choices) input$grid3_batch else "all"
      updateSelectInput(session, "grid3_batch", choices = choices, selected = sel)
    })

    # Grid 1: By Sample Size (fix innov_dist + phi, vary n)
    grid1_data <- reactive({
      req(input$grid1_innov, input$grid1_phi)
      grid_query(con,
                 innov_dist = input$grid1_innov,
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
      grid_query(con,
                 innov_dist = input$grid2_innov,
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
      grid_query(con,
                 phi = input$grid3_phi,
                 n = input$grid3_n,
                 order_col = "innov_dist",
                 batch_id = input$grid3_batch)
    })

    # Grid 3 grouped bar chart
    output$grid3_barchart <- renderPlot(bg = "transparent", {
      df <- grid3_data()
      if (nrow(df) == 0) {
        plot.new()
        text(0.5, 0.5, "No data for this combination", cex = 1.2, col = "grey50")
        return()
      }

      # Build long-format data for bar chart
      bar_rows <- list()
      for (i in seq_len(nrow(df))) {
        if ("reject_asymp_05" %in% names(df)) {
          bar_rows <- c(bar_rows, list(data.frame(
            distribution = df$innov_dist[i], method = "Asymptotic (CO)",
            rate = df$reject_asymp_05[i], stringsAsFactors = FALSE)))
        }
        if ("reject_05" %in% names(df)) {
          bar_rows <- c(bar_rows, list(data.frame(
            distribution = df$innov_dist[i], method = "Bootstrap (COB)",
            rate = df$reject_05[i], stringsAsFactors = FALSE)))
        }
        if ("reject_adj_05" %in% names(df)) {
          bar_rows <- c(bar_rows, list(data.frame(
            distribution = df$innov_dist[i], method = "COBA",
            rate = df$reject_adj_05[i], stringsAsFactors = FALSE)))
        }
      }
      bar_df <- do.call(rbind, bar_rows)
      bar_df$method <- factor(bar_df$method,
                              levels = c("Asymptotic (CO)", "Bootstrap (COB)", "COBA"))

      p <- ggplot(bar_df, aes(x = rate, y = distribution, fill = method)) +
        geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
        geom_vline(xintercept = 0.05, linetype = "dashed", linewidth = 0.8) +
        scale_fill_manual(values = c("Asymptotic (CO)" = "#e41a1c",
                                     "Bootstrap (COB)" = "#377eb8",
                                     "COBA" = "#4daf4a")) +
        labs(x = "Rejection Rate", y = NULL, fill = "Method",
             title = "Rejection Rates by Innovation Distribution") +
        viewer_plot_theme(base_size = 13) +
        theme(legend.position = "bottom")

      print(p)
    })

    output$grid3_table <- renderDT({
      format_grid_dt(grid3_data(), row_label_col = "innov_dist",
                     row_label_name = "Distribution")
    })
  })
}
