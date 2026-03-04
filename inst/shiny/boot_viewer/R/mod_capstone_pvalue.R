# Capstone Sub-Tab: P-value Diagnostics (Tab 6)
# Histograms, QQ plots, scatter comparing p-value methods

mod_capstone_pvalue_ui <- function(ns) {
  tabPanel("P-value Diagnostics",
    br(),
    h4("P-value Diagnostics"),
    p(class = "text-body-secondary",
      "Validates p-value calibration across methods. ",
      "Well-calibrated p-values should follow a Uniform(0,1) distribution under the null."),
    div(class = "plot-controls",
      fluidRow(
        column(4,
          selectInput(ns("pval_cell"), "Grid Cell:", choices = c("All Cells" = "all"))
        ),
        column(4,
          div(style = "margin-top: 25px;",
            textOutput(ns("pval_status")))
        )
      )
    ),
    wellPanel(
      h5("P-value Histograms"),
      plotOutput(ns("pval_hist"), height = "300px")
    ),
    wellPanel(
      h5("P-value QQ Plot (vs Uniform)"),
      plotOutput(ns("pval_qq"), height = "350px")
    ),
    wellPanel(
      h5("P-value Method Comparison"),
      plotOutput(ns("pval_scatter"), height = "350px")
    )
  )
}

mod_capstone_pvalue_server <- function(input, output, session,
                                       cap_sim_data, cap_cell_choices) {

  observe({
    ch <- cap_cell_choices()
    all_choices <- c("All Cells" = "all", ch)
    updateSelectInput(session, "pval_cell", choices = all_choices)
  })

  # Collect p-values from selected cell(s)
  pval_df <- reactive({
    req(input$pval_cell)
    ch <- cap_cell_choices()
    if (length(ch) == 0) return(data.frame())

    keys <- if (input$pval_cell == "all") unname(ch) else input$pval_cell
    rows <- list()
    for (key in keys) {
      sims <- cap_sim_data(key)
      if (is.null(sims) || length(sims) == 0) next
      rows <- c(rows, list(data.frame(
        pvalue = vapply(sims, function(r) r$pvalue, numeric(1)),
        pvalue_asymp = vapply(sims, function(r) r$pvalue_asymp, numeric(1)),
        pvalue_adj = vapply(sims, function(r) {
          if (is.null(r$pvalue_adj)) NA_real_ else r$pvalue_adj
        }, numeric(1)),
        stringsAsFactors = FALSE
      )))
    }
    if (length(rows) == 0) return(data.frame())
    do.call(rbind, rows)
  })

  output$pval_status <- renderText({
    df <- pval_df()
    if (nrow(df) == 0) return("No data")
    paste0(nrow(df), " simulations")
  })

  output$pval_hist <- renderPlot(bg = "transparent", {
    df <- pval_df()
    if (nrow(df) == 0) {
      plot.new(); text(0.5, 0.5, "No data", cex = 1.2, col = viewer_plot_fg()); return()
    }
    p <- plot_pvalue_histogram(df)
    if (!is.null(p)) print(p)
  })

  output$pval_qq <- renderPlot(bg = "transparent", {
    df <- pval_df()
    if (nrow(df) == 0) {
      plot.new(); text(0.5, 0.5, "No data", cex = 1.2, col = viewer_plot_fg()); return()
    }
    p <- plot_pvalue_qq(df)
    if (!is.null(p)) print(p)
  })

  output$pval_scatter <- renderPlot(bg = "transparent", {
    df <- pval_df()
    if (nrow(df) == 0) {
      plot.new(); text(0.5, 0.5, "No data", cex = 1.2, col = viewer_plot_fg()); return()
    }
    p <- plot_pvalue_scatter(df)
    if (!is.null(p)) print(p)
  })
}
