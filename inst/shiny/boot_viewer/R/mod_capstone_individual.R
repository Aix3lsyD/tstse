# Capstone Sub-Tab: Individual Results (Tab 3)
# Per-simulation result table for a selected grid cell

mod_capstone_individual_ui <- function(ns) {
  tabPanel("Individual Results",
    br(),
    h4("Per-Simulation Results"),
    p(class = "text-body-secondary",
      "Select a grid cell to view individual simulation results. ",
      "Click a row to view its bootstrap distribution on the Bootstrap tab."),
    div(class = "plot-controls",
      fluidRow(
        column(4,
          selectInput(ns("indiv_cell"), "Grid Cell:", choices = NULL)
        ),
        column(4,
          div(style = "margin-top: 25px;",
            textOutput(ns("indiv_status")))
        )
      )
    ),
    DTOutput(ns("indiv_table")),
    hr(),
    h5("Rejection Rate Summary"),
    wellPanel(
      fluidRow(
        column(4,
          selectizeInput(ns("rr_n_filter"), "Sample Size (n):",
                         choices = NULL, multiple = TRUE,
                         options = list(placeholder = "All..."))
        ),
        column(4,
          selectizeInput(ns("rr_phi_filter"), "Phi:",
                         choices = NULL, multiple = TRUE,
                         options = list(placeholder = "All..."))
        ),
        column(4,
          selectizeInput(ns("rr_innov_filter"), "Innovation:",
                         choices = NULL, multiple = TRUE,
                         options = list(placeholder = "All..."))
        )
      ),
      DTOutput(ns("rr_table"))
    )
  )
}

mod_capstone_individual_server <- function(input, output, session,
                                           cap_sim_data, cap_cell_choices,
                                           cap_combined_results) {
  # Populate cell selector
  observe({
    ch <- cap_cell_choices()
    updateSelectInput(session, "indiv_cell", choices = ch)
  })

  # Get sim results for selected cell
  indiv_results <- reactive({
    req(input$indiv_cell)
    cap_sim_data(input$indiv_cell)
  })

  output$indiv_status <- renderText({
    res <- indiv_results()
    if (is.null(res) || length(res) == 0) return("No data for this cell")
    paste0(length(res), " simulations")
  })

  output$indiv_table <- renderDT({
    res <- indiv_results()
    if (is.null(res) || length(res) == 0) {
      return(datatable(data.frame(Message = "No data for this cell"),
                       rownames = FALSE, options = list(dom = "t")))
    }

    df <- data.frame(
      Sim = seq_along(res),
      `Obs. Stat` = vapply(res, function(r) r$tco_obs %||% r$obs_stat %||% NA_real_, numeric(1)),
      `P-Value (Boot)` = vapply(res, function(r) r$pvalue, numeric(1)),
      `P-Value (Asymp)` = vapply(res, function(r) r$pvalue_asymp, numeric(1)),
      `P-Value (COBA)` = vapply(res, function(r) {
        if (is.null(r$pvalue_adj)) NA_real_ else r$pvalue_adj
      }, numeric(1)),
      `AR Order` = vapply(res, function(r) r$p, integer(1)),
      `Innov. Var.` = vapply(res, function(r) r$vara, numeric(1)),
      check.names = FALSE
    )

    datatable(df, rownames = FALSE, selection = "single",
              options = list(pageLength = 50, scrollX = TRUE,
                             dom = "ltip")) |>
      formatRound(columns = c("Obs. Stat", "P-Value (Boot)", "P-Value (Asymp)",
                               "P-Value (COBA)", "Innov. Var."), digits = 4)
  })

  # --- Rejection Rates Summary Table ---
  observe({
    df <- cap_combined_results()
    if (nrow(df) == 0) return()
    updateSelectizeInput(session, "rr_n_filter",
                         choices = as.character(sort(unique(df$n))), server = FALSE)
    updateSelectizeInput(session, "rr_phi_filter",
                         choices = as.character(sort(unique(df$phi))), server = FALSE)
    updateSelectizeInput(session, "rr_innov_filter",
                         choices = sort(unique(df$innov_dist)), server = FALSE)
  })

  rr_filtered <- reactive({
    df <- cap_combined_results()
    if (nrow(df) == 0) return(df)
    if (length(input$rr_n_filter) > 0)
      df <- df[df$n %in% as.integer(input$rr_n_filter), , drop = FALSE]
    if (length(input$rr_phi_filter) > 0)
      df <- df[df$phi %in% as.numeric(input$rr_phi_filter), , drop = FALSE]
    if (length(input$rr_innov_filter) > 0)
      df <- df[df$innov_dist %in% input$rr_innov_filter, , drop = FALSE]
    df[order(df$innov_dist, df$n, df$phi), , drop = FALSE]
  })

  output$rr_table <- renderDT({
    df <- rr_filtered()
    if (nrow(df) == 0) {
      return(datatable(data.frame(Message = "No data"),
                       rownames = FALSE, options = list(dom = "t")))
    }
    display <- data.frame(
      Innovation = df$innov_dist, n = df$n, phi = df$phi,
      Sims = df$n_sims,
      `CO Rate` = df$reject_asymp_05, `CO SE` = df$reject_asymp_05_se,
      `COB Rate` = df$reject_05, `COB SE` = df$reject_05_se,
      `COBA Rate` = df$reject_adj_05, `COBA SE` = df$reject_adj_05_se,
      check.names = FALSE
    )
    dt <- datatable(display, rownames = FALSE,
                    options = list(dom = "tip", pageLength = 50, scrollX = TRUE))
    dt <- formatRound(dt, columns = c("CO Rate", "CO SE", "COB Rate", "COB SE",
                                       "COBA Rate", "COBA SE"), digits = 4)
    for (col in c("CO Rate", "COB Rate", "COBA Rate")) {
      dt <- formatStyle(dt, col,
        backgroundColor = styleInterval(c(0.03, 0.07),
                                         c("#fff3cd", "#d4edda", "#f8d7da")))
    }
    dt
  })

  # Return selected sim index for Bootstrap tab
  reactive(input$indiv_table_rows_selected)
}
