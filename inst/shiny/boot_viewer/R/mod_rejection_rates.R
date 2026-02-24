# Module: Rejection Rates (Tab 2)
# Extracted from the monolithic app.R

mod_rejection_rates_ui <- function(id) {
  ns <- NS(id)

  tabPanel(
    "Rejection Rates",
    br(),
    wellPanel(class = "plot-controls",
      fluidRow(
        column(2,
          selectizeInput(ns("rr_n_filter"), "Sample Size (n)",
                         choices = NULL, multiple = TRUE,
                         options = list(placeholder = "All..."))
        ),
        column(2,
          selectizeInput(ns("rr_phi_filter"), "Phi Value",
                         choices = NULL, multiple = TRUE,
                         options = list(placeholder = "All..."))
        ),
        column(3,
          selectizeInput(ns("rr_innov_filter"), "Innovation Distribution",
                         choices = NULL, multiple = TRUE,
                         options = list(placeholder = "All..."))
        ),
        column(2,
          div(style = "margin-top: 25px;",
            checkboxInput(ns("rr_compare_batches"), "Compare Batches", value = FALSE)
          )
        ),
        column(1,
          div(style = "margin-top: 25px;",
            actionButton(ns("rr_clear"), "Clear",
                         icon = icon("times-circle"),
                         class = "btn-sm btn-outline-secondary")
          )
        ),
        column(2,
          div(style = "margin-top: 25px;",
            downloadButton(ns("export_csv"), "Export CSV",
                           class = "btn-sm btn-outline-secondary")
          )
        )
      )
    ),
    h4("Rejection Rate Summary"),
    DTOutput(ns("rejection_table"))
  )
}

mod_rejection_rates_server <- function(id, con, init_choices) {
  moduleServer(id, function(input, output, session) {

    # Populate filter dropdowns from initial choices
    updateSelectizeInput(session, "rr_n_filter",
                         choices = init_choices$n, server = FALSE)
    updateSelectizeInput(session, "rr_phi_filter",
                         choices = init_choices$phi, server = FALSE)
    updateSelectizeInput(session, "rr_innov_filter",
                         choices = init_choices$innov, server = FALSE)

    # Rejection rates data reactive
    rr_data <- reactive({
      view_name <- if (isTRUE(input$rr_compare_batches)) {
        "v_rejection_rates_by_batch"
      } else {
        "v_rejection_rates"
      }
      q <- build_query(view_name,
                       n_vals = input$rr_n_filter,
                       phi_vals = input$rr_phi_filter,
                       innov_vals = input$rr_innov_filter)
      tryCatch(
        dbGetQuery(con, q$sql, params = q$params),
        error = function(e) data.frame()
      )
    })

    # Clear filters
    observeEvent(input$rr_clear, {
      updateSelectizeInput(session, "rr_n_filter", selected = character(0))
      updateSelectizeInput(session, "rr_phi_filter", selected = character(0))
      updateSelectizeInput(session, "rr_innov_filter", selected = character(0))
      updateCheckboxInput(session, "rr_compare_batches", value = FALSE)
    })

    # Rejection rates table
    output$rejection_table <- renderDT({
      df <- rr_data()
      if (nrow(df) == 0) return(datatable(data.frame(Message = "No data")))

      # Select display columns based on mode
      if (isTRUE(input$rr_compare_batches)) {
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
        write.csv(rr_data(), file, row.names = FALSE)
      }
    )
  })
}
