# Capstone Sub-Tab: Null Model Diagnostics (Tab 9)
# AR order distribution, variance distribution, rejection by order,
# MC convergence, batch consistency, test stat distribution

mod_capstone_null_diag_ui <- function(ns) {
  tabPanel("Null Model Diagnostics",
    br(),
    h4("Null Model & MC Diagnostics"),
    p(class = "text-body-secondary",
      "Diagnostics for the fitted null models and Monte Carlo machinery. ",
      "Top section uses per-simulation data; bottom section uses aggregate DB data."),

    div(class = "plot-controls",
      fluidRow(
        column(4,
          selectInput(ns("ndiag_cell"), "Grid Cell:",
                      choices = c("All Cells" = "all"))
        ),
        column(4,
          div(style = "margin-top: 25px;",
            textOutput(ns("ndiag_status")))
        )
      )
    ),

    # Per-simulation diagnostics (from raw results)
    wellPanel(
      h5("Null Model Diagnostics (per-simulation)"),
      plotOutput(ns("ndiag_null_model"), height = "500px")
    ),

    # ggplot-based diagnostics
    wellPanel(
      h5("AR Order Distribution"),
      plotOutput(ns("ndiag_ar_order"), height = "300px")
    ),
    wellPanel(
      h5("Estimated AR(1) Coefficient Distribution"),
      plotOutput(ns("ndiag_phi_dist"), height = "300px")
    ),
    wellPanel(
      h5("Test Statistic Distribution"),
      plotOutput(ns("ndiag_tstat"), height = "300px")
    ),
    wellPanel(
      h5("Monte Carlo Convergence"),
      fluidRow(
        column(4,
          selectInput(ns("ndiag_conv_method"), "Method:",
                      choices = c("Bootstrap (COB)" = "pvalue",
                                  "Asymptotic (CO)" = "pvalue_asymp",
                                  "COBA" = "pvalue_adj"))
        )
      ),
      plotOutput(ns("ndiag_convergence"), height = "350px")
    ),
    wellPanel(
      h5("Batch-to-Batch Consistency"),
      p(class = "text-body-secondary", style = "margin-bottom: 10px;",
        "Requires DB with multiple batches."),
      fluidRow(
        column(4,
          selectInput(ns("ndiag_batch_method"), "Method:",
                      choices = c("Bootstrap (COB)" = "reject_05",
                                  "Asymptotic (CO)" = "reject_asymp_05",
                                  "COBA" = "reject_adj_05"))
        )
      ),
      plotOutput(ns("ndiag_batch_plot"), height = "350px")
    )
  )
}

mod_capstone_null_diag_server <- function(input, output, session,
                                          cap_sim_data, cap_cell_choices,
                                          con, db_refresh_trigger) {

  observe({
    ch <- cap_cell_choices()
    all_choices <- c("All Cells" = "all", ch)
    updateSelectInput(session, "ndiag_cell", choices = all_choices)
  })

  # Collect raw results across selected cell(s)
  ndiag_sims <- reactive({
    req(input$ndiag_cell)
    ch <- cap_cell_choices()
    if (length(ch) == 0) return(NULL)
    keys <- if (input$ndiag_cell == "all") unname(ch) else input$ndiag_cell
    all_sims <- list()
    for (key in keys) {
      sims <- cap_sim_data(key)
      if (!is.null(sims) && length(sims) > 0) {
        all_sims <- c(all_sims, sims)
      }
    }
    if (length(all_sims) == 0) return(NULL)
    all_sims
  })

  output$ndiag_status <- renderText({
    sims <- ndiag_sims()
    if (is.null(sims)) return("No data")
    paste0(length(sims), " simulations")
  })

  # Null model diagnostics (base graphics, 3 panels)
  output$ndiag_null_model <- renderPlot(bg = "transparent", {
    sims <- ndiag_sims()
    if (is.null(sims) || length(sims) == 0) {
      plot.new(); text(0.5, 0.5, "No data", cex = 1.2, col = viewer_plot_fg()); return()
    }
    fg <- viewer_plot_fg()
    # Need maxp - guess from max observed order
    maxp <- max(vapply(sims, function(r) r$p, integer(1)), na.rm = TRUE)
    plot_null_model_diagnostics(sims, nsims = length(sims), maxp = maxp,
                                 min_p = 1L, fg = fg)
  })

  # Build diag data frame for ggplot-based plots
  ndiag_df <- reactive({
    sims <- ndiag_sims()
    if (is.null(sims)) return(data.frame())

    # Parse cell key(s) for config info
    ch <- cap_cell_choices()
    keys <- if (input$ndiag_cell == "all") unname(ch) else input$ndiag_cell

    rows <- list()
    for (key in keys) {
      cell_sims <- cap_sim_data(key)
      if (is.null(cell_sims)) next
      parts <- strsplit(key, "\\|")[[1]]
      for (sim in cell_sims) {
        phi_hat <- if (is.null(sim$phi) || length(sim$phi) == 0) NA_real_ else sim$phi[1]
        rows <- c(rows, list(data.frame(
          null_ar_order = sim$p,
          null_phi1 = phi_hat,
          obs_stat = sim$tco_obs %||% sim$obs_stat %||% NA_real_,
          pvalue = sim$pvalue,
          pvalue_asymp = sim$pvalue_asymp,
          pvalue_adj = if (is.null(sim$pvalue_adj)) NA_real_ else sim$pvalue_adj,
          n = as.integer(parts[1]),
          phi = as.numeric(parts[2]),
          innov_dist = parts[3],
          stringsAsFactors = FALSE
        )))
      }
    }
    if (length(rows) == 0) return(data.frame())
    do.call(rbind, rows)
  })

  # AR order distribution (ggplot)
  output$ndiag_ar_order <- renderPlot(bg = "transparent", {
    df <- ndiag_df()
    if (nrow(df) == 0) {
      plot.new(); text(0.5, 0.5, "No data", cex = 1.2, col = viewer_plot_fg()); return()
    }
    p <- plot_ar_order_distribution(df)
    if (!is.null(p)) print(p)
  })

  # AR(1) coefficient distribution
  output$ndiag_phi_dist <- renderPlot(bg = "transparent", {
    df <- ndiag_df()
    if (nrow(df) == 0) {
      plot.new(); text(0.5, 0.5, "No data", cex = 1.2, col = viewer_plot_fg()); return()
    }
    p <- plot_ar_coefficient_distribution(df)
    if (!is.null(p)) print(p) else {
      plot.new(); text(0.5, 0.5, "No AR coefficient data", cex = 1.2, col = viewer_plot_fg())
    }
  })

  # Test statistic distribution
  output$ndiag_tstat <- renderPlot(bg = "transparent", {
    df <- ndiag_df()
    if (nrow(df) == 0 || all(is.na(df$obs_stat))) {
      plot.new(); text(0.5, 0.5, "No data", cex = 1.2, col = viewer_plot_fg()); return()
    }
    p <- plot_test_statistic_distribution(df[!is.na(df$obs_stat), ])
    if (!is.null(p)) print(p)
  })

  # MC convergence
  output$ndiag_convergence <- renderPlot(bg = "transparent", {
    df <- ndiag_df()
    pval_col <- input$ndiag_conv_method %||% "pvalue"
    if (nrow(df) == 0 || !pval_col %in% names(df)) {
      plot.new(); text(0.5, 0.5, "No data", cex = 1.2, col = viewer_plot_fg()); return()
    }
    method_labels <- c(pvalue = "Bootstrap (COB)",
                       pvalue_asymp = "Asymptotic (CO)",
                       pvalue_adj = "COBA")
    p <- plot_mc_convergence(df, pval_col = pval_col,
                              method_label = method_labels[pval_col])
    if (!is.null(p)) print(p)
  })

  # Batch consistency (requires DB)
  output$ndiag_batch_plot <- renderPlot(bg = "transparent", {
    if (is.null(con)) {
      plot.new(); text(0.5, 0.5, "Requires database connection", cex = 1.2,
                       col = viewer_plot_fg()); return()
    }
    db_refresh_trigger()
    rate_col <- input$ndiag_batch_method %||% "reject_05"
    df <- tryCatch(
      DBI::dbGetQuery(con, "SELECT * FROM v_rejection_rates_by_batch ORDER BY innov_dist, n, phi"),
      error = function(e) data.frame())
    if (nrow(df) == 0 || !rate_col %in% names(df)) {
      plot.new(); text(0.5, 0.5, "No per-batch data available", cex = 1.2,
                       col = viewer_plot_fg()); return()
    }
    p <- plot_batch_consistency(df, rate_col = rate_col)
    if (!is.null(p)) print(p)
  })
}
