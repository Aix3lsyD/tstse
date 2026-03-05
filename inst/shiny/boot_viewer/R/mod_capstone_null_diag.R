# Capstone Sub-Tab: Null Model Diagnostics (Tab 9)
# AR order distribution, variance distribution, rejection by order,
# MC convergence, test stat distribution

mod_capstone_null_diag_ui <- function(ns) {
  tabPanel("Null Model Diagnostics",
    br(),
    h4("Null Model & MC Diagnostics"),
    p(class = "text-body-secondary",
      "Diagnostics for fitted null models and Monte Carlo behavior in the selected scenario."),

    div(class = "plot-controls",
      fluidRow(
        column(4,
          selectInput(ns("ndiag_cell"), "Grid Cell:",
                      choices = character(0))
        ),
        column(4,
          selectInput(ns("ndiag_reject_mode"), "Rejection rates shown:",
                      choices = c("All available" = "all",
                                  "COB" = "cob",
                                  "CO" = "co",
                                  "COBA" = "coba"),
                      selected = "all")
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
      plotOutput(ns("ndiag_null_model"), height = "900px")
    ),

    fluidRow(
      column(12,
        wellPanel(
          h5("Estimated AR(1) Coefficient Distribution"),
          plotOutput(ns("ndiag_phi_dist"), height = "420px")
        )
      )
    ),

    fluidRow(
      column(6,
        wellPanel(
          h5("Test Statistic Distribution"),
          plotOutput(ns("ndiag_tstat"), height = "420px")
        )
      ),
      column(6,
        wellPanel(
          h5("Monte Carlo Convergence"),
          fluidRow(
            column(12,
              selectInput(ns("ndiag_conv_method"), "Method:",
                          choices = c("Bootstrap (COB)" = "pvalue",
                                      "Asymptotic (CO)" = "pvalue_asymp",
                                      "COBA" = "pvalue_adj"))
            )
          ),
          plotOutput(ns("ndiag_convergence"), height = "420px")
        )
      )
    )
  )
}

mod_capstone_null_diag_server <- function(input, output, session,
                                          cap_sim_data, cap_cell_choices) {

  observe({
    ch <- cap_cell_choices()
    current <- isolate(input$ndiag_cell)
    selected <- if (!is.null(current) && current %in% unname(ch)) current else {
      if (length(ch) > 0) unname(ch)[1] else character(0)
    }
    updateSelectInput(session, "ndiag_cell", choices = ch,
                      selected = selected)
  })

  # Collect raw results for selected scenario
  ndiag_sims <- reactive({
    req(input$ndiag_cell)
    ch <- cap_cell_choices()
    if (length(ch) == 0) return(NULL)
    keys <- input$ndiag_cell
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
    paste0(input$ndiag_cell, " (", length(sims), " simulations)")
  })

  # Null model diagnostics (base graphics, 3 panels)
  output$ndiag_null_model <- renderPlot(bg = "transparent", {
    sims <- ndiag_sims()
    if (is.null(sims) || length(sims) == 0) {
      op <- viewer_base_par(scale = 1.2)
      on.exit(par(op), add = TRUE)
      plot.new(); text(0.5, 0.5, "No data", cex = 1.2, col = viewer_plot_fg(session = session)); return()
    }
    op <- viewer_base_par(scale = 1.2)
    on.exit(par(op), add = TRUE)
    fg <- viewer_plot_fg(session = session)
    # Need maxp - guess from max observed order
    maxp <- max(vapply(sims, function(r) r$p, integer(1)), na.rm = TRUE)
    reject_mode <- input$ndiag_reject_mode %||% "all"
    plot_null_model_diagnostics(sims, nsims = length(sims), maxp = maxp,
                                  min_p = 1L, fg = fg, reject_mode = reject_mode)
  })

  # Build diag data frame for ggplot-based plots
  ndiag_df <- reactive({
    sims <- ndiag_sims()
    if (is.null(sims)) return(data.frame())

    # Parse cell key(s) for config info
    ch <- cap_cell_choices()
    keys <- input$ndiag_cell

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

  # AR(1) coefficient distribution
  output$ndiag_phi_dist <- renderPlot(bg = "transparent", {
    df <- ndiag_df()
    if (nrow(df) == 0) {
      plot.new(); text(0.5, 0.5, "No data", cex = 1.2, col = viewer_plot_fg()); return()
    }
    p <- plot_ar_coefficient_distribution(df)
    if (!is.null(p)) print(p) else {
      plot.new(); text(0.5, 0.5, "No AR(1)-selected fits for this scenario", cex = 1.2, col = viewer_plot_fg())
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

}
