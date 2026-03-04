# Capstone Sub-Tab: Bootstrap Distribution (Tab 7)
# Per-simulation bootstrap t-statistic histogram

mod_capstone_bootdist_ui <- function(ns) {
  tabPanel("Bootstrap Distribution",
    br(),
    h4("Bootstrap Distribution"),
    p(class = "text-body-secondary",
      "View the bootstrap t-statistic distribution for a specific simulation. ",
      "The red line shows the observed test statistic."),
    div(class = "plot-controls",
      fluidRow(
        column(4,
          selectInput(ns("bd_cell"), "Grid Cell:", choices = NULL)
        ),
        column(4,
          selectizeInput(ns("bd_sim"), "Simulation:",
                         choices = NULL, options = list(maxOptions = 5000))
        ),
        column(4,
          div(style = "margin-top: 25px;",
            textOutput(ns("bd_status")))
        )
      )
    ),
    plotOutput(ns("bd_histogram"), height = "400px"),
    hr(),
    verbatimTextOutput(ns("bd_summary"))
  )
}

mod_capstone_bootdist_server <- function(input, output, session,
                                          cap_sim_data, cap_cell_choices,
                                          selected_sim_from_individual = NULL) {

  observe({
    ch <- cap_cell_choices()
    updateSelectInput(session, "bd_cell", choices = ch)
  })

  # Populate simulation dropdown when cell changes
  observe({
    req(input$bd_cell)
    sims <- cap_sim_data(input$bd_cell)
    if (is.null(sims) || length(sims) == 0) {
      updateSelectizeInput(session, "bd_sim", choices = character(0))
      return()
    }
    n_sims <- length(sims)
    labels <- vapply(seq_len(n_sims), function(i) {
      t_val <- sims[[i]]$tco_obs %||% sims[[i]]$obs_stat %||% NA_real_
      sprintf("#%d (t=%.3f)", i, t_val)
    }, character(1))
    choices <- setNames(as.character(seq_len(n_sims)), labels)
    updateSelectizeInput(session, "bd_sim", choices = choices, server = FALSE)
  })

  # If user clicked a row in Individual Results, jump to that sim
  observe({
    if (!is.null(selected_sim_from_individual)) {
      sel <- selected_sim_from_individual()
      if (!is.null(sel) && length(sel) > 0) {
        updateSelectizeInput(session, "bd_sim", selected = as.character(sel[1]))
      }
    }
  })

  output$bd_status <- renderText({
    sims <- cap_sim_data(input$bd_cell)
    if (is.null(sims)) return("No data")
    paste0(length(sims), " simulations available")
  })

  output$bd_histogram <- renderPlot(bg = "transparent", {
    req(input$bd_cell, input$bd_sim)
    sims <- cap_sim_data(input$bd_cell)
    idx <- as.integer(input$bd_sim)
    if (is.null(sims) || idx > length(sims) || idx < 1) {
      plot.new(); text(0.5, 0.5, "No data", cex = 1.2, col = viewer_plot_fg()); return()
    }
    sim <- sims[[idx]]
    boot_dist <- sim$boot_tstats
    obs_stat <- sim$tco_obs %||% sim$obs_stat
    pval <- sim$pvalue

    if (is.null(boot_dist) || length(boot_dist) == 0) {
      plot.new(); text(0.5, 0.5, "No bootstrap distribution stored", cex = 1.2,
                       col = viewer_plot_fg()); return()
    }

    parts <- strsplit(input$bd_cell, "\\|")[[1]]
    title <- sprintf("Bootstrap Distribution (sim #%d, n=%s, phi=%s, %s)",
                     idx, parts[1], parts[2], parts[3])
    p <- plot_bootstrap_distribution(boot_dist, obs_stat, pval, title = title)
    print(p)
  })

  output$bd_summary <- renderText({
    req(input$bd_cell, input$bd_sim)
    sims <- cap_sim_data(input$bd_cell)
    idx <- as.integer(input$bd_sim)
    if (is.null(sims) || idx > length(sims)) return("")
    sim <- sims[[idx]]
    boot_dist <- sim$boot_tstats
    obs_stat <- sim$tco_obs %||% sim$obs_stat

    lines <- c(
      sprintf("Simulation: #%d", idx),
      sprintf("Observed t-statistic: %.4f", obs_stat),
      sprintf("P-value (bootstrap): %.4f", sim$pvalue),
      sprintf("P-value (asymptotic): %.4f", sim$pvalue_asymp)
    )
    if (!is.null(sim$pvalue_adj))
      lines <- c(lines, sprintf("P-value (COBA): %.4f", sim$pvalue_adj))
    if (!is.null(boot_dist) && length(boot_dist) > 0) {
      lines <- c(lines, "",
                  sprintf("Bootstrap distribution: n=%d, mean=%.4f, sd=%.4f",
                          length(boot_dist), mean(boot_dist), sd(boot_dist)))
    }
    paste(lines, collapse = "\n")
  })
}
