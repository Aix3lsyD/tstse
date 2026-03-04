# Capstone Sub-Tab: Coverage / Overview (Tab 2)
# Shows in-session coverage KPIs + heatmap + coverage table

mod_capstone_overview_ui <- function(ns) {
  tabPanel("Coverage / Overview",
    br(),
    uiOutput(ns("overview_kpis")),
    p(class = "text-body-secondary",
      "Shows how many simulations have been generated in this session for each grid cell."),
    plotOutput(ns("coverage_heatmap"), height = "450px"),
    br(),
    DTOutput(ns("coverage_table"))
  )
}

mod_capstone_overview_server <- function(input, output, session,
                                         grid_rv, raw_results_rv) {

  # Coverage data from in-memory results
  coverage_data <- reactive({
    grid <- grid_rv()
    if (nrow(grid) == 0) return(data.frame())

    n_existing <- integer(nrow(grid))
    raw <- raw_results_rv()
    for (i in seq_len(nrow(grid))) {
      cell_key <- paste(grid$n[i], grid$phi[i], grid$innov_dist_str[i], sep = "|")
      if (cell_key %in% names(raw)) {
        n_existing[i] <- length(raw[[cell_key]]$results)
      }
    }

    data.frame(
      phi            = grid$phi,
      n              = grid$n,
      innov_label    = grid$innov_label,
      innov_dist_str = grid$innov_dist_str,
      n_existing     = n_existing,
      stringsAsFactors = FALSE
    )
  })

  # KPI cards (DB summary)
  output$overview_kpis <- renderUI({
    raw <- raw_results_rv()
    n_cells <- length(raw)
    n_sims <- sum(vapply(raw, function(x) length(x$results), integer(1)))

    fluidRow(
      column(3, wellPanel(h5("Total Sims"), h3(format(n_sims, big.mark = ",")))),
      column(3, wellPanel(h5("Configs"), h3(n_cells))),
      column(3, wellPanel(h5("Source"), h3("Session"))),
      column(3, wellPanel(h5("Persistence"), h3("Disabled")))
    )
  })

  # Coverage heatmap
  output$coverage_heatmap <- renderPlot(bg = "transparent", {
    df <- coverage_data()
    if (nrow(df) == 0) {
      plot.new()
      text(0.5, 0.5, "No grid cells defined.", cex = 1.2,
           col = viewer_plot_fg())
      return()
    }
    df$phi_f <- factor(df$phi)
    df$n_f <- factor(df$n)
    df$label <- as.character(df$n_existing)

    ggplot2::ggplot(df, ggplot2::aes(x = phi_f, y = n_f, fill = n_existing)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.5) +
      ggplot2::geom_text(ggplot2::aes(label = label), size = 4.5, fontface = "bold") +
      ggplot2::scale_fill_gradient(low = "#f8f9fa", high = "#08519c",
                                   na.value = "#f8f9fa") +
      ggplot2::facet_wrap(~ innov_label, scales = "free") +
      viewer_plot_theme() +
      ggplot2::labs(x = "Phi", y = "n", fill = "Existing\nSims",
                    title = "Coverage by Grid Cell")
  })

  # Coverage table
  output$coverage_table <- renderDT({
    df <- coverage_data()
    if (nrow(df) == 0) {
      return(datatable(data.frame(Message = "No grid cells defined."),
                       rownames = FALSE, options = list(dom = "t")))
    }
    display <- data.frame(
      Phi = df$phi, n = df$n, Innovation = df$innov_label,
      `DB Key` = df$innov_dist_str, `Existing Sims` = df$n_existing,
      check.names = FALSE)
    datatable(display, rownames = FALSE,
              options = list(dom = "tip", pageLength = 25, ordering = TRUE))
  }, server = FALSE)
}
