# Capstone Sub-Tab: Plots (Tab 5)
# Power curves, heatmaps, deviation, rate vs n, forest plot, parallel coords

mod_capstone_plots_ui <- function(ns) {
  tabPanel("Plots",
    br(),
    h4("Visualization"),
    div(class = "plot-controls",
      fluidRow(
        column(3,
          selectInput(ns("plot_type"), "Plot Type:",
                      choices = c("Power Curve" = "power",
                                  "Heatmap" = "heatmap",
                                  "Deviation from Nominal" = "deviation",
                                  "Rate vs Sample Size" = "rate_n",
                                  "Method Forest Plot" = "forest",
                                  "Parallel Coordinates: Rate Profile" = "pc_rate",
                                  "Parallel Coordinates: DGP Scatter" = "pc_dgp",
                                  "Parallel Coordinates: Method Comparison" = "pc_method"))
        ),
        column(3,
          conditionalPanel(
            condition = "input.plot_type == 'power' || input.plot_type == 'heatmap' || input.plot_type == 'rate_n'",
            ns = ns,
            selectInput(ns("plot_rate_col"), "Rate Type:",
                        choices = c("Bootstrap (COB)" = "reject_05",
                                    "Asymptotic (CO)" = "reject_asymp_05",
                                    "COBA" = "reject_adj_05"))
          )
        ),
        column(3,
          conditionalPanel(
            condition = "input.plot_type == 'power'", ns = ns,
            selectInput(ns("plot_facet"), "Facet By:",
                        choices = c("None" = "none", "n" = "n",
                                    "Innovation" = "innov_dist"))
          ),
          conditionalPanel(
            condition = "input.plot_type == 'rate_n'", ns = ns,
            selectInput(ns("plot_rate_n_facet"), "Facet By:",
                        choices = c("None" = "none",
                                    "Innovation" = "innov_dist",
                                    "Phi" = "phi"))
          ),
          conditionalPanel(
            condition = "input.plot_type == 'pc_method'", ns = ns,
            selectInput(ns("pc_color_by"), "Color By:",
                        choices = c("Innovation" = "innov_dist",
                                    "Sample Size" = "n",
                                    "Phi" = "phi"))
          )
        ),
        column(3,
          conditionalPanel(
            condition = "input.plot_type == 'pc_method'", ns = ns,
            selectInput(ns("pc_method"), "Method:",
                        choices = c("Bootstrap (COB)" = "reject_05",
                                    "Asymptotic (CO)" = "reject_asymp_05",
                                    "COBA" = "reject_adj_05"))
          )
        )
      )
    ),
    plotOutput(ns("cap_plot"), height = "500px"),
    br(),
    downloadButton(ns("download_plot_png"), "Download PNG",
                   class = "btn-sm btn-outline-secondary")
  )
}

mod_capstone_plots_server <- function(input, output, session,
                                      cap_combined_results) {

  cap_plot_obj <- reactive({
    df <- cap_combined_results()
    if (nrow(df) == 0) return(NULL)

    switch(input$plot_type,
      "power" = plot_power_curve(df, rate_col = input$plot_rate_col,
                                  facet_by = input$plot_facet %||% "none"),
      "heatmap" = plot_heatmap_rejection(df, rate_col = input$plot_rate_col),
      "deviation" = plot_deviation_from_nominal(df),
      "rate_n" = plot_rate_vs_n(df, rate_col = input$plot_rate_col,
                                 facet_by = input$plot_rate_n_facet %||% "none"),
      "forest" = plot_forest_methods(df),
      "pc_rate" = build_parcoord(
        df,
        cols = c(reject_asymp_05 = "CO Rate", reject_05 = "COB Rate",
                 reject_adj_05 = "COBA Rate"),
        color_col = "innov_dist", color_label = "Innovation",
        title = "Rate Profile: CO \u2192 COB \u2192 COBA"),
      "pc_dgp" = plot_pc_dgp_scatter(df),
      "pc_method" = {
        color_col <- input$pc_color_by %||% "innov_dist"
        df[[color_col]] <- as.character(df[[color_col]])
        color_labels <- c(innov_dist = "Innovation", n = "Sample Size", phi = "Phi")

        rate_col <- input$pc_method %||% "reject_05"
        rate_labels_pc <- c(reject_05 = "COB Rate",
                            reject_asymp_05 = "CO Rate",
                            reject_adj_05 = "COBA Rate")
        named_cols <- c(n = "n", phi = "phi")
        named_cols[[rate_col]] <- rate_labels_pc[rate_col]

        build_parcoord(df, cols = named_cols,
                       color_col = color_col,
                       color_label = color_labels[color_col],
                       title = paste0("Sensitivity: n \u2192 phi \u2192 ",
                                      rate_labels_pc[rate_col]))
      },
      NULL
    )
  })

  output$cap_plot <- renderPlot(bg = "transparent", {
    p <- cap_plot_obj()
    if (is.null(p)) {
      plot.new()
      text(0.5, 0.5, "No data available", cex = 1.2, col = viewer_plot_fg())
      return()
    }
    print(p)
  })

  output$download_plot_png <- downloadHandler(
    filename = function() {
      paste0("capstone_plot_", input$plot_type, "_",
             format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      p <- cap_plot_obj()
      if (!is.null(p)) {
        ggplot2::ggsave(file, plot = p, width = 10, height = 7, dpi = 150,
                        bg = "white")
      }
    }
  )
}
