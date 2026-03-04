# Capstone Sub-Tab: Innovation Diagnostics (Tab 8)
# Regenerates innovation sample for a selected cell and shows diagnostics

mod_capstone_innov_diag_ui <- function(ns) {
  tabPanel("Innovation Diagnostics",
    br(),
    h4("Innovation Diagnostics"),
    p(class = "text-body-secondary",
      "Regenerates a diagnostic innovation sample for the selected cell. ",
      "For IID distributions: histogram with theoretical overlay + QQ plot. ",
      "For time-dependent (GARCH/Hetero): time series + histogram + squared innovations."),
    div(class = "plot-controls",
      fluidRow(
        column(4,
          selectInput(ns("idiag_cell"), "Grid Cell:", choices = NULL)
        ),
        column(3,
          numericInput(ns("idiag_n"), "Sample Size", value = 1000,
                       min = 100, max = 10000, step = 100)
        ),
        column(3,
          numericInput(ns("idiag_seed"), "Seed", value = 42)
        ),
        column(2,
          div(style = "margin-top: 25px;",
            actionButton(ns("idiag_gen"), "Generate",
                         icon = icon("refresh"), class = "btn-sm btn-primary"))
        )
      )
    ),
    plotOutput(ns("idiag_plot"), height = "500px")
  )
}

mod_capstone_innov_diag_server <- function(input, output, session,
                                            cap_cell_choices, cap_cell_meta) {

  observe({
    ch <- cap_cell_choices()
    updateSelectInput(session, "idiag_cell", choices = ch)
  })

  innov_sample <- reactiveVal(NULL)

  observeEvent(input$idiag_gen, {
    req(input$idiag_cell)
    meta <- cap_cell_meta(input$idiag_cell)
    if (is.null(meta)) {
      showNotification("No metadata for this cell", type = "warning")
      return()
    }

    n_sample <- as.integer(input$idiag_n %||% 1000)
    seed_val <- as.integer(input$idiag_seed %||% 42)

    # Build generator
    gen <- tryCatch(
      build_innov_gen(meta$innov_label, meta$innov_params),
      error = function(e) NULL)
    if (is.null(gen)) {
      showNotification("Could not build innovation generator", type = "warning")
      return()
    }

    set.seed(seed_val)
    sample <- tryCatch(gen(n_sample), error = function(e) NULL)
    if (is.null(sample)) {
      showNotification("Generator returned NULL", type = "warning")
      return()
    }

    is_td <- meta$innov_label %in% c("GARCH", "Heteroscedastic")
    innov_sample(list(sample = sample, is_td = is_td,
                      label = meta$innov_label,
                      dist = meta$innov_label,
                      params = meta$innov_params))
  })

  output$idiag_plot <- renderPlot(bg = "transparent", {
    data <- innov_sample()
    if (is.null(data)) {
      plot.new()
      text(0.5, 0.5, "Select a cell and click Generate", cex = 1.2,
           col = viewer_plot_fg())
      return()
    }

    fg <- viewer_plot_fg()
    # Map params to innov_theoretical_density format
    theo_params <- list()
    p <- data$params
    if (data$dist == "Normal") {
      theo_params <- list(sd = p$norm_sd %||% 1)
    } else if (data$dist == "Student's t") {
      theo_params <- list(df = p$t_df %||% 3, scale = isTRUE(p$t_scale))
    } else if (data$dist == "Skew-t") {
      theo_params <- list(df = p$skt_df %||% 5, alpha = p$skt_alpha %||% 0,
                          scale = isTRUE(p$skt_scale))
    } else if (data$dist == "GED") {
      theo_params <- list(nu = p$ged_nu %||% 2, sd = p$ged_sd %||% 1)
    } else if (data$dist == "Laplace") {
      theo_params <- list(scale = p$lap_scale %||% (1 / sqrt(2)))
    } else if (data$dist == "Uniform") {
      theo_params <- list(half_width = p$unif_hw %||% sqrt(3))
    } else if (data$dist == "Mixture Normal") {
      theo_params <- list(sd1 = p$mix_sd1 %||% 1, sd2 = p$mix_sd2 %||% 3,
                          prob1 = p$mix_prob1 %||% 0.9)
    }

    plot_innovation_diagnostics(
      innov = data$sample,
      is_time_dependent = data$is_td,
      dist_label = data$label,
      innov_dist = data$dist,
      innov_params = theo_params,
      fg = fg
    )
  })
}
