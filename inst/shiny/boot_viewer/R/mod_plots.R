# Module: Plots tab (Tab 3)
# Extracted from monolithic app.R

mod_plots_ui <- function(id) {
  ns <- NS(id)

  tabPanel(
    "Plots",
    br(),
    wellPanel(class = "plot-controls",
      fluidRow(
        column(2,
          selectizeInput(ns("plot_n_filter"), "Sample Size (n)",
                         choices = NULL, multiple = TRUE,
                         options = list(placeholder = "All..."))
        ),
        column(2,
          selectizeInput(ns("plot_phi_filter"), "Phi Value",
                         choices = NULL, multiple = TRUE,
                         options = list(placeholder = "All..."))
        ),
        column(2,
          selectizeInput(ns("plot_innov_filter"), "Innovation Distribution",
                         choices = NULL, multiple = TRUE,
                         options = list(placeholder = "All..."))
        ),
        column(1,
          div(style = "margin-top: 25px;",
            actionButton(ns("plot_clear"), "Clear",
                         icon = icon("times-circle"),
                         class = "btn-sm btn-outline-secondary")
          )
        ),
        column(2,
          selectInput(ns("plot_type"), "Plot Type",
                      choices = c("Power Curve" = "power_curve",
                                  "Heatmap" = "heatmap",
                                  "Deviation from Nominal" = "deviation",
                                  "Rate vs Sample Size" = "rate_vs_n",
                                  "Method Forest Plot" = "forest"))
        ),
        column(2,
          conditionalPanel(
            condition = "input.plot_type == 'power_curve' || input.plot_type == 'heatmap' || input.plot_type == 'rate_vs_n'",
            ns = ns,
            selectInput(ns("rate_type"), "Rate Type",
                        choices = c("Bootstrap" = "reject_05",
                                    "Asymptotic" = "reject_asymp_05",
                                    "COBA" = "reject_adj_05"))
          )
        ),
        column(1,
          conditionalPanel(
            condition = "input.plot_type == 'power_curve'",
            ns = ns,
            selectInput(ns("facet_by"), "Facet By",
                        choices = c("n" = "n",
                                    "Innov" = "innov_dist",
                                    "None" = "none"))
          ),
          conditionalPanel(
            condition = "input.plot_type == 'rate_vs_n'",
            ns = ns,
            selectInput(ns("rate_vs_n_facet"), "Facet By",
                        choices = c("Innov" = "innov_dist",
                                    "Phi" = "phi",
                                    "None" = "none"))
          )
        )
      )
    ),
    plotOutput(ns("main_plot"), height = "600px"),
    br(),
    downloadButton(ns("export_plot"), "Download Plot (PNG)",
                   class = "btn-sm btn-outline-secondary")
  )
}

mod_plots_server <- function(id, con, init_choices) {
  moduleServer(id, function(input, output, session) {

    # Populate filter dropdowns
    updateSelectizeInput(session, "plot_n_filter", choices = init_choices$n, server = FALSE)
    updateSelectizeInput(session, "plot_phi_filter", choices = init_choices$phi, server = FALSE)
    updateSelectizeInput(session, "plot_innov_filter", choices = init_choices$innov, server = FALSE)

    # Clear filters
    observeEvent(input$plot_clear, {
      updateSelectizeInput(session, "plot_n_filter", selected = character(0))
      updateSelectizeInput(session, "plot_phi_filter", selected = character(0))
      updateSelectizeInput(session, "plot_innov_filter", selected = character(0))
    })

    # Query data for plots
    plot_data <- reactive({
      q <- build_query("v_rejection_rates",
                       n_vals = input$plot_n_filter,
                       phi_vals = input$plot_phi_filter,
                       innov_vals = input$plot_innov_filter)
      tryCatch(
        dbGetQuery(con, q$sql, params = q$params),
        error = function(e) data.frame()
      )
    })

    # Build the current plot object
    current_plot <- reactive({
      df <- plot_data()
      if (nrow(df) == 0) return(NULL)

      rate_col <- input$rate_type
      if (!rate_col %in% names(df)) return(NULL)

      rate_labels <- c(
        reject_05 = "Bootstrap Rejection Rate",
        reject_asymp_05 = "Asymptotic Rejection Rate",
        reject_adj_05 = "COBA Rejection Rate"
      )

      if (input$plot_type == "power_curve") {
        # Power curve: rejection rate vs phi
        df$phi <- as.numeric(df$phi)

        p <- ggplot(df, aes(x = phi, y = .data[[rate_col]],
                            color = factor(innov_dist),
                            group = factor(innov_dist))) +
          geom_line(linewidth = 0.8) +
          geom_point(size = 2) +
          geom_hline(yintercept = 0.05, linetype = "dashed") +
          labs(
            x = expression(phi),
            y = "Rejection Rate",
            color = "Innovation",
            title = rate_labels[rate_col]
          ) +
          viewer_plot_theme(base_size = 14) +
          theme(
            legend.position = "bottom"
          )

        # SE error bars
        se_col <- paste0(rate_col, "_se")
        if (se_col %in% names(df)) {
          p <- p + geom_errorbar(
            aes(ymin = .data[[rate_col]] - .data[[se_col]],
                ymax = .data[[rate_col]] + .data[[se_col]]),
            width = 0.02, alpha = 0.5
          )
        }

        # Faceting
        if (input$facet_by == "n") {
          p <- p + facet_wrap(~ n, scales = "free_x",
                              labeller = labeller(n = function(x) paste0("n = ", x)))
        } else if (input$facet_by == "innov_dist") {
          p <- p + facet_wrap(~ innov_dist)
        }

        p

      } else if (input$plot_type == "heatmap") {
        # Heatmap: n vs phi, fill = rejection rate
        df$n <- factor(df$n)
        df$phi <- factor(df$phi)
        df$label <- sprintf("%.3f", df[[rate_col]])
        df$label_col <- ifelse(
          is.na(df[[rate_col]]),
          "#212529",
          ifelse(df[[rate_col]] <= 0.02 | df[[rate_col]] >= 0.08, "#f8f9fa", "#212529")
        )

        p <- ggplot(df, aes(x = phi, y = n, fill = .data[[rate_col]])) +
          geom_tile(color = grDevices::adjustcolor(viewer_plot_fg(session), alpha.f = 0.20),
                    linewidth = 0.25) +
          scale_fill_gradient2(
            low = "#2b8cbe", mid = "#fee8c8", high = "#d7301f",
            midpoint = 0.05, name = "Rejection\nRate"
          ) +
          geom_text(aes(label = label, colour = label_col),
                    size = 4.2, fontface = "bold") +
          scale_colour_identity() +
          labs(
            x = expression(phi),
            y = "Sample Size (n)",
            title = rate_labels[rate_col]
          ) +
          viewer_plot_theme(base_size = 14) +
          theme(
            panel.grid = element_blank(),
            legend.position = "right"
          ) +
          facet_wrap(~ innov_dist)

        p

      } else if (input$plot_type == "deviation") {
        # Deviation from nominal: dot plot of (rate - 0.05) for all 3 methods
        # Build long-format data with one row per (config, method)
        dev_rows <- list()
        for (i in seq_len(nrow(df))) {
          config_label <- sprintf("%s / n=%s / phi=%s", df$innov_dist[i], df$n[i], df$phi[i])
          if ("reject_asymp_05" %in% names(df)) {
            dev_rows <- c(dev_rows, list(data.frame(
              config = config_label, method = "CO",
              deviation = df$reject_asymp_05[i] - 0.05,
              se = if ("reject_asymp_05_se" %in% names(df)) df$reject_asymp_05_se[i] else NA_real_,
              stringsAsFactors = FALSE
            )))
          }
          if ("reject_05" %in% names(df)) {
            dev_rows <- c(dev_rows, list(data.frame(
              config = config_label, method = "COB",
              deviation = df$reject_05[i] - 0.05,
              se = if ("reject_05_se" %in% names(df)) df$reject_05_se[i] else NA_real_,
              stringsAsFactors = FALSE
            )))
          }
          if ("reject_adj_05" %in% names(df)) {
            dev_rows <- c(dev_rows, list(data.frame(
              config = config_label, method = "COBA",
              deviation = df$reject_adj_05[i] - 0.05,
              se = if ("reject_adj_05_se" %in% names(df)) df$reject_adj_05_se[i] else NA_real_,
              stringsAsFactors = FALSE
            )))
          }
        }
        dev_df <- do.call(rbind, dev_rows)
        dev_df$method <- factor(dev_df$method, levels = c("CO", "COB", "COBA"))

        p <- ggplot(dev_df, aes(x = deviation, y = config, color = method)) +
          geom_vline(xintercept = 0, linetype = "dashed") +
          geom_vline(xintercept = c(-0.02, 0.02), linetype = "dotted") +
          geom_point(size = 2.5, position = position_dodge(width = 0.6)) +
          {
            if (!all(is.na(dev_df$se))) {
              geom_errorbarh(
                aes(xmin = deviation - 1.96 * se, xmax = deviation + 1.96 * se),
                height = 0.2, position = position_dodge(width = 0.6), alpha = 0.5
              )
            }
          } +
          scale_color_manual(values = c(CO = "#e41a1c", COB = "#377eb8", COBA = "#4daf4a")) +
          labs(
            x = "Deviation from Nominal (Rate - 0.05)",
            y = NULL,
            color = "Method",
            title = "Deviation from Nominal 5% Rate"
          ) +
          viewer_plot_theme(base_size = 13) +
          theme(
            legend.position = "bottom"
          )

        p

      } else if (input$plot_type == "rate_vs_n") {
        # Rate vs sample size: line plot of rate vs n
        df$n <- as.numeric(df$n)

        p <- ggplot(df, aes(x = n, y = .data[[rate_col]],
                            color = factor(innov_dist),
                            group = factor(innov_dist))) +
          geom_line(linewidth = 0.8) +
          geom_point(size = 2) +
          geom_hline(yintercept = 0.05, linetype = "dashed") +
          labs(
            x = "Sample Size (n)",
            y = "Rejection Rate",
            color = "Innovation",
            title = paste(rate_labels[rate_col], "vs Sample Size")
          ) +
          viewer_plot_theme(base_size = 14) +
          theme(
            legend.position = "bottom"
          )

        # SE error bars
        se_col <- paste0(rate_col, "_se")
        if (se_col %in% names(df)) {
          p <- p + geom_errorbar(
            aes(ymin = .data[[rate_col]] - .data[[se_col]],
                ymax = .data[[rate_col]] + .data[[se_col]]),
            width = 10, alpha = 0.5
          )
        }

        # Faceting
        facet_choice <- input$rate_vs_n_facet
        if (!is.null(facet_choice) && facet_choice == "innov_dist") {
          p <- p + facet_wrap(~ innov_dist)
        } else if (!is.null(facet_choice) && facet_choice == "phi") {
          p <- p + facet_wrap(~ phi,
                              labeller = labeller(phi = function(x) paste0("phi = ", x)))
        }

        p

      } else if (input$plot_type == "forest") {
        # Method comparison forest plot: 3 methods per config, horizontal
        forest_rows <- list()
        for (i in seq_len(nrow(df))) {
          config_label <- sprintf("%s / n=%s / phi=%s", df$innov_dist[i], df$n[i], df$phi[i])
          if ("reject_asymp_05" %in% names(df)) {
            forest_rows <- c(forest_rows, list(data.frame(
              config = config_label, method = "CO",
              rate = df$reject_asymp_05[i],
              se = if ("reject_asymp_05_se" %in% names(df)) df$reject_asymp_05_se[i] else NA_real_,
              stringsAsFactors = FALSE
            )))
          }
          if ("reject_05" %in% names(df)) {
            forest_rows <- c(forest_rows, list(data.frame(
              config = config_label, method = "COB",
              rate = df$reject_05[i],
              se = if ("reject_05_se" %in% names(df)) df$reject_05_se[i] else NA_real_,
              stringsAsFactors = FALSE
            )))
          }
          if ("reject_adj_05" %in% names(df)) {
            forest_rows <- c(forest_rows, list(data.frame(
              config = config_label, method = "COBA",
              rate = df$reject_adj_05[i],
              se = if ("reject_adj_05_se" %in% names(df)) df$reject_adj_05_se[i] else NA_real_,
              stringsAsFactors = FALSE
            )))
          }
        }
        forest_df <- do.call(rbind, forest_rows)
        forest_df$method <- factor(forest_df$method, levels = c("CO", "COB", "COBA"))

        p <- ggplot(forest_df, aes(x = rate, y = config, color = method)) +
          geom_vline(xintercept = 0.05, linetype = "dashed") +
          geom_point(size = 2.5, position = position_dodge(width = 0.6)) +
          {
            if (!all(is.na(forest_df$se))) {
              geom_errorbarh(
                aes(xmin = rate - 1.96 * se, xmax = rate + 1.96 * se),
                height = 0.2, position = position_dodge(width = 0.6), alpha = 0.5
              )
            }
          } +
          scale_color_manual(values = c(CO = "#e41a1c", COB = "#377eb8", COBA = "#4daf4a")) +
          labs(
            x = "Rejection Rate",
            y = NULL,
            color = "Method",
            title = "Method Comparison: Rejection Rate with 95% CI"
          ) +
          viewer_plot_theme(base_size = 13) +
          theme(
            legend.position = "bottom"
          )

        p
      }
    })

    # Render the plot
    output$main_plot <- renderPlot(bg = "transparent", {
      p <- current_plot()
      if (is.null(p)) {
        plot.new()
        text(0.5, 0.5, "No data matching filters", cex = 1.2, col = "grey50")
        return()
      }
      print(p)
    })

    # PNG download
    output$export_plot <- downloadHandler(
      filename = function() {
        paste0("mc_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
      },
      content = function(file) {
        p <- current_plot()
        if (!is.null(p)) {
          ggsave(file, plot = p, width = 10, height = 6, dpi = 300)
        }
      }
    )

  })
}
