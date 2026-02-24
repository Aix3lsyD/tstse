# Module: Parallel Coordinates tab (Tab 7)
# Extracted from monolithic app.R

mod_parallel_coords_ui <- function(id) {
  ns <- NS(id)

  tabPanel(
    "Parallel Coordinates",
    br(),
    h4("Parallel Coordinates Plots"),
    p(class = "text-body-secondary",
      "Each line represents one (n, phi, innov_dist) configuration. ",
      "Rate axes show raw rejection rate values. ",
      "Use the filters below to reduce clutter."
    ),

    # Shared filters for this tab
    div(class = "plot-controls",
      fluidRow(
        column(3,
          selectizeInput(ns("pc_n_filter"), "Sample Size (n)",
                         choices = NULL, multiple = TRUE,
                         options = list(placeholder = "All..."))
        ),
        column(3,
          selectizeInput(ns("pc_phi_filter"), "Phi Value",
                         choices = NULL, multiple = TRUE,
                         options = list(placeholder = "All..."))
        ),
        column(3,
          selectizeInput(ns("pc_innov_filter"), "Innovation Distribution",
                         choices = NULL, multiple = TRUE,
                         options = list(placeholder = "All..."))
        ),
        column(1,
          div(style = "margin-top: 25px;",
            actionButton(ns("pc_clear"), "Clear",
                         icon = icon("times-circle"),
                         class = "btn-sm btn-outline-secondary")
          )
        )
      )
    ),

    # Plot 1a: Rate Profile (CO -> COB -> COBA)
    wellPanel(
      h5("Rate Profile"),
      p(class = "text-body-secondary", style = "margin-bottom: 10px;",
        "Axes: CO \u2192 COB \u2192 COBA (raw rejection rates). ",
        "Color by innovation distribution."),
      plotOutput(ns("pc_rate_profile_plot"), height = "350px")
    ),

    # Plot 1b: DGP Parameter Relationships
    wellPanel(
      h5("DGP Parameter Relationships"),
      p(class = "text-body-secondary", style = "margin-bottom: 10px;",
        "Rejection rates vs. DGP parameters (n and phi). ",
        "Faceted scatter view for parameter relationships."),
      plotOutput(ns("pc_dgp_scatter_plot"), height = "350px")
    ),

    # Plot 2: Method Comparison
    wellPanel(
      h5("Method Comparison"),
      p(class = "text-body-secondary", style = "margin-bottom: 10px;",
        "Axes: CO \u2192 COB \u2192 COBA. ",
        "Focused view comparing how the three methods relate."),
      fluidRow(
        column(3,
          selectInput(ns("pc_color_by"), "Color by:",
                      choices = c("Innovation" = "innov_dist",
                                  "Sample Size (n)" = "n",
                                  "Phi" = "phi"))
        )
      ),
      plotOutput(ns("pc_method_plot"), height = "350px")
    ),

    # Plot 3: Single-Method Sensitivity
    wellPanel(
      h5("Single-Method Sensitivity"),
      p(class = "text-body-secondary", style = "margin-bottom: 10px;",
        "Axes: n \u2192 phi \u2192 Rate. ",
        "How one method responds to DGP parameters."),
      fluidRow(
        column(3,
          selectInput(ns("pc_method"), "Method:",
                      choices = c("Bootstrap (COB)" = "reject_05",
                                  "Asymptotic (CO)" = "reject_asymp_05",
                                  "COBA" = "reject_adj_05"))
        )
      ),
      plotOutput(ns("pc_sensitivity_plot"), height = "350px")
    )
  )
}

mod_parallel_coords_server <- function(id, con, init_choices) {
  moduleServer(id, function(input, output, session) {

    # --- Filter initialization -----------------------------------------------
    updateSelectizeInput(session, "pc_n_filter",
                         choices = init_choices$n, server = FALSE)
    updateSelectizeInput(session, "pc_phi_filter",
                         choices = init_choices$phi, server = FALSE)
    updateSelectizeInput(session, "pc_innov_filter",
                         choices = init_choices$innov, server = FALSE)

    # --- Clear filters -------------------------------------------------------
    observeEvent(input$pc_clear, {
      updateSelectizeInput(session, "pc_n_filter", selected = character(0))
      updateSelectizeInput(session, "pc_phi_filter", selected = character(0))
      updateSelectizeInput(session, "pc_innov_filter", selected = character(0))
    })

    # --- Filtered data (independent of sidebar) ------------------------------
    pc_data <- reactive({
      q <- build_query("v_rejection_rates",
                        n_vals     = input$pc_n_filter,
                        phi_vals   = input$pc_phi_filter,
                        innov_vals = input$pc_innov_filter)
      tryCatch(
        dbGetQuery(con, q$sql, params = q$params),
        error = function(e) data.frame()
      )
    })

    # --- Local helper: build a parallel coordinates ggplot -------------------
    # cols      = named vector of column_name -> axis_label
    # color_col = column to use for colour aesthetic
    # color_label = legend title
    # rate_cols = column names that are rejection rates (shared raw scale)
    build_parcoord <- function(df, cols, color_col, color_label, title = "",
                               rate_cols = c("reject_asymp_05", "reject_05",
                                             "reject_adj_05")) {
      if (nrow(df) == 0) return(NULL)

      col_names  <- names(cols)
      axis_labels <- unname(cols)
      n_axes     <- length(col_names)

      # Assign config ID to each row
      df$config_id <- seq_len(nrow(df))

      # Determine fixed rate scale: [0, max(max_rate * 1.1, 0.5)]
      rate_in_cols <- intersect(col_names, rate_cols)
      if (length(rate_in_cols) > 0) {
        all_rates <- unlist(lapply(rate_in_cols, function(c) as.numeric(df[[c]])))
        rate_max <- max(max(all_rates, na.rm = TRUE) * 1.1, 0.5)
      } else {
        rate_max <- 1
      }

      # Scale each column: rate columns use shared [0, rate_max] scale,
      # other columns use their own [min, max] range
      ranges <- list()
      scaled <- data.frame(config_id = df$config_id)

      for (col in col_names) {
        vals <- as.numeric(df[[col]])
        if (col %in% rate_cols) {
          # Raw rate scale: map [0, rate_max] to [0, 1] for plotting
          ranges[[col]] <- c(0, rate_max)
          scaled[[col]] <- vals / rate_max
        } else {
          rng <- range(vals, na.rm = TRUE)
          ranges[[col]] <- rng
          if (rng[1] == rng[2]) {
            scaled[[col]] <- 0.5
          } else {
            scaled[[col]] <- (vals - rng[1]) / (rng[2] - rng[1])
          }
        }
      }

      # Build long-format data for geom_line
      long_rows <- vector("list", nrow(df) * n_axes)
      idx <- 1
      for (i in seq_len(nrow(df))) {
        for (j in seq_len(n_axes)) {
          long_rows[[idx]] <- data.frame(
            config_id  = i,
            axis_num   = j,
            axis_label = axis_labels[j],
            value      = scaled[[col_names[j]]][i],
            color_var  = as.character(df[[color_col]][i]),
            stringsAsFactors = FALSE
          )
          idx <- idx + 1
        }
      }
      long_df <- do.call(rbind, long_rows)
      long_df$axis_num <- as.integer(long_df$axis_num)

      # Build annotation data for axis tick labels (min/max at each axis)
      tick_labels <- unlist(lapply(col_names, function(col) {
        rng <- ranges[[col]]
        c(format(rng[1], digits = 3), format(rng[2], digits = 3))
      }))
      tick_df <- data.frame(
        axis_num = rep(seq_len(n_axes), each = 2),
        y        = rep(c(0, 1), times = n_axes),
        label    = tick_labels,
        stringsAsFactors = FALSE
      )

      # Add 0.05 reference line for rate axes (scaled position)
      ref_lines <- data.frame()
      for (j in seq_along(col_names)) {
        if (col_names[j] %in% rate_cols) {
          ref_y <- 0.05 / rate_max
          ref_lines <- rbind(ref_lines, data.frame(
            axis_num = j, y = ref_y, stringsAsFactors = FALSE))
        }
      }

      p <- ggplot(long_df, aes(x = axis_num, y = value,
                                group = config_id, color = color_var)) +
        geom_line(alpha = 0.6, linewidth = 0.7) +
        geom_point(size = 1.5, alpha = 0.7) +
        scale_x_continuous(breaks = seq_len(n_axes), labels = axis_labels,
                           limits = c(0.5, n_axes + 0.5)) +
        scale_y_continuous(limits = c(-0.05, 1.05),
                           breaks = seq(0, 1, by = 0.25)) +
        # Axis tick labels showing original values
        geom_text(data = tick_df,
                  aes(x = axis_num, y = y, label = label),
                  inherit.aes = FALSE, size = 2.8,
                  hjust = 1.2, fontface = "italic") +
        # Vertical axis lines
        geom_vline(xintercept = seq_len(n_axes),
                   linewidth = 0.3) +
        labs(x = NULL, y = NULL, color = color_label, title = title) +
        viewer_plot_theme(base_size = 13) +
        theme(
          panel.grid.major.x = element_blank(),
          panel.grid.minor   = element_blank(),
          legend.position     = "bottom"
        )

      # Add 0.05 reference markers on rate axes
      if (nrow(ref_lines) > 0) {
        p <- p + geom_point(data = ref_lines,
                            aes(x = axis_num, y = y),
                            inherit.aes = FALSE, shape = 4, size = 3,
                            stroke = 1.2)
      }

      p
    }

    # --- Plot 1a: Rate Profile (CO -> COB -> COBA, raw scales) ---------------
    output$pc_rate_profile_plot <- renderPlot(bg = "transparent", {
      df <- pc_data()
      if (nrow(df) == 0) {
        plot.new()
        text(0.5, 0.5, "No data matching filters",
             cex = 1.2, col = "grey50")
        return()
      }

      p <- build_parcoord(
        df,
        cols = c(reject_asymp_05 = "CO Rate",
                 reject_05       = "COB Rate",
                 reject_adj_05   = "COBA Rate"),
        color_col   = "innov_dist",
        color_label = "Innovation",
        title = "Rate Profile: CO \u2192 COB \u2192 COBA (raw rejection rates)"
      )

      if (!is.null(p)) print(p)
    })

    # --- Plot 1b: DGP Parameter Relationships (scatter / faceted view) -------
    output$pc_dgp_scatter_plot <- renderPlot(bg = "transparent", {
      df <- pc_data()
      if (nrow(df) == 0) {
        plot.new()
        text(0.5, 0.5, "No data matching filters",
             cex = 1.2, col = "grey50")
        return()
      }

      # Build long-format: one row per (config, method)
      rate_map <- c(reject_asymp_05 = "CO",
                    reject_05       = "COB",
                    reject_adj_05   = "COBA")
      long_rows <- list()
      for (i in seq_len(nrow(df))) {
        for (rate_col in names(rate_map)) {
          if (rate_col %in% names(df)) {
            long_rows <- c(long_rows, list(data.frame(
              n          = df$n[i],
              phi        = df$phi[i],
              innov_dist = df$innov_dist[i],
              method     = rate_map[rate_col],
              rate       = df[[rate_col]][i],
              stringsAsFactors = FALSE
            )))
          }
        }
      }
      long_df <- do.call(rbind, long_rows)
      long_df$method <- factor(long_df$method, levels = c("CO", "COB", "COBA"))

      p <- ggplot(long_df, aes(x = phi, y = rate, color = factor(n),
                                shape = innov_dist)) +
        geom_point(size = 2.5, alpha = 0.7) +
        geom_hline(yintercept = 0.05, linetype = "dashed") +
        facet_wrap(~ method) +
        labs(x     = expression(phi),
             y     = "Rejection Rate",
             color = "n",
             shape = "Innovation",
             title = "DGP Parameter Relationships: Rate vs phi by n") +
        viewer_plot_theme(base_size = 13) +
        theme(legend.position = "bottom")

      print(p)
    })

    # --- Plot 2: Method Comparison (CO -> COB -> COBA, user-selected color) --
    output$pc_method_plot <- renderPlot(bg = "transparent", {
      df <- pc_data()
      if (nrow(df) == 0) {
        plot.new()
        text(0.5, 0.5, "No data matching filters",
             cex = 1.2, col = "grey50")
        return()
      }

      color_col    <- input$pc_color_by
      color_labels <- c(innov_dist = "Innovation",
                        n          = "Sample Size",
                        phi        = "Phi")
      # Ensure color column is character for consistent legend
      df[[color_col]] <- as.character(df[[color_col]])

      p <- build_parcoord(
        df,
        cols = c(reject_asymp_05 = "CO Rate",
                 reject_05       = "COB Rate",
                 reject_adj_05   = "COBA Rate"),
        color_col   = color_col,
        color_label = color_labels[color_col],
        title = "Method Comparison: CO \u2192 COB \u2192 COBA"
      )

      if (!is.null(p)) print(p)
    })

    # --- Plot 3: Single-Method Sensitivity (n -> phi -> rate) ----------------
    output$pc_sensitivity_plot <- renderPlot(bg = "transparent", {
      df <- pc_data()
      if (nrow(df) == 0) {
        plot.new()
        text(0.5, 0.5, "No data matching filters",
             cex = 1.2, col = "grey50")
        return()
      }

      rate_col    <- input$pc_method
      rate_labels <- c(reject_05       = "COB Rate",
                       reject_asymp_05 = "CO Rate",
                       reject_adj_05   = "COBA Rate")

      cols <- c("n", "phi", rate_col)
      names(cols) <- NULL
      named_cols <- setNames(cols, cols)
      named_cols <- c(n = "n", phi = "phi")
      named_cols[[rate_col]] <- rate_labels[rate_col]

      p <- build_parcoord(
        df,
        cols        = named_cols,
        color_col   = "innov_dist",
        color_label = "Innovation",
        title = paste0("Sensitivity: n \u2192 phi \u2192 ",
                       rate_labels[rate_col])
      )

      if (!is.null(p)) print(p)
    })
  })
}
