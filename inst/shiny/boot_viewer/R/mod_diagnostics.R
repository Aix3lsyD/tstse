# Module: Diagnostics (Tab 8)
# Extracted from app.R -- simulation diagnostics: null model fits,
# convergence of rejection rates, batch-to-batch consistency,
# and test statistic distributions.

mod_diagnostics_ui <- function(id) {
  ns <- NS(id)

  tabPanel(
    "Diagnostics",
    br(),
    h4("Simulation Diagnostics"),
    p(class = "text-body-secondary",
      "Validates the Monte Carlo machinery: null model fits, ",
      "convergence of rejection rates, batch-to-batch consistency, ",
      "and test statistic distributions."
    ),

    div(class = "plot-controls",
      fluidRow(
        column(3,
          selectizeInput(ns("diag_n_filter"), "Sample Size (n)",
                         choices = NULL, multiple = TRUE,
                         options = list(placeholder = "All..."))
        ),
        column(3,
          selectizeInput(ns("diag_phi_filter"), "Phi Value",
                         choices = NULL, multiple = TRUE,
                         options = list(placeholder = "All..."))
        ),
        column(3,
          selectizeInput(ns("diag_innov_filter"), "Innovation Distribution",
                         choices = NULL, multiple = TRUE,
                         options = list(placeholder = "All..."))
        ),
        column(3,
          div(style = "margin-top: 25px;",
            textOutput(ns("diag_status"))
          )
        )
      )
    ),

    # AR order distribution
    wellPanel(
      h5("Fitted AR Order Distribution"),
      p(class = "text-body-secondary", style = "margin-bottom: 10px;",
        "How often each AR order is selected under the null. ",
        "For AR(1) DGPs, order 1 should dominate."),
      plotOutput(ns("diag_ar_order_plot"), height = "300px")
    ),

    # Estimated phi distribution
    wellPanel(
      h5("Estimated Null AR Coefficient Distribution"),
      p(class = "text-body-secondary", style = "margin-bottom: 10px;",
        "Distribution of the first fitted AR coefficient across simulations. ",
        "Should cluster near the true phi value."),
      plotOutput(ns("diag_phi_plot"), height = "300px")
    ),

    # MC Convergence
    wellPanel(
      h5("Monte Carlo Convergence"),
      p(class = "text-body-secondary", style = "margin-bottom: 10px;",
        "Running rejection rate as simulations accumulate. ",
        "The rate should stabilize, confirming sufficient MC sample size."),
      fluidRow(
        column(4,
          selectInput(ns("diag_conv_method"), "Method:",
                      choices = c("Bootstrap (COB)" = "pvalue",
                                  "Asymptotic (CO)" = "pvalue_asymp",
                                  "COBA" = "pvalue_adj"))
        )
      ),
      plotOutput(ns("diag_convergence_plot"), height = "350px")
    ),

    # Batch Consistency
    wellPanel(
      h5("Batch-to-Batch Consistency"),
      p(class = "text-body-secondary", style = "margin-bottom: 10px;",
        "Per-batch rejection rates for each configuration. ",
        "Low spread indicates stable Monte Carlo estimates."),
      fluidRow(
        column(4,
          selectInput(ns("diag_batch_method"), "Method:",
                      choices = c("Bootstrap (COB)" = "reject_05",
                                  "Asymptotic (CO)" = "reject_asymp_05",
                                  "COBA" = "reject_adj_05"))
        )
      ),
      plotOutput(ns("diag_batch_plot"), height = "350px")
    ),

    # Test statistic distribution
    wellPanel(
      h5("Test Statistic Distribution"),
      p(class = "text-body-secondary", style = "margin-bottom: 10px;",
        "Distribution of the observed test statistic across simulations. ",
        "Helps identify outliers or unexpected distributional shapes."),
      plotOutput(ns("diag_tstat_plot"), height = "300px")
    )
  )
}

mod_diagnostics_server <- function(id, con, init_choices) {
  moduleServer(id, function(input, output, session) {

    # --- Populate filters ---
    updateSelectizeInput(session, "diag_n_filter",
                         choices = init_choices$n, server = FALSE)
    updateSelectizeInput(session, "diag_phi_filter",
                         choices = init_choices$phi, server = FALSE)
    updateSelectizeInput(session, "diag_innov_filter",
                         choices = init_choices$innov, server = FALSE)

    # --- Helper: build filter clauses used by several queries ---
    diag_filter_clauses <- reactive({
      sql_extra <- ""
      params <- list()

      if (length(input$diag_n_filter) > 0) {
        ph <- paste(rep("?", length(input$diag_n_filter)), collapse = ", ")
        sql_extra <- paste0(sql_extra, " AND n IN (", ph, ")")
        params <- c(params, as.list(as.integer(input$diag_n_filter)))
      }
      if (length(input$diag_phi_filter) > 0) {
        ph <- paste(rep("?", length(input$diag_phi_filter)), collapse = ", ")
        sql_extra <- paste0(sql_extra, " AND phi IN (", ph, ")")
        params <- c(params, as.list(as.numeric(input$diag_phi_filter)))
      }
      if (length(input$diag_innov_filter) > 0) {
        ph <- paste(rep("?", length(input$diag_innov_filter)), collapse = ", ")
        sql_extra <- paste0(sql_extra, " AND innov_dist IN (", ph, ")")
        params <- c(params, as.list(input$diag_innov_filter))
      }

      list(sql = sql_extra, params = params)
    })

    # --- Filtered diagnostics data (independent of sidebar) ---
    diag_data <- reactive({
      fc <- diag_filter_clauses()
      sql <- paste0(
        "SELECT null_ar_order, null_ar_phi[1] AS null_phi1, n, phi, innov_dist ",
        "FROM simulations WHERE null_ar_order IS NOT NULL",
        fc$sql
      )

      tryCatch(
        dbGetQuery(con, sql, params = fc$params),
        error = function(e) data.frame()
      )
    })

    # --- Status text ---
    output$diag_status <- renderText({
      df <- diag_data()
      if (nrow(df) == 0) return("No data")
      n_configs <- nrow(unique(df[, c("n", "phi", "innov_dist"), drop = FALSE]))
      paste0(format(nrow(df), big.mark = ","), " sims, ", n_configs, " configs")
    })

    # --- AR order distribution bar chart ---
    output$diag_ar_order_plot <- renderPlot(bg = "transparent", {
      df <- diag_data()
      if (nrow(df) == 0) {
        plot.new()
        text(0.5, 0.5, "No data matching filters", cex = 1.2, col = "grey50")
        return()
      }

      df$null_ar_order <- factor(df$null_ar_order)

      # Build config label for faceting
      df$config <- sprintf("n=%s, phi=%s, %s", df$n, df$phi, df$innov_dist)

      # If many configs, facet; otherwise single plot
      n_configs <- length(unique(df$config))

      p <- ggplot(df, aes(x = null_ar_order)) +
        geom_bar(fill = "#4292c6", color = NA, alpha = 0.8) +
        labs(
          x = "Fitted AR Order",
          y = "Count",
          title = "Distribution of Fitted AR Orders Under the Null"
        ) +
        viewer_plot_theme(base_size = 14)

      if (n_configs > 1 && n_configs <= 12) {
        p <- p + facet_wrap(~ config, scales = "free_y")
      } else if (n_configs > 12) {
        p <- p + facet_wrap(~ config, scales = "free_y", ncol = 4)
      }

      print(p)
    })

    # --- Estimated phi distribution ---
    output$diag_phi_plot <- renderPlot(bg = "transparent", {
      df <- diag_data()
      if (nrow(df) == 0 || all(is.na(df$null_phi1))) {
        plot.new()
        text(0.5, 0.5, "No data matching filters", cex = 1.2, col = "grey50")
        return()
      }

      df <- df[!is.na(df$null_phi1), ]
      df$true_phi_label <- paste0("True phi = ", df$phi)

      p <- ggplot(df, aes(x = null_phi1)) +
        geom_histogram(bins = 50, fill = "#4292c6", color = NA, alpha = 0.8) +
        geom_vline(aes(xintercept = phi), color = "red", linewidth = 1, linetype = "dashed") +
        labs(
          x = "Estimated First AR Coefficient",
          y = "Count",
          title = "Distribution of Estimated Null AR(1) Coefficient"
        ) +
        viewer_plot_theme(base_size = 14)

      n_configs <- length(unique(df$true_phi_label))
      if (n_configs > 1) {
        p <- p + facet_wrap(~ true_phi_label, scales = "free")
      }

      print(p)
    })

    # --- MC Convergence: running rejection rate as sims accumulate ---
    output$diag_convergence_plot <- renderPlot(bg = "transparent", {
      pval_col <- input$diag_conv_method
      method_labels <- c(pvalue = "Bootstrap (COB)",
                         pvalue_asymp = "Asymptotic (CO)",
                         pvalue_adj = "COBA")

      fc <- diag_filter_clauses()

      # Query raw sims ordered by sim_id, with diag filters
      sql <- paste0("SELECT sim_id, ", pval_col, " AS pval, n, phi, innov_dist ",
                    "FROM simulations WHERE ", pval_col, " IS NOT NULL",
                    fc$sql,
                    " ORDER BY n, phi, innov_dist, sim_id")
      df <- tryCatch(dbGetQuery(con, sql, params = fc$params),
                     error = function(e) data.frame())

      if (nrow(df) == 0) {
        plot.new()
        text(0.5, 0.5, "No data matching filters", cex = 1.2, col = "grey50")
        return()
      }

      # Compute running rejection rate per config
      df$config <- sprintf("%s / n=%s / phi=%s", df$innov_dist, df$n, df$phi)
      configs <- unique(df$config)

      conv_rows <- list()
      for (cfg in configs) {
        sub <- df[df$config == cfg, ]
        rejected <- sub$pval < 0.05
        cum_rate <- cumsum(rejected) / seq_along(rejected)
        conv_rows <- c(conv_rows, list(data.frame(
          sim_num = seq_along(cum_rate),
          cum_reject_rate = cum_rate,
          config = cfg,
          stringsAsFactors = FALSE
        )))
      }
      conv_df <- do.call(rbind, conv_rows)

      p <- ggplot(conv_df, aes(x = sim_num, y = cum_reject_rate, color = config)) +
        geom_line(alpha = 0.7, linewidth = 0.6) +
        geom_hline(yintercept = 0.05, linetype = "dashed") +
        labs(
          x = "Simulation Number",
          y = "Cumulative Rejection Rate",
          color = "Configuration",
          title = paste0("MC Convergence: ", method_labels[pval_col])
        ) +
        viewer_plot_theme(base_size = 13) +
        theme(legend.position = "bottom",
              legend.text = element_text(size = 8, colour = viewer_plot_fg(session)))

      print(p)
    })

    # --- Batch consistency: boxplot of per-batch rejection rates ---
    output$diag_batch_plot <- renderPlot(bg = "transparent", {
      rate_col <- input$diag_batch_method
      rate_labels <- c(reject_05 = "Bootstrap (COB)",
                       reject_asymp_05 = "Asymptotic (CO)",
                       reject_adj_05 = "COBA")

      fc <- diag_filter_clauses()

      # Query per-batch view with diag filters
      sql <- paste0("SELECT * FROM v_rejection_rates_by_batch WHERE 1=1",
                    fc$sql,
                    " ORDER BY innov_dist, n, phi")
      df <- tryCatch(dbGetQuery(con, sql, params = fc$params),
                     error = function(e) data.frame())

      if (nrow(df) == 0 || !rate_col %in% names(df)) {
        plot.new()
        text(0.5, 0.5, "No per-batch data available", cex = 1.2, col = "grey50")
        return()
      }

      df$config <- sprintf("%s / n=%s / phi=%s", df$innov_dist, df$n, df$phi)

      p <- ggplot(df, aes(x = config, y = .data[[rate_col]])) +
        geom_hline(yintercept = 0.05, linetype = "dashed") +
        geom_boxplot(fill = "#5fad8c", alpha = 0.6, outlier.shape = NA) +
        geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "#2171b5") +
        labs(
          x = NULL,
          y = "Rejection Rate",
          title = paste0("Batch-to-Batch Consistency: ", rate_labels[rate_col])
        ) +
        viewer_plot_theme(base_size = 13) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9,
                                         colour = viewer_plot_fg(session)))

      print(p)
    })

    # --- Test statistic distribution: histogram of obs_stat per config ---
    output$diag_tstat_plot <- renderPlot(bg = "transparent", {
      fc <- diag_filter_clauses()

      sql <- paste0("SELECT obs_stat, n, phi, innov_dist FROM simulations ",
                    "WHERE obs_stat IS NOT NULL",
                    fc$sql)
      df <- tryCatch(dbGetQuery(con, sql, params = fc$params),
                     error = function(e) data.frame())

      if (nrow(df) == 0) {
        plot.new()
        text(0.5, 0.5, "No data matching filters", cex = 1.2, col = "grey50")
        return()
      }

      df$config <- sprintf("%s / n=%s / phi=%s", df$innov_dist, df$n, df$phi)
      n_configs <- length(unique(df$config))

      p <- ggplot(df, aes(x = obs_stat)) +
        geom_histogram(bins = 50, fill = "#4292c6", color = NA, alpha = 0.8) +
        labs(
          x = "Observed Test Statistic",
          y = "Count",
          title = "Distribution of Observed Test Statistics"
        ) +
        viewer_plot_theme(base_size = 14)

      if (n_configs > 1 && n_configs <= 12) {
        p <- p + facet_wrap(~ config, scales = "free")
      } else if (n_configs > 12) {
        p <- p + facet_wrap(~ config, scales = "free", ncol = 4)
      }

      print(p)
    })

  })
}
