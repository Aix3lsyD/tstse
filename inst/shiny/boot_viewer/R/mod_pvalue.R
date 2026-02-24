# Module: P-value Diagnostics (Tab 4)
# Extracted from monolithic app.R

mod_pvalue_ui <- function(id) {
  ns <- NS(id)

  tabPanel(
    "P-value Diagnostics",
    br(),
    h4("P-value Calibration Diagnostics"),
    p(class = "text-body-secondary",
      "Under H0 (phi = 0), p-values should be uniformly distributed on [0, 1]. ",
      "Non-uniformity indicates size distortion. Filter to phi = 0 configs to check calibration."
    ),

    wellPanel(class = "plot-controls",
      fluidRow(
        column(3,
          selectizeInput(ns("pval_n_filter"), "Sample Size (n)",
                         choices = NULL, multiple = TRUE,
                         options = list(placeholder = "All..."))
        ),
        column(3,
          selectizeInput(ns("pval_phi_filter"), "Phi Value",
                         choices = NULL, multiple = TRUE,
                         options = list(placeholder = "All..."))
        ),
        column(3,
          selectizeInput(ns("pval_innov_filter"), "Innovation Distribution",
                         choices = NULL, multiple = TRUE,
                         options = list(placeholder = "All..."))
        ),
        column(1,
          div(style = "margin-top: 25px;",
            actionButton(ns("pval_clear"), "Clear",
                         icon = icon("times-circle"),
                         class = "btn-sm btn-outline-secondary")
          )
        ),
        column(2,
          div(style = "margin-top: 25px;",
            textOutput(ns("pval_status"))
          )
        )
      )
    ),

    # P-value histogram
    wellPanel(
      h5("P-value Histograms"),
      p(class = "text-body-secondary", style = "margin-bottom: 10px;",
        "Histograms of p-values for each method. ",
        "Under H0, these should be approximately uniform (flat)."),
      plotOutput(ns("pval_hist"), height = "300px")
    ),

    # QQ plot
    wellPanel(
      h5("P-value QQ Plot"),
      p(class = "text-body-secondary", style = "margin-bottom: 10px;",
        "Quantile-quantile plot against Uniform(0, 1). ",
        "Points on the diagonal indicate well-calibrated p-values."),
      plotOutput(ns("pval_qq"), height = "350px")
    ),

    # Method comparison scatter
    wellPanel(
      h5("P-value Method Comparison"),
      p(class = "text-body-secondary", style = "margin-bottom: 10px;",
        "Scatter plots comparing p-values between methods. ",
        "Points on the diagonal indicate agreement."),
      plotOutput(ns("pval_scatter"), height = "350px")
    )
  )
}

mod_pvalue_server <- function(id, con, init_choices) {
  moduleServer(id, function(input, output, session) {

    # --- Filter initialization ------------------------------------------------
    updateSelectizeInput(session, "pval_n_filter",
                         choices = init_choices$n, server = FALSE)
    updateSelectizeInput(session, "pval_phi_filter",
                         choices = init_choices$phi, server = FALSE)
    updateSelectizeInput(session, "pval_innov_filter",
                         choices = init_choices$innov, server = FALSE)

    # --- Clear button ---------------------------------------------------------
    observeEvent(input$pval_clear, {
      updateSelectizeInput(session, "pval_n_filter", selected = character(0))
      updateSelectizeInput(session, "pval_phi_filter", selected = character(0))
      updateSelectizeInput(session, "pval_innov_filter", selected = character(0))
    })

    # --- Filtered p-value data (independent of sidebar) -----------------------
    pval_data <- reactive({
      sql <- "SELECT pvalue, pvalue_asymp, pvalue_adj, n, phi, innov_dist FROM simulations WHERE 1=1"
      params <- list()

      if (length(input$pval_n_filter) > 0) {
        ph <- paste(rep("?", length(input$pval_n_filter)), collapse = ", ")
        sql <- paste0(sql, " AND n IN (", ph, ")")
        params <- c(params, as.list(as.integer(input$pval_n_filter)))
      }
      if (length(input$pval_phi_filter) > 0) {
        ph <- paste(rep("?", length(input$pval_phi_filter)), collapse = ", ")
        sql <- paste0(sql, " AND phi IN (", ph, ")")
        params <- c(params, as.list(as.numeric(input$pval_phi_filter)))
      }
      if (length(input$pval_innov_filter) > 0) {
        ph <- paste(rep("?", length(input$pval_innov_filter)), collapse = ", ")
        sql <- paste0(sql, " AND innov_dist IN (", ph, ")")
        params <- c(params, as.list(input$pval_innov_filter))
      }

      tryCatch(
        dbGetQuery(con, sql, params = params),
        error = function(e) data.frame()
      )
    })

    # --- Status text ----------------------------------------------------------
    output$pval_status <- renderText({
      df <- pval_data()
      if (nrow(df) == 0) return("No data")
      n_configs <- nrow(unique(df[, c("n", "phi", "innov_dist"), drop = FALSE]))
      paste0(format(nrow(df), big.mark = ","), " sims, ", n_configs, " configs")
    })

    # --- P-value histogram: side-by-side for 3 methods ------------------------
    output$pval_hist <- renderPlot(bg = "transparent", {
      df <- pval_data()
      if (nrow(df) == 0) {
        plot.new()
        text(0.5, 0.5, "No data matching filters", cex = 1.2, col = "grey50")
        return()
      }

      # Build long-format: one row per (sim, method)
      long_rows <- list()
      if ("pvalue" %in% names(df) && !all(is.na(df$pvalue))) {
        long_rows <- c(long_rows, list(data.frame(
          pval = df$pvalue[!is.na(df$pvalue)], method = "COB (Bootstrap)",
          stringsAsFactors = FALSE
        )))
      }
      if ("pvalue_asymp" %in% names(df) && !all(is.na(df$pvalue_asymp))) {
        long_rows <- c(long_rows, list(data.frame(
          pval = df$pvalue_asymp[!is.na(df$pvalue_asymp)], method = "CO (Asymptotic)",
          stringsAsFactors = FALSE
        )))
      }
      if ("pvalue_adj" %in% names(df) && !all(is.na(df$pvalue_adj))) {
        long_rows <- c(long_rows, list(data.frame(
          pval = df$pvalue_adj[!is.na(df$pvalue_adj)], method = "COBA (Adjusted)",
          stringsAsFactors = FALSE
        )))
      }

      if (length(long_rows) == 0) {
        plot.new()
        text(0.5, 0.5, "No p-value data available", cex = 1.2, col = "grey50")
        return()
      }

      pval_long <- do.call(rbind, long_rows)

      p <- ggplot(pval_long, aes(x = pval)) +
        geom_histogram(bins = 20, fill = "#4292c6", color = NA, alpha = 0.8,
                       boundary = 0) +
        facet_wrap(~ method, scales = "free_y") +
        labs(x = "P-value", y = "Count",
             title = "P-value Distributions (uniform = well-calibrated)") +
        viewer_plot_theme(base_size = 13) +
        theme(strip.text = element_text(face = "bold", colour = viewer_plot_fg(session)))

      # Add uniform reference line (expected count per bin)
      n_per_method <- table(pval_long$method)
      ref_df <- data.frame(
        method = names(n_per_method),
        yint = as.numeric(n_per_method) / 20,
        stringsAsFactors = FALSE
      )
      p <- p + geom_hline(data = ref_df, aes(yintercept = yint),
                          color = "red", linetype = "dashed", linewidth = 0.6)

      print(p)
    })

    # --- P-value QQ plot against U(0,1) ---------------------------------------
    output$pval_qq <- renderPlot(bg = "transparent", {
      df <- pval_data()
      if (nrow(df) == 0) {
        plot.new()
        text(0.5, 0.5, "No data matching filters", cex = 1.2, col = "grey50")
        return()
      }

      # Build long-format for QQ
      qq_rows <- list()
      methods <- c(pvalue = "COB (Bootstrap)", pvalue_asymp = "CO (Asymptotic)",
                   pvalue_adj = "COBA (Adjusted)")

      for (col in names(methods)) {
        if (col %in% names(df)) {
          vals <- df[[col]][!is.na(df[[col]])]
          if (length(vals) > 0) {
            n_vals <- length(vals)
            theoretical <- (seq_len(n_vals) - 0.5) / n_vals
            qq_rows <- c(qq_rows, list(data.frame(
              theoretical = theoretical,
              observed = sort(vals),
              method = methods[col],
              stringsAsFactors = FALSE
            )))
          }
        }
      }

      if (length(qq_rows) == 0) {
        plot.new()
        text(0.5, 0.5, "No p-value data available", cex = 1.2, col = "grey50")
        return()
      }

      qq_df <- do.call(rbind, qq_rows)

      p <- ggplot(qq_df, aes(x = theoretical, y = observed)) +
        geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
        geom_point(size = 0.8, alpha = 0.5, color = "#2171b5") +
        facet_wrap(~ method) +
        labs(x = "Theoretical Quantiles (Uniform)",
             y = "Observed P-value Quantiles",
             title = "QQ Plot: P-values vs Uniform(0, 1)") +
        coord_equal() +
        viewer_plot_theme(base_size = 13) +
        theme(strip.text = element_text(face = "bold", colour = viewer_plot_fg(session)))

      print(p)
    })

    # --- P-value scatter: bootstrap vs asymptotic, bootstrap vs COBA ----------
    output$pval_scatter <- renderPlot(bg = "transparent", {
      df <- pval_data()
      if (nrow(df) == 0) {
        plot.new()
        text(0.5, 0.5, "No data matching filters", cex = 1.2, col = "grey50")
        return()
      }

      has_asymp <- "pvalue_asymp" %in% names(df) && !all(is.na(df$pvalue_asymp))
      has_adj <- "pvalue_adj" %in% names(df) && !all(is.na(df$pvalue_adj))
      has_boot <- "pvalue" %in% names(df) && !all(is.na(df$pvalue))

      if (!has_boot || (!has_asymp && !has_adj)) {
        plot.new()
        text(0.5, 0.5, "Need at least 2 methods with p-values", cex = 1.2, col = "grey50")
        return()
      }

      scatter_rows <- list()
      if (has_asymp) {
        idx <- !is.na(df$pvalue) & !is.na(df$pvalue_asymp)
        scatter_rows <- c(scatter_rows, list(data.frame(
          boot_pval = df$pvalue[idx],
          other_pval = df$pvalue_asymp[idx],
          comparison = "COB vs CO (Asymptotic)",
          stringsAsFactors = FALSE
        )))
      }
      if (has_adj) {
        idx <- !is.na(df$pvalue) & !is.na(df$pvalue_adj)
        scatter_rows <- c(scatter_rows, list(data.frame(
          boot_pval = df$pvalue[idx],
          other_pval = df$pvalue_adj[idx],
          comparison = "COB vs COBA (Adjusted)",
          stringsAsFactors = FALSE
        )))
      }

      scatter_df <- do.call(rbind, scatter_rows)

      p <- ggplot(scatter_df, aes(x = boot_pval, y = other_pval)) +
        geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
        geom_point(size = 0.6, alpha = 0.2, color = "#2171b5") +
        facet_wrap(~ comparison) +
        labs(x = "Bootstrap P-value (COB)",
             y = "Comparison Method P-value",
             title = "P-value Agreement Between Methods") +
        coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
        viewer_plot_theme(base_size = 13) +
        theme(strip.text = element_text(face = "bold", colour = viewer_plot_fg(session)))

      print(p)
    })

  })
}
