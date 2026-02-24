# Module: Bootstrap Distribution (Tab 5)
# Cascading selectors: Batch -> Config -> Simulation
# Displays histogram of bootstrap t-statistics with observed statistic overlay

mod_bootstrap_dist_ui <- function(id) {
  ns <- NS(id)
  tabPanel(
    "Bootstrap Distribution",
    br(),
    wellPanel(class = "plot-controls",
      fluidRow(
        column(3,
          selectInput(ns("dist_batch"), "Batch",
                      choices = c("All" = "all"))
        ),
        column(4,
          selectizeInput(
            ns("dist_config"),
            "Configuration (n, phi, innov_dist)",
            choices = NULL,
            options = list(placeholder = "Select configuration...")
          )
        ),
        column(5,
          selectizeInput(
            ns("dist_sim"),
            "Simulation",
            choices = NULL,
            options = list(placeholder = "Select simulation...")
          )
        )
      )
    ),
    plotOutput(ns("dist_plot"), height = "400px"),
    verbatimTextOutput(ns("dist_summary"))
  )
}

mod_bootstrap_dist_server <- function(id, con) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # -----------------------------------------------------------------------
    # Populate batch selector
    # -----------------------------------------------------------------------
    observe({
      batches <- tryCatch(
        dbGetQuery(con, "SELECT batch_id, label FROM batches ORDER BY batch_id"),
        error = function(e) data.frame()
      )
      if (nrow(batches) == 0) {
        choices <- c("All" = "all")
      } else {
        batch_labels <- ifelse(
          is.na(batches$label) | batches$label == "",
          paste("Batch", batches$batch_id),
          paste0("Batch ", batches$batch_id, ": ", batches$label)
        )
        choices <- c("All" = "all",
                      setNames(as.character(batches$batch_id), batch_labels))
      }
      updateSelectInput(session, "dist_batch", choices = choices)
    })

    # -----------------------------------------------------------------------
    # Cascade: batch -> config
    # -----------------------------------------------------------------------
    observe({
      batch_sel <- input$dist_batch
      req(batch_sel)

      if (batch_sel == "all") {
        configs <- tryCatch(
          dbGetQuery(con,
            "SELECT DISTINCT n, phi, innov_dist FROM simulations ORDER BY innov_dist, n, phi"),
          error = function(e) data.frame()
        )
      } else {
        configs <- tryCatch(
          dbGetQuery(con,
            "SELECT DISTINCT n, phi, innov_dist FROM simulations WHERE batch_id = ? ORDER BY innov_dist, n, phi",
            params = list(as.integer(batch_sel))),
          error = function(e) data.frame()
        )
      }

      if (nrow(configs) == 0) {
        updateSelectizeInput(
          session, "dist_config",
          choices = c("No configs" = ""),
          selected = "",
          server = TRUE
        )
        return()
      }

      config_labels <- sprintf("n=%d, phi=%.2f, %s", configs$n, configs$phi, configs$innov_dist)
      config_vals <- sprintf("%d|%.4f|%s", configs$n, configs$phi, configs$innov_dist)
      choices <- setNames(config_vals, config_labels)
      updateSelectizeInput(
        session, "dist_config",
        choices = choices,
        selected = unname(config_vals[1]),
        server = TRUE
      )
    })

    # -----------------------------------------------------------------------
    # Cascade: config -> simulation
    # -----------------------------------------------------------------------
    observe({
      config_sel <- input$dist_config
      req(config_sel, nchar(config_sel) > 0)

      parts <- strsplit(config_sel, "\\|")[[1]]
      if (length(parts) != 3) return()
      sel_n <- as.integer(parts[1])
      sel_phi <- as.numeric(parts[2])
      sel_innov <- parts[3]

      batch_sel <- input$dist_batch

      if (batch_sel == "all") {
        sims <- tryCatch(
          dbGetQuery(con,
            "SELECT sim_id, iteration, obs_stat FROM simulations
             WHERE n = ? AND phi = ? AND innov_dist = ?
             ORDER BY iteration",
            params = list(sel_n, sel_phi, sel_innov)),
          error = function(e) data.frame()
        )
      } else {
        sims <- tryCatch(
          dbGetQuery(con,
            "SELECT sim_id, iteration, obs_stat FROM simulations
             WHERE batch_id = ? AND n = ? AND phi = ? AND innov_dist = ?
             ORDER BY iteration",
            params = list(as.integer(batch_sel), sel_n, sel_phi, sel_innov)),
          error = function(e) data.frame()
        )
      }

      if (nrow(sims) == 0) {
        updateSelectizeInput(
          session, "dist_sim",
          choices = c("No simulations" = ""),
          selected = "",
          server = TRUE
        )
        return()
      }

      sim_labels <- sprintf("Iter %d (t = %.3f)", sims$iteration, sims$obs_stat)
      choices <- setNames(as.character(sims$sim_id), sim_labels)
      updateSelectizeInput(
        session, "dist_sim",
        choices = choices,
        selected = as.character(sims$sim_id[1]),
        server = TRUE
      )
    })

    # -----------------------------------------------------------------------
    # Reactive: load selected simulation data
    # -----------------------------------------------------------------------
    dist_data <- reactive({
      sim_id <- input$dist_sim
      req(sim_id, nchar(sim_id) > 0)
      tryCatch({
        result <- dbGetQuery(con,
          "SELECT sim_id, obs_stat, boot_dist, pvalue, n, phi, innov_dist
           FROM simulations WHERE sim_id = ?",
          params = list(as.integer(sim_id)))
        if (nrow(result) == 0) return(NULL)
        result
      }, error = function(e) NULL)
    })

    # -----------------------------------------------------------------------
    # Plot: bootstrap distribution histogram
    # -----------------------------------------------------------------------
    output$dist_plot <- renderPlot(bg = "transparent", {
      d <- dist_data()
      if (is.null(d)) {
        plot.new()
        text(0.5, 0.5, "Select a batch, config, and simulation above",
             cex = 1.2, col = "grey50")
        return()
      }

      boot_dist <- d$boot_dist[[1]]
      obs_stat <- d$obs_stat[1]
      pval <- d$pvalue[1]

      hist_data <- data.frame(t_stat = boot_dist)

      p <- ggplot(hist_data, aes(x = t_stat)) +
        geom_histogram(bins = 40, fill = "#4292c6", color = NA, alpha = 0.8) +
        geom_vline(xintercept = obs_stat, color = "red", linewidth = 1.2, linetype = "solid") +
        annotate("text", x = obs_stat, y = Inf, vjust = 2, hjust = -0.1,
                 label = sprintf("T_obs = %.3f\np = %.4f", obs_stat, pval),
                 color = "red", size = 4, fontface = "bold") +
        labs(
          x = "Bootstrap t-statistic",
          y = "Count",
          title = sprintf("Bootstrap Distribution (sim_id = %d, n=%d, phi=%.2f, %s)",
                          d$sim_id[1], d$n[1], d$phi[1], d$innov_dist[1])
        ) +
        viewer_plot_theme(base_size = 14)

      print(p)
    })

    # -----------------------------------------------------------------------
    # Text summary: simulation details and p-value recalculation
    # -----------------------------------------------------------------------
    output$dist_summary <- renderPrint({
      d <- dist_data()
      if (is.null(d)) return(cat("Select a simulation using the selectors above"))

      boot_dist <- d$boot_dist[[1]]
      cat(sprintf("Simulation ID: %d\n", d$sim_id[1]))
      cat(sprintf("Config: n=%d, phi=%.2f, innov_dist=%s\n", d$n[1], d$phi[1], d$innov_dist[1]))
      cat(sprintf("Observed t-stat: %.4f\n", d$obs_stat[1]))
      cat(sprintf("Bootstrap replicates: %d\n", length(boot_dist)))
      cat(sprintf("Bootstrap mean: %.4f, sd: %.4f\n", mean(boot_dist), sd(boot_dist)))
      cat(sprintf("P-value (stored): %.4f\n", d$pvalue[1]))

      # Recalculate
      nb <- length(boot_dist)
      recalc_p <- (sum(abs(boot_dist) >= abs(d$obs_stat[1])) + 1) / (nb + 1)
      cat(sprintf("P-value (recalculated): %.4f\n", recalc_p))
    })
  })
}
