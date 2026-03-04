# Capstone Sub-Tab: Coverage / Overview (Tab 2)
# Shows DB coverage KPIs + heatmap + coverage table

mod_capstone_overview_ui <- function(ns) {
  tabPanel("Coverage / Overview",
    br(),
    uiOutput(ns("overview_kpis")),
    p(class = "text-body-secondary",
      "Shows how many simulations already exist for each grid cell."),
    plotOutput(ns("coverage_heatmap"), height = "450px"),
    br(),
    DTOutput(ns("coverage_table"))
  )
}

mod_capstone_overview_server <- function(input, output, session,
                                         con, grid_rv, raw_results_rv,
                                         db_refresh_trigger, db_save_path_input) {

  # Coverage data: merge DB + in-memory counts
  coverage_data <- reactive({
    grid <- grid_rv()
    if (nrow(grid) == 0) return(data.frame())
    db_refresh_trigger()

    read_con <- con
    own_con <- FALSE
    if (is.null(read_con)) {
      save_path <- trimws(db_save_path_input() %||% "")
      if (nzchar(save_path) && file.exists(save_path)) {
        read_con <- tryCatch(mc_db_connect(save_path, read_only = TRUE),
                             error = function(e) NULL)
        own_con <- TRUE
      }
    }

    n_existing <- integer(nrow(grid))
    if (!is.null(read_con)) {
      if (own_con) on.exit(DBI::dbDisconnect(read_con, shutdown = TRUE))
      for (i in seq_len(nrow(grid))) {
        cnt <- tryCatch(
          DBI::dbGetQuery(read_con,
            "SELECT COUNT(*) AS cnt FROM simulations WHERE n = ? AND phi = ? AND innov_dist = ?",
            params = list(as.integer(grid$n[i]),
                          as.numeric(grid$phi[i]),
                          grid$innov_dist_str[i]))$cnt,
          error = function(e) 0L)
        n_existing[i] <- cnt
      }
    }

    # Also check in-memory raw results
    raw <- raw_results_rv()
    for (i in seq_len(nrow(grid))) {
      cell_key <- paste(grid$n[i], grid$phi[i], grid$innov_dist_str[i], sep = "|")
      if (cell_key %in% names(raw)) {
        mem_count <- length(raw[[cell_key]]$results)
        if (n_existing[i] == 0) n_existing[i] <- mem_count
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
    ns <- session$ns
    db_refresh_trigger()

    read_con <- con
    own_con <- FALSE
    if (is.null(read_con)) {
      save_path <- trimws(db_save_path_input() %||% "")
      if (nzchar(save_path) && file.exists(save_path)) {
        read_con <- tryCatch(mc_db_connect(save_path, read_only = TRUE),
                             error = function(e) NULL)
        own_con <- TRUE
      }
    }

    if (is.null(read_con)) {
      # Show in-memory summary
      raw <- raw_results_rv()
      n_cells <- length(raw)
      n_sims <- sum(vapply(raw, function(x) length(x$results), integer(1)))
      return(fluidRow(
        column(3, wellPanel(h5("Total Sims"), h3(format(n_sims, big.mark = ",")))),
        column(3, wellPanel(h5("Configs"), h3(n_cells))),
        column(3, wellPanel(h5("Source"), h3("In-Memory"))),
        column(3, wellPanel(h5("Batches"), h3("N/A")))
      ))
    }
    if (own_con) on.exit(DBI::dbDisconnect(read_con, shutdown = TRUE))

    total_sims <- tryCatch(
      DBI::dbGetQuery(read_con, "SELECT COUNT(*) AS cnt FROM simulations")$cnt,
      error = function(e) 0)
    n_configs <- tryCatch(
      DBI::dbGetQuery(read_con,
        "SELECT COUNT(DISTINCT (n, phi, innov_dist)) AS cnt FROM simulations")$cnt,
      error = function(e) 0)
    n_batches <- tryCatch(
      DBI::dbGetQuery(read_con, "SELECT COUNT(*) AS cnt FROM batches")$cnt,
      error = function(e) 0)
    date_range <- tryCatch(
      DBI::dbGetQuery(read_con,
        "SELECT MIN(created_at) AS mn, MAX(created_at) AS mx FROM batches"),
      error = function(e) data.frame(mn = NA, mx = NA))

    fluidRow(
      column(3, wellPanel(h5("Total Sims"), h3(format(total_sims, big.mark = ",")))),
      column(3, wellPanel(h5("Configs"), h3(n_configs))),
      column(3, wellPanel(h5("Batches"), h3(n_batches))),
      column(3, wellPanel(h5("Date Range"),
                           h6(if (is.na(date_range$mn)) "N/A"
                              else paste(substr(date_range$mn, 1, 10), "to",
                                         substr(date_range$mx, 1, 10)))))
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
