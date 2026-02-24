# Module: Database Overview (Tab 1)
# Extracted from app.R -- provides summary stats and coverage matrix

mod_overview_ui <- function(id) {
  ns <- NS(id)

  tabPanel(
    "Database Overview",
    br(),
    h4("Database Overview"),
    fluidRow(
      column(3,
        wellPanel(style = "text-align: center;",
          h6(class = "text-body-secondary", style = "margin-bottom: 5px;", "Total Simulations"),
          h3(style = "margin: 0;", textOutput(ns("overview_total_sims"), inline = TRUE))
        )
      ),
      column(3,
        wellPanel(style = "text-align: center;",
          h6(class = "text-body-secondary", style = "margin-bottom: 5px;", "Configurations"),
          h3(style = "margin: 0;", textOutput(ns("overview_n_configs"), inline = TRUE))
        )
      ),
      column(3,
        wellPanel(style = "text-align: center;",
          h6(class = "text-body-secondary", style = "margin-bottom: 5px;", "Batches"),
          h3(style = "margin: 0;", textOutput(ns("overview_n_batches"), inline = TRUE))
        )
      ),
      column(3,
        wellPanel(style = "text-align: center;",
          h6(class = "text-body-secondary", style = "margin-bottom: 5px;", "Date Range"),
          div(style = "font-size: 0.95em;",
            textOutput(ns("overview_date_range"), inline = TRUE)
          )
        )
      )
    ),
    hr(),
    h5("Coverage Matrix: Simulations per Configuration"),
    p(class = "text-body-secondary",
      "Each cell shows the number of simulations for that (n, phi) combination. ",
      "Faceted by innovation distribution. Empty/missing cells indicate gaps in coverage."
    ),
    plotOutput(ns("coverage_matrix"), height = "450px")
  )
}

mod_overview_server <- function(id, con) {
  moduleServer(id, function(input, output, session) {

    overview_stats <- reactive({
      tryCatch({
        dbGetQuery(con, "
          SELECT COUNT(*) as total_sims,
                 COUNT(DISTINCT batch_id) as n_batches,
                 MIN(created_at) as first_sim,
                 MAX(created_at) as last_sim
          FROM simulations
        ")
      }, error = function(e) data.frame(total_sims = 0, n_batches = 0,
                                         first_sim = NA, last_sim = NA))
    })

    overview_coverage <- reactive({
      tryCatch({
        dbGetQuery(con, "
          SELECT n, phi, innov_dist, COUNT(*) as n_sims
          FROM simulations
          GROUP BY n, phi, innov_dist
          ORDER BY innov_dist, n, phi
        ")
      }, error = function(e) data.frame())
    })

    output$overview_total_sims <- renderText({
      stats <- overview_stats()
      format(stats$total_sims[1], big.mark = ",")
    })

    output$overview_n_configs <- renderText({
      df <- overview_coverage()
      as.character(nrow(df))
    })

    output$overview_n_batches <- renderText({
      stats <- overview_stats()
      as.character(stats$n_batches[1])
    })

    output$overview_date_range <- renderText({
      stats <- overview_stats()
      if (is.na(stats$first_sim[1])) return("No data")
      paste0(format(stats$first_sim[1], "%Y-%m-%d"), " to ",
             format(stats$last_sim[1], "%Y-%m-%d"))
    })

    output$coverage_matrix <- renderPlot(bg = "transparent", {
      df <- overview_coverage()
      if (nrow(df) == 0) {
        plot.new()
        text(0.5, 0.5, "No data", cex = 1.2, col = "grey50")
        return()
      }

      df$n <- factor(df$n)
      df$phi <- factor(df$phi)
      n_mid <- stats::median(df$n_sims, na.rm = TRUE)
      df$label_col <- ifelse(df$n_sims >= n_mid, "#f8f9fa", "#212529")

      p <- ggplot(df, aes(x = phi, y = n, fill = n_sims)) +
        geom_tile(color = grDevices::adjustcolor(viewer_plot_fg(session), alpha.f = 0.20),
                  linewidth = 0.25) +
        geom_text(aes(label = format(n_sims, big.mark = ","), colour = label_col),
                  size = 4.2, fontface = "bold") +
        scale_colour_identity() +
        scale_fill_gradient(low = "#4a90c4", high = "#08519c",
                            name = "Simulations") +
        labs(x = expression(phi), y = "Sample Size (n)") +
        viewer_plot_theme(base_size = 14) +
        theme(
          panel.grid = element_blank(),
          legend.position = "right",
          strip.text = element_text(face = "bold", colour = viewer_plot_fg(session))
        ) +
        facet_wrap(~ innov_dist)

      print(p)
    })

  })
}
