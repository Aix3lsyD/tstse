# Module: Innovation Comparison (Tab 9)
# Extracted from app.R -- compares how different innovation distributions
# affect AR(1) realizations using the same seed, phi, and n.

mod_innovation_comp_ui <- function(id) {
  ns <- NS(id)

  tabPanel(
    "Innovation Comparison",
    br(),
    h4("Innovation Distribution Comparison"),
    p(class = "text-body-secondary", style = "margin-bottom: 15px;",
      "Compare how different innovation distributions affect AR(1) realizations. ",
      "All plots use the same seed, phi, and n so differences are purely due to the innovation distribution."),

    fluidRow(
      # --- Controls sidebar ---
      column(3,
        wellPanel(
          h5("Parameters"),
          numericInput(ns("innov_phi"), "AR(1) Coefficient (phi):",
                       value = 0.95, min = 0.01, max = 0.999, step = 0.01),
          numericInput(ns("innov_n"), "Sample Size (n):",
                       value = 200, min = 10, max = 2000, step = 10),
          numericInput(ns("innov_seed"), "Seed:", value = 42, min = 1, step = 1),
          checkboxInput(ns("innov_compact_ref"), "Compact reference plot", value = FALSE),
          actionButton(ns("innov_generate"), "Generate",
                       icon = icon("play"), class = "btn-primary btn-block"),

          hr(),
          h5("Distribution Settings"),

          # Slot 1: Student's t
          tags$strong("Student's t"),
          numericInput(ns("innov_t_df"), "df:", value = 3, min = 1, step = 1),

          # Slot 2: ARCH
          tags$strong("ARCH"),
          textAreaInput(ns("innov_arch_alpha"), "alpha (comma-separated):",
                        value = "0.2, 0.175, 0.15, 0.125, 0.1, 0.075, 0.05, 0.025",
                        rows = 2),
          numericInput(ns("innov_arch_omega"), "omega:", value = 0.1, min = 0.001, step = 0.01),

          # Slot 3: Heteroscedastic
          tags$strong("Heteroscedastic"),
          numericInput(ns("innov_hetero_sd"), "sd:", value = 1, min = 0.01, step = 0.1),

          # Slot 4: Laplace
          tags$strong("Laplace"),
          numericInput(ns("innov_lap_scale"), "scale:", value = 0.707, min = 0.01, step = 0.1)
        )
      ),

      # --- Plots area ---
      column(9,
        # Reference normal plot: full-width or compact (half-width)
        conditionalPanel(
          condition = "!input.innov_compact_ref",
          ns = ns,
          plotOutput(ns("innov_ref_plot_full"), height = "350px")
        ),
        conditionalPanel(
          condition = "input.innov_compact_ref",
          ns = ns,
          fluidRow(
            column(6, plotOutput(ns("innov_ref_plot_compact"), height = "300px"))
          )
        ),

        # 2x2 comparison grid
        fluidRow(
          column(6, plotOutput(ns("innov_comp_plot_1"), height = "300px")),
          column(6, plotOutput(ns("innov_comp_plot_2"), height = "300px"))
        ),
        fluidRow(
          column(6, plotOutput(ns("innov_comp_plot_3"), height = "300px")),
          column(6, plotOutput(ns("innov_comp_plot_4"), height = "300px"))
        )
      )
    )
  )
}

mod_innovation_comp_server <- function(id, con, init_choices) {
  moduleServer(id, function(input, output, session) {

    # Reactive: generate all five realizations on button click
    innov_realizations <- eventReactive(input$innov_generate, {
      phi <- input$innov_phi
      n <- as.integer(input$innov_n)
      if (is.na(n) || n < 10) n <- 10L
      if (n > 2000) n <- 2000L
      seed <- input$innov_seed

      # Clamp phi to valid AR(1) range
      if (is.null(phi) || is.na(phi)) phi <- 0.95
      phi <- max(0.01, min(0.999, phi))

      colors <- c("#e41a1c", "#ff7f00", "#984ea3", "#4daf4a")

      withProgress(message = "Generating realizations...", value = 0, {
        # Reference: normal
        incProgress(1/5, detail = "Normal (reference)")
        ref_y <- gen_aruma_flex(n, phi = phi, innov_gen = make_gen_norm(sd = 1),
                                seed = seed, plot = FALSE)$y

        # Slot 1: Student's t
        incProgress(1/5, detail = "Student's t")
        t_df <- input$innov_t_df
        if (is.null(t_df) || is.na(t_df) || t_df < 1) t_df <- 3
        t_y <- gen_aruma_flex(n, phi = phi,
                              innov_gen = make_gen_t(df = t_df, scale = FALSE),
                              seed = seed, plot = FALSE)$y

        # Slot 2: ARCH
        incProgress(1/5, detail = "ARCH")
        arch_slot <- if (!requireNamespace("rugarch", quietly = TRUE)) {
          list(y = NULL, label = "ARCH", color = colors[2],
               error = "Install the 'rugarch' package to use ARCH innovations")
        } else {
          omega <- input$innov_arch_omega
          if (is.null(omega) || is.na(omega)) omega <- 0.1
          alpha_str <- input$innov_arch_alpha
          alpha <- tryCatch(
            as.numeric(trimws(strsplit(alpha_str, ",")[[1]])),
            warning = function(w) NULL,
            error = function(e) NULL
          )
          if (is.null(alpha) || any(is.na(alpha)) || length(alpha) == 0) {
            alpha <- c(0.2, 0.175, 0.15, 0.125, 0.1, 0.075, 0.05, 0.025)
          }
          arch_gen <- make_gen_garch(omega = omega, alpha = alpha)
          arch_y <- gen_aruma_flex(n, phi = phi, innov_gen = arch_gen,
                                   seed = seed, plot = FALSE)$y
          list(y = arch_y,
               label = sprintf("ARCH(%d)", length(alpha)),
               color = colors[2], error = NULL)
        }

        # Slot 3: Heteroscedastic
        incProgress(1/5, detail = "Heteroscedastic")
        h_sd <- input$innov_hetero_sd
        if (is.null(h_sd) || is.na(h_sd) || h_sd <= 0) h_sd <- 1
        hetero_y <- gen_aruma_flex(n, phi = phi,
                                   innov_gen = make_gen_hetero(sd = h_sd),
                                   seed = seed, plot = FALSE)$y

        # Slot 4: Laplace
        incProgress(1/5, detail = "Laplace")
        lap_sc <- input$innov_lap_scale
        if (is.null(lap_sc) || is.na(lap_sc) || lap_sc <= 0) lap_sc <- 1 / sqrt(2)
        laplace_y <- gen_aruma_flex(n, phi = phi,
                                    innov_gen = make_gen_laplace(scale = lap_sc),
                                    seed = seed, plot = FALSE)$y
      })

      slots <- list(
        list(y = t_y, label = sprintf("Student's t(df=%d)", t_df),
             color = colors[1], error = NULL),
        arch_slot,
        list(y = hetero_y, label = sprintf("Heteroscedastic(sd=%.2f)", h_sd),
             color = colors[3], error = NULL),
        list(y = laplace_y, label = sprintf("Laplace(scale=%.3f)", lap_sc),
             color = colors[4], error = NULL)
      )

      list(ref_y = ref_y, slots = slots,
           phi = phi, n = n, seed = seed)
    })

    # Helper to render the reference normal plot
    .render_innov_ref <- function() {
      op <- viewer_base_par()
      on.exit(par(op), add = TRUE)
      data <- tryCatch(innov_realizations(), error = function(e) NULL)

      if (is.null(data)) {
        plot.new()
        text(0.5, 0.5, "Click Generate to create realizations",
             cex = 1.4, col = "grey50")
        return()
      }

      y <- data$ref_y
      plot(y, type = "l", col = "steelblue", lwd = 1.2,
           xlab = "Time", ylab = "Value", ylim = range(y),
           main = sprintf("Normal (Reference)  |  phi=%.2f, n=%d, seed=%d",
                          data$phi, data$n, data$seed))
      abline(h = 0, lty = 3, col = "grey60")
    }

    # Full-width reference plot (default)
    output$innov_ref_plot_full <- renderPlot(bg = "transparent", { .render_innov_ref() })

    # Compact reference plot (same size as comparison plots)
    output$innov_ref_plot_compact <- renderPlot(bg = "transparent", { .render_innov_ref() })

    # Comparison plots (slots 1-4)
    local({
      colors <- c("#e41a1c", "#ff7f00", "#984ea3", "#4daf4a")
      for (idx in 1:4) {
        local({
          i <- idx
          output[[paste0("innov_comp_plot_", i)]] <- renderPlot(bg = "transparent", {
            op <- viewer_base_par()
            on.exit(par(op), add = TRUE)
            data <- tryCatch(innov_realizations(), error = function(e) NULL)

            if (is.null(data)) {
              plot.new()
              text(0.5, 0.5, "Click Generate to create realizations",
                   cex = 1.2, col = "grey50")
              return()
            }

            slot <- data$slots[[i]]

            if (!is.null(slot$error)) {
              plot.new()
              text(0.5, 0.5, slot$error, cex = 1.1, col = "grey50")
              return()
            }

            y <- slot$y
            plot(y, type = "l", col = slot$color, lwd = 1.2,
                 xlab = "Time", ylab = "Value", ylim = range(y),
                 main = slot$label)
            abline(h = 0, lty = 3, col = "grey60")
          })
        })
      }
    })

  })
}
