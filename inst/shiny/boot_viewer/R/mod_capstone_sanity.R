# Top-level Tab: Sanity Checks
# Fast vs standard innovation parity + phi=0 generator diagnostics

.SANITY_INNOV_CHOICES <- c(
  "Normal", "Student's t", "Skew-t", "GED",
  "Laplace", "Uniform", "Mixture Normal",
  "GARCH", "Heteroscedastic"
)

.sanity_picker_input <- function(inputId, label, choices, selected = NULL, options = list()) {
  if (requireNamespace("shinyWidgets", quietly = TRUE)) {
    shinyWidgets::pickerInput(
      inputId = inputId,
      label = label,
      choices = choices,
      selected = selected,
      options = options
    )
  } else {
    selectInput(
      inputId = inputId,
      label = label,
      choices = choices,
      selected = selected
    )
  }
}

.sanity_toggle_input <- function(inputId, label, value = FALSE) {
  if (requireNamespace("shinyWidgets", quietly = TRUE)) {
    shinyWidgets::materialSwitch(
      inputId = inputId,
      label = label,
      value = value,
      status = "primary",
      right = FALSE
    )
  } else {
    checkboxInput(inputId = inputId, label = label, value = value)
  }
}

.sanity_action_button <- function(inputId, label, icon_name = NULL) {
  btn_icon <- if (is.null(icon_name)) NULL else icon(icon_name)
  if (requireNamespace("shinyWidgets", quietly = TRUE)) {
    shinyWidgets::actionBttn(
      inputId = inputId,
      label = label,
      icon = btn_icon,
      style = "material-flat",
      color = "primary",
      size = "sm"
    )
  } else {
    actionButton(inputId = inputId, label = label, icon = btn_icon, class = "btn-primary")
  }
}

mod_capstone_sanity_ui <- function(id) {
  ns <- NS(id)
  tabPanel("Sanity Checks",
    br(),
    h4("Sanity Checks"),
    p(class = "text-body-secondary",
      "Compare fast vs standard innovation generators and inspect phi=0 ARUMA output behavior."),

    div(class = "plot-controls",
      fluidRow(
        column(3,
          .sanity_picker_input(
            ns("san_dist"), "Innovation Distribution",
            choices = .SANITY_INNOV_CHOICES,
            selected = "Normal",
            options = list(`live-search` = FALSE, size = 10)
          )
        ),
        column(2,
          numericInput(ns("san_n"), "Sample Size (n)", value = 250, min = 200, max = 100000, step = 50)
        ),
        column(2,
          numericInput(ns("san_reps"), "Replications", value = 3, min = 1, max = 50, step = 1)
        ),
        column(2,
          numericInput(ns("san_seed"), "Seed", value = 2026, min = 1, step = 1)
        ),
        column(3,
          div(style = "margin-top: 30px;",
            p(class = "text-body-secondary", style = "margin-bottom: 0;",
              "Relevant metrics are shown by distribution type.")
          )
        )
      ),
      div(style = "margin-top: 6px;", uiOutput(ns("san_dist_params_ui")))
    ),

    div(id = ns("sanity_content"),
      tabsetPanel(
        id = ns("san_subtabs"),

        tabPanel("Fast vs Standard",
          br(),
          .sanity_action_button(ns("san_run_equiv"), "Run Equivalence Check", "flask"),
          br(), br(),
          DT::dataTableOutput(ns("san_equiv_table")),
          br(),
          plotOutput(ns("san_equiv_plot"), height = "420px")
        ),

        tabPanel("ARUMA phi=0",
          br(),
          .sanity_action_button(ns("san_run_phi0"), "Run phi=0 Diagnostic", "chart-line"),
          br(), br(),
          DT::dataTableOutput(ns("san_phi0_table")),
          br(),
          plotOutput(ns("san_phi0_plot"), height = "700px")
        )
      )
    )
  )
}

mod_capstone_sanity_server <- function(id) {
  moduleServer(id, function(input, output, session) {
  ns <- session$ns
  has_waiter <- requireNamespace("waiter", quietly = TRUE)
  has_shinyvalidate <- requireNamespace("shinyvalidate", quietly = TRUE)

  .with_waiter <- function(expr) {
    if (!has_waiter) return(force(expr))
    w <- waiter::Waiter$new(
      id = ns("sanity_content"),
      html = waiter::spin_fading_circles(),
      color = "rgba(0, 0, 0, 0.25)"
    )
    w$show()
    on.exit(w$hide(), add = TRUE)
    force(expr)
  }

  .sample_stats <- function(x) {
    m <- mean(x)
    s <- stats::sd(x)
    c(
      mean = m,
      sd = s,
      skew = if (s > 0) mean(((x - m) / s)^3) else NA_real_,
      kurt = if (s > 0) mean(((x - m) / s)^4) else NA_real_
    )
  }

  output$san_dist_params_ui <- renderUI({
    dist <- input$san_dist %||% "Normal"

    switch(dist,
      "Normal" = fluidRow(
        column(3, numericInput(ns("san_norm_sd"), "Normal SD", value = 1, min = 0.001, step = 0.1))
      ),

      "Student's t" = fluidRow(
        column(3, numericInput(ns("san_t_df"), "t df", value = 5, min = 1, step = 1)),
        column(3, div(style = "margin-top: 24px;", .sanity_toggle_input(ns("san_t_scale"), "Scale to unit variance", value = TRUE)))
      ),

      "Skew-t" = fluidRow(
        column(3, numericInput(ns("san_skt_df"), "Skew-t df", value = 6, min = 1, step = 1)),
        column(3, numericInput(ns("san_skt_alpha"), "Skew alpha", value = 1, step = 0.1)),
        column(3, div(style = "margin-top: 24px;", .sanity_toggle_input(ns("san_skt_scale"), "Scale to unit variance", value = TRUE)))
      ),

      "GED" = fluidRow(
        column(3, numericInput(ns("san_ged_nu"), "GED nu", value = 2, min = 0.05, step = 0.1)),
        column(3, numericInput(ns("san_ged_sd"), "GED SD", value = 1, min = 0.001, step = 0.1))
      ),

      "Laplace" = fluidRow(
        column(3, numericInput(ns("san_lap_scale"), "Laplace scale", value = 1, min = 0.001, step = 0.05))
      ),

      "Uniform" = fluidRow(
        column(3, numericInput(ns("san_unif_hw"), "Uniform half-width", value = sqrt(3), min = 0.001, step = 0.1))
      ),

      "Mixture Normal" = fluidRow(
        column(3, numericInput(ns("san_mix_sd1"), "Mixture sd1", value = 1, min = 0.001, step = 0.1)),
        column(3, numericInput(ns("san_mix_sd2"), "Mixture sd2", value = 3, min = 0.001, step = 0.1)),
        column(3, numericInput(ns("san_mix_prob1"), "Mixture prob1", value = 0.9, min = 0.01, max = 0.99, step = 0.01))
      ),

      "GARCH" = fluidRow(
        column(3, numericInput(ns("san_garch_omega"), "GARCH omega", value = 0.1, min = 0.001, step = 0.01)),
        column(3, textInput(ns("san_garch_alpha"), "GARCH alpha list (comma-separated)", value = "0.15",
                            placeholder = "e.g., 0.15, 0.05")),
        column(3, textInput(ns("san_garch_beta"), "GARCH beta list (comma-separated)", value = "0.8",
                            placeholder = "e.g., 0.8, 0.1")),
        column(3, div(style = "margin-top: 30px;",
                      p(class = "text-body-secondary", style = "margin-bottom: 0;",
                        "Length of alpha and beta lists sets GARCH(q, p) order.")))
      ),

      "Heteroscedastic" = fluidRow(
        column(3,
          .sanity_picker_input(
            ns("san_hetero_shape"), "Hetero shape",
            choices = c("linear", "sqrt", "log", "exp", "power", "step", "periodic"),
            selected = "linear"
          )
        ),
        column(3, numericInput(ns("san_hetero_from"), "From", value = 1, min = 0.01, step = 0.1)),
        column(3, numericInput(ns("san_hetero_to"), "To", value = 10, min = 0.01, step = 0.1)),
        column(3, numericInput(ns("san_hetero_sd"), "Hetero SD", value = 1, min = 0.001, step = 0.1))
      ),

      fluidRow(column(12, p(class = "text-body-secondary", "No additional parameters.")))
    )
  })

  .params_from_input <- function(input) {
    list(
      norm_sd = input$san_norm_sd,
      t_df = input$san_t_df,
      t_scale = isTRUE(input$san_t_scale),
      skt_df = input$san_skt_df,
      skt_alpha = input$san_skt_alpha,
      skt_scale = isTRUE(input$san_skt_scale),
      ged_nu = input$san_ged_nu,
      ged_sd = input$san_ged_sd,
      lap_scale = input$san_lap_scale,
      unif_hw = input$san_unif_hw,
      mix_sd1 = input$san_mix_sd1,
      mix_sd2 = input$san_mix_sd2,
      mix_prob1 = input$san_mix_prob1,
      garch_omega = input$san_garch_omega,
      garch_alpha = input$san_garch_alpha,
      garch_beta = input$san_garch_beta,
      hetero_shape = input$san_hetero_shape,
      hetero_from = input$san_hetero_from,
      hetero_to = input$san_hetero_to,
      hetero_sd = input$san_hetero_sd
    )
  }

  .validate_params <- function(dist, p, n_val, reps_val, seed_val) {
    errs <- character(0)

    if (is.na(n_val) || n_val < 200) errs <- c(errs, "Sample size (n) must be at least 200.")
    if (is.na(reps_val) || reps_val < 1) errs <- c(errs, "Replications must be at least 1.")
    if (is.na(seed_val) || seed_val < 1) errs <- c(errs, "Seed must be a positive integer.")

    if (dist == "Normal") {
      if (is.na(p$norm_sd) || p$norm_sd <= 0) errs <- c(errs, "Normal SD must be positive.")
    } else if (dist == "Student's t") {
      if (is.na(p$t_df) || p$t_df <= 0) errs <- c(errs, "t df must be positive.")
    } else if (dist == "Skew-t") {
      if (is.na(p$skt_df) || p$skt_df <= 0) errs <- c(errs, "Skew-t df must be positive.")
    } else if (dist == "GED") {
      if (is.na(p$ged_nu) || p$ged_nu <= 0) errs <- c(errs, "GED nu must be positive.")
      if (is.na(p$ged_sd) || p$ged_sd <= 0) errs <- c(errs, "GED SD must be positive.")
    } else if (dist == "Laplace") {
      if (is.na(p$lap_scale) || p$lap_scale <= 0) errs <- c(errs, "Laplace scale must be positive.")
    } else if (dist == "Uniform") {
      if (is.na(p$unif_hw) || p$unif_hw <= 0) errs <- c(errs, "Uniform half-width must be positive.")
    } else if (dist == "Mixture Normal") {
      if (is.na(p$mix_sd1) || p$mix_sd1 <= 0 || is.na(p$mix_sd2) || p$mix_sd2 <= 0) {
        errs <- c(errs, "Mixture sd1 and sd2 must be positive.")
      }
      if (is.na(p$mix_prob1) || p$mix_prob1 <= 0 || p$mix_prob1 >= 1) {
        errs <- c(errs, "Mixture prob1 must be strictly between 0 and 1.")
      }
    } else if (dist == "GARCH") {
      alpha <- suppressWarnings(as.numeric(strsplit(p$garch_alpha %||% "", "\\s*,\\s*")[[1]]))
      beta <- suppressWarnings(as.numeric(strsplit(p$garch_beta %||% "", "\\s*,\\s*")[[1]]))
      if (is.na(p$garch_omega) || p$garch_omega <= 0) errs <- c(errs, "GARCH omega must be positive.")
      if (length(alpha) == 0 || any(is.na(alpha)) || any(alpha < 0)) errs <- c(errs, "GARCH alpha must be non-negative numeric values.")
      if (length(beta) == 0 || any(is.na(beta)) || any(beta < 0)) errs <- c(errs, "GARCH beta must be non-negative numeric values.")
    } else if (dist == "Heteroscedastic") {
      if (is.na(p$hetero_sd) || p$hetero_sd <= 0) errs <- c(errs, "Hetero SD must be positive.")
      if (is.na(p$hetero_from) || p$hetero_from <= 0) errs <- c(errs, "Hetero 'From' must be positive.")
      if (is.na(p$hetero_to) || p$hetero_to <= 0) errs <- c(errs, "Hetero 'To' must be positive.")
      if (!is.na(p$hetero_from) && !is.na(p$hetero_to) && p$hetero_to < p$hetero_from) {
        errs <- c(errs, "Hetero 'To' should be greater than or equal to 'From'.")
      }
    }

    errs
  }

  if (has_shinyvalidate) {
    sv <- shinyvalidate::InputValidator$new()
    sv$add_rule("san_n", function(value) {
      v <- suppressWarnings(as.numeric(value)); if (is.na(v) || v < 200) "Must be >= 200" else NULL
    })
    sv$add_rule("san_reps", function(value) {
      v <- suppressWarnings(as.numeric(value)); if (is.na(v) || v < 1) "Must be >= 1" else NULL
    })
    sv$add_rule("san_seed", function(value) {
      v <- suppressWarnings(as.numeric(value)); if (is.na(v) || v < 1) "Must be positive" else NULL
    })
    sv$enable()
  }

  .theo_params_from_input <- function(dist, p) {
    if (dist == "Normal") {
      list(sd = p$norm_sd %||% 1)
    } else if (dist == "Student's t") {
      list(df = p$t_df %||% 3, scale = isTRUE(p$t_scale))
    } else if (dist == "Skew-t") {
      list(df = p$skt_df %||% 5, alpha = p$skt_alpha %||% 0,
           scale = isTRUE(p$skt_scale))
    } else if (dist == "GED") {
      list(nu = p$ged_nu %||% 2, sd = p$ged_sd %||% 1)
    } else if (dist == "Laplace") {
      list(scale = p$lap_scale %||% 1)
    } else if (dist == "Uniform") {
      list(half_width = p$unif_hw %||% sqrt(3))
    } else if (dist == "Mixture Normal") {
      list(sd1 = p$mix_sd1 %||% 1, sd2 = p$mix_sd2 %||% 3,
           prob1 = p$mix_prob1 %||% 0.9)
    } else {
      list()
    }
  }

  .acf1 <- function(x) {
    out <- tryCatch(stats::acf(x, lag.max = 1, plot = FALSE)$acf[2],
                    error = function(e) NA_real_)
    as.numeric(out)
  }

  .equiv_metrics <- function(x, dist) {
    base <- .sample_stats(x)
    if (dist %in% c("GARCH", "Heteroscedastic")) {
      c(
        mean = base[["mean"]],
        sd = base[["sd"]],
        acf1 = .acf1(x),
        acf1_sq = .acf1(x^2)
      )
    } else {
      c(
        mean = base[["mean"]],
        sd = base[["sd"]],
        skew = base[["skew"]],
        kurt = base[["kurt"]]
      )
    }
  }

  equiv_data <- eventReactive(input$san_run_equiv, {
    .with_waiter({
      n <- as.integer(input$san_n %||% 250L)
      reps <- as.integer(input$san_reps %||% 3L)
      seed <- as.integer(input$san_seed %||% 2026L)
      dist <- input$san_dist
      params <- .params_from_input(input)

      errs <- .validate_params(dist, params, n, reps, seed)
      validate(need(length(errs) == 0, paste(errs, collapse = " ")))

      gen_std <- tryCatch(build_innov_gen(dist, params, use_fast = FALSE), error = function(e) NULL)
      gen_fast <- tryCatch(build_innov_gen(dist, params, use_fast = TRUE), error = function(e) NULL)

      validate(
        need(!is.null(gen_std), "Could not build standard innovation generator."),
        need(!is.null(gen_fast), "Could not build fast innovation generator.")
      )

      set.seed(seed)
      x_std <- numeric(0)
      x_fast <- numeric(0)
      for (i in seq_len(reps)) {
        x_std <- c(x_std, gen_std(n))
        x_fast <- c(x_fast, gen_fast(n))
      }

      s_std <- .equiv_metrics(x_std, dist)
      s_fast <- .equiv_metrics(x_fast, dist)
      delta <- abs(s_fast - s_std)

      ks_p <- if (dist %in% c("GARCH", "Heteroscedastic")) NA_real_ else {
        tryCatch(stats::ks.test(x_fast, x_std)$p.value, error = function(e) NA_real_)
      }

      list(
        x_std = x_std,
        x_fast = x_fast,
        table = data.frame(
          metric = names(s_std),
          standard = as.numeric(s_std),
          fast = as.numeric(s_fast),
          abs_delta = as.numeric(delta),
          stringsAsFactors = FALSE
        ),
        ks_p = ks_p,
        dist = dist,
        theo_params = .theo_params_from_input(dist, params)
      )
    })
  })

  output$san_equiv_table <- DT::renderDataTable({
    d <- equiv_data()
    tbl <- d$table
    if (!is.na(d$ks_p)) {
      footer <- data.frame(
        metric = "KS p-value (fast vs standard)",
        standard = NA_real_, fast = NA_real_, abs_delta = d$ks_p,
        stringsAsFactors = FALSE
      )
      tbl <- rbind(tbl, footer)
    }

    DT::datatable(tbl, rownames = FALSE, options = list(pageLength = 10, dom = "t"))
  })

  output$san_equiv_plot <- renderPlot(bg = "transparent", {
    d <- equiv_data()

    dens_std <- stats::density(d$x_std)
    dens_fast <- stats::density(d$x_fast)
    plot(dens_std, main = sprintf("Fast vs Standard Density Overlay (%s)", d$dist),
         xlab = "Innovation value", ylab = "Density",
         col = "#1f77b4", lwd = 2)
    lines(dens_fast, col = "#d62728", lwd = 2, lty = 2)

    lgd <- c("Standard", "Fast")
    lgd_col <- c("#1f77b4", "#d62728")
    lgd_lty <- c(1, 2)

    theo <- innov_theoretical_density(d$dist, d$theo_params)
    if (!is.null(theo)) {
      xlim <- range(c(dens_std$x, dens_fast$x))
      x_seq <- seq(xlim[1], xlim[2], length.out = 600)
      y_theo <- tryCatch(theo$dfun(x_seq), error = function(e) NULL)
      if (!is.null(y_theo) && all(is.finite(y_theo))) {
        lines(x_seq, y_theo, col = "#2ca02c", lwd = 2, lty = 3)
        lgd <- c(lgd, "Theoretical")
        lgd_col <- c(lgd_col, "#2ca02c")
        lgd_lty <- c(lgd_lty, 3)
      }
    }

    legend("topright", legend = lgd,
           col = lgd_col, lwd = 2, lty = lgd_lty, bty = "n")
  })

  phi0_data <- eventReactive(input$san_run_phi0, {
    .with_waiter({
      n <- as.integer(input$san_n %||% 250L)
      seed <- as.integer(input$san_seed %||% 2026L)
      dist <- input$san_dist
      params <- .params_from_input(input)

      errs <- .validate_params(dist, params, n, 1L, seed)
      validate(need(length(errs) == 0, paste(errs, collapse = " ")))

      gen_std <- tryCatch(build_innov_gen(dist, params, use_fast = FALSE), error = function(e) NULL)
      gen_fast <- tryCatch(build_innov_gen(dist, params, use_fast = TRUE), error = function(e) NULL)

      validate(
        need(!is.null(gen_std), "Could not build standard innovation generator."),
        need(!is.null(gen_fast), "Could not build fast innovation generator.")
      )

      set.seed(seed)
      y_std <- gen_aruma_flex(n = n, phi = 0, innov_gen = gen_std, plot = FALSE)$y
      set.seed(seed)
      y_fast <- gen_aruma_flex(n = n, phi = 0, innov_gen = gen_fast, plot = FALSE)$y

      data.frame(
        series = c("ARUMA phi=0 standard", "ARUMA phi=0 fast"),
        mean = c(mean(y_std), mean(y_fast)),
        sd = c(stats::sd(y_std), stats::sd(y_fast)),
        lag1_acf = c(stats::acf(y_std, plot = FALSE, lag.max = 1)$acf[2],
                     stats::acf(y_fast, plot = FALSE, lag.max = 1)$acf[2]),
        stringsAsFactors = FALSE
      ) -> summary_tbl

      list(
        y_std = y_std,
        y_fast = y_fast,
        table = summary_tbl,
        dist = dist,
        theo_params = .theo_params_from_input(dist, params)
      )
    })
  })

  output$san_phi0_table <- DT::renderDataTable({
    DT::datatable(phi0_data()$table, rownames = FALSE, options = list(dom = "t"))
  })

  output$san_phi0_plot <- renderPlot(bg = "transparent", {
    d <- phi0_data()
    y_std <- as.numeric(d$y_std)
    y_fast <- as.numeric(d$y_fast)
    n_plot <- min(length(y_std), length(y_fast))
    validate(need(n_plot > 1, "Not enough data to plot."))
    y_std <- y_std[seq_len(n_plot)]
    y_fast <- y_fast[seq_len(n_plot)]
    fg <- viewer_plot_fg()

    op <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(op), add = TRUE)
    graphics::layout(matrix(c(1, 2, 3, 4, 5, 5), nrow = 3, byrow = TRUE),
                     heights = c(1, 1, 1.1))
    graphics::par(mar = c(4, 4, 2.5, 1.5))

    idx <- seq_len(n_plot)
    yr <- range(c(y_std, y_fast), na.rm = TRUE)

    graphics::plot(idx, y_std, type = "l",
                   col = "#1f77b4", lwd = 1,
                   main = sprintf("ARUMA phi=0: Standard (%s)", d$dist),
                   ylab = "Value", xlab = "Index", ylim = yr,
                   col.axis = fg, col.lab = fg, col.main = fg)

    graphics::plot(idx, y_fast, type = "l",
                   col = "#d62728", lwd = 1,
                   main = sprintf("ARUMA phi=0: Fast (%s)", d$dist),
                   ylab = "Value", xlab = "Index", ylim = yr,
                   col.axis = fg, col.lab = fg, col.main = fg)

    acf_std <- stats::acf(y_std, plot = FALSE, lag.max = 25)
    graphics::plot(acf_std$lag[-1], acf_std$acf[-1], type = "h", lwd = 2,
                   col = "#1f77b4", xlab = "Lag", ylab = "ACF",
                   main = "ACF (Standard)",
                   col.axis = fg, col.lab = fg, col.main = fg)
    graphics::abline(h = 0, col = fg)

    acf_fast <- stats::acf(y_fast, plot = FALSE, lag.max = 25)
    graphics::plot(acf_fast$lag[-1], acf_fast$acf[-1], type = "h", lwd = 2,
                   col = "#d62728", xlab = "Lag", ylab = "ACF",
                   main = "ACF (Fast)",
                   col.axis = fg, col.lab = fg, col.main = fg)
    graphics::abline(h = 0, col = fg)

    h1 <- graphics::hist(y_std, breaks = 40, plot = FALSE)
    h2 <- graphics::hist(y_fast, breaks = 40, plot = FALSE)
    ymax <- max(h1$density, h2$density)
    graphics::hist(y_std, breaks = 40,
                   col = grDevices::adjustcolor("#1f77b4", alpha.f = 0.4),
                   border = NA, freq = FALSE,
                   main = "Density Overlay", xlab = "Value", ylim = c(0, ymax * 1.1),
                   col.axis = fg, col.lab = fg, col.main = fg)
    graphics::hist(y_fast, breaks = 40,
                   col = grDevices::adjustcolor("#d62728", alpha.f = 0.35),
                   border = NA, freq = FALSE, add = TRUE)

    theo <- innov_theoretical_density(d$dist, d$theo_params)
    if (!is.null(theo)) {
      xr <- graphics::par("usr")[1:2]
      x_seq <- seq(xr[1], xr[2], length.out = 600)
      y_theo <- tryCatch(theo$dfun(x_seq), error = function(e) NULL)
      if (!is.null(y_theo) && all(is.finite(y_theo))) {
        graphics::lines(x_seq, y_theo, col = "#2ca02c", lwd = 2, lty = 3)
      }
    }
    graphics::legend("topright", legend = c("Standard", "Fast", "Theoretical"),
                     pch = c(15, 15, NA),
                     pt.cex = c(1.4, 1.4, NA),
                     col = c("#1f77b4", "#d62728", "#2ca02c"),
                     lty = c(NA, NA, 3), lwd = c(NA, NA, 2), bty = "n",
                     text.col = fg)
  })
  })
}
