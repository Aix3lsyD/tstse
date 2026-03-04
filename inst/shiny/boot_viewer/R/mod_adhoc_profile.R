# Ad-Hoc Sub-Tab: Performance Profile
# Times each component of the simulation pipeline to identify bottlenecks
# Includes parent-child breakdown of the bootstrap kernel sub-components

mod_adhoc_profile_ui <- function(ns) {
  tabPanel("Performance Profile",
    br(),
    h4("Performance Profiling"),
    p(class = "text-body-secondary",
      "Benchmarks each component of a single simulation iteration. ",
      "Identifies bottlenecks in series generation, null model fitting, ",
      "and bootstrap kernel sub-operations."),
    div(class = "plot-controls",
      fluidRow(
        column(3,
          checkboxInput(ns("prof_use_current"), "Use Current Ad-Hoc Settings",
                        value = TRUE)
        ),
        column(2,
          conditionalPanel(
            condition = sprintf("!input['%s']", ns("prof_use_current")),
            numericInput(ns("prof_n"), "n", value = 200, min = 20, max = 2000)
          )
        ),
        column(2,
          conditionalPanel(
            condition = sprintf("!input['%s']", ns("prof_use_current")),
            numericInput(ns("prof_nb"), "nb", value = 399, min = 49, max = 4999)
          )
        ),
        column(2,
          conditionalPanel(
            condition = sprintf("!input['%s']", ns("prof_use_current")),
            numericInput(ns("prof_maxp"), "maxp", value = 5, min = 1, max = 20)
          )
        ),
        column(3,
          div(style = "margin-top: 25px;",
            actionButton(ns("prof_run"), "Run Profile",
                         icon = icon("stopwatch"), class = "btn-sm btn-info"))
        )
      )
    ),
    wellPanel(
      h5("Component Breakdown"),
      plotOutput(ns("prof_bar"), height = "420px")
    ),
    wellPanel(
      h5("Timing Details"),
      DTOutput(ns("prof_table"))
    ),
    wellPanel(
      h5("Scaling Analysis"),
      p(class = "text-body-secondary",
        "How time scales with nb, n, or both. ",
        "The n x nb grid can be slow for large ranges."),
      fluidRow(
        column(3,
          selectInput(ns("prof_scale_mode"), "Scaling Dimension:",
                      choices = c("nb (bootstrap replicates)" = "nb",
                                  "n (sample size)" = "n",
                                  "n x nb (grid)" = "nxnb"))
        ),
        column(3,
          div(style = "margin-top: 25px;",
            actionButton(ns("prof_run_scaling"), "Run Scaling",
                         icon = icon("chart-line"), class = "btn-sm btn-outline-info"))
        ),
        column(6,
          conditionalPanel(
            condition = sprintf("input['%s'] == 'nxnb'", ns("prof_scale_mode")),
            div(class = "alert alert-warning", style = "margin-top: 10px; padding: 8px;",
              "Grid mode runs n x nb combinations (30 points). May take a minute.")
          )
        )
      ),
      plotOutput(ns("prof_scaling"), height = "350px")
    )
  )
}

mod_adhoc_profile_server <- function(input, output, session, adhoc_params) {

  ns <- session$ns

  profile_results_rv <- reactiveVal(NULL)
  scaling_results_rv <- reactiveVal(NULL)

  # Resolve params: either from ad-hoc inputs or local overrides
  .get_profile_params <- function() {
    if (isTRUE(input$prof_use_current)) {
      adhoc_params()
    } else {
      list(
        n         = as.integer(input$prof_n %||% 200),
        phi       = adhoc_params()$phi,
        innov_gen = adhoc_params()$innov_gen,
        nb        = as.integer(input$prof_nb %||% 399),
        maxp      = as.integer(input$prof_maxp %||% 5),
        criterion = adhoc_params()$criterion,
        bootadj   = adhoc_params()$bootadj,
        min_p     = adhoc_params()$min_p
      )
    }
  }

  # Profile a single iteration, timing each component + kernel sub-components
  .profile_once <- function(params) {
    n_val     <- params$n
    phi_val   <- params$phi
    innov_gen <- params$innov_gen
    nb_val    <- params$nb
    maxp_val  <- params$maxp
    criterion <- params$criterion
    min_p     <- params$min_p
    bootadj   <- params$bootadj

    nreps <- 5L  # repetitions for stable timing

    # 1) Series generation
    t_gen <- system.time(for (r in seq_len(nreps)) {
      x <- gen_aruma_flex(n_val, phi = phi_val, innov_gen = innov_gen,
                          seed = NULL, plot = FALSE)$y
    })[["elapsed"]] / nreps

    # Use last generated series for remaining components
    # 2) Observed t-statistic
    t_tstat <- system.time(for (r in seq_len(nreps)) {
      tco <- co_tstat_pure_cpp(x, maxp_val, criterion, as.integer(min_p))
    })[["elapsed"]] / nreps

    # 3) Null AR model fit
    t_burg <- system.time(for (r in seq_len(nreps)) {
      fit <- burg_aic_select_cpp(x, maxp_val, criterion, as.integer(min_p))
    })[["elapsed"]] / nreps

    phi_null <- as.numeric(fit$phi)
    vara_null <- fit$vara
    if (fit$p == 0) phi_null <- numeric(0)

    # 4) Bootstrap kernel (parallel, actual wall time)
    num_threads <- RcppParallel::defaultNumThreads()
    grain_size <- max(16L, as.integer(nb_val / (4L * num_threads)))
    boot_seeds <- as.numeric(sample.int(.Machine$integer.max, nb_val))

    if (isTRUE(bootadj)) {
      t_boot <- system.time(for (r in seq_len(nreps)) {
        wbg_bootstrap_coba_kernel_grain_cpp(
          n = n_val, phi = phi_null, vara = vara_null,
          seeds = boot_seeds, maxp = maxp_val, criterion = criterion,
          grain_size = grain_size, min_p = as.integer(min_p))
      })[["elapsed"]] / nreps
      boot_label <- "Bootstrap (COBA)"
    } else {
      t_boot <- system.time(for (r in seq_len(nreps)) {
        wbg_bootstrap_kernel_grain_cpp(
          n = n_val, phi = phi_null, vara = vara_null,
          seeds = boot_seeds, maxp = maxp_val, criterion = criterion,
          grain_size = grain_size, min_p = as.integer(min_p))
      })[["elapsed"]] / nreps
      boot_label <- "Bootstrap Kernel"
    }

    # 5) Full wbg_boot_fast (for overhead calculation)
    t_full <- system.time(for (r in seq_len(nreps)) {
      wbg_boot_fast(x, nb = nb_val, maxp = maxp_val, criterion = criterion,
                    bootadj = bootadj, min_p = min_p)
    })[["elapsed"]] / nreps

    components <- t_gen + t_tstat + t_burg + t_boot
    t_overhead <- max(0, t_full - components)

    # 6) Kernel sub-component profiling (sequential, for ratio extraction)
    n_profile <- min(50L, nb_val)
    profile_seeds <- boot_seeds[seq_len(n_profile)]
    kp <- wbg_profile_kernel_components_cpp(
      n = n_val, phi = phi_null, vara = vara_null,
      seeds = profile_seeds, maxp = maxp_val, criterion = criterion,
      min_p = as.integer(min_p), coba = isTRUE(bootadj)
    )

    # Compute sub-component ratios from sequential profiling
    sub_us <- c(kp$gen_ar_us, kp$ols_us, kp$burg_us, kp$co_fused_us)
    if (isTRUE(bootadj)) sub_us <- c(sub_us, kp$burg_coba_us)
    total_sub_us <- sum(sub_us)

    # Scale ratios to actual parallel kernel time
    if (total_sub_us > 0) {
      ratios <- sub_us / total_sub_us
    } else {
      ratios <- rep(1 / length(sub_us), length(sub_us))
    }
    kernel_ms <- t_boot * 1000
    sub_ms <- ratios * kernel_ms
    kernel_overhead_ms <- max(0, kernel_ms - sum(sub_ms))

    # Build sub-component labels and values
    sub_labels <- c("  AR Generation", "  OLS Detrend",
                    "  Burg AR Fit", "  CO Transform")
    if (isTRUE(bootadj)) sub_labels <- c(sub_labels, "  Burg (COBA)")

    # Build combined data.frame: top-level + kernel children
    top_components <- c("gen_aruma_flex", "co_tstat_pure_cpp",
                        "burg_aic_select_cpp", boot_label)
    top_ms <- c(t_gen, t_tstat, t_burg, t_boot) * 1000

    # Insert children after the kernel parent
    all_labels <- c(top_components, sub_labels, "  Kernel Overhead",
                    "R Overhead", "TOTAL (wbg_boot_fast)")
    all_ms <- c(top_ms, sub_ms, kernel_overhead_ms,
                t_overhead * 1000, t_full * 1000)
    all_pct <- all_ms / (t_full * 1000) * 100

    # Mark which rows are children (for chart/table formatting)
    n_subs <- length(sub_labels) + 1L  # +1 for Kernel Overhead
    is_child <- c(rep(FALSE, length(top_components)),
                  rep(TRUE, n_subs),
                  FALSE, FALSE)

    data.frame(
      Component = all_labels,
      Time_ms   = all_ms,
      Pct       = all_pct,
      Reps      = nreps,
      IsChild   = is_child,
      stringsAsFactors = FALSE
    )
  }

  # Time a single wbg_boot_fast call with given n and nb
  .time_boot <- function(params, n_val, nb_val, nreps = 3L) {
    x <- gen_aruma_flex(n_val, phi = params$phi,
                        innov_gen = params$innov_gen,
                        seed = NULL, plot = FALSE)$y
    elapsed <- system.time(for (r in seq_len(nreps)) {
      wbg_boot_fast(x, nb = nb_val, maxp = params$maxp,
                    criterion = params$criterion,
                    bootadj = params$bootadj, min_p = params$min_p)
    })[["elapsed"]] / nreps
    elapsed * 1000  # return ms
  }

  # Run profile
  observeEvent(input$prof_run, {
    params <- .get_profile_params()
    if (is.null(params$phi)) {
      showNotification("Set phi in the ad-hoc parameters first.", type = "warning")
      return()
    }

    withProgress(message = "Profiling components...", value = 0.5, {
      result <- tryCatch(.profile_once(params),
                         error = function(e) {
                           showNotification(paste("Profile error:", e$message),
                                            type = "error")
                           NULL
                         })
    })
    if (!is.null(result)) {
      profile_results_rv(result)
      showNotification("Profile complete", type = "message", duration = 4)
    }
  })

  # Run scaling analysis
  observeEvent(input$prof_run_scaling, {
    params <- .get_profile_params()
    if (is.null(params$phi)) {
      showNotification("Set phi in the ad-hoc parameters first.", type = "warning")
      return()
    }

    mode <- input$prof_scale_mode %||% "nb"

    if (mode == "nb") {
      # --- Scale nb, hold n fixed ---
      nb_values <- c(49L, 99L, 199L, 399L, 799L, 1599L)
      n_fixed <- params$n

      withProgress(message = "Scaling analysis (nb)...", value = 0, {
        rows <- vector("list", length(nb_values))
        for (j in seq_along(nb_values)) {
          incProgress(1 / length(nb_values),
                      detail = sprintf("nb = %d", nb_values[j]))
          rows[[j]] <- data.frame(
            n = n_fixed, nb = nb_values[j],
            time_ms = .time_boot(params, n_fixed, nb_values[j]),
            stringsAsFactors = FALSE)
        }
      })
      result <- do.call(rbind, rows)
      result$.mode <- "nb"

    } else if (mode == "n") {
      # --- Scale n, hold nb fixed ---
      n_values <- c(50L, 100L, 200L, 500L, 1000L, 2000L)
      nb_fixed <- params$nb

      withProgress(message = "Scaling analysis (n)...", value = 0, {
        rows <- vector("list", length(n_values))
        for (j in seq_along(n_values)) {
          incProgress(1 / length(n_values),
                      detail = sprintf("n = %d", n_values[j]))
          rows[[j]] <- data.frame(
            n = n_values[j], nb = nb_fixed,
            time_ms = .time_boot(params, n_values[j], nb_fixed),
            stringsAsFactors = FALSE)
        }
      })
      result <- do.call(rbind, rows)
      result$.mode <- "n"

    } else {
      # --- n x nb grid ---
      n_values  <- c(50L, 100L, 250L, 500L, 1000L)
      nb_values <- c(99L, 199L, 399L, 799L, 1599L)
      # Reduce reps to 2 for the grid (30 combos)
      combos <- expand.grid(n = n_values, nb = nb_values,
                            stringsAsFactors = FALSE)
      n_combos <- nrow(combos)

      withProgress(message = "Scaling analysis (n x nb grid)...", value = 0, {
        rows <- vector("list", n_combos)
        for (j in seq_len(n_combos)) {
          incProgress(1 / n_combos,
                      detail = sprintf("n=%d, nb=%d",
                                       combos$n[j], combos$nb[j]))
          rows[[j]] <- data.frame(
            n = combos$n[j], nb = combos$nb[j],
            time_ms = .time_boot(params, combos$n[j], combos$nb[j], nreps = 2L),
            stringsAsFactors = FALSE)
        }
      })
      result <- do.call(rbind, rows)
      result$.mode <- "nxnb"
    }

    scaling_results_rv(result)
    showNotification("Scaling analysis complete", type = "message", duration = 4)
  })

  # Component breakdown bar chart (with parent-child kernel breakdown)
  output$prof_bar <- renderPlot(bg = "transparent", {
    df <- profile_results_rv()
    if (is.null(df)) {
      plot.new()
      text(0.5, 0.5, "Click 'Run Profile' to benchmark components",
           cex = 1.2, col = viewer_plot_fg())
      return()
    }

    # Exclude the TOTAL row for the bar chart
    bar_df <- df[df$Component != "TOTAL (wbg_boot_fast)", ]
    bar_df$Component <- factor(bar_df$Component,
                               levels = rev(bar_df$Component))

    # Color palette: parents get distinct colors, children get kernel-family shades
    fill_colors <- c(
      "gen_aruma_flex"       = "#4e79a7",
      "co_tstat_pure_cpp"    = "#f28e2b",
      "burg_aic_select_cpp"  = "#e15759",
      "Bootstrap Kernel"     = "#76b7b2",
      "Bootstrap (COBA)"     = "#76b7b2",
      "  AR Generation"      = "#5a9e97",
      "  OLS Detrend"        = "#8ec4bf",
      "  Burg AR Fit"        = "#479089",
      "  CO Transform"       = "#a8d8d3",
      "  Burg (COBA)"        = "#3d7a74",
      "  Kernel Overhead"    = "#c8e8e5",
      "R Overhead"           = "#bab0ac"
    )

    p <- ggplot2::ggplot(bar_df,
           ggplot2::aes(x = Time_ms, y = Component, fill = Component)) +
      ggplot2::geom_col(show.legend = FALSE, width = 0.6) +
      ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%.1f ms (%.0f%%)", Time_ms, Pct)),
        hjust = -0.05, size = 3.5) +
      ggplot2::scale_x_continuous(
        expand = ggplot2::expansion(mult = c(0, 0.35))) +
      ggplot2::scale_fill_manual(values = fill_colors) +
      ggplot2::labs(x = "Time (ms)", y = NULL,
                    title = sprintf("Single Iteration Breakdown (avg of %d reps)",
                                    bar_df$Reps[1])) +
      viewer_plot_theme() +
      ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())

    print(p)
  })

  # Timing details table (with parent-child indentation)
  output$prof_table <- renderDT({
    df <- profile_results_rv()
    if (is.null(df)) {
      return(datatable(data.frame(Message = "No profile data yet."),
                       rownames = FALSE, options = list(dom = "t")))
    }

    # Find the kernel parent row to compute % of kernel for children
    kernel_row <- which(df$Component %in% c("Bootstrap Kernel", "Bootstrap (COBA)"))
    kernel_ms <- if (length(kernel_row) == 1) df$Time_ms[kernel_row] else NA_real_

    # Build display component names with tree prefixes for children
    display_names <- df$Component
    child_indices <- which(df$IsChild)
    if (length(child_indices) > 0) {
      for (i in child_indices) {
        base_name <- trimws(display_names[i])
        if (i == max(child_indices)) {
          display_names[i] <- paste0("\u2514\u2500 ", base_name)
        } else {
          display_names[i] <- paste0("\u251c\u2500 ", base_name)
        }
      }
    }

    # Build % of Kernel column
    pct_kernel <- ifelse(
      df$IsChild & !is.na(kernel_ms) & kernel_ms > 0,
      sprintf("%.0f%%", df$Time_ms / kernel_ms * 100),
      ""
    )

    display <- data.frame(
      Component       = display_names,
      `Time (ms)`     = sprintf("%.2f", df$Time_ms),
      `% of Total`    = sprintf("%.1f%%", df$Pct),
      `% of Kernel`   = pct_kernel,
      `Avg of N Reps` = df$Reps,
      check.names = FALSE, stringsAsFactors = FALSE
    )

    dt <- datatable(display, rownames = FALSE,
                    options = list(dom = "t", pageLength = 20))
    # Bold the TOTAL and kernel parent rows
    bold_rows <- c("TOTAL (wbg_boot_fast)", "Bootstrap Kernel", "Bootstrap (COBA)")
    dt <- formatStyle(dt, "Component",
      fontWeight = styleEqual(bold_rows, rep("bold", length(bold_rows))))
    dt
  })

  # Scaling analysis chart (adapts to mode) — unchanged
  output$prof_scaling <- renderPlot(bg = "transparent", {
    df <- scaling_results_rv()
    if (is.null(df)) {
      plot.new()
      text(0.5, 0.5, "Select a scaling dimension and click 'Run Scaling'",
           cex = 1.2, col = viewer_plot_fg())
      return()
    }

    mode <- df$.mode[1]

    if (mode == "nb") {
      # Line chart: nb vs time
      p <- ggplot2::ggplot(df, ggplot2::aes(x = nb, y = time_ms)) +
        ggplot2::geom_line(linewidth = 0.8, color = "#4e79a7") +
        ggplot2::geom_point(size = 3, color = "#4e79a7") +
        ggplot2::geom_text(
          ggplot2::aes(label = sprintf("%.0f ms", time_ms)),
          vjust = -1, size = 3.5) +
        ggplot2::scale_x_continuous(breaks = df$nb) +
        ggplot2::scale_y_continuous(
          expand = ggplot2::expansion(mult = c(0.05, 0.15))) +
        ggplot2::labs(x = "Bootstrap Replicates (nb)",
                      y = "Time per Iteration (ms)",
                      title = sprintf("Scaling with nb (n = %d)", df$n[1])) +
        viewer_plot_theme()

    } else if (mode == "n") {
      # Line chart: n vs time
      p <- ggplot2::ggplot(df, ggplot2::aes(x = n, y = time_ms)) +
        ggplot2::geom_line(linewidth = 0.8, color = "#e15759") +
        ggplot2::geom_point(size = 3, color = "#e15759") +
        ggplot2::geom_text(
          ggplot2::aes(label = sprintf("%.0f ms", time_ms)),
          vjust = -1, size = 3.5) +
        ggplot2::scale_x_continuous(breaks = df$n) +
        ggplot2::scale_y_continuous(
          expand = ggplot2::expansion(mult = c(0.05, 0.15))) +
        ggplot2::labs(x = "Sample Size (n)",
                      y = "Time per Iteration (ms)",
                      title = sprintf("Scaling with n (nb = %d)", df$nb[1])) +
        viewer_plot_theme()

    } else {
      # Heatmap: n x nb
      df$n_label  <- factor(df$n, levels = sort(unique(df$n)))
      df$nb_label <- factor(df$nb, levels = sort(unique(df$nb)))

      p <- ggplot2::ggplot(df,
             ggplot2::aes(x = nb_label, y = n_label, fill = time_ms)) +
        ggplot2::geom_tile(color = "white", linewidth = 0.5) +
        ggplot2::geom_text(
          ggplot2::aes(label = sprintf("%.0f", time_ms)),
          size = 3.5) +
        ggplot2::scale_fill_gradient(
          low = "#d4edda", high = "#e15759",
          name = "Time (ms)") +
        ggplot2::labs(x = "Bootstrap Replicates (nb)",
                      y = "Sample Size (n)",
                      title = "Time per Iteration (ms): n x nb Grid") +
        viewer_plot_theme() +
        ggplot2::theme(
          panel.grid = ggplot2::element_blank(),
          legend.position = "right")
    }

    print(p)
  })
}
