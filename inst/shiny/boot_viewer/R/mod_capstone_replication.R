# Capstone Sub-Tab: Replication Comparison
# Two modes:
#   1. Published mode: compare against Woodward Table 3 / TSDS Table 8.3
#   2. Seed comparison mode: run a second simulation with a different seed

# Published reference data -----------------------------------------------

# Woodward, Bottone & Gray (1997) Table 3 — COB and COBA (percentage)
.WOODWARD_TABLE3 <- data.frame(
  phi = rep(c(0.8, 0.9, 0.95), each = 5),
  n   = rep(c(50L, 100L, 250L, 500L, 1000L), 3),
  COB_ref  = c( 5.9,  6.4,  3.8, 4.6, 4.1,
                8.1,  5.3,  5.1, 5.7, 3.8,
               10.0,  7.6,  6.4, 6.4, 6.1),
  COBA_ref = c( 5.3,  5.9,  4.1, 5.1, 4.0,
                4.8,  4.1,  5.0, 5.6, 3.9,
                6.2,  4.3,  5.7, 6.2, 6.2),
  stringsAsFactors = FALSE
)

# TSDS (Woodward, Sadler & Robertson) Table 8.3 — COB only (percentage)
.TSDS_TABLE83 <- data.frame(
  n   = rep(c(100L, 200L, 500L), each = 4),
  phi = rep(c(0.8, 0.9, 0.95, 0.99), 3),
  COB_ref = c(4.7, 5.2, 7.4, 11.9,
              4.2, 5.2, 5.5, 10.2,
              4.6, 5.3, 5.5,  8.0),
  stringsAsFactors = FALSE
)

# Map preset names to their reference data
.REPLICATION_PRESETS <- list(
  "Woodward Table 3" = list(
    ref       = .WOODWARD_TABLE3,
    label     = "Woodward, Bottone & Gray (1997) Table 3",
    has_coba  = TRUE
  ),
  "TSDS Table 8.3" = list(
    ref       = .TSDS_TABLE83,
    label     = "TSDS (Woodward, Sadler & Robertson) Table 8.3",
    has_coba  = FALSE
  )
)

# Seed comparison helper --------------------------------------------------

# Run the grid with a given seed and return aggregated rejection rates.
# Simplified version of the main loop in mod_capstone.R (no DB, no log).
.run_seed_comparison <- function(grid, params, base_seed, session = NULL) {
  nsims     <- params$nsims
  nb        <- params$nb
  maxp      <- params$maxp
  criterion <- params$criterion
  bootadj   <- params$bootadj
  min_p     <- params$min_p
  n_cells   <- nrow(grid)

  all_rates <- vector("list", n_cells)

  for (row_i in seq_len(n_cells)) {
    phi_val        <- grid$phi[row_i]
    n_val          <- as.integer(grid$n[row_i])
    innov_label    <- grid$innov_label[row_i]
    innov_dist_str <- grid$innov_dist_str[row_i]
    innov_params   <- if (!is.null(grid$innov_params))
                        grid$innov_params[[row_i]] else list()

    if (!is.null(session)) {
      shiny::incProgress(0, detail = sprintf("Cell %d/%d: phi=%.2f, n=%d, %s",
                                              row_i, n_cells, phi_val, n_val,
                                              innov_dist_str))
    }

    # Build innovation generator (same logic as main loop)
    use_vara <- FALSE
    innov_gen <- NULL

    if (innov_label == "Normal" && length(innov_params) == 0) {
      use_vara <- TRUE
      vara <- 1 - phi_val^2
    } else {
      innov_gen <- tryCatch(
        build_innov_gen(innov_label, innov_params),
        error = function(e) NULL)
      if (is.null(innov_gen)) {
        if (!is.null(session)) shiny::incProgress(1 / n_cells)
        next
      }
    }

    # Cell-specific seed (same formula as main loop)
    if (!is.null(base_seed) && !is.na(base_seed)) {
      set.seed(as.integer(base_seed) * 10000L + row_i)
    }

    results_list <- vector("list", nsims)

    for (i in seq_len(nsims)) {
      y <- tryCatch({
        if (use_vara) {
          gen_aruma_flex(n_val, phi = phi_val, vara = vara, plot = FALSE)$y
        } else {
          gen_aruma_flex(n_val, phi = phi_val, innov_gen = innov_gen, plot = FALSE)$y
        }
      }, error = function(e) NULL)
      if (is.null(y)) next

      results_list[[i]] <- tryCatch(
        wbg_boot_fast(y, nb = nb, maxp = maxp, criterion = criterion,
                      bootadj = bootadj, min_p = min_p),
        error = function(e) NULL)
    }

    results_list <- Filter(Negate(is.null), results_list)

    if (length(results_list) > 0) {
      pvals      <- vapply(results_list, function(r) r$pvalue, numeric(1))
      pvals_asym <- vapply(results_list, function(r) r$pvalue_asymp, numeric(1))
      pvals_adj  <- vapply(results_list, function(r) {
        if (is.null(r$pvalue_adj)) NA_real_ else r$pvalue_adj
      }, numeric(1))
      n_ok <- length(results_list)
      .bse <- function(p, nn) sqrt(p * (1 - p) / nn)

      rej_05      <- mean(pvals < 0.05, na.rm = TRUE)
      rej_asym_05 <- mean(pvals_asym < 0.05, na.rm = TRUE)
      rej_adj_05  <- mean(pvals_adj < 0.05, na.rm = TRUE)

      all_rates[[row_i]] <- data.frame(
        n = n_val, phi = phi_val,
        innov_dist = innov_dist_str, innov_label = innov_label,
        n_sims = n_ok,
        reject_05      = rej_05,      reject_05_se      = .bse(rej_05, n_ok),
        reject_asymp_05 = rej_asym_05, reject_asymp_05_se = .bse(rej_asym_05, n_ok),
        reject_adj_05  = rej_adj_05,  reject_adj_05_se  = .bse(rej_adj_05, n_ok),
        stringsAsFactors = FALSE
      )
    }

    if (!is.null(session)) shiny::incProgress(1 / n_cells)
  }

  do.call(rbind, Filter(Negate(is.null), all_rates))
}

# UI ---------------------------------------------------------------------

mod_capstone_replication_ui <- function(ns) {
  tabPanel("Replication Comparison",
    br(),
    h4("Replication Comparison"),
    uiOutput(ns("repl_status")),
    # Seed comparison controls (only shown when no published preset)
    uiOutput(ns("seed_comp_controls")),
    wellPanel(
      h5("Forest Plot"),
      plotOutput(ns("repl_forest"), height = "450px")
    ),
    wellPanel(
      h5("Comparison Table"),
      uiOutput(ns("repl_table_caption")),
      DTOutput(ns("repl_table"))
    )
  )
}

# Server -----------------------------------------------------------------

mod_capstone_replication_server <- function(input, output, session,
                                            preset_name,
                                            cap_combined_results,
                                            grid_rv,
                                            test_params) {

  ns <- session$ns

  # Which replication preset (if any) is active?
  repl_info <- reactive({
    pn <- preset_name()
    if (is.null(pn) || !pn %in% names(.REPLICATION_PRESETS)) return(NULL)
    .REPLICATION_PRESETS[[pn]]
  })

  # --- Published mode reactives (unchanged) --------------------------------

  repl_merged <- reactive({
    info <- repl_info()
    if (is.null(info)) return(NULL)
    ours <- cap_combined_results()
    if (nrow(ours) == 0) return(NULL)

    ref <- info$ref
    merged <- merge(ref, ours, by = c("n", "phi"), all.x = TRUE)
    if (nrow(merged) == 0) return(NULL)
    merged
  })

  # --- Seed comparison mode reactives --------------------------------------

  comp_results_rv <- reactiveVal(NULL)
  comp_seed_used  <- reactiveVal(NULL)

  # Reset comparison results when preset changes
  observeEvent(preset_name(), { comp_results_rv(NULL); comp_seed_used(NULL) })

  seed_comp_merged <- reactive({
    ours   <- cap_combined_results()
    theirs <- comp_results_rv()
    if (is.null(ours) || nrow(ours) == 0) return(NULL)
    if (is.null(theirs) || nrow(theirs) == 0) return(NULL)
    merged <- merge(ours, theirs, by = c("n", "phi", "innov_dist"),
                    suffixes = c("_a", "_b"))
    if (nrow(merged) == 0) return(NULL)
    merged
  })

  # Seed comparison controls
  output$seed_comp_controls <- renderUI({
    info <- repl_info()
    if (!is.null(info)) return(NULL)  # published mode — no controls

    ours <- cap_combined_results()
    if (is.null(ours) || nrow(ours) == 0) return(NULL)  # no primary results yet

    tp <- test_params()
    default_seed <- as.integer(tp$base_seed %||% 3024) + 1L

    wellPanel(class = "plot-controls",
      fluidRow(
        column(3,
          numericInput(ns("comp_seed"), "Comparison Seed", value = default_seed)
        ),
        column(3,
          div(style = "margin-top: 25px;",
            actionButton(ns("run_comp"), "Run Comparison",
                         icon = icon("shuffle"), class = "btn-sm btn-warning"))
        ),
        column(6,
          div(style = "margin-top: 25px;",
            uiOutput(ns("comp_status")))
        )
      )
    )
  })

  # Run comparison simulation
  observeEvent(input$run_comp, {
    grid <- grid_rv()
    if (is.null(grid) || nrow(grid) == 0) {
      showNotification("No grid defined. Load a preset first.", type = "warning")
      return()
    }
    comp_seed <- input$comp_seed
    if (is.null(comp_seed) || is.na(comp_seed)) {
      showNotification("Enter a valid comparison seed.", type = "warning")
      return()
    }
    tp <- test_params()

    withProgress(message = "Running seed comparison...", value = 0, {
      result <- .run_seed_comparison(grid, tp, base_seed = comp_seed,
                                     session = session)
    })

    if (is.null(result) || nrow(result) == 0) {
      showNotification("Comparison simulation produced no results.", type = "error")
      return()
    }

    comp_results_rv(result)
    comp_seed_used(comp_seed)
    showNotification(
      sprintf("Comparison complete: %d cells with seed %d", nrow(result), comp_seed),
      type = "message", duration = 6)
  })

  output$comp_status <- renderUI({
    comp <- comp_results_rv()
    if (is.null(comp)) return(NULL)
    seed_b <- comp_seed_used()
    tags$span(class = "text-success fw-bold",
      sprintf("Seed %s comparison ready (%d cells)", seed_b, nrow(comp)))
  })

  # --- Status banner (both modes) -----------------------------------------

  output$repl_status <- renderUI({
    info <- repl_info()
    if (!is.null(info)) {
      # Published mode
      merged <- repl_merged()
      if (is.null(merged)) {
        return(div(class = "alert alert-info",
          sprintf("Preset: %s. Run the simulation to see comparison.", info$label)))
      }
      n_matched <- sum(!is.na(merged$reject_05))
      n_total <- nrow(merged)
      return(div(class = "alert alert-success",
        sprintf("Comparing against %s (%d / %d cells have results).",
                info$label, n_matched, n_total)))
    }

    # Seed comparison mode
    ours <- cap_combined_results()
    if (is.null(ours) || nrow(ours) == 0) {
      return(div(class = "alert alert-secondary",
        "Run the main simulation first to establish baseline results (Seed A). ",
        "Then use this tab to run a comparison with a different seed."))
    }

    tp <- test_params()
    comp <- comp_results_rv()
    if (is.null(comp)) {
      return(div(class = "alert alert-info",
        sprintf("Primary results available (%d cells, seed = %s). ",
                nrow(ours), tp$base_seed),
        "Enter a comparison seed below and click 'Run Comparison'."))
    }

    merged <- seed_comp_merged()
    if (is.null(merged)) {
      return(div(class = "alert alert-warning",
        "Comparison results do not match primary grid cells."))
    }

    # Count agreements (COB)
    .check_agree <- function(ra, se_a, rb, se_b) {
      pooled_se <- sqrt(se_a^2 + se_b^2)
      ifelse(pooled_se < 1e-10, TRUE, abs(ra - rb) <= 2 * pooled_se)
    }
    agrees <- .check_agree(merged$reject_05_a, merged$reject_05_se_a,
                           merged$reject_05_b, merged$reject_05_se_b)
    n_agree <- sum(agrees, na.rm = TRUE)
    seed_b <- comp_seed_used()

    div(class = "alert alert-success",
      sprintf("Seed comparison: %s vs %s (%d / %d COB cells agree).",
              tp$base_seed, seed_b, n_agree, nrow(merged)))
  })

  # --- Table caption -------------------------------------------------------

  output$repl_table_caption <- renderUI({
    info <- repl_info()
    if (!is.null(info)) {
      return(p(class = "text-body-secondary",
        "Pass: published value falls within our result \u00B1 2 SE. ",
        "Warn: outside \u00B1 2 SE, but published values have no reported SE ",
        "so comparison is one-sided. Rates shown as percentages."))
    }
    p(class = "text-body-secondary",
      "Agree: |rate_A \u2212 rate_B| \u2264 2 \u00D7 sqrt(SE_A\u00B2 + SE_B\u00B2). ",
      "Rates shown as percentages.")
  })

  # --- Comparison table (both modes) ---------------------------------------

  output$repl_table <- renderDT({
    info <- repl_info()

    if (!is.null(info)) {
      # ---- Published mode (unchanged) ----
      merged <- repl_merged()
      if (is.null(merged)) {
        return(datatable(data.frame(Message = "No comparison data available."),
                         rownames = FALSE, options = list(dom = "t")))
      }

      mc_se_pct <- function(rate_pct, n_sims) {
        p <- rate_pct / 100
        sqrt(p * (1 - p) / n_sims) * 100
      }
      check_match <- function(ours_pct, ref_pct, n_sims) {
        se <- mc_se_pct(ours_pct, n_sims)
        ifelse(is.na(ours_pct), NA_character_,
               ifelse(abs(ours_pct - ref_pct) <= 2 * se, "Pass", "Warn"))
      }
      pm <- "\u00B1"

      if (info$has_coba) {
        cob_ours <- merged$reject_05 * 100
        coba_ours <- merged$reject_adj_05 * 100
        n_sims <- merged$n_sims

        display <- data.frame(
          n   = merged$n,
          phi = merged$phi,
          `COB (Published)` = sprintf("%.1f", merged$COB_ref),
          `COB (Ours)` = ifelse(is.na(cob_ours), "\u2014",
            sprintf("%.1f %s %.1f", cob_ours, pm, mc_se_pct(cob_ours, n_sims))),
          `COB` = check_match(cob_ours, merged$COB_ref, n_sims),
          `COBA (Published)` = sprintf("%.1f", merged$COBA_ref),
          `COBA (Ours)` = ifelse(is.na(coba_ours), "\u2014",
            sprintf("%.1f %s %.1f", coba_ours, pm, mc_se_pct(coba_ours, n_sims))),
          `COBA` = check_match(coba_ours, merged$COBA_ref, n_sims),
          check.names = FALSE, stringsAsFactors = FALSE
        )
        display <- display[order(display$phi, display$n), ]
        dt <- datatable(display, rownames = FALSE,
                        options = list(dom = "t", pageLength = 50, ordering = TRUE))
        for (col in c("COB", "COBA")) {
          dt <- formatStyle(dt, col,
            backgroundColor = styleEqual(c("Pass", "Warn"),
                                          c("#d4edda", "#fff3cd")),
            fontWeight = "bold")
        }
      } else {
        cob_ours <- merged$reject_05 * 100
        n_sims <- merged$n_sims

        display <- data.frame(
          n   = merged$n,
          phi = merged$phi,
          `COB (Published)` = sprintf("%.1f", merged$COB_ref),
          `COB (Ours)` = ifelse(is.na(cob_ours), "\u2014",
            sprintf("%.1f %s %.1f", cob_ours, pm, mc_se_pct(cob_ours, n_sims))),
          Match = check_match(cob_ours, merged$COB_ref, n_sims),
          check.names = FALSE, stringsAsFactors = FALSE
        )
        display <- display[order(display$phi, display$n), ]
        dt <- datatable(display, rownames = FALSE,
                        options = list(dom = "t", pageLength = 50, ordering = TRUE))
        dt <- formatStyle(dt, "Match",
          backgroundColor = styleEqual(c("Pass", "Warn"),
                                        c("#d4edda", "#fff3cd")),
          fontWeight = "bold")
      }
      return(dt)
    }

    # ---- Seed comparison mode ----
    merged <- seed_comp_merged()
    if (is.null(merged)) {
      return(datatable(data.frame(Message = "No comparison data available."),
                       rownames = FALSE, options = list(dom = "t")))
    }

    tp <- test_params()
    seed_a <- tp$base_seed
    seed_b <- comp_seed_used()
    if (identical(as.integer(seed_a), as.integer(seed_b))) {
      seed_a_label <- paste0("Seed ", seed_a, " (primary)")
      seed_b_label <- paste0("Seed ", seed_b, " (comparison)")
    } else {
      seed_a_label <- paste0("Seed ", seed_a)
      seed_b_label <- paste0("Seed ", seed_b)
    }
    pm <- "\u00B1"

    mc_se_pct <- function(rate, se) {
      # Already on proportion scale; convert to pct
      rate_pct <- rate * 100
      se_pct <- se * 100
      ifelse(is.na(rate_pct), "\u2014",
             sprintf("%.1f %s %.1f", rate_pct, pm, se_pct))
    }
    check_agree <- function(ra, se_a, rb, se_b) {
      pooled_se <- sqrt(se_a^2 + se_b^2)
      ifelse(is.na(ra) | is.na(rb), NA_character_,
             ifelse(pooled_se < 1e-10 | abs(ra - rb) <= 2 * pooled_se,
                    "Agree", "Disagree"))
    }

    has_coba <- any(!is.na(merged$reject_adj_05_a) & !is.na(merged$reject_adj_05_b))

    # Build display columns
    innov_col <- if ("innov_label_a" %in% names(merged)) merged$innov_label_a
                 else merged$innov_dist
    display <- data.frame(
      n     = merged$n,
      phi   = merged$phi,
      Innov = innov_col,
      check.names = FALSE, stringsAsFactors = FALSE
    )
    display[[paste0("COB (", seed_a_label, ")")]] <-
      mc_se_pct(merged$reject_05_a, merged$reject_05_se_a)
    display[[paste0("COB (", seed_b_label, ")")]] <-
      mc_se_pct(merged$reject_05_b, merged$reject_05_se_b)
    display[["COB"]] <- check_agree(
      merged$reject_05_a, merged$reject_05_se_a,
      merged$reject_05_b, merged$reject_05_se_b)

    if (has_coba) {
      display[[paste0("COBA (", seed_a_label, ")")]] <-
        mc_se_pct(merged$reject_adj_05_a, merged$reject_adj_05_se_a)
      display[[paste0("COBA (", seed_b_label, ")")]] <-
        mc_se_pct(merged$reject_adj_05_b, merged$reject_adj_05_se_b)
      display[["COBA"]] <- check_agree(
        merged$reject_adj_05_a, merged$reject_adj_05_se_a,
        merged$reject_adj_05_b, merged$reject_adj_05_se_b)
    }

    display <- display[order(display$phi, display$n), ]
    dt <- datatable(display, rownames = FALSE,
                    options = list(dom = "t", pageLength = 50, ordering = TRUE))

    match_cols <- "COB"
    if (has_coba) match_cols <- c(match_cols, "COBA")
    for (col in match_cols) {
      dt <- formatStyle(dt, col,
        backgroundColor = styleEqual(c("Agree", "Disagree"),
                                      c("#d4edda", "#f8d7da")),
        fontWeight = "bold")
    }
    dt
  })

  # --- Forest plot (both modes) --------------------------------------------

  output$repl_forest <- renderPlot(bg = "transparent", {
    info <- repl_info()

    if (!is.null(info)) {
      # ---- Published mode (unchanged) ----
      merged <- repl_merged()
      if (is.null(merged) || all(is.na(merged$reject_05))) {
        plot.new()
        text(0.5, 0.5, "Run the simulation to see comparison",
             cex = 1.2, col = viewer_plot_fg())
        return()
      }

      mc_se_pct <- function(rate_pct, n_sims) {
        p <- rate_pct / 100
        sqrt(p * (1 - p) / n_sims) * 100
      }

      source_label <- info$label
      forest_rows <- list()

      for (i in seq_len(nrow(merged))) {
        r <- merged[i, ]
        if (is.na(r$reject_05)) next
        config <- sprintf("n = %d", r$n)
        cob_ours_pct <- r$reject_05 * 100
        cob_se <- mc_se_pct(cob_ours_pct, r$n_sims)

        forest_rows <- c(forest_rows, list(
          data.frame(phi = r$phi, config = config, method = "COB",
                     source = source_label, rate = r$COB_ref,
                     se = NA_real_, stringsAsFactors = FALSE),
          data.frame(phi = r$phi, config = config, method = "COB",
                     source = "tstse", rate = cob_ours_pct,
                     se = cob_se, stringsAsFactors = FALSE)
        ))

        if (info$has_coba && !is.na(r$reject_adj_05)) {
          coba_ours_pct <- r$reject_adj_05 * 100
          coba_se <- mc_se_pct(coba_ours_pct, r$n_sims)
          forest_rows <- c(forest_rows, list(
            data.frame(phi = r$phi, config = config, method = "COBA",
                       source = source_label, rate = r$COBA_ref,
                       se = NA_real_, stringsAsFactors = FALSE),
            data.frame(phi = r$phi, config = config, method = "COBA",
                       source = "tstse", rate = coba_ours_pct,
                       se = coba_se, stringsAsFactors = FALSE)
          ))
        }
      }

      if (length(forest_rows) == 0) {
        plot.new()
        text(0.5, 0.5, "No matched results", cex = 1.2, col = viewer_plot_fg())
        return()
      }

      forest_df <- do.call(rbind, forest_rows)
      forest_df$source <- factor(forest_df$source,
                                 levels = c("tstse", source_label))
      n_vals <- sort(unique(merged$n), decreasing = TRUE)
      config_levels <- sprintf("n = %d", n_vals)
      forest_df$config <- factor(forest_df$config, levels = config_levels)
      forest_df$phi_label <- sprintf("phi[1] == %.2f", forest_df$phi)

      p <- ggplot2::ggplot(forest_df,
             ggplot2::aes(x = rate, y = config, shape = source)) +
        ggplot2::geom_vline(xintercept = 5, linetype = "dashed", linewidth = 0.4) +
        ggplot2::geom_errorbarh(
          ggplot2::aes(xmin = rate - 2 * se, xmax = rate + 2 * se),
          height = 0.25, linewidth = 0.4, na.rm = TRUE,
          position = ggplot2::position_dodge(width = 0.45)) +
        ggplot2::geom_point(
          ggplot2::aes(fill = source), size = 2.5,
          position = ggplot2::position_dodge(width = 0.45), stroke = 0.4) +
        ggplot2::scale_shape_manual(
          values = stats::setNames(c(21, 24), c("tstse", source_label))) +
        ggplot2::scale_fill_manual(
          values = stats::setNames(c("black", "white"), c("tstse", source_label))) +
        ggplot2::labs(x = "Rejection Rate (%)", y = NULL,
                      shape = NULL, fill = NULL) +
        ggplot2::guides(
          shape = ggplot2::guide_legend(override.aes = list(size = 2.5)),
          fill  = ggplot2::guide_legend(override.aes = list(size = 2.5)))

      if (info$has_coba) {
        p <- p + ggplot2::facet_grid(method ~ phi_label,
                                      labeller = ggplot2::label_parsed,
                                      scales = "free_x")
      } else {
        p <- p + ggplot2::facet_wrap(~ phi_label,
                                      labeller = ggplot2::label_parsed,
                                      nrow = 1, scales = "free_x")
      }

      p <- p + viewer_plot_theme() +
        ggplot2::theme(
          legend.position = "bottom",
          legend.margin = ggplot2::margin(t = 0),
          strip.text = ggplot2::element_text(face = "bold"),
          panel.grid.major.y = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank())

      print(p)
      return()
    }

    # ---- Seed comparison mode ----
    merged <- seed_comp_merged()
    if (is.null(merged)) {
      plot.new()
      msg <- if (is.null(cap_combined_results()) || nrow(cap_combined_results()) == 0)
               "Run the main simulation first"
             else "Run a comparison to see the forest plot"
      text(0.5, 0.5, msg, cex = 1.2, col = viewer_plot_fg())
      return()
    }

    tp <- test_params()
    seed_a <- tp$base_seed
    seed_b <- comp_seed_used()
    if (identical(as.integer(seed_a), as.integer(seed_b))) {
      label_a <- paste0("Seed ", seed_a, " (primary)")
      label_b <- paste0("Seed ", seed_b, " (comparison)")
    } else {
      label_a <- paste0("Seed ", seed_a)
      label_b <- paste0("Seed ", seed_b)
    }

    mc_se_pct_val <- function(rate, se) {
      # rate and se are on proportion scale; return pct
      list(rate_pct = rate * 100, se_pct = se * 100)
    }

    has_coba <- any(!is.na(merged$reject_adj_05_a) & !is.na(merged$reject_adj_05_b))

    forest_rows <- list()
    for (i in seq_len(nrow(merged))) {
      r <- merged[i, ]
      if (is.na(r$reject_05_a) || is.na(r$reject_05_b)) next
      config <- sprintf("n = %d", r$n)

      # COB
      forest_rows <- c(forest_rows, list(
        data.frame(phi = r$phi, config = config, method = "COB",
                   source = label_a, rate = r$reject_05_a * 100,
                   se = r$reject_05_se_a * 100, stringsAsFactors = FALSE),
        data.frame(phi = r$phi, config = config, method = "COB",
                   source = label_b, rate = r$reject_05_b * 100,
                   se = r$reject_05_se_b * 100, stringsAsFactors = FALSE)
      ))

      # COBA (if available)
      if (has_coba && !is.na(r$reject_adj_05_a) && !is.na(r$reject_adj_05_b)) {
        forest_rows <- c(forest_rows, list(
          data.frame(phi = r$phi, config = config, method = "COBA",
                     source = label_a, rate = r$reject_adj_05_a * 100,
                     se = r$reject_adj_05_se_a * 100, stringsAsFactors = FALSE),
          data.frame(phi = r$phi, config = config, method = "COBA",
                     source = label_b, rate = r$reject_adj_05_b * 100,
                     se = r$reject_adj_05_se_b * 100, stringsAsFactors = FALSE)
        ))
      }
    }

    if (length(forest_rows) == 0) {
      plot.new()
      text(0.5, 0.5, "No matched results", cex = 1.2, col = viewer_plot_fg())
      return()
    }

    forest_df <- do.call(rbind, forest_rows)
    forest_df$source <- factor(forest_df$source, levels = c(label_a, label_b))

    n_vals <- sort(unique(merged$n), decreasing = TRUE)
    config_levels <- sprintf("n = %d", n_vals)
    forest_df$config <- factor(forest_df$config, levels = config_levels)
    forest_df$phi_label <- sprintf("phi[1] == %.2f", forest_df$phi)

    p <- ggplot2::ggplot(forest_df,
           ggplot2::aes(x = rate, y = config, shape = source)) +
      ggplot2::geom_vline(xintercept = 5, linetype = "dashed", linewidth = 0.4) +
      ggplot2::geom_errorbarh(
        ggplot2::aes(xmin = rate - 2 * se, xmax = rate + 2 * se),
        height = 0.25, linewidth = 0.4, na.rm = TRUE,
        position = ggplot2::position_dodge(width = 0.45)) +
      ggplot2::geom_point(
        ggplot2::aes(fill = source), size = 2.5,
        position = ggplot2::position_dodge(width = 0.45), stroke = 0.4) +
      ggplot2::scale_shape_manual(
        values = stats::setNames(c(21, 23), c(label_a, label_b))) +
      ggplot2::scale_fill_manual(
        values = stats::setNames(c("black", "white"), c(label_a, label_b))) +
      ggplot2::labs(x = "Rejection Rate (%)", y = NULL,
                    shape = NULL, fill = NULL) +
      ggplot2::guides(
        shape = ggplot2::guide_legend(override.aes = list(size = 2.5)),
        fill  = ggplot2::guide_legend(override.aes = list(size = 2.5)))

    if (has_coba) {
      p <- p + ggplot2::facet_grid(method ~ phi_label,
                                    labeller = ggplot2::label_parsed,
                                    scales = "free_x")
    } else {
      p <- p + ggplot2::facet_wrap(~ phi_label,
                                    labeller = ggplot2::label_parsed,
                                    nrow = 1, scales = "free_x")
    }

    p <- p + viewer_plot_theme() +
      ggplot2::theme(
        legend.position = "bottom",
        legend.margin = ggplot2::margin(t = 0),
        strip.text = ggplot2::element_text(face = "bold"),
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank())

    print(p)
  })
}
