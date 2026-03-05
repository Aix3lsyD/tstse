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

.repl_action_button <- function(inputId, label, icon_name = NULL,
                                class = "btn-primary", full_width = FALSE) {
  btn_icon <- if (is.null(icon_name)) NULL else icon(icon_name)
  if (requireNamespace("shinyWidgets", quietly = TRUE)) {
    btn_color <- if (grepl("danger", class, fixed = TRUE)) {
      "danger"
    } else if (grepl("success", class, fixed = TRUE)) {
      "success"
    } else if (grepl("secondary", class, fixed = TRUE)) {
      "default"
    } else if (grepl("warning", class, fixed = TRUE)) {
      "warning"
    } else {
      "primary"
    }
    shinyWidgets::actionBttn(
      inputId = inputId,
      label = label,
      icon = btn_icon,
      style = "material-flat",
      color = btn_color,
      size = "sm",
      block = isTRUE(full_width)
    )
  } else {
    actionButton(
      inputId = inputId,
      label = label,
      icon = btn_icon,
      class = class,
      style = if (isTRUE(full_width)) "width: 100%;" else NULL
    )
  }
}

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
  use_fast  <- isTRUE(params$use_fast)
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
        build_innov_gen(innov_label, innov_params, use_fast = use_fast),
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
    div(id = ns("repl_content"),
      br(),
      h4("Replication Comparison"),
      uiOutput(ns("repl_status")),
      # Seed comparison controls (only shown when no published preset)
      uiOutput(ns("seed_comp_controls")),
      uiOutput(ns("repl_innov_filter_ui")),
      wellPanel(
        h5("Forest Plot"),
        plotOutput(ns("repl_forest"), height = "550px")
      ),
      wellPanel(
        h5("Comparison Table"),
        uiOutput(ns("repl_table_caption")),
        DTOutput(ns("repl_table"))
      )
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
  has_waiter <- requireNamespace("waiter", quietly = TRUE)
  has_shinyvalidate <- requireNamespace("shinyvalidate", quietly = TRUE)
  repl_validator <- NULL

  if (has_shinyvalidate) {
    v <- shinyvalidate::InputValidator$new()
    v$add_rule("comp_seed", function(value) {
      num <- suppressWarnings(as.numeric(value))
      int_val <- suppressWarnings(as.integer(value))
      if (is.na(num) || is.na(int_val) || num != int_val || int_val < 1) {
        "Must be an integer >= 1"
      } else {
        NULL
      }
    })
    v$enable()
    repl_validator <- v
  }

  .with_repl_waiter <- function(expr) {
    if (!has_waiter) return(force(expr))
    w <- waiter::Waiter$new(
      id = ns("repl_content"),
      html = waiter::spin_fading_circles(),
      color = "rgba(0, 0, 0, 0.25)"
    )
    w$show()
    on.exit(w$hide(), add = TRUE)
    force(expr)
  }

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
  comp_state <- new.env(parent = emptyenv())
  comp_state$active <- FALSE
  comp_state$cancel_requested <- FALSE
  comp_state$ctx <- NULL
  comp_state$toast_id <- NULL

  .repl_toast_ui <- function(ctx, detail = NULL) {
    n_cells <- max(1L, as.integer(ctx$n_cells))
    done <- max(0L, min(n_cells, as.integer(ctx$cell_idx) - 1L))
    pct <- round(100 * done / n_cells)
    status_line <- sprintf("Running seed comparison... Cell %d/%d", min(as.integer(ctx$cell_idx), n_cells), n_cells)
    if (isTRUE(comp_state$cancel_requested)) {
      status_line <- paste0(status_line, " (stop requested)")
    }
    if (is.null(detail) || !nzchar(detail)) detail <- ""

    tags$div(
      style = "min-width: 320px;",
      tags$div(class = "progress", style = "height: 8px; margin-bottom: 8px;",
               tags$div(class = "progress-bar bg-info",
                        role = "progressbar",
                        style = sprintf("width: %d%%;", pct),
                        `aria-valuenow` = pct, `aria-valuemin` = "0", `aria-valuemax` = "100")),
      tags$div(tags$strong(status_line)),
      if (nzchar(detail)) tags$div(class = "text-muted", style = "font-size: 0.9em;", detail)
    )
  }

  .show_repl_toast <- function(detail = NULL) {
    if (!isTRUE(comp_state$active) || is.null(comp_state$ctx)) return(invisible(NULL))
    if (is.null(comp_state$toast_id)) comp_state$toast_id <- paste0(ns("comp_progress"), "_toast")

    showNotification(
      ui = .repl_toast_ui(comp_state$ctx, detail = detail),
      id = comp_state$toast_id,
      duration = NULL,
      closeButton = FALSE,
      type = "message",
      action = actionLink(
        inputId = ns("cancel_comp_toast"),
        label = icon("ban"),
        style = "padding: 0; color: var(--bs-warning); text-decoration: none;"
      ),
      session = session
    )
    invisible(NULL)
  }

  .clear_repl_toast <- function() {
    if (!is.null(comp_state$toast_id)) {
      try(removeNotification(comp_state$toast_id, session = session), silent = TRUE)
    }
    comp_state$toast_id <- NULL
  }

  .finalize_comp_run <- function(cancelled = FALSE) {
    if (!is.null(comp_state$ctx)) {
      results_df <- do.call(rbind, Filter(Negate(is.null), comp_state$ctx$all_rates))
      if (!is.null(results_df) && nrow(results_df) > 0) {
        comp_results_rv(results_df)
        comp_seed_used(comp_state$ctx$comp_seed)
      }
      if (isTRUE(cancelled)) {
        showNotification(
          sprintf("Comparison cancelled after %d/%d cells", comp_state$ctx$cells_done, comp_state$ctx$n_cells),
          type = "warning", duration = 5
        )
      } else {
        n_out <- if (is.null(results_df) || nrow(results_df) == 0) 0L else nrow(results_df)
        showNotification(
          sprintf("Comparison complete: %d cells with seed %d", n_out, comp_state$ctx$comp_seed),
          type = "message", duration = 6
        )
      }
    }

    .clear_repl_toast()
    comp_state$active <- FALSE
    comp_state$cancel_requested <- FALSE
    comp_state$ctx <- NULL
    invisible(NULL)
  }

  .schedule_comp_next <- function() {
    later::later(function() {
      shiny::withReactiveDomain(session, {
        tryCatch(
          .run_comp_next_cell(),
          error = function(e) {
            .clear_repl_toast()
            comp_state$active <- FALSE
            comp_state$cancel_requested <- FALSE
            comp_state$ctx <- NULL
            showNotification(paste("Comparison run failed:", e$message), type = "error", duration = 8)
          }
        )
      })
    }, delay = 0)
  }

  .run_comp_next_cell <- function() {
    if (!isTRUE(comp_state$active) || is.null(comp_state$ctx)) return(invisible(NULL))
    ctx <- comp_state$ctx

    if (isTRUE(comp_state$cancel_requested)) {
      .finalize_comp_run(cancelled = TRUE)
      return(invisible(NULL))
    }

    if (ctx$cell_idx > ctx$n_cells) {
      .finalize_comp_run(cancelled = FALSE)
      return(invisible(NULL))
    }

    row_i <- ctx$cell_idx
    phi_val        <- ctx$grid$phi[row_i]
    n_val          <- as.integer(ctx$grid$n[row_i])
    innov_label    <- ctx$grid$innov_label[row_i]
    innov_dist_str <- ctx$grid$innov_dist_str[row_i]
    innov_params   <- if (!is.null(ctx$grid$innov_params)) ctx$grid$innov_params[[row_i]] else list()

    .show_repl_toast(detail = sprintf("phi=%.2f, n=%d, %s", phi_val, n_val, innov_dist_str))

    use_vara <- FALSE
    vara <- NULL
    innov_gen <- NULL

    if (innov_label == "Normal" && length(innov_params) == 0) {
      use_vara <- TRUE
      vara <- 1 - phi_val^2
    } else {
      innov_gen <- tryCatch(
        build_innov_gen(innov_label, innov_params, use_fast = ctx$use_fast),
        error = function(e) NULL
      )
      if (is.null(innov_gen)) {
        ctx$cell_idx <- ctx$cell_idx + 1L
        ctx$cells_done <- ctx$cells_done + 1L
        comp_state$ctx <- ctx
        .show_repl_toast()
        .schedule_comp_next()
        return(invisible(NULL))
      }
    }

    if (!is.null(ctx$comp_seed) && !is.na(ctx$comp_seed)) {
      set.seed(as.integer(ctx$comp_seed) * 10000L + row_i)
    }

    results_list <- vector("list", ctx$nsims)
    for (i in seq_len(ctx$nsims)) {
      y <- tryCatch({
        if (use_vara) {
          gen_aruma_flex(n_val, phi = phi_val, vara = vara, plot = FALSE)$y
        } else {
          gen_aruma_flex(n_val, phi = phi_val, innov_gen = innov_gen, plot = FALSE)$y
        }
      }, error = function(e) NULL)
      if (is.null(y)) next

      results_list[[i]] <- tryCatch(
        wbg_boot_fast(y, nb = ctx$nb, maxp = ctx$maxp, criterion = ctx$criterion,
                      bootadj = ctx$bootadj, min_p = ctx$min_p),
        error = function(e) NULL
      )
    }

    results_list <- Filter(Negate(is.null), results_list)
    if (length(results_list) > 0) {
      pvals      <- vapply(results_list, function(r) r$pvalue, numeric(1))
      pvals_asym <- vapply(results_list, function(r) r$pvalue_asymp, numeric(1))
      pvals_adj  <- vapply(results_list, function(r) if (is.null(r$pvalue_adj)) NA_real_ else r$pvalue_adj, numeric(1))
      n_ok <- length(results_list)
      .bse <- function(p, nn) sqrt(p * (1 - p) / nn)

      rej_05      <- mean(pvals < 0.05, na.rm = TRUE)
      rej_asym_05 <- mean(pvals_asym < 0.05, na.rm = TRUE)
      rej_adj_05  <- mean(pvals_adj < 0.05, na.rm = TRUE)

      ctx$all_rates[[row_i]] <- data.frame(
        n = n_val, phi = phi_val,
        innov_dist = innov_dist_str, innov_label = innov_label,
        n_sims = n_ok,
        reject_05 = rej_05, reject_05_se = .bse(rej_05, n_ok),
        reject_asymp_05 = rej_asym_05, reject_asymp_05_se = .bse(rej_asym_05, n_ok),
        reject_adj_05 = rej_adj_05, reject_adj_05_se = .bse(rej_adj_05, n_ok),
        stringsAsFactors = FALSE
      )
    }

    ctx$cell_idx <- ctx$cell_idx + 1L
    ctx$cells_done <- ctx$cells_done + 1L
    comp_state$ctx <- ctx
    .show_repl_toast()
    .schedule_comp_next()
    invisible(NULL)
  }

  # Reset comparison results when preset changes
  observeEvent(preset_name(), {
    if (isTRUE(comp_state$active)) {
      comp_state$cancel_requested <- TRUE
    }
    .clear_repl_toast()
    comp_results_rv(NULL)
    comp_seed_used(NULL)
  })

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

  .innov_col <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    if ("innov_label_a" %in% names(df)) return("innov_label_a")
    if ("innov_label" %in% names(df)) return("innov_label")
    if ("innov_dist" %in% names(df)) return("innov_dist")
    NULL
  }

  .innov_values <- function(df) {
    col <- .innov_col(df)
    if (is.null(col)) return(character(0))
    vals <- unique(as.character(df[[col]]))
    vals <- vals[!is.na(vals) & nzchar(vals)]
    sort(vals)
  }

  output$repl_innov_filter_ui <- renderUI({
    info <- repl_info()
    data_for_choices <- if (!is.null(info)) repl_merged() else seed_comp_merged()
    innov_vals <- .innov_values(data_for_choices)
    if (length(innov_vals) == 0) return(NULL)

    wellPanel(class = "plot-controls",
      fluidRow(
        column(8,
          selectizeInput(
            ns("repl_innov_filter"), "Innovation Filter",
            choices = character(0),
            selected = NULL,
            multiple = FALSE,
            options = list(placeholder = "Select an innovation type")
          )
        )
      )
    )
  })

  observeEvent({
    info <- repl_info()
    data_for_choices <- if (!is.null(info)) repl_merged() else seed_comp_merged()
    paste(.innov_values(data_for_choices), collapse = "|")
  }, {
    info <- repl_info()
    data_for_choices <- if (!is.null(info)) repl_merged() else seed_comp_merged()
    innov_vals <- .innov_values(data_for_choices)
    if (length(innov_vals) == 0) return()

    current <- as.character(input$repl_innov_filter %||% "")
    selected <- if (nzchar(current) && current %in% innov_vals) current else innov_vals[1]

    freezeReactiveValue(input, "repl_innov_filter")
    updateSelectizeInput(
      session, "repl_innov_filter",
      choices = innov_vals,
      selected = selected,
      server = TRUE
    )
  }, ignoreInit = FALSE)

  .filter_by_innov <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(df)
    col <- .innov_col(df)
    if (is.null(col)) return(df)

    sel <- as.character(input$repl_innov_filter %||% "")
    if (!nzchar(sel)) return(df)

    out <- df[as.character(df[[col]]) == sel, , drop = FALSE]
    if (nrow(out) == 0) return(df)
    out
  }

  repl_merged_filtered <- reactive({
    .filter_by_innov(repl_merged())
  })

  seed_comp_merged_filtered <- reactive({
    .filter_by_innov(seed_comp_merged())
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
            .repl_action_button(ns("run_comp"), "Run Comparison",
                                icon_name = "shuffle",
                                class = "btn-sm btn-warning"))
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
    if (!is.null(repl_validator) && !isTRUE(repl_validator$is_valid())) {
      showNotification("Please fix highlighted input errors before running.",
                       type = "warning", duration = 5)
      return()
    }

    if (isTRUE(comp_state$active)) {
      showNotification("Comparison is already running.", type = "warning", duration = 3)
      return()
    }
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

    comp_state$ctx <- list(
      grid = grid,
      n_cells = nrow(grid),
      cell_idx = 1L,
      cells_done = 0L,
      all_rates = vector("list", nrow(grid)),
      comp_seed = as.integer(comp_seed),
      nsims = as.integer(tp$nsims),
      nb = as.integer(tp$nb),
      maxp = as.integer(tp$maxp),
      criterion = tp$criterion,
      bootadj = isTRUE(tp$bootadj),
      min_p = as.integer(tp$min_p),
      use_fast = isTRUE(tp$use_fast)
    )
    comp_results_rv(NULL)
    comp_seed_used(NULL)
    comp_state$cancel_requested <- FALSE
    comp_state$active <- TRUE
    .show_repl_toast()
    .schedule_comp_next()
  })

  observeEvent(input$cancel_comp_toast, {
    if (!isTRUE(comp_state$active)) return()
    comp_state$cancel_requested <- TRUE
    .show_repl_toast()
    showNotification("Stop requested. Current scenario will finish first.",
                     type = "warning", duration = 4)
  }, ignoreInit = TRUE)

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
      merged <- repl_merged_filtered()
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
    merged <- seed_comp_merged_filtered()
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
      merged <- repl_merged_filtered()
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
        innov <- if ("innov_label" %in% names(r)) as.character(r$innov_label) else
                 if ("innov_dist" %in% names(r)) as.character(r$innov_dist) else ""
        cob_ours_pct <- r$reject_05 * 100
        cob_se <- mc_se_pct(cob_ours_pct, r$n_sims)

        forest_rows <- c(forest_rows, list(
          data.frame(phi = r$phi, config = config, method = "COB",
                     source = source_label, rate = r$COB_ref,
                     se = NA_real_, innov = innov, stringsAsFactors = FALSE),
          data.frame(phi = r$phi, config = config, method = "COB",
                     source = "tstse", rate = cob_ours_pct,
                     se = cob_se, innov = innov, stringsAsFactors = FALSE)
        ))

        if (info$has_coba && !is.na(r$reject_adj_05)) {
          coba_ours_pct <- r$reject_adj_05 * 100
          coba_se <- mc_se_pct(coba_ours_pct, r$n_sims)
          forest_rows <- c(forest_rows, list(
            data.frame(phi = r$phi, config = config, method = "COBA",
                       source = source_label, rate = r$COBA_ref,
                       se = NA_real_, innov = innov, stringsAsFactors = FALSE),
            data.frame(phi = r$phi, config = config, method = "COBA",
                       source = "tstse", rate = coba_ours_pct,
                       se = coba_se, innov = innov, stringsAsFactors = FALSE)
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
      forest_df$innov <- factor(forest_df$innov)
      n_innov <- length(unique(forest_df$innov))

      p <- ggplot2::ggplot(forest_df,
             ggplot2::aes(x = rate, y = config, shape = source)) +
        ggplot2::geom_vline(xintercept = 5, linetype = "dashed", linewidth = 0.4) +
        ggplot2::geom_errorbar(
          ggplot2::aes(xmin = rate - 2 * se, xmax = rate + 2 * se),
          width = 0.25, linewidth = 0.4, na.rm = TRUE, orientation = "y",
          position = ggplot2::position_dodge(width = 0.45)) +
        ggplot2::geom_point(
          ggplot2::aes(fill = source), size = 2.5,
          position = ggplot2::position_dodge(width = 0.45), stroke = 0.4) +
        ggplot2::scale_shape_manual(
          values = stats::setNames(c(21, 24), c("tstse", source_label))) +
        ggplot2::scale_fill_manual(
          values = stats::setNames(c("black", "white"), c("tstse", source_label))) +
        ggplot2::labs(x = "Rejection Rate (%)", y = NULL,
                      shape = NULL, fill = NULL)

      if (info$has_coba) {
        if (n_innov > 1) {
          p <- p + ggplot2::facet_grid(method + innov ~ phi_label,
                                        labeller = ggplot2::labeller(phi_label = ggplot2::label_parsed),
                                        scales = "free_x")
        } else {
          p <- p + ggplot2::facet_grid(method ~ phi_label,
                                        labeller = ggplot2::labeller(phi_label = ggplot2::label_parsed),
                                        scales = "free_x")
        }
      } else {
        if (n_innov > 1) {
          p <- p + ggplot2::facet_grid(innov ~ phi_label,
                                        labeller = ggplot2::labeller(phi_label = ggplot2::label_parsed),
                                        scales = "free_x")
        } else {
          p <- p + ggplot2::facet_wrap(~ phi_label,
                                        labeller = ggplot2::labeller(phi_label = ggplot2::label_parsed),
                                        nrow = 1, scales = "free_x")
        }
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
    merged <- seed_comp_merged_filtered()
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

    has_coba <- any(!is.na(merged$reject_adj_05_a) & !is.na(merged$reject_adj_05_b))

    forest_rows <- list()
    for (i in seq_len(nrow(merged))) {
      r <- merged[i, ]
      if (is.na(r$reject_05_a) || is.na(r$reject_05_b)) next
      config <- sprintf("n = %d", r$n)
      innov <- if ("innov_label_a" %in% names(r)) as.character(r$innov_label_a)
               else if ("innov_dist" %in% names(r)) as.character(r$innov_dist)
               else ""

      # COB
      forest_rows <- c(forest_rows, list(
        data.frame(phi = r$phi, config = config, method = "COB",
                   source = label_a, rate = r$reject_05_a * 100,
                   se = r$reject_05_se_a * 100, innov = innov, stringsAsFactors = FALSE),
        data.frame(phi = r$phi, config = config, method = "COB",
                   source = label_b, rate = r$reject_05_b * 100,
                   se = r$reject_05_se_b * 100, innov = innov, stringsAsFactors = FALSE)
      ))

      # COBA (if available)
      if (has_coba && !is.na(r$reject_adj_05_a) && !is.na(r$reject_adj_05_b)) {
        forest_rows <- c(forest_rows, list(
          data.frame(phi = r$phi, config = config, method = "COBA",
                     source = label_a, rate = r$reject_adj_05_a * 100,
                     se = r$reject_adj_05_se_a * 100, innov = innov, stringsAsFactors = FALSE),
          data.frame(phi = r$phi, config = config, method = "COBA",
                     source = label_b, rate = r$reject_adj_05_b * 100,
                     se = r$reject_adj_05_se_b * 100, innov = innov, stringsAsFactors = FALSE)
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
    forest_df$innov <- factor(forest_df$innov)
    n_innov <- length(unique(forest_df$innov))

    p <- ggplot2::ggplot(forest_df,
           ggplot2::aes(x = rate, y = config, shape = source)) +
      ggplot2::geom_vline(xintercept = 5, linetype = "dashed", linewidth = 0.4) +
      ggplot2::geom_errorbar(
        ggplot2::aes(xmin = rate - 2 * se, xmax = rate + 2 * se),
        width = 0.25, linewidth = 0.4, na.rm = TRUE, orientation = "y",
        position = ggplot2::position_dodge(width = 0.45)) +
      ggplot2::geom_point(
        ggplot2::aes(fill = source), size = 2.5,
        position = ggplot2::position_dodge(width = 0.45), stroke = 0.4) +
      ggplot2::scale_shape_manual(
        values = stats::setNames(c(21, 23), c(label_a, label_b))) +
      ggplot2::scale_fill_manual(
        values = stats::setNames(c("black", "white"), c(label_a, label_b))) +
      ggplot2::labs(x = "Rejection Rate (%)", y = NULL,
                    shape = NULL, fill = NULL)

    if (has_coba) {
      if (n_innov > 1) {
        p <- p + ggplot2::facet_grid(method + innov ~ phi_label,
                                      labeller = ggplot2::labeller(phi_label = ggplot2::label_parsed),
                                      scales = "free_x")
      } else {
        p <- p + ggplot2::facet_grid(method ~ phi_label,
                                      labeller = ggplot2::labeller(phi_label = ggplot2::label_parsed),
                                      scales = "free_x")
      }
    } else {
      if (n_innov > 1) {
        p <- p + ggplot2::facet_grid(innov ~ phi_label,
                                      labeller = ggplot2::labeller(phi_label = ggplot2::label_parsed),
                                      scales = "free_x")
      } else {
        p <- p + ggplot2::facet_wrap(~ phi_label,
                                      labeller = ggplot2::labeller(phi_label = ggplot2::label_parsed),
                                      nrow = 1, scales = "free_x")
      }
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
