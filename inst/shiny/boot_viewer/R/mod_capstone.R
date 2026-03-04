# =============================================================================
# Module: Capstone Simulation (Grid-based Monte Carlo runner)
# =============================================================================

# Preset configurations matching CapstoneSimulations_V1.Rmd
.CAPSTONE_PRESETS <- list(
  "Capstone V1 (Full)" = list(
    description = "Full grid from CapstoneSimulations_V1.Rmd + Laplace scenario",
    n     = c(50L, 100L, 250L, 500L),
    phi   = c(0.80, 0.95, 0.99),
    innov = list(
      list(label = "Normal",          dist_str = "norm",
           params = list()),
      list(label = "GARCH",           dist_str = "arch(8)",
           params = list(garch_omega = 0.1,
                         garch_alpha = c(0.2, 0.175, 0.15, 0.125, 0.1, 0.075, 0.05, 0.025))),
      list(label = "Heteroscedastic", dist_str = "hetero",
           params = list(hetero_w = function(n) seq(1, 10, length.out = n))),
      list(label = "Student's t",     dist_str = "t(3)",
           params = list(t_df = 3, t_scale = FALSE)),
      list(label = "Laplace",         dist_str = "laplace",
           params = list(lap_scale = 1))
    ),
    test_defaults = list(nsims = 1000L, nb = 399L, maxp = 5L,
                         criterion = "aic", bootadj = TRUE, min_p = 1L,
                         seed = 3024L)
  ),
  "Size Study Only" = list(
    description = "Null hypothesis only (phi=0) for size calibration",
    n     = c(50L, 100L, 250L, 500L),
    phi   = c(0.0),
    innov = list(
      list(label = "Normal",      dist_str = "norm",    params = list()),
      list(label = "Student's t", dist_str = "t(3)",    params = list(t_df = 3, t_scale = FALSE)),
      list(label = "Laplace",     dist_str = "laplace", params = list(lap_scale = 1))
    ),
    test_defaults = list(nsims = 1000L, nb = 399L, maxp = 5L,
                         criterion = "aic", bootadj = TRUE, min_p = 1L,
                         seed = 4001L)
  ),
  "Quick Smoke Test" = list(
    description = "Small grid for testing (2 cells, 50 sims each)",
    n     = c(100L),
    phi   = c(0.80, 0.95),
    innov = list(
      list(label = "Normal", dist_str = "norm", params = list())
    ),
    test_defaults = list(nsims = 50L, nb = 199L, maxp = 5L,
                         criterion = "aic", bootadj = TRUE, min_p = 1L,
                         seed = 9999L)
  ),
  "Woodward Table 3" = list(
    description = "Replication of Woodward, Bottone & Gray (1997) Table 3",
    n     = c(50L, 100L, 250L, 500L, 1000L),
    phi   = c(0.80, 0.90, 0.95),
    innov = list(
      list(label = "Normal", dist_str = "norm", params = list())
    ),
    test_defaults = list(nsims = 1000L, nb = 399L, maxp = 1L,
                         criterion = "aic", bootadj = TRUE, min_p = 1L,
                         seed = 1997L)
  ),
  "TSDS Table 8.3" = list(
    description = "Replication of TSDS (Woodward, Sadler & Robertson) Table 8.3",
    n     = c(100L, 200L, 500L),
    phi   = c(0.80, 0.90, 0.95, 0.99),
    innov = list(
      list(label = "Normal", dist_str = "norm", params = list())
    ),
    test_defaults = list(nsims = 1000L, nb = 199L, maxp = 5L,
                         criterion = "aic", bootadj = FALSE, min_p = 1L,
                         seed = 2017L)
  )
)

# Innovation distribution choices for the "Add Row" form
.CAP_INNOV_CHOICES <- c("Normal", "Student's t", "Skew-t", "GED",
                         "Laplace", "Uniform", "Mixture Normal",
                         "GARCH", "Heteroscedastic")

# Expand a preset into a grid data.frame
.expand_preset_grid <- function(preset) {
  rows <- list()
  idx <- 1L
  for (inv in preset$innov) {
    for (phi_val in preset$phi) {
      for (n_val in preset$n) {
        rows[[idx]] <- list(
          phi            = phi_val,
          n              = n_val,
          innov_label    = inv$label,
          innov_dist_str = inv$dist_str,
          innov_params   = list(inv$params)
        )
        idx <- idx + 1L
      }
    }
  }
  df <- data.frame(
    phi            = vapply(rows, `[[`, numeric(1), "phi"),
    n              = vapply(rows, `[[`, integer(1), "n"),
    innov_label    = vapply(rows, `[[`, character(1), "innov_label"),
    innov_dist_str = vapply(rows, `[[`, character(1), "innov_dist_str"),
    stringsAsFactors = FALSE
  )
  df$innov_params <- lapply(rows, function(r) r$innov_params[[1]])
  df
}

# =============================================================================
# UI
# =============================================================================

mod_capstone_ui <- function(id) {
  ns <- NS(id)

  tabPanel(
    "Capstone Simulation",
    br(),
    fluidRow(
      # --- Sidebar (col-3) ---
      column(3,
        # Panel 1: Run Controls (top for visibility)
        wellPanel(
          h4("Run"),
          actionButton(ns("run_grid"), "Run Grid Simulation",
                       icon = icon("play"), class = "btn-primary",
                       style = "width: 100%;"),
          uiOutput(ns("run_status_ui")),
          hr(),
          textInput(ns("db_save_path"), "DuckDB Path (optional)",
                    value = "", placeholder = "Leave blank to skip saving"),
          p(class = "text-body-secondary", style = "font-size: 0.8em; margin-top: -5px;",
            "Creates DB if it doesn't exist. Leave empty to run without saving."),
          textInput(ns("batch_label"), "Batch Label", value = "Capstone Run"),
          radioButtons(ns("run_mode"), NULL,
                       choices = c("All cells" = "all",
                                   "Selected rows" = "selected"),
                       selected = "all", inline = TRUE)
        ),

        # Panel 2: Preset & Grid
        wellPanel(
          h4("Grid Preset"),
          selectInput(ns("preset"), "Load Preset",
                      choices = c("Custom", names(.CAPSTONE_PRESETS))),
          actionButton(ns("load_preset"), "Load Preset",
                       icon = icon("upload"), class = "btn-sm btn-outline-secondary"),
          hr(),
          h5("Add Config Row"),
          fluidRow(
            column(6, numericInput(ns("new_phi"), "Phi", 0.95,
                                   min = 0, max = 0.999, step = 0.01)),
            column(6, numericInput(ns("new_n"), "n", 200,
                                   min = 10, max = 5000, step = 10))
          ),
          selectInput(ns("new_innov"), "Innovation",
                      choices = .CAP_INNOV_CHOICES, selected = "Normal"),

          # Distribution-specific params (minimal set for custom rows)
          conditionalPanel(
            condition = "input.new_innov == \"Student's t\"", ns = ns,
            numericInput(ns("new_t_df"), "df", value = 3, min = 1, step = 1)
          ),
          conditionalPanel(
            condition = "input.new_innov == 'Laplace'", ns = ns,
            numericInput(ns("new_lap_scale"), "Scale", value = 1, min = 0.001, step = 0.1)
          ),
          conditionalPanel(
            condition = "input.new_innov == 'GARCH'", ns = ns,
            numericInput(ns("new_garch_omega"), "Omega", value = 0.1, min = 0.001, step = 0.01),
            textInput(ns("new_garch_alpha"), "Alpha (comma-separated)",
                      value = "0.2, 0.175, 0.15, 0.125, 0.1, 0.075, 0.05, 0.025")
          ),
          conditionalPanel(
            condition = "input.new_innov == 'Heteroscedastic'", ns = ns,
            selectInput(ns("new_hetero_shape"), "Shape",
                        choices = c("linear", "sqrt", "log", "exp", "power",
                                    "step", "periodic"),
                        selected = "linear"),
            numericInput(ns("new_hetero_from"), "From", value = 1, min = 0.01, step = 0.5),
            numericInput(ns("new_hetero_to"), "To", value = 10, min = 0.01, step = 0.5)
          ),
          conditionalPanel(
            condition = "input.new_innov == 'Mixture Normal'", ns = ns,
            numericInput(ns("new_mix_sd1"), "sd1", value = 1, min = 0.001, step = 0.1),
            numericInput(ns("new_mix_sd2"), "sd2", value = 3, min = 0.001, step = 0.1),
            numericInput(ns("new_mix_prob1"), "prob1", value = 0.9, min = 0.01, max = 0.99, step = 0.05)
          ),
          conditionalPanel(
            condition = "input.new_innov == 'Skew-t'", ns = ns,
            numericInput(ns("new_skt_df"), "df", value = 5, min = 3, step = 1),
            numericInput(ns("new_skt_alpha"), "alpha", value = 0, step = 0.1)
          ),
          conditionalPanel(
            condition = "input.new_innov == 'GED'", ns = ns,
            numericInput(ns("new_ged_nu"), "nu", value = 2, min = 0.1, step = 0.1)
          ),
          conditionalPanel(
            condition = "input.new_innov == 'Uniform'", ns = ns,
            numericInput(ns("new_unif_hw"), "Half-width", value = 1.732, min = 0.001, step = 0.1)
          ),

          fluidRow(
            column(6, actionButton(ns("add_row"), "Add Row",
                                   icon = icon("plus"), class = "btn-sm btn-success")),
            column(6, actionButton(ns("clear_grid"), "Clear All",
                                   icon = icon("trash"), class = "btn-sm btn-outline-danger"))
          )
        ),

        # Panel 3: Test Parameters
        wellPanel(
          h4("Test Parameters"),
          numericInput(ns("nsims"), "Simulations per cell", value = 1000,
                       min = 1, max = 10000, step = 100),
          numericInput(ns("nb"), "Bootstrap replicates (nb)", value = 399,
                       min = 1, max = 9999, step = 100),
          numericInput(ns("maxp"), "Max AR order", value = 5,
                       min = 1, max = 20, step = 1),
          selectInput(ns("criterion"), "Criterion",
                      choices = c("aic", "aicc", "bic"), selected = "aic"),
          checkboxInput(ns("bootadj"), "COBA adjustment", value = TRUE),
          checkboxInput(ns("minp1"), "Min AR order = 1", value = TRUE),
          numericInput(ns("seed"), "Base seed", value = 3024)
        )
      ),

      # --- Main area (col-9) ---
      column(9,
        tabsetPanel(
          id = ns("cap_tabs"),

          # Tab 1: Setup (Grid Definition)
          tabPanel("Setup",
            br(),
            p(class = "text-body-secondary",
              "Define the simulation grid. Each row is one (n, phi, innov_dist) cell.",
              "Select rows to run a subset. Load a preset for published configurations."),
            DTOutput(ns("grid_table")),
            br(),
            fluidRow(
              column(6, actionButton(ns("delete_selected"), "Delete Selected Rows",
                                     icon = icon("minus"), class = "btn-sm btn-outline-danger"))
            )
          ),

          # Tab 2: Coverage / Overview
          mod_capstone_overview_ui(ns),
          # Tab 3: Individual Results
          mod_capstone_individual_ui(ns),
          # Tab 4: Analysis Grids
          mod_capstone_grids_ui(ns),
          # Tab 5: Replication Comparison
          mod_capstone_replication_ui(ns),
          # Tab 6: Plots
          mod_capstone_plots_ui(ns),
          # Tab 7: P-value Diagnostics
          mod_capstone_pvalue_ui(ns),
          # Tab 8: Bootstrap Distribution
          mod_capstone_bootdist_ui(ns),
          # Tab 9: Innovation Diagnostics
          mod_capstone_innov_diag_ui(ns),
          # Tab 10: Null Model Diagnostics
          mod_capstone_null_diag_ui(ns),

          # Tab 11: Run Log
          tabPanel("Run Log",
            br(),
            verbatimTextOutput(ns("run_log"))
          )
        )
      )
    )
  )
}

# =============================================================================
# Server
# =============================================================================

mod_capstone_server <- function(id, con, db_path, db_refresh_trigger,
                                init_choices) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # When a DB is loaded at startup, it's read-only for viewing.
    # Save path defaults empty so new saves go to a fresh DB.
    # Store the loaded path to prevent overwriting.
    loaded_db_path <- if (!is.null(db_path) && nzchar(db_path)) db_path else NULL

    # --- Grid state ---
    grid_rv <- reactiveVal(data.frame(
      phi = numeric(0), n = integer(0),
      innov_label = character(0), innov_dist_str = character(0),
      stringsAsFactors = FALSE
    ))

    # In-memory accumulated rejection rates (persists across runs within session)
    run_results_rv <- reactiveVal(data.frame(
      n = integer(0), phi = numeric(0), innov_dist = character(0),
      innov_label = character(0), n_sims = integer(0),
      reject_05 = numeric(0), reject_05_se = numeric(0),
      reject_asymp_05 = numeric(0), reject_asymp_05_se = numeric(0),
      reject_adj_05 = numeric(0), reject_adj_05_se = numeric(0),
      stringsAsFactors = FALSE
    ))

    # Raw per-cell simulation results (for detail sub-tabs)
    # Keyed by "n|phi|innov_dist" -> list(results, n, phi, innov_dist,
    #   innov_label, innov_params, timestamp)
    raw_results_rv <- reactiveVal(list())

    # Run log
    run_log_rv <- reactiveVal("")

    # =========================================================================
    # Preset loading
    # =========================================================================

    observeEvent(input$load_preset, {
      preset_name <- input$preset
      if (preset_name == "Custom" || is.null(preset_name)) {
        grid_rv(data.frame(
          phi = numeric(0), n = integer(0),
          innov_label = character(0), innov_dist_str = character(0),
          stringsAsFactors = FALSE
        ))
        return()
      }
      preset <- .CAPSTONE_PRESETS[[preset_name]]
      if (is.null(preset)) return()

      grid_df <- .expand_preset_grid(preset)
      grid_rv(grid_df)

      # Also update test params from preset defaults
      defs <- preset$test_defaults
      if (!is.null(defs$nsims))     updateNumericInput(session, "nsims", value = defs$nsims)
      if (!is.null(defs$nb))        updateNumericInput(session, "nb", value = defs$nb)
      if (!is.null(defs$maxp))      updateNumericInput(session, "maxp", value = defs$maxp)
      if (!is.null(defs$criterion)) updateSelectInput(session, "criterion", selected = defs$criterion)
      if (!is.null(defs$bootadj))   updateCheckboxInput(session, "bootadj", value = defs$bootadj)
      if (!is.null(defs$min_p))     updateCheckboxInput(session, "minp1", value = defs$min_p == 1L)
      if (!is.null(defs$seed))      updateNumericInput(session, "seed", value = defs$seed)

      showNotification(
        sprintf("Loaded '%s': %d cells", preset_name, nrow(grid_df)),
        type = "message", duration = 4)
    })

    # =========================================================================
    # Add / Delete / Clear rows
    # =========================================================================

    observeEvent(input$add_row, {
      phi_val <- input$new_phi
      n_val <- as.integer(input$new_n)
      innov_label <- input$new_innov

      if (is.null(phi_val) || is.na(phi_val)) phi_val <- 0.95
      if (is.null(n_val) || is.na(n_val)) n_val <- 200L

      # Build params from the add-row form
      params <- switch(innov_label,
        "Normal"      = list(),
        "Student's t" = list(t_df = input$new_t_df %||% 3, t_scale = FALSE),
        "Skew-t"      = list(skt_df = input$new_skt_df %||% 5,
                             skt_alpha = input$new_skt_alpha %||% 0),
        "GED"         = list(ged_nu = input$new_ged_nu %||% 2),
        "Laplace"     = list(lap_scale = input$new_lap_scale %||% 1),
        "Uniform"     = list(unif_hw = input$new_unif_hw %||% sqrt(3)),
        "Mixture Normal" = list(mix_sd1 = input$new_mix_sd1 %||% 1,
                                mix_sd2 = input$new_mix_sd2 %||% 3,
                                mix_prob1 = input$new_mix_prob1 %||% 0.9),
        "GARCH"       = list(garch_omega = input$new_garch_omega %||% 0.1,
                             garch_alpha = input$new_garch_alpha %||%
                               "0.2,0.175,0.15,0.125,0.1,0.075,0.05,0.025"),
        "Heteroscedastic" = list(hetero_shape = input$new_hetero_shape %||% "linear",
                                  hetero_from = input$new_hetero_from %||% 1,
                                  hetero_to = input$new_hetero_to %||% 10),
        list()
      )

      dist_str <- build_innov_dist_str(innov_label, params)

      new_row <- data.frame(
        phi = phi_val, n = n_val,
        innov_label = innov_label, innov_dist_str = dist_str,
        stringsAsFactors = FALSE
      )
      new_row$innov_params <- list(params)

      current <- grid_rv()
      if (nrow(current) > 0 && !is.null(current$innov_params)) {
        grid_rv(rbind(current, new_row))
      } else {
        grid_rv(new_row)
      }
    })

    observeEvent(input$delete_selected, {
      sel <- input$grid_table_rows_selected
      if (length(sel) == 0) {
        showNotification("No rows selected.", type = "warning", duration = 3)
        return()
      }
      current <- grid_rv()
      grid_rv(current[-sel, , drop = FALSE])
    })

    observeEvent(input$clear_grid, {
      grid_rv(data.frame(
        phi = numeric(0), n = integer(0),
        innov_label = character(0), innov_dist_str = character(0),
        stringsAsFactors = FALSE
      ))
    })

    # =========================================================================
    # Grid table
    # =========================================================================

    output$grid_table <- renderDT({
      grid <- grid_rv()
      if (nrow(grid) == 0) {
        return(datatable(
          data.frame(Message = "No grid cells defined. Load a preset or add rows."),
          rownames = FALSE, options = list(dom = "t"),
          selection = "none"
        ))
      }

      display <- data.frame(
        Phi     = grid$phi,
        n       = grid$n,
        Innovation = grid$innov_label,
        `DB Key`   = grid$innov_dist_str,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )

      datatable(
        display,
        rownames = FALSE,
        selection = "multiple",
        options = list(
          dom = "tip",
          pageLength = 25,
          ordering = TRUE,
          columnDefs = list(
            list(targets = c(0, 1), className = "dt-body-right")
          )
        )
      )
    }, server = FALSE)

    # =========================================================================
    # Simulation execution
    # =========================================================================

    observeEvent(input$run_grid, {
      grid <- grid_rv()
      if (nrow(grid) == 0) {
        showNotification("Grid is empty. Load a preset or add rows.", type = "warning")
        return()
      }

      # Determine rows to run
      if (input$run_mode == "selected") {
        rows_to_run <- input$grid_table_rows_selected
        if (length(rows_to_run) == 0) {
          showNotification("No rows selected in the grid table.", type = "warning")
          return()
        }
      } else {
        rows_to_run <- seq_len(nrow(grid))
      }

      # Snapshot test params
      nsims     <- as.integer(input$nsims %||% 1000)
      nb        <- as.integer(input$nb %||% 399)
      maxp      <- as.integer(input$maxp %||% 5)
      maxp      <- min(maxp, 20L)
      criterion <- match.arg(input$criterion, c("aic", "aicc", "bic"))
      bootadj   <- isTRUE(input$bootadj)
      min_p     <- if (isTRUE(input$minp1)) 1L else 0L
      base_seed <- input$seed
      label     <- trimws(input$batch_label %||% "")
      if (!nzchar(label)) label <- NULL

      # Initialize log
      log_lines <- character(0)
      .log <- function(msg) {
        log_lines <<- c(log_lines,
          sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), msg))
        run_log_rv(paste(log_lines, collapse = "\n"))
      }
      run_log_rv("")
      .log(sprintf("Starting grid simulation: %d cells, %d sims each",
                    length(rows_to_run), nsims))

      # Determine DB save path (may be empty → run without saving)
      save_path <- trimws(input$db_save_path %||% "")
      save_to_db <- nzchar(save_path)

      # Prevent overwriting the loaded (read-only) database
      if (save_to_db && !is.null(loaded_db_path) &&
          normalizePath(save_path, mustWork = FALSE) ==
          normalizePath(loaded_db_path, mustWork = FALSE)) {
        showNotification(
          "Cannot save to the loaded database (read-only). Enter a different path.",
          type = "error", duration = 8)
        return()
      }

      write_con <- NULL
      batch_id <- NULL

      if (save_to_db) {
        write_con <- tryCatch(
          mc_db_connect(save_path, read_only = FALSE),
          error = function(e) {
            showNotification(paste("DB connection error:", e$message),
                             type = "error", duration = 8)
            NULL
          }
        )
        if (is.null(write_con)) return()
        on.exit(DBI::dbDisconnect(write_con, shutdown = TRUE), add = TRUE)

        mc_db_init(write_con)
        batch_id <- .mc_create_batch(write_con, label = label)
        .log(sprintf("Created batch #%d ('%s') in %s", batch_id,
                      label %||% "", save_path))
      } else {
        .log("Running without DB save (no path specified)")
      }

      n_cells <- length(rows_to_run)
      t_total_start <- proc.time()[["elapsed"]]
      cells_ok <- 0L
      cells_err <- 0L

      withProgress(message = "Running capstone grid...", value = 0, {

        for (cell_idx in seq_along(rows_to_run)) {
          row_i <- rows_to_run[cell_idx]
          phi_val        <- grid$phi[row_i]
          n_val          <- as.integer(grid$n[row_i])
          innov_label    <- grid$innov_label[row_i]
          innov_dist_str <- grid$innov_dist_str[row_i]
          innov_params   <- if (!is.null(grid$innov_params))
                              grid$innov_params[[row_i]] else list()

          incProgress(0, detail = sprintf("Cell %d/%d: phi=%.2f, n=%d, %s",
                                           cell_idx, n_cells, phi_val, n_val,
                                           innov_dist_str))
          .log(sprintf("Cell %d/%d: phi=%.2f n=%d %s",
                       cell_idx, n_cells, phi_val, n_val, innov_dist_str))

          # Build innovation generator
          use_vara <- FALSE
          innov_gen <- NULL

          if (innov_label == "Normal" && length(innov_params) == 0) {
            # Match Rmd: innov_gen=NULL, vara=1-phi^2 for Normal scenario
            use_vara <- TRUE
            vara <- 1 - phi_val^2
          } else {
            innov_gen <- tryCatch(
              build_innov_gen(innov_label, innov_params),
              error = function(e) {
                .log(sprintf("  Generator error: %s", e$message))
                NULL
              }
            )
            if (is.null(innov_gen) && !use_vara) {
              .log("  SKIPPED (generator failed)")
              cells_err <- cells_err + 1L
              incProgress(1 / n_cells)
              next
            }
          }

          # Cell-specific seed for reproducibility
          if (!is.null(base_seed) && !is.na(base_seed)) {
            set.seed(as.integer(base_seed) * 10000L + row_i)
          }

          # Inner simulation loop
          results_list <- vector("list", nsims)
          failed_count <- 0L
          t_cell_start <- proc.time()[["elapsed"]]

          for (i in seq_len(nsims)) {
            # Generate series
            y <- tryCatch({
              if (use_vara) {
                gen_aruma_flex(n_val, phi = phi_val, vara = vara,
                               plot = FALSE)$y
              } else {
                gen_aruma_flex(n_val, phi = phi_val, innov_gen = innov_gen,
                               plot = FALSE)$y
              }
            }, error = function(e) NULL)
            if (is.null(y)) { failed_count <- failed_count + 1L; next }

            # Bootstrap test
            results_list[[i]] <- tryCatch(
              wbg_boot_fast(y, nb = nb, maxp = maxp, criterion = criterion,
                            bootadj = bootadj, min_p = min_p),
              error = function(e) { failed_count <<- failed_count + 1L; NULL }
            )
          }

          # Remove NULLs from failed iterations
          results_list <- Filter(Negate(is.null), results_list)
          elapsed_cell <- proc.time()[["elapsed"]] - t_cell_start

          # Compute and store in-memory rejection rates
          if (length(results_list) > 0) {
            pvals      <- vapply(results_list, function(r) r$pvalue, numeric(1))
            pvals_asym <- vapply(results_list, function(r) r$pvalue_asymp, numeric(1))
            pvals_adj  <- vapply(results_list, function(r) {
              if (is.null(r$pvalue_adj)) NA_real_ else r$pvalue_adj
            }, numeric(1))
            n_ok <- length(results_list)
            .bse <- function(p, nn) sqrt(p * (1 - p) / nn)

            cell_rates <- data.frame(
              n = n_val, phi = phi_val, innov_dist = innov_dist_str,
              innov_label = innov_label, n_sims = n_ok,
              reject_05          = mean(pvals < 0.05, na.rm = TRUE),
              reject_05_se       = .bse(mean(pvals < 0.05, na.rm = TRUE), n_ok),
              reject_asymp_05    = mean(pvals_asym < 0.05, na.rm = TRUE),
              reject_asymp_05_se = .bse(mean(pvals_asym < 0.05, na.rm = TRUE), n_ok),
              reject_adj_05      = mean(pvals_adj < 0.05, na.rm = TRUE),
              reject_adj_05_se   = .bse(mean(pvals_adj < 0.05, na.rm = TRUE), n_ok),
              stringsAsFactors = FALSE
            )
            # Dedup: replace prior run of same cell
            existing <- run_results_rv()
            if (nrow(existing) > 0) {
              dup <- existing$n == n_val &
                     abs(existing$phi - phi_val) < 1e-9 &
                     existing$innov_dist == innov_dist_str
              existing <- existing[!dup, , drop = FALSE]
            }
            run_results_rv(rbind(existing, cell_rates))

            # Store raw results for detail sub-tabs
            cell_key <- paste(n_val, phi_val, innov_dist_str, sep = "|")
            raw <- raw_results_rv()
            raw[[cell_key]] <- list(
              results      = results_list,
              n            = n_val,
              phi          = phi_val,
              innov_dist   = innov_dist_str,
              innov_label  = innov_label,
              innov_params = innov_params,
              timestamp    = Sys.time()
            )
            raw_results_rv(raw)
          }

          # Write cell to DB (if saving)
          if (!is.null(write_con)) {
            tryCatch({
              mc_db_write_batch(write_con, results_list, batch_id,
                                n = n_val, phi = phi_val,
                                innov_dist = innov_dist_str)
              .log(sprintf("  Saved %d/%d sims in %.1f s%s",
                           length(results_list), nsims, elapsed_cell,
                           if (failed_count > 0)
                             sprintf(" (%d failed)", failed_count) else ""))
              cells_ok <- cells_ok + 1L
            }, error = function(e) {
              .log(sprintf("  WRITE ERROR: %s", e$message))
              cells_err <<- cells_err + 1L
            })
          } else {
            .log(sprintf("  Completed %d/%d sims in %.1f s (not saved)%s",
                         length(results_list), nsims, elapsed_cell,
                         if (failed_count > 0)
                           sprintf(" (%d failed)", failed_count) else ""))
            cells_ok <- cells_ok + 1L
          }

          incProgress(1 / n_cells,
                      detail = sprintf("Completed %d/%d", cell_idx, n_cells))
        }
      })

      total_elapsed <- proc.time()[["elapsed"]] - t_total_start

      if (save_to_db) {
        .log(sprintf("Grid complete: %d OK, %d errors, %.1f s total (batch #%d)",
                      cells_ok, cells_err, total_elapsed, batch_id))
        # Trigger DB refresh for all tabs
        db_refresh_trigger(db_refresh_trigger() + 1L)
        showNotification(
          sprintf("Completed %d/%d cells. Batch #%d saved. (%.0f s)",
                  cells_ok, n_cells, batch_id, total_elapsed),
          type = "message", duration = 8)
      } else {
        .log(sprintf("Grid complete: %d OK, %d errors, %.1f s total (not saved)",
                      cells_ok, cells_err, total_elapsed))
        showNotification(
          sprintf("Completed %d/%d cells in %.0f s (results not saved to DB)",
                  cells_ok, n_cells, total_elapsed),
          type = "message", duration = 8)
      }
    })

    # =========================================================================
    # Run status & log
    # =========================================================================

    output$run_status_ui <- renderUI({
      grid <- grid_rv()
      n_cells <- nrow(grid)
      if (n_cells == 0) {
        return(p(class = "text-body-secondary", "No cells in grid."))
      }
      n_sel <- length(input$grid_table_rows_selected)
      to_run <- if (input$run_mode == "selected") n_sel else n_cells
      nsims <- input$nsims %||% 1000
      tagList(
        p(class = "text-body-secondary",
          sprintf("%d cells in grid, %d to run (%s sims each)",
                  n_cells, to_run, format(nsims, big.mark = ",")))
      )
    })

    output$run_log <- renderText({
      log_text <- run_log_rv()
      if (!nzchar(log_text)) "No simulation runs yet." else log_text
    })

    # =========================================================================
    # Combined results: in-memory + DB (in-memory takes priority)
    # =========================================================================

    cap_combined_results <- reactive({
      mem <- run_results_rv()
      db_refresh_trigger()

      # Try to read DB rates (shared con -> save_path fallback)
      db_rates <- data.frame()
      read_con <- con
      own_con <- FALSE
      if (is.null(read_con)) {
        save_path <- trimws(input$db_save_path %||% "")
        if (nzchar(save_path) && file.exists(save_path)) {
          read_con <- tryCatch(mc_db_connect(save_path, read_only = TRUE),
                               error = function(e) NULL)
          own_con <- TRUE
        }
      }
      if (!is.null(read_con)) {
        if (own_con) on.exit(DBI::dbDisconnect(read_con, shutdown = TRUE))
        db_rates <- tryCatch(
          DBI::dbGetQuery(read_con, "SELECT * FROM v_rejection_rates"),
          error = function(e) data.frame())
      }

      if (nrow(db_rates) == 0) return(mem)

      # Align DB columns to in-memory schema
      db_rates$innov_label <- db_rates$innov_dist
      target_cols <- c("n", "phi", "innov_dist", "innov_label", "n_sims",
                       "reject_05", "reject_05_se",
                       "reject_asymp_05", "reject_asymp_05_se",
                       "reject_adj_05", "reject_adj_05_se")
      db_rates <- db_rates[, intersect(target_cols, names(db_rates)), drop = FALSE]

      if (nrow(mem) == 0) return(db_rates)

      # In-memory takes priority: remove DB rows that overlap
      keep_db <- !vapply(seq_len(nrow(db_rates)), function(i) {
        any(mem$n == db_rates$n[i] &
            abs(mem$phi - db_rates$phi[i]) < 1e-9 &
            mem$innov_dist == db_rates$innov_dist[i])
      }, logical(1))
      db_extra <- db_rates[keep_db, , drop = FALSE]
      if (nrow(db_extra) > 0) rbind(mem, db_extra) else mem
    })

    # =========================================================================
    # Unified data accessors for detail sub-tabs
    # =========================================================================

    # Return list of wbg_boot_fast results for a cell_key "n|phi|innov_dist".
    # Checks in-memory first, falls back to DB query.
    cap_sim_data <- function(cell_key) {
      if (is.null(cell_key) || !nzchar(cell_key)) return(NULL)
      # Check in-memory first
      raw <- raw_results_rv()
      if (cell_key %in% names(raw)) return(raw[[cell_key]]$results)

      # Fall back to DB
      parts <- strsplit(cell_key, "\\|")[[1]]
      if (length(parts) != 3) return(NULL)
      read_con <- con
      own_con <- FALSE
      if (is.null(read_con)) {
        save_path <- trimws(input$db_save_path %||% "")
        if (nzchar(save_path) && file.exists(save_path)) {
          read_con <- tryCatch(mc_db_connect(save_path, read_only = TRUE),
                               error = function(e) NULL)
          own_con <- TRUE
        }
      }
      if (is.null(read_con)) return(NULL)
      if (own_con) on.exit(DBI::dbDisconnect(read_con, shutdown = TRUE))

      db_sims <- tryCatch(
        DBI::dbGetQuery(read_con,
          "SELECT obs_stat, boot_dist, pvalue, pvalue_upper, pvalue_lower,
                  pvalue_asymp, pvalue_adj, null_ar_order, null_ar_phi, null_vara
           FROM simulations WHERE n = ? AND phi = ? AND innov_dist = ?",
          params = list(as.integer(parts[1]), as.numeric(parts[2]), parts[3])),
        error = function(e) data.frame())
      if (nrow(db_sims) == 0) return(NULL)

      # Convert DB rows to wbg_boot_fast-like lists
      lapply(seq_len(nrow(db_sims)), function(i) {
        list(
          tco_obs       = db_sims$obs_stat[i],
          boot_tstats   = db_sims$boot_dist[[i]],
          pvalue        = db_sims$pvalue[i],
          pvalue_upper  = db_sims$pvalue_upper[i],
          pvalue_lower  = db_sims$pvalue_lower[i],
          pvalue_asymp  = db_sims$pvalue_asymp[i],
          pvalue_adj    = db_sims$pvalue_adj[i],
          p             = as.integer(db_sims$null_ar_order[i]),
          phi           = db_sims$null_ar_phi[[i]],
          vara          = db_sims$null_vara[i]
        )
      })
    }

    # Available cell keys from both memory and DB
    cap_cell_choices <- reactive({
      raw <- raw_results_rv()
      mem_keys <- names(raw)

      db_keys <- character(0)
      read_con <- con
      own_con <- FALSE
      if (is.null(read_con)) {
        save_path <- trimws(input$db_save_path %||% "")
        if (nzchar(save_path) && file.exists(save_path)) {
          read_con <- tryCatch(mc_db_connect(save_path, read_only = TRUE),
                               error = function(e) NULL)
          own_con <- TRUE
        }
      }
      if (!is.null(read_con)) {
        if (own_con) on.exit(DBI::dbDisconnect(read_con, shutdown = TRUE))
        db_configs <- tryCatch(
          DBI::dbGetQuery(read_con,
            "SELECT DISTINCT n, phi, innov_dist FROM simulations
             ORDER BY innov_dist, n, phi"),
          error = function(e) data.frame())
        if (nrow(db_configs) > 0) {
          db_keys <- paste(db_configs$n, db_configs$phi,
                           db_configs$innov_dist, sep = "|")
        }
      }

      all_keys <- unique(c(mem_keys, db_keys))
      if (length(all_keys) == 0) return(character(0))
      labels <- gsub("\\|", ", ", all_keys)
      setNames(all_keys, labels)
    })

    # Return raw entry metadata for a cell (innov_label, innov_params, etc.)
    cap_cell_meta <- function(cell_key) {
      raw <- raw_results_rv()
      if (cell_key %in% names(raw)) {
        entry <- raw[[cell_key]]
        return(list(innov_label = entry$innov_label,
                    innov_params = entry$innov_params,
                    n = entry$n, phi = entry$phi,
                    innov_dist = entry$innov_dist))
      }
      # DB fallback: parse key
      parts <- strsplit(cell_key, "\\|")[[1]]
      if (length(parts) == 3) {
        return(list(innov_label = parts[3], innov_params = list(),
                    n = as.integer(parts[1]), phi = as.numeric(parts[2]),
                    innov_dist = parts[3]))
      }
      NULL
    }

    # =========================================================================
    # Sub-tab server fragments
    # =========================================================================

    # Reactive wrapper for db_save_path input (needed by overview module)
    db_save_path_reactive <- reactive({ input$db_save_path })

    # Tab 2: Coverage / Overview
    mod_capstone_overview_server(input, output, session,
      con = con, grid_rv = grid_rv, raw_results_rv = raw_results_rv,
      db_refresh_trigger = db_refresh_trigger,
      db_save_path_input = db_save_path_reactive)

    # Tab 3: Individual Results (returns selected row for Bootstrap tab)
    selected_sim_rv <- mod_capstone_individual_server(input, output, session,
      cap_sim_data = cap_sim_data, cap_cell_choices = cap_cell_choices,
      cap_combined_results = cap_combined_results)

    # Tab 4: Analysis Grids
    mod_capstone_grids_server(input, output, session,
      cap_combined_results = cap_combined_results,
      run_results_rv = run_results_rv)

    # Tab 5: Replication Comparison
    mod_capstone_replication_server(input, output, session,
      preset_name          = reactive(input$preset),
      cap_combined_results = cap_combined_results,
      grid_rv              = grid_rv,
      test_params          = reactive(list(
        nsims     = as.integer(input$nsims %||% 1000),
        nb        = as.integer(input$nb %||% 399),
        maxp      = min(as.integer(input$maxp %||% 5), 20L),
        criterion = input$criterion %||% "aic",
        bootadj   = isTRUE(input$bootadj),
        min_p     = if (isTRUE(input$minp1)) 1L else 0L,
        base_seed = input$seed
      )))

    # Tab 6: Plots
    mod_capstone_plots_server(input, output, session,
      cap_combined_results = cap_combined_results)

    # Tab 6: P-value Diagnostics
    mod_capstone_pvalue_server(input, output, session,
      cap_sim_data = cap_sim_data, cap_cell_choices = cap_cell_choices)

    # Tab 7: Bootstrap Distribution
    mod_capstone_bootdist_server(input, output, session,
      cap_sim_data = cap_sim_data, cap_cell_choices = cap_cell_choices,
      selected_sim_from_individual = selected_sim_rv)

    # Tab 8: Innovation Diagnostics
    mod_capstone_innov_diag_server(input, output, session,
      cap_cell_choices = cap_cell_choices, cap_cell_meta = cap_cell_meta)

    # Tab 9: Null Model Diagnostics
    mod_capstone_null_diag_server(input, output, session,
      cap_sim_data = cap_sim_data, cap_cell_choices = cap_cell_choices,
      con = con, db_refresh_trigger = db_refresh_trigger)

  })
}
