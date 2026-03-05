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
           params = list(hetero_shape = "linear", hetero_from = 1,
                         hetero_to = 10, hetero_sd = 1)),
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

.cap_picker_input <- function(inputId, label, choices, selected = NULL, options = list()) {
  if (requireNamespace("shinyWidgets", quietly = TRUE)) {
    shinyWidgets::pickerInput(
      inputId = inputId,
      label = label,
      choices = choices,
      selected = selected,
      options = options
    )
  } else {
    selectInput(inputId = inputId, label = label, choices = choices, selected = selected)
  }
}

.cap_toggle_input <- function(inputId, label, value = FALSE) {
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

.cap_action_button <- function(inputId, label, icon_name = NULL,
                               class = "btn-primary", full_width = FALSE) {
  btn_icon <- if (is.null(icon_name)) NULL else icon(icon_name)
  if (requireNamespace("shinyWidgets", quietly = TRUE)) {
    btn_color <- if (grepl("danger", class, fixed = TRUE)) {
      "danger"
    } else if (grepl("success", class, fixed = TRUE)) {
      "success"
    } else if (grepl("secondary", class, fixed = TRUE)) {
      "default"
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

.cap_radio_buttons <- function(inputId, label, choices, selected = NULL, inline = TRUE) {
  if (requireNamespace("shinyWidgets", quietly = TRUE)) {
    shinyWidgets::radioGroupButtons(
      inputId = inputId,
      label = label,
      choices = choices,
      selected = selected,
      justified = !isTRUE(inline),
      direction = if (isTRUE(inline)) "horizontal" else "vertical",
      size = "sm",
      checkIcon = list(yes = icon("check"))
    )
  } else {
    radioButtons(inputId = inputId, label = label, choices = choices,
                 selected = selected, inline = isTRUE(inline))
  }
}

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
    div(id = ns("cap_content"), fluidRow(
      # --- Sidebar (col-3) ---
      column(3,
        # Panel 1: Run Controls (top for visibility)
        wellPanel(
          h4("Run"),
          .cap_action_button(ns("run_grid"), "Run Grid Simulation",
                             icon_name = "play", class = "btn-primary",
                             full_width = TRUE),
          hr(),
          p(class = "text-body-secondary", style = "font-size: 0.85em; margin-bottom: 8px;",
            "Use the toast cancel icon to stop at scenario boundaries."),
          .cap_radio_buttons(ns("run_mode"), NULL,
                             choices = c("All cells" = "all",
                                         "Selected rows" = "selected"),
                             selected = "all", inline = TRUE)
        ),

        # Panel 2: Preset & Grid
        wellPanel(
          h4("Grid Preset"),
          .cap_picker_input(ns("preset"), "Load Preset",
                            choices = c("Custom", names(.CAPSTONE_PRESETS)),
                            selected = "Custom"),
          .cap_action_button(ns("load_preset"), "Load Preset",
                             icon_name = "upload", class = "btn-sm btn-secondary",
                             full_width = TRUE),
          hr(),
          h5("Add Config Row"),
          fluidRow(
            column(6, numericInput(ns("new_phi"), "Phi", 0.95,
                                   min = 0, max = 0.999, step = 0.01)),
            column(6, numericInput(ns("new_n"), "n", 200,
                                   min = 10, max = 5000, step = 10))
          ),
          .cap_picker_input(ns("new_innov"), "Innovation",
                            choices = .CAP_INNOV_CHOICES, selected = "Normal"),

          # Distribution-specific params (minimal set for custom rows)
          conditionalPanel(
            condition = "input.new_innov == 'Normal'", ns = ns,
            numericInput(ns("new_norm_sd"), "sd", value = 1, min = 0.001, step = 0.1)
          ),
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
                      value = "0.2, 0.175, 0.15, 0.125, 0.1, 0.075, 0.05, 0.025"),
            textInput(ns("new_garch_beta"), "Beta (comma-separated, optional)",
                      value = "")
          ),
          conditionalPanel(
            condition = "input.new_innov == 'Heteroscedastic'", ns = ns,
            .cap_picker_input(ns("new_hetero_shape"), "Shape",
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
            column(6, .cap_action_button(ns("add_row"), "Add Row",
                                         icon_name = "plus", class = "btn-sm btn-success")),
            column(6, .cap_action_button(ns("clear_grid"), "Clear All",
                                         icon_name = "trash", class = "btn-sm btn-outline-danger"))
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
          .cap_picker_input(ns("criterion"), "Criterion",
                            choices = c("aic", "aicc", "bic"), selected = "aic"),
          .cap_toggle_input(ns("bootadj"), "COBA adjustment", value = TRUE),
          .cap_toggle_input(ns("minp1"), "Min AR order = 1", value = TRUE),
          .cap_toggle_input(ns("use_fast"), "Use fast innovation generators", value = FALSE),
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
              column(6, .cap_action_button(ns("delete_selected"), "Delete Selected Rows",
                                           icon_name = "minus", class = "btn-sm btn-outline-danger"))
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
          # Tab 6: P-value Diagnostics
          mod_capstone_pvalue_ui(ns),
          # Tab 7: Bootstrap Distribution
          mod_capstone_bootdist_ui(ns),
          # Tab 8: Null Model Diagnostics
          mod_capstone_null_diag_ui(ns),

          # Tab 9: Run Log
          tabPanel("Run Log",
            br(),
            verbatimTextOutput(ns("run_log"))
          )
        )
      )
    ))
  )
}

# =============================================================================
# Server
# =============================================================================

mod_capstone_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    has_waiter <- requireNamespace("waiter", quietly = TRUE)
    has_shinyvalidate <- requireNamespace("shinyvalidate", quietly = TRUE)
    has_shinywidgets <- requireNamespace("shinyWidgets", quietly = TRUE)

    .update_picker <- function(input_id, selected = NULL, choices = NULL) {
      if (has_shinywidgets) {
        shinyWidgets::updatePickerInput(session, input_id,
          selected = selected, choices = choices)
      } else {
        updateSelectInput(session, input_id,
          selected = selected, choices = choices)
      }
    }

    .update_toggle <- function(input_id, value) {
      if (has_shinywidgets) {
        shinyWidgets::updateMaterialSwitch(session, input_id, value = isTRUE(value))
      } else {
        updateCheckboxInput(session, input_id, value = isTRUE(value))
      }
    }

    .with_cap_waiter <- function(expr) {
      if (!has_waiter) return(force(expr))
      w <- waiter::Waiter$new(
        id = ns("cap_content"),
        html = waiter::spin_fading_circles(),
        color = "rgba(0, 0, 0, 0.25)"
      )
      w$show()
      on.exit(w$hide(), add = TRUE)
      force(expr)
    }

    .has_nonblank_value <- function(x) {
      if (is.null(x)) return(FALSE)
      any(nzchar(trimws(as.character(x))))
    }

    .parse_num_list_or_null <- function(x) {
      if (is.null(x)) return(NULL)

      txt <- trimws(as.character(x))
      txt <- txt[nzchar(txt)]
      if (length(txt) == 0) return(NULL)

      if (length(txt) == 1 && grepl(",", txt, fixed = TRUE)) {
        txt <- trimws(strsplit(txt, "\\s*,\\s*")[[1]])
      }

      vals <- suppressWarnings(as.numeric(txt))
      vals <- vals[!is.na(vals)]
      if (length(vals) == 0) return(NULL)
      vals
    }

    .validate_add_row <- function(innov_label, params, phi_val, n_val) {
      errs <- character(0)
      if (is.na(phi_val) || phi_val < 0 || phi_val >= 1) {
        errs <- c(errs, "Phi must be in [0, 1).")
      }
      if (is.na(n_val) || n_val < 10) {
        errs <- c(errs, "n must be at least 10.")
      }
      if (innov_label == "Student's t" && (is.na(params$t_df) || params$t_df <= 0)) {
        errs <- c(errs, "Student's t df must be positive.")
      }
      if (innov_label == "Normal" &&
          (!is.null(params$norm_sd)) &&
          (is.na(params$norm_sd) || params$norm_sd <= 0)) {
        errs <- c(errs, "Normal sd must be positive.")
      }
      if (innov_label == "Laplace" && (is.na(params$lap_scale) || params$lap_scale <= 0)) {
        errs <- c(errs, "Laplace scale must be positive.")
      }
      if (innov_label == "GARCH") {
        if (is.na(params$garch_omega) || params$garch_omega <= 0) {
          errs <- c(errs, "GARCH omega must be positive.")
        }
        alpha <- .parse_num_list(params$garch_alpha)
        if (length(alpha) == 0 || any(alpha < 0)) {
          errs <- c(errs, "GARCH alpha must be a comma-separated list of non-negative numbers.")
        }
        if (.has_nonblank_value(params$garch_beta)) {
          beta <- .parse_num_list_or_null(params$garch_beta)
          if (is.null(beta) || any(beta < 0)) {
            errs <- c(errs, "GARCH beta must be a comma-separated list of non-negative numbers.")
          }
        }
      }
      if (innov_label == "Heteroscedastic") {
        if (is.na(params$hetero_from) || params$hetero_from <= 0 ||
            is.na(params$hetero_to) || params$hetero_to <= 0) {
          errs <- c(errs, "Heteroscedastic from/to must be positive.")
        }
        if (!is.na(params$hetero_from) && !is.na(params$hetero_to) &&
            params$hetero_to < params$hetero_from) {
          errs <- c(errs, "Heteroscedastic 'To' must be >= 'From'.")
        }
      }
      if (innov_label == "Mixture Normal") {
        if (is.na(params$mix_sd1) || params$mix_sd1 <= 0 ||
            is.na(params$mix_sd2) || params$mix_sd2 <= 0) {
          errs <- c(errs, "Mixture sd1/sd2 must be positive.")
        }
        if (is.na(params$mix_prob1) || params$mix_prob1 <= 0 || params$mix_prob1 >= 1) {
          errs <- c(errs, "Mixture prob1 must be strictly between 0 and 1.")
        }
      }
      if (innov_label == "Skew-t" && (is.na(params$skt_df) || params$skt_df <= 0)) {
        errs <- c(errs, "Skew-t df must be positive.")
      }
      if (innov_label == "GED" && (is.na(params$ged_nu) || params$ged_nu <= 0)) {
        errs <- c(errs, "GED nu must be positive.")
      }
      if (innov_label == "Uniform" && (is.na(params$unif_hw) || params$unif_hw <= 0)) {
        errs <- c(errs, "Uniform half-width must be positive.")
      }
      errs
    }

    .validate_run_params <- function(nsims, nb, maxp, seed) {
      errs <- character(0)
      if (is.na(nsims) || nsims < 1) errs <- c(errs, "Simulations per cell must be at least 1.")
      if (is.na(nb) || nb < 1) errs <- c(errs, "Bootstrap replicates must be at least 1.")
      if (is.na(maxp) || maxp < 1 || maxp > 20) errs <- c(errs, "Max AR order must be between 1 and 20.")
      if (!is.null(seed) && !is.na(seed) && seed < 1) errs <- c(errs, "Base seed must be positive when provided.")
      errs
    }

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
    run_status_tick <- reactiveVal(0L)
    run_state <- new.env(parent = emptyenv())
    run_state$active <- FALSE
    run_state$cancel_requested <- FALSE
    run_state$ctx <- NULL
    run_state$waiter <- NULL
    run_state$toast_id <- NULL
    run_state$status_seq <- 0L

    .bump_run_status <- function() {
      run_state$status_seq <- run_state$status_seq + 1L
      run_status_tick(run_state$status_seq)
    }

    v <- NULL
    if (has_shinyvalidate) {
      v <- shinyvalidate::InputValidator$new()
      v$add_rule("new_phi", function(value) {
        vv <- suppressWarnings(as.numeric(value))
        if (is.na(vv) || vv < 0 || vv >= 1) "Must be in [0, 1)" else NULL
      })
      v$add_rule("new_n", function(value) {
        vv <- suppressWarnings(as.integer(value))
        if (is.na(vv) || vv < 10) "Must be >= 10" else NULL
      })
      v$add_rule("nsims", function(value) {
        vv <- suppressWarnings(as.integer(value))
        if (is.na(vv) || vv < 1) "Must be >= 1" else NULL
      })
      v$add_rule("nb", function(value) {
        vv <- suppressWarnings(as.integer(value))
        if (is.na(vv) || vv < 1) "Must be >= 1" else NULL
      })
      v$add_rule("maxp", function(value) {
        vv <- suppressWarnings(as.integer(value))
        if (is.na(vv) || vv < 1 || vv > 20) "Must be 1..20" else NULL
      })
      v$add_rule("seed", function(value) {
        if (is.null(value) || is.na(value) || !nzchar(as.character(value))) return(NULL)
        vv <- suppressWarnings(as.integer(value))
        if (is.na(vv) || vv < 1) "Must be positive" else NULL
      })

      # Distribution-specific add-row rules
      v$add_rule("new_norm_sd", function(value) {
        if (!identical(input$new_innov, "Normal")) return(NULL)
        vv <- suppressWarnings(as.numeric(value))
        if (is.na(vv) || vv <= 0) "sd must be positive" else NULL
      })
      v$add_rule("new_t_df", function(value) {
        if (!identical(input$new_innov, "Student's t")) return(NULL)
        vv <- suppressWarnings(as.numeric(value))
        if (is.na(vv) || vv <= 0) "df must be positive" else NULL
      })
      v$add_rule("new_lap_scale", function(value) {
        if (!identical(input$new_innov, "Laplace")) return(NULL)
        vv <- suppressWarnings(as.numeric(value))
        if (is.na(vv) || vv <= 0) "Scale must be positive" else NULL
      })
      v$add_rule("new_garch_omega", function(value) {
        if (!identical(input$new_innov, "GARCH")) return(NULL)
        vv <- suppressWarnings(as.numeric(value))
        if (is.na(vv) || vv <= 0) "Omega must be positive" else NULL
      })
      v$add_rule("new_garch_alpha", function(value) {
        if (!identical(input$new_innov, "GARCH")) return(NULL)
        vals <- .parse_num_list(value)
        if (length(vals) == 0) return("Enter at least one alpha value")
        if (any(vals < 0)) return("Alpha values must be non-negative")
        NULL
      })
      v$add_rule("new_garch_beta", function(value) {
        if (!identical(input$new_innov, "GARCH")) return(NULL)
        txt <- trimws(as.character(value %||% ""))
        if (!nzchar(txt)) return(NULL)
        vals <- .parse_num_list(txt)
        if (length(vals) == 0) return("Enter comma-separated beta values")
        if (any(vals < 0)) return("Beta values must be non-negative")
        NULL
      })
      v$add_rule("new_hetero_from", function(value) {
        if (!identical(input$new_innov, "Heteroscedastic")) return(NULL)
        vv <- suppressWarnings(as.numeric(value))
        if (is.na(vv) || vv <= 0) "From must be positive" else NULL
      })
      v$add_rule("new_hetero_to", function(value) {
        if (!identical(input$new_innov, "Heteroscedastic")) return(NULL)
        to_v <- suppressWarnings(as.numeric(value))
        from_v <- suppressWarnings(as.numeric(input$new_hetero_from))
        if (is.na(to_v) || to_v <= 0) return("To must be positive")
        if (!is.na(from_v) && to_v < from_v) return("To must be >= From")
        NULL
      })
      v$add_rule("new_mix_sd1", function(value) {
        if (!identical(input$new_innov, "Mixture Normal")) return(NULL)
        vv <- suppressWarnings(as.numeric(value))
        if (is.na(vv) || vv <= 0) "sd1 must be positive" else NULL
      })
      v$add_rule("new_mix_sd2", function(value) {
        if (!identical(input$new_innov, "Mixture Normal")) return(NULL)
        vv <- suppressWarnings(as.numeric(value))
        if (is.na(vv) || vv <= 0) "sd2 must be positive" else NULL
      })
      v$add_rule("new_mix_prob1", function(value) {
        if (!identical(input$new_innov, "Mixture Normal")) return(NULL)
        vv <- suppressWarnings(as.numeric(value))
        if (is.na(vv) || vv <= 0 || vv >= 1) "prob1 must be in (0,1)" else NULL
      })
      v$add_rule("new_skt_df", function(value) {
        if (!identical(input$new_innov, "Skew-t")) return(NULL)
        vv <- suppressWarnings(as.numeric(value))
        if (is.na(vv) || vv <= 0) "df must be positive" else NULL
      })
      v$add_rule("new_ged_nu", function(value) {
        if (!identical(input$new_innov, "GED")) return(NULL)
        vv <- suppressWarnings(as.numeric(value))
        if (is.na(vv) || vv <= 0) "nu must be positive" else NULL
      })
      v$add_rule("new_unif_hw", function(value) {
        if (!identical(input$new_innov, "Uniform")) return(NULL)
        vv <- suppressWarnings(as.numeric(value))
        if (is.na(vv) || vv <= 0) "Half-width must be positive" else NULL
      })

      v$enable()
    }

    # =========================================================================
    # Preset loading
    # =========================================================================

    observeEvent(input$load_preset, {
      .with_cap_waiter({
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
      if (!is.null(defs$criterion)) .update_picker("criterion", selected = defs$criterion)
      if (!is.null(defs$bootadj))   .update_toggle("bootadj", value = defs$bootadj)
      if (!is.null(defs$min_p))     .update_toggle("minp1", value = defs$min_p == 1L)
      if (!is.null(defs$seed))      updateNumericInput(session, "seed", value = defs$seed)

      showNotification(
        sprintf("Loaded '%s': %d cells", preset_name, nrow(grid_df)),
        type = "message", duration = 4)
      })
    })

    # =========================================================================
    # Add / Delete / Clear rows
    # =========================================================================

    observeEvent(input$add_row, {
      .with_cap_waiter({
      if (!is.null(v) && !isTRUE(v$is_valid())) {
        showNotification("Please fix highlighted input errors before adding a row.",
                         type = "warning", duration = 4)
        return()
      }
      phi_val <- input$new_phi
      n_val <- as.integer(input$new_n)
      innov_label <- input$new_innov

      if (is.null(phi_val) || is.na(phi_val)) phi_val <- 0.95
      if (is.null(n_val) || is.na(n_val)) n_val <- 200L

      # Build params from the add-row form
      params <- switch(innov_label,
        "Normal"      = list(norm_sd = input$new_norm_sd %||% 1),
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
                               "0.2,0.175,0.15,0.125,0.1,0.075,0.05,0.025",
                             garch_beta = input$new_garch_beta %||% ""),
        "Heteroscedastic" = list(hetero_shape = input$new_hetero_shape %||% "linear",
                                   hetero_from = input$new_hetero_from %||% 1,
                                   hetero_to = input$new_hetero_to %||% 10),
        list()
      )

      errs <- .validate_add_row(innov_label, params, phi_val, n_val)
      if (length(errs) > 0) {
        showNotification(paste(errs, collapse = " "), type = "error", duration = 6)
        return()
      }

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

      showNotification("Row added to grid.", type = "message", duration = 2)
      })
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

      .compact_vec <- function(x, max_items = 4L) {
        if (is.function(x)) return("<function>")
        if (is.null(x)) return("NULL")
        if (is.character(x) && length(x) == 1L) {
          txt <- trimws(x)
          if (grepl(",", txt, fixed = TRUE)) {
            vals <- suppressWarnings(as.numeric(trimws(strsplit(txt, ",")[[1]])))
            if (!any(is.na(vals))) x <- vals
          }
        }
        if (is.numeric(x)) {
          vals <- signif(as.numeric(x), 4)
          if (length(vals) <= max_items) {
            return(paste0("[", paste(vals, collapse = ","), "]"))
          }
          return(paste0("[", paste(vals[seq_len(max_items)], collapse = ","),
                        ",+", length(vals) - max_items, "]"))
        }
        txt <- paste(deparse(x, width.cutoff = 120L), collapse = "")
        if (nchar(txt) > 60) txt <- paste0(substr(txt, 1, 57), "...")
        txt
      }

      .scenario_params_full <- function(label, params) {
        p <- params
        if (is.null(p)) p <- list()
        switch(label,
          "Normal" = sprintf("sd=%s", p$norm_sd %||% 1),
          "Student's t" = sprintf("df=%s, scale=%s", p$t_df %||% 3, isTRUE(p$t_scale)),
          "Skew-t" = sprintf("df=%s, alpha=%s, scale=%s",
                              p$skt_df %||% 5, p$skt_alpha %||% 0, isTRUE(p$skt_scale)),
          "GED" = sprintf("nu=%s, sd=%s", p$ged_nu %||% 2, p$ged_sd %||% 1),
          "Laplace" = sprintf("scale=%s", p$lap_scale %||% (1 / sqrt(2))),
          "Uniform" = sprintf("half_width=%s", p$unif_hw %||% sqrt(3)),
          "Mixture Normal" = sprintf("sd1=%s, sd2=%s, prob1=%s",
                                      p$mix_sd1 %||% 1, p$mix_sd2 %||% 3, p$mix_prob1 %||% 0.9),
          "GARCH" = {
            alpha <- .compact_vec(p$garch_alpha %||% c(0.2, 0.175, 0.15, 0.125, 0.1, 0.075, 0.05, 0.025), max_items = 8L)
            beta <- if (is.null(p$garch_beta) || identical(trimws(as.character(p$garch_beta)), "")) {
              "NULL"
            } else {
              .compact_vec(p$garch_beta, max_items = 8L)
            }
            sprintf("omega=%s, alpha=%s, beta=%s", p$garch_omega %||% 0.1, alpha, beta)
          },
          "Heteroscedastic" = {
            if (!is.null(p$hetero_w)) {
              sprintf("w=%s, sd=%s", .compact_vec(p$hetero_w), p$hetero_sd %||% 1)
            } else {
              shape <- p$hetero_shape %||% "linear"
              if (shape %in% c("linear", "sqrt", "log", "exp")) {
                sprintf("shape=%s, from=%s, to=%s, sd=%s",
                        shape, p$hetero_from %||% 1, p$hetero_to %||% 10, p$hetero_sd %||% 1)
              } else if (shape == "power") {
                sprintf("shape=power, from=%s, to=%s, power=%s, sd=%s",
                        p$hetero_from %||% 1, p$hetero_to %||% 10, p$hetero_power %||% 2, p$hetero_sd %||% 1)
              } else if (shape == "step") {
                sprintf("shape=step, breaks=%s, levels=%s, sd=%s",
                        .compact_vec(p$hetero_breaks %||% "0.5"),
                        .compact_vec(p$hetero_levels %||% "1,5"),
                        p$hetero_sd %||% 1)
              } else {
                sprintf("shape=periodic, base_w=%s, amplitude=%s, period=%s, sd=%s",
                        p$hetero_base_w %||% 1, p$hetero_amplitude %||% 0.5,
                        p$hetero_period %||% 12, p$hetero_sd %||% 1)
              }
            }
          },
          {
            if (length(p) == 0) return("default")
            paste(sprintf("%s=%s", names(p), vapply(p, .compact_vec, character(1))), collapse = ", ")
          }
        )
      }

      .scenario_params_short <- function(full_txt, max_chars = 72L) {
        out <- full_txt
        out <- gsub("omega", "om", out, fixed = TRUE)
        out <- gsub("alpha", "a", out, fixed = TRUE)
        out <- gsub("beta", "b", out, fixed = TRUE)
        out <- gsub("amplitude", "amp", out, fixed = TRUE)
        out <- gsub("period", "per", out, fixed = TRUE)
        out <- gsub("half_width", "hw", out, fixed = TRUE)
        out <- gsub("shape", "sh", out, fixed = TRUE)
        out <- gsub("scale", "sc", out, fixed = TRUE)
        if (nchar(out) > max_chars) {
          out <- paste0(substr(out, 1, max_chars - 3L), "...")
        }
        out
      }

      params_full <- mapply(
        .scenario_params_full,
        label = grid$innov_label,
        params = grid$innov_params,
        SIMPLIFY = TRUE,
        USE.NAMES = FALSE
      )
      params_short <- vapply(params_full, .scenario_params_short, character(1))
      params_cell <- sprintf(
        "<span title=\"%s\">%s</span>",
        htmltools::htmlEscape(params_full, attribute = TRUE),
        htmltools::htmlEscape(params_short)
      )

      display <- data.frame(
        Phi     = grid$phi,
        n       = grid$n,
        Innovation = grid$innov_label,
        `Scenario Key` = grid$innov_dist_str,
        `Scenario Params` = params_cell,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )

      datatable(
        display,
        rownames = FALSE,
        escape = c(1, 2, 3, 4),
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

    .log_line <- function(msg) {
      ctx <- run_state$ctx
      if (is.null(ctx)) return()
      ctx$log_lines <- c(ctx$log_lines,
                         sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), msg))
      run_state$ctx <- ctx
      run_log_rv(paste(ctx$log_lines, collapse = "\n"))
      .bump_run_status()
    }

    .ctx_with_latest_logs <- function(ctx) {
      current <- run_state$ctx
      if (!is.null(current) && !is.null(current$log_lines)) {
        ctx$log_lines <- current$log_lines
      }
      ctx
    }

    .run_toast_ui <- function(ctx, detail = NULL) {
      n_cells <- max(1L, as.integer(ctx$n_cells))
      done <- max(0L, min(n_cells, as.integer(ctx$cell_idx) - 1L))
      pct <- round(100 * done / n_cells)
      status_line <- sprintf("Running grid simulation... Cell %d/%d", min(as.integer(ctx$cell_idx), n_cells), n_cells)
      if (isTRUE(run_state$cancel_requested)) {
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

    .show_run_toast <- function(detail = NULL) {
      if (!isTRUE(run_state$active) || is.null(run_state$ctx)) return(invisible(NULL))
      if (is.null(run_state$toast_id)) run_state$toast_id <- paste0(ns("run_progress"), "_toast")

      showNotification(
        ui = .run_toast_ui(run_state$ctx, detail = detail),
        id = run_state$toast_id,
        duration = NULL,
        closeButton = FALSE,
        type = "message",
        action = actionLink(
          inputId = ns("cancel_run_toast"),
          label = icon("ban"),
          style = "padding: 0; color: var(--bs-warning); text-decoration: none;"
        ),
        session = session
      )
      invisible(NULL)
    }

    .fmt_log_value <- function(x) {
      if (is.function(x)) return("<function>")
      txt <- paste(deparse(x, width.cutoff = 500L), collapse = "")
      if (nchar(txt) > 240) {
        txt <- paste0(substr(txt, 1, 237), "...")
      }
      txt
    }

    .parse_num_list <- function(x) {
      if (is.null(x)) return(NULL)
      if (is.character(x)) {
        vals <- suppressWarnings(as.numeric(trimws(strsplit(x, ",")[[1]])))
      } else {
        vals <- suppressWarnings(as.numeric(x))
      }
      if (length(vals) == 0 || any(is.na(vals))) return(NULL)
      vals
    }

    .format_make_gen_call <- function(innov_label, params, use_fast) {
      suffix <- if (isTRUE(use_fast)) "_fast" else ""
      switch(innov_label,
        "Normal" = {
          sd_val <- params$norm_sd
          if (is.null(sd_val) || is.na(sd_val) || sd_val <= 0) sd_val <- 1
          sprintf("make_gen_norm%s(sd = %s)", suffix, .fmt_log_value(sd_val))
        },
        "Student's t" = {
          df_val <- params$t_df
          if (is.null(df_val) || is.na(df_val) || df_val < 1) df_val <- 3
          sc_val <- isTRUE(params$t_scale)
          sprintf("make_gen_t%s(df = %s, scale = %s)",
                  suffix, .fmt_log_value(df_val), .fmt_log_value(sc_val))
        },
        "Skew-t" = {
          df_val <- params$skt_df
          if (is.null(df_val) || is.na(df_val) || df_val < 3) df_val <- 5
          alpha_val <- params$skt_alpha
          if (is.null(alpha_val) || is.na(alpha_val)) alpha_val <- 0
          sc_val <- isTRUE(params$skt_scale)
          sprintf("make_gen_skt%s(df = %s, alpha = %s, scale = %s)",
                  suffix, .fmt_log_value(df_val), .fmt_log_value(alpha_val), .fmt_log_value(sc_val))
        },
        "GED" = {
          nu_val <- params$ged_nu
          if (is.null(nu_val) || is.na(nu_val) || nu_val <= 0) nu_val <- 2
          sd_val <- params$ged_sd
          if (is.null(sd_val) || is.na(sd_val) || sd_val <= 0) sd_val <- 1
          sprintf("make_gen_ged%s(nu = %s, sd = %s)",
                  suffix, .fmt_log_value(nu_val), .fmt_log_value(sd_val))
        },
        "Laplace" = {
          sc <- params$lap_scale
          if (is.null(sc) || is.na(sc) || sc <= 0) sc <- 1 / sqrt(2)
          sprintf("make_gen_laplace%s(scale = %s)", suffix, .fmt_log_value(sc))
        },
        "Uniform" = {
          hw <- params$unif_hw
          if (is.null(hw) || is.na(hw) || hw <= 0) hw <- sqrt(3)
          sprintf("make_gen_unif%s(half_width = %s)", suffix, .fmt_log_value(hw))
        },
        "Mixture Normal" = {
          sd1 <- params$mix_sd1
          sd2 <- params$mix_sd2
          p1 <- params$mix_prob1
          if (is.null(sd1) || is.na(sd1) || sd1 <= 0) sd1 <- 1
          if (is.null(sd2) || is.na(sd2) || sd2 <= 0) sd2 <- 3
          if (is.null(p1) || is.na(p1) || p1 <= 0 || p1 >= 1) p1 <- 0.9
          sprintf("make_gen_mixnorm%s(sd1 = %s, sd2 = %s, prob1 = %s)",
                  suffix, .fmt_log_value(sd1), .fmt_log_value(sd2), .fmt_log_value(p1))
        },
        "GARCH" = {
          omega <- params$garch_omega
          if (is.null(omega) || is.na(omega) || omega <= 0) omega <- 0.1
          alpha <- .parse_num_list_or_null(params$garch_alpha)
          if (is.null(alpha) || length(alpha) == 0) {
            alpha <- c(0.2, 0.175, 0.15, 0.125, 0.1, 0.075, 0.05, 0.025)
          }
          beta <- .parse_num_list_or_null(params$garch_beta)
          if (is.null(beta)) {
            sprintf("make_gen_garch%s(omega = %s, alpha = %s, beta = NULL)",
                    suffix, .fmt_log_value(omega), .fmt_log_value(alpha))
          } else {
            sprintf("make_gen_garch%s(omega = %s, alpha = %s, beta = %s)",
                    suffix, .fmt_log_value(omega), .fmt_log_value(alpha), .fmt_log_value(beta))
          }
        },
        "Heteroscedastic" = {
          if (!is.null(params$hetero_w)) {
            sd_val <- params$hetero_sd
            if (is.null(sd_val) || is.na(sd_val) || sd_val <= 0) sd_val <- 1
            sprintf("make_gen_hetero%s(w = %s, sd = %s)",
                    suffix, .fmt_log_value(params$hetero_w), .fmt_log_value(sd_val))
          } else {
            sd_val <- params$hetero_sd
            if (is.null(sd_val) || is.na(sd_val) || sd_val <= 0) sd_val <- 1
            shape <- params$hetero_shape %||% "linear"
            if (shape %in% c("linear", "sqrt", "log", "exp")) {
              from_val <- params$hetero_from %||% 1
              to_val <- params$hetero_to %||% 10
              sprintf("make_gen_hetero%s(shape = %s, from = %s, to = %s, sd = %s)",
                      suffix, .fmt_log_value(shape), .fmt_log_value(from_val),
                      .fmt_log_value(to_val), .fmt_log_value(sd_val))
            } else if (shape == "power") {
              from_val <- params$hetero_from %||% 1
              to_val <- params$hetero_to %||% 10
              p_val <- params$hetero_power %||% 2
              sprintf("make_gen_hetero%s(shape = \"power\", from = %s, to = %s, power = %s, sd = %s)",
                      suffix, .fmt_log_value(from_val), .fmt_log_value(to_val),
                      .fmt_log_value(p_val), .fmt_log_value(sd_val))
            } else if (shape == "step") {
              brk <- .parse_num_list_or_null(params$hetero_breaks)
              lvl <- .parse_num_list_or_null(params$hetero_levels)
              if (is.null(brk)) brk <- 0.5
              if (is.null(lvl)) lvl <- c(1, 5)
              sprintf("make_gen_hetero%s(shape = \"step\", breaks = %s, levels = %s, sd = %s)",
                      suffix, .fmt_log_value(brk), .fmt_log_value(lvl), .fmt_log_value(sd_val))
            } else {
              bw <- params$hetero_base_w %||% 1
              amp <- params$hetero_amplitude %||% 0.5
              per <- params$hetero_period %||% 12
              sprintf("make_gen_hetero%s(shape = \"periodic\", base_w = %s, amplitude = %s, period = %s, sd = %s)",
                      suffix, .fmt_log_value(bw), .fmt_log_value(amp),
                      .fmt_log_value(per), .fmt_log_value(sd_val))
            }
          }
        },
        sprintf("make_gen_<unknown>(label = %s, params = %s)",
                .fmt_log_value(innov_label), .fmt_log_value(params))
      )
    }

    .log_cell_calls <- function(ctx, phi_val, n_val, innov_label, innov_params, use_vara, vara) {
      .log_line(sprintf("  Calls for Cell %d/%d:", ctx$cell_idx, ctx$n_cells))
      if (isTRUE(use_vara)) {
        .log_line("    make_gen_* call: (not used; Normal with default variance path)")
        .log_line(sprintf("    gen_aruma_flex(n = %s, phi = %s, vara = %s, plot = FALSE)$y",
                          .fmt_log_value(as.integer(n_val)), .fmt_log_value(phi_val), .fmt_log_value(vara)))
      } else {
        .log_line(sprintf("    %s", .format_make_gen_call(innov_label, innov_params, use_fast = ctx$use_fast)))
        .log_line(sprintf("    gen_aruma_flex(n = %s, phi = %s, innov_gen = <make_gen_*>, plot = FALSE)$y",
                          .fmt_log_value(as.integer(n_val)), .fmt_log_value(phi_val)))
      }
      .log_line(sprintf("    wbg_boot_fast(y, nb = %s, maxp = %s, criterion = %s, bootadj = %s, min_p = %s)",
                        .fmt_log_value(as.integer(ctx$nb)), .fmt_log_value(as.integer(ctx$maxp)),
                        .fmt_log_value(ctx$criterion), .fmt_log_value(isTRUE(ctx$bootadj)),
                        .fmt_log_value(as.integer(ctx$min_p))))
    }

    .finalize_run <- function(cancelled = FALSE) {
      ctx <- run_state$ctx
      if (!is.null(ctx)) {
        # Publish accumulated results at run boundary (completion or cancel)
        if (!is.null(ctx$accum_results)) run_results_rv(ctx$accum_results)
        if (!is.null(ctx$accum_raw)) raw_results_rv(ctx$accum_raw)

        total_elapsed <- proc.time()[["elapsed"]] - ctx$t_total_start
        if (isTRUE(cancelled)) {
          .log_line(sprintf("Run cancelled at scenario boundary: %d/%d completed (%.1f s total)",
                            ctx$cells_ok, ctx$n_cells, total_elapsed))
          showNotification(
            sprintf("Run cancelled after %d/%d cells (%.0f s)",
                    ctx$cells_ok, ctx$n_cells, total_elapsed),
            type = "warning", duration = 8, session = session)
        } else {
          .log_line(sprintf("Grid complete: %d OK, %d errors, %.1f s total",
                            ctx$cells_ok, ctx$cells_err, total_elapsed))
          showNotification(
            sprintf("Completed %d/%d cells in %.0f s",
                    ctx$cells_ok, ctx$n_cells, total_elapsed),
            type = "message", duration = 8, session = session)
        }
      }

      if (!is.null(run_state$waiter)) {
        try(run_state$waiter$hide(), silent = TRUE)
      }
      if (!is.null(run_state$toast_id)) {
        try(removeNotification(run_state$toast_id, session = session), silent = TRUE)
      }
      run_state$waiter <- NULL
      run_state$toast_id <- NULL
      run_state$active <- FALSE
      run_state$cancel_requested <- FALSE
      run_state$ctx <- NULL
      .bump_run_status()
    }

    .run_next_cell <- function() {
      if (!isTRUE(run_state$active)) return(invisible(NULL))
      ctx <- run_state$ctx
      if (is.null(ctx)) return(invisible(NULL))

      # Boundary check: cancel before starting next scenario
      if (isTRUE(run_state$cancel_requested)) {
        .finalize_run(cancelled = TRUE)
        return(invisible(NULL))
      }

      if (ctx$cell_idx > ctx$n_cells) {
        .finalize_run(cancelled = FALSE)
        return(invisible(NULL))
      }

      row_i <- ctx$rows_to_run[ctx$cell_idx]
      phi_val        <- ctx$grid$phi[row_i]
      n_val          <- as.integer(ctx$grid$n[row_i])
      innov_label    <- ctx$grid$innov_label[row_i]
      innov_dist_str <- ctx$grid$innov_dist_str[row_i]
      innov_params   <- if (!is.null(ctx$grid$innov_params)) ctx$grid$innov_params[[row_i]] else list()

      .show_run_toast(detail = sprintf("phi=%.2f, n=%d, %s", phi_val, n_val, innov_dist_str))

      .log_line(sprintf("Cell %d/%d: phi=%.2f n=%d %s",
                        ctx$cell_idx, ctx$n_cells, phi_val, n_val, innov_dist_str))

      # Build innovation generator
      use_vara <- FALSE
      vara <- NULL
      innov_gen <- NULL

      if (innov_label == "Normal" && length(innov_params) == 0) {
        use_vara <- TRUE
        vara <- 1 - phi_val^2
      } else {
        innov_gen <- tryCatch(
          build_innov_gen(innov_label, innov_params, use_fast = ctx$use_fast),
          error = function(e) {
            .log_line(sprintf("  Generator error: %s", e$message))
            NULL
          }
        )
        if (is.null(innov_gen) && !use_vara) {
          .log_line("  SKIPPED (generator failed)")
          ctx$cells_err <- ctx$cells_err + 1L
          ctx$cell_idx <- ctx$cell_idx + 1L
          ctx <- .ctx_with_latest_logs(ctx)
          run_state$ctx <- ctx
          .show_run_toast()
          .bump_run_status()
          later::later(function() {
            shiny::withReactiveDomain(session, .run_next_cell())
          }, delay = 0)
          return(invisible(NULL))
        }
      }

      .log_cell_calls(ctx, phi_val, n_val, innov_label, innov_params, use_vara, vara)

      if (!is.null(ctx$base_seed) && !is.na(ctx$base_seed)) {
        set.seed(as.integer(ctx$base_seed) * 10000L + row_i)
      }

      results_list <- vector("list", ctx$nsims)
      failed_count <- 0L
      t_cell_start <- proc.time()[["elapsed"]]

      for (i in seq_len(ctx$nsims)) {
        y <- tryCatch({
          if (use_vara) {
            gen_aruma_flex(n_val, phi = phi_val, vara = vara, plot = FALSE)$y
          } else {
            gen_aruma_flex(n_val, phi = phi_val, innov_gen = innov_gen, plot = FALSE)$y
          }
        }, error = function(e) NULL)
        if (is.null(y)) { failed_count <- failed_count + 1L; next }

        results_list[[i]] <- tryCatch(
          wbg_boot_fast(y, nb = ctx$nb, maxp = ctx$maxp, criterion = ctx$criterion,
                        bootadj = ctx$bootadj, min_p = ctx$min_p),
          error = function(e) { failed_count <<- failed_count + 1L; NULL }
        )
      }

      results_list <- Filter(Negate(is.null), results_list)
      elapsed_cell <- proc.time()[["elapsed"]] - t_cell_start

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
        existing <- ctx$accum_results
        if (nrow(existing) > 0) {
          dup <- existing$n == n_val &
            abs(existing$phi - phi_val) < 1e-9 &
            existing$innov_dist == innov_dist_str
          existing <- existing[!dup, , drop = FALSE]
        }
        ctx$accum_results <- rbind(existing, cell_rates)

        cell_key <- paste(n_val, phi_val, innov_dist_str, sep = "|")
        raw <- ctx$accum_raw
        raw[[cell_key]] <- list(
          results      = results_list,
          n            = n_val,
          phi          = phi_val,
          innov_dist   = innov_dist_str,
          innov_label  = innov_label,
          innov_params = innov_params,
          use_fast     = ctx$use_fast,
          timestamp    = Sys.time()
        )
        ctx$accum_raw <- raw
      }

      .log_line(sprintf("  Completed %d/%d sims in %.1f s%s",
                        length(results_list), ctx$nsims, elapsed_cell,
                        if (failed_count > 0) sprintf(" (%d failed)", failed_count) else ""))

      ctx$cells_ok <- ctx$cells_ok + 1L
      ctx$cell_idx <- ctx$cell_idx + 1L
      ctx <- .ctx_with_latest_logs(ctx)
      run_state$ctx <- ctx
      .show_run_toast()
      .bump_run_status()

      later::later(function() {
        shiny::withReactiveDomain(session, .run_next_cell())
      }, delay = 0)
      invisible(NULL)
    }

    observeEvent(input$cancel_run_toast, {
      if (!isTRUE(run_state$active)) {
        return()
      }
      run_state$cancel_requested <- TRUE
      .show_run_toast()
      .bump_run_status()
      showNotification("Stop requested. Current scenario will finish first.",
                       type = "warning", duration = 4)
      .log_line("Stop requested: will cancel at next scenario boundary.")
    }, ignoreInit = TRUE)

    observeEvent(input$run_grid, {
      if (isTRUE(run_state$active)) {
        showNotification("A run is already in progress.", type = "warning", duration = 3)
        return()
      }

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
      use_fast  <- isTRUE(input$use_fast)
      base_seed <- input$seed

      run_errs <- .validate_run_params(nsims = nsims, nb = nb, maxp = maxp, seed = base_seed)
      if (length(run_errs) > 0) {
        showNotification(paste(run_errs, collapse = " "), type = "error", duration = 6)
        return()
      }

      run_state$waiter <- NULL

      run_log_rv("")
      run_state$ctx <- list(
        grid = grid,
        rows_to_run = rows_to_run,
        n_cells = length(rows_to_run),
        cell_idx = 1L,
        nsims = nsims,
        nb = nb,
        maxp = maxp,
        criterion = criterion,
        bootadj = bootadj,
        min_p = min_p,
        use_fast = use_fast,
        base_seed = base_seed,
        cells_ok = 0L,
        cells_err = 0L,
        t_total_start = proc.time()[["elapsed"]],
        log_lines = character(0),
        accum_results = run_results_rv(),
        accum_raw = raw_results_rv()
      )
      run_state$cancel_requested <- FALSE
      run_state$active <- TRUE
      .show_run_toast()
      .bump_run_status()

      .log_line(sprintf("Starting grid simulation: %d cells, %d sims each",
                        length(rows_to_run), nsims))
      later::later(function() {
        shiny::withReactiveDomain(session, .run_next_cell())
      }, delay = 0)
    })

    # =========================================================================
    # Run status & log
    # =========================================================================

    output$run_status_ui <- renderUI({
      run_status_tick()
      grid <- grid_rv()
      n_cells <- nrow(grid)
      if (n_cells == 0) {
        return(p(class = "text-body-secondary", "No cells in grid."))
      }

      if (isTRUE(run_state$active) && !is.null(run_state$ctx)) {
        ctx <- run_state$ctx
        current <- min(ctx$cell_idx, ctx$n_cells)
        status_msg <- sprintf("Running %d/%d scenarios (%s sims each)",
                              current, ctx$n_cells,
                              format(ctx$nsims, big.mark = ","))
        if (isTRUE(run_state$cancel_requested)) {
          status_msg <- paste0(status_msg, " - stop requested")
        }
        return(tagList(p(class = "text-warning", status_msg)))
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
    # Combined results: in-memory only
    # =========================================================================

    cap_combined_results <- reactive({
      run_results_rv()
    })

    # =========================================================================
    # Unified data accessors for detail sub-tabs
    # =========================================================================

    # Return list of wbg_boot_fast results for a cell_key "n|phi|innov_dist".
    # In-memory only.
    cap_sim_data <- function(cell_key) {
      if (is.null(cell_key) || !nzchar(cell_key)) return(NULL)
      raw <- raw_results_rv()
      if (cell_key %in% names(raw)) return(raw[[cell_key]]$results)
      NULL
    }

    # Available cell keys from in-memory results
    cap_cell_choices <- reactive({
      raw <- raw_results_rv()
      all_keys <- names(raw)
      if (length(all_keys) == 0) return(character(0))
      labels <- gsub("\\|", ", ", all_keys)
      setNames(all_keys, labels)
    })

    # =========================================================================
    # Sub-tab server fragments
    # =========================================================================

    # Tab 2: Coverage / Overview
    mod_capstone_overview_server(input, output, session,
      grid_rv = grid_rv, raw_results_rv = raw_results_rv)

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
        use_fast  = isTRUE(input$use_fast),
        base_seed = input$seed
      )))

    # Tab 6: P-value Diagnostics
    mod_capstone_pvalue_server(input, output, session,
      cap_sim_data = cap_sim_data, cap_cell_choices = cap_cell_choices)

    # Tab 7: Bootstrap Distribution
    mod_capstone_bootdist_server(input, output, session,
      cap_sim_data = cap_sim_data, cap_cell_choices = cap_cell_choices,
      selected_sim_from_individual = selected_sim_rv)

    # Tab 8: Null Model Diagnostics
    mod_capstone_null_diag_server(input, output, session,
      cap_sim_data = cap_sim_data, cap_cell_choices = cap_cell_choices)

  })
}
