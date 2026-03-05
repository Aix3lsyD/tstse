# =============================================================================
# Module: Ad-Hoc Simulation (Tab 10)
# =============================================================================

#' @description UI function for the Ad-Hoc Simulation module.
#' @param id Character string module namespace ID.
.adhoc_picker_input <- function(inputId, label, choices, selected = NULL, options = list()) {
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

.adhoc_toggle_input <- function(inputId, label, value = FALSE) {
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

.adhoc_action_button <- function(inputId, label, icon_name = NULL,
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

mod_adhoc_sim_ui <- function(id) {
  ns <- NS(id)

  tabPanel(
    "Ad-Hoc Simulation",
    br(),
    fluidRow(
      # --- Sidebar controls (col-3) ---
      column(3,
        wellPanel(
          h4("Run"),
          .adhoc_action_button(ns("sim_run"), "Run Simulation",
                               icon_name = "play", class = "btn-primary",
                               full_width = TRUE)
        ),

        wellPanel(
          h4("DGP Parameters"),
          numericInput(ns("sim_phi"), "AR(1) Phi", value = 0.95,
                       min = 0, max = 0.999, step = 0.01),
          numericInput(ns("sim_n"), "Sample Size (n)", value = 200,
                       min = 10, max = 2000, step = 10),
          .adhoc_picker_input(ns("sim_innov_dist"), "Innovation Distribution",
                              choices = c("Normal", "Student's t", "Skew-t", "GED",
                                          "Laplace", "Uniform", "Mixture Normal",
                                          "GARCH", "Heteroscedastic"),
                              selected = "Normal"),

          # --- Distribution-specific parameter panels ---
          conditionalPanel(
            condition = "input.sim_innov_dist == 'Normal'",
            ns = ns,
            numericInput(ns("sim_norm_sd"), "sd", value = 1, min = 0.001, step = 0.1)
          ),
          conditionalPanel(
            condition = "input.sim_innov_dist == \"Student's t\"",
            ns = ns,
            numericInput(ns("sim_t_df"), "df", value = 3, min = 1, step = 1),
            checkboxInput(ns("sim_t_scale"), "Scale to unit variance", value = FALSE)
          ),
          conditionalPanel(
            condition = "input.sim_innov_dist == 'Skew-t'",
            ns = ns,
            numericInput(ns("sim_skt_df"), "df", value = 5, min = 3, step = 1),
            numericInput(ns("sim_skt_alpha"), "Skewness (alpha)", value = 0, step = 0.1),
            checkboxInput(ns("sim_skt_scale"), "Scale to unit variance", value = FALSE)
          ),
          conditionalPanel(
            condition = "input.sim_innov_dist == 'GED'",
            ns = ns,
            numericInput(ns("sim_ged_nu"), "Shape (nu)", value = 2, min = 0.1, step = 0.1),
            numericInput(ns("sim_ged_sd"), "sd", value = 1, min = 0.001, step = 0.1)
          ),
          conditionalPanel(
            condition = "input.sim_innov_dist == 'Laplace'",
            ns = ns,
            numericInput(ns("sim_lap_scale"), "Scale", value = 0.707, min = 0.001, step = 0.01)
          ),
          conditionalPanel(
            condition = "input.sim_innov_dist == 'Uniform'",
            ns = ns,
            numericInput(ns("sim_unif_hw"), "Half-width", value = 1.732, min = 0.001, step = 0.1)
          ),
          conditionalPanel(
            condition = "input.sim_innov_dist == 'Mixture Normal'",
            ns = ns,
            numericInput(ns("sim_mix_sd1"), "sd1", value = 1, min = 0.001, step = 0.1),
            numericInput(ns("sim_mix_sd2"), "sd2", value = 3, min = 0.001, step = 0.1),
            numericInput(ns("sim_mix_prob1"), "prob1 (weight of component 1)",
                         value = 0.9, min = 0.01, max = 0.99, step = 0.05)
          ),
          conditionalPanel(
            condition = "input.sim_innov_dist == 'GARCH'",
            ns = ns,
            numericInput(ns("sim_garch_omega"), "omega", value = 0.1, min = 0.001, step = 0.01),
            textAreaInput(ns("sim_garch_alpha"), "alpha (comma-separated)",
                          value = "0.2, 0.175, 0.15, 0.125, 0.1, 0.075, 0.05, 0.025",
                          rows = 2),
            textInput(ns("sim_garch_beta"), "beta (comma-separated, optional)", value = "")
          ),
          conditionalPanel(
            condition = "input.sim_innov_dist == 'Heteroscedastic'",
            ns = ns,
            .adhoc_picker_input(ns("sim_hetero_shape"), "Weight Shape",
                                choices = c("linear", "sqrt", "log", "power", "exp",
                                            "step", "periodic"),
                                selected = "linear"),
            conditionalPanel(
              condition = paste0(
                "input.sim_hetero_shape == 'linear' || ",
                "input.sim_hetero_shape == 'sqrt' || ",
                "input.sim_hetero_shape == 'log' || ",
                "input.sim_hetero_shape == 'power' || ",
                "input.sim_hetero_shape == 'exp'"),
              ns = ns,
              numericInput(ns("sim_hetero_from"), "From (start weight)",
                           value = 1, min = 0.01, step = 0.1),
              numericInput(ns("sim_hetero_to"), "To (end weight)",
                           value = 10, min = 0.01, step = 0.1)
            ),
            conditionalPanel(
              condition = "input.sim_hetero_shape == 'power'",
              ns = ns,
              numericInput(ns("sim_hetero_power"), "Power", value = 2, min = 0.1, step = 0.1)
            ),
            conditionalPanel(
              condition = "input.sim_hetero_shape == 'step'",
              ns = ns,
              textInput(ns("sim_hetero_breaks"), "Breaks (comma-separated, 0-1)", value = "0.5"),
              textInput(ns("sim_hetero_levels"), "Levels (comma-separated, SD weights)", value = "1, 5")
            ),
            conditionalPanel(
              condition = "input.sim_hetero_shape == 'periodic'",
              ns = ns,
              numericInput(ns("sim_hetero_base_w"), "Base Weight",
                           value = 1, min = 0.01, step = 0.1),
              numericInput(ns("sim_hetero_amplitude"), "Amplitude", value = 0.5, step = 0.1),
              numericInput(ns("sim_hetero_period"), "Period (observations)",
                           value = 12, min = 1, step = 1)
            ),
            numericInput(ns("sim_hetero_sd"), "sd (base normal)", value = 1, min = 0.001, step = 0.1)
          )
        ),

        wellPanel(
          h4("Test Parameters"),
          numericInput(ns("sim_nsims"), "Number of Simulations", value = 1000,
                       min = 1, max = 10000, step = 10),
          numericInput(ns("sim_nb"), "Bootstrap Replicates (nb)", value = 399,
                       min = 1, max = 9999, step = 100),
          numericInput(ns("sim_maxp"), "Max AR Order (maxp)", value = 5,
                       min = 1, max = 20, step = 1),
          .adhoc_picker_input(ns("sim_criterion"), "Information Criterion",
                              choices = c("aic", "aicc", "bic"), selected = "aic"),
          .adhoc_toggle_input(ns("sim_bootadj"), "COBA Adjustment", value = FALSE),
          .adhoc_toggle_input(ns("sim_minp1"), "Minimum AR order = 1 (exclude AR(0))",
                              value = TRUE),
          .adhoc_toggle_input(ns("sim_use_fast"), "Use fast innovation generators",
                              value = FALSE),
          numericInput(ns("sim_seed"), "Random Seed (empty = random)",
                        value = NA, min = 1, step = 1)
        )
      ),

      # --- Results area (col-9) ---
      column(9,
        tabsetPanel(
          id = ns("sim_result_tabs"),
          tabPanel("Summary",
            br(),
            plotOutput(ns("sim_realization_plot"), height = "300px"),
            hr(),
            uiOutput(ns("sim_summary_ui"))
          ),
          tabPanel("Individual Results",
            br(),
            DT::dataTableOutput(ns("sim_results_table")),
            helpText("Click a row to view its bootstrap distribution in the Bootstrap Detail tab.")
          ),
          tabPanel("P-Value Distribution",
            br(),
            plotOutput(ns("sim_pval_hist"), height = "450px")
          ),
          tabPanel("Bootstrap Detail",
            br(),
            plotOutput(ns("sim_boot_detail"), height = "450px")
          ),
          tabPanel("Innovation Diagnostics",
            br(),
            plotOutput(ns("sim_innov_diag"), height = "900px")
          ),
          tabPanel("Null Model Diagnostics",
            br(),
            fluidRow(
              column(4,
                selectInput(ns("sim_null_reject_mode"), "Rejection rates shown:",
                            choices = c("All available" = "all",
                                        "COB" = "cob",
                                        "CO" = "co",
                                        "COBA" = "coba"),
                            selected = "all")
              )
            ),
            wellPanel(
              h5("Null Model Diagnostics (per-simulation)"),
              plotOutput(ns("sim_null_diag"), height = "900px")
            ),
            fluidRow(
              column(12,
                wellPanel(
                  h5("Estimated AR(1) Coefficient Distribution"),
                  plotOutput(ns("sim_null_phi_dist"), height = "420px")
                )
              )
            ),
            fluidRow(
              column(6,
                wellPanel(
                  h5("Test Statistic Distribution"),
                  plotOutput(ns("sim_null_tstat"), height = "420px")
                )
              ),
              column(6,
                wellPanel(
                  h5("Monte Carlo Convergence"),
                  selectInput(ns("sim_null_conv_method"), "Method:",
                              choices = c("Bootstrap (COB)" = "pvalue",
                                          "Asymptotic (CO)" = "pvalue_asymp",
                                          "COBA" = "pvalue_adj"),
                              selected = "pvalue"),
                  plotOutput(ns("sim_null_convergence"), height = "420px")
                )
              )
            )
          ),
          tabPanel("Run Log",
            br(),
            verbatimTextOutput(ns("sim_run_log"))
          ),
          mod_adhoc_profile_ui(ns)
        )
      )
    )
  )
}

#' @description Server function for the Ad-Hoc Simulation module.
#' @param id Character string module namespace ID.
mod_adhoc_sim_server <- function(id) {
  moduleServer(id, function(input, output, session) {

    # Storage for simulation results
    sim_results_rv <- reactiveVal(NULL)
    run_log_rv <- reactiveVal("")

    # Helper: build innov_dist string from selection and params
    # Delegates to shared helpers in utils.R via .adhoc_params_from_input()
    .build_innov_dist_str <- function(dist, input) {
      build_innov_dist_str(dist, .adhoc_params_from_input(input))
    }

    .build_innov_gen <- function(dist, input, use_fast = FALSE) {
      build_innov_gen(dist, .adhoc_params_from_input(input), use_fast = use_fast)
    }

    .log_line <- function(msg) {
      old <- run_log_rv()
      line <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), msg)
      run_log_rv(paste(c(old, line), collapse = if (nzchar(old)) "\n" else ""))
    }

    .fmt_log_value <- function(x) {
      if (is.function(x)) return("<function>")
      txt <- paste(deparse(x, width.cutoff = 500L), collapse = "")
      if (nchar(txt) > 240) txt <- paste0(substr(txt, 1, 237), "...")
      txt
    }

    .parse_num_list_or_null <- function(x) {
      if (is.null(x)) return(NULL)
      vals <- if (is.character(x)) {
        suppressWarnings(as.numeric(trimws(strsplit(x, ",")[[1]])))
      } else {
        suppressWarnings(as.numeric(x))
      }
      if (length(vals) == 0 || any(is.na(vals))) return(NULL)
      vals
    }

    .parse_num_list <- function(x) {
      if (is.null(x)) return(numeric(0))
      txt <- trimws(as.character(x))
      if (!nzchar(txt)) return(numeric(0))
      vals <- suppressWarnings(as.numeric(trimws(strsplit(txt, ",")[[1]])))
      if (length(vals) == 0 || any(is.na(vals))) return(NULL)
      vals
    }

    has_shinyvalidate <- requireNamespace("shinyvalidate", quietly = TRUE)
    adhoc_validator <- NULL

    if (has_shinyvalidate) {
      v <- shinyvalidate::InputValidator$new()

      .rule_number <- function(expr, msg) {
        force(expr)
        function(value) {
          if (is.null(value) || is.na(value) || !isTRUE(expr(value))) msg else NULL
        }
      }

      .rule_optional_integer <- function(min_val = 1L) {
        force(min_val)
        function(value) {
          if (is.null(value) || is.na(value) || !nzchar(trimws(as.character(value)))) return(NULL)
          int_val <- suppressWarnings(as.integer(value))
          if (is.na(int_val) || int_val < min_val) {
            sprintf("Must be an integer >= %d", as.integer(min_val))
          } else {
            NULL
          }
        }
      }

      .rule_when <- function(cond_fn, rule_fn) {
        force(cond_fn)
        force(rule_fn)
        function(value) {
          if (!isTRUE(cond_fn())) return(NULL)
          rule_fn(value)
        }
      }

      .rule_csv_nonneg <- function(required = TRUE) {
        force(required)
        function(value) {
          vals <- .parse_num_list(value)
          if (is.null(vals)) return("Must be a comma-separated numeric list")
          if (length(vals) == 0) {
            if (isTRUE(required)) return("Provide at least one value")
            return(NULL)
          }
          if (any(vals < 0)) return("Values must be >= 0")
          NULL
        }
      }

      .rule_hetero_breaks <- function(value) {
        brk <- .parse_num_list(value)
        if (is.null(brk) || length(brk) == 0) return("Provide comma-separated numeric breaks")
        if (any(brk <= 0 | brk >= 1)) return("Breaks must be between 0 and 1")
        if (is.unsorted(brk, strictly = TRUE)) return("Breaks must be strictly increasing")
        NULL
      }

      .rule_hetero_levels <- function(value) {
        lvl <- .parse_num_list(value)
        if (is.null(lvl) || length(lvl) == 0) return("Provide comma-separated numeric levels")
        if (any(lvl <= 0)) return("Levels must be > 0")
        brk <- .parse_num_list(input$sim_hetero_breaks)
        if (!is.null(brk) && length(brk) > 0 && length(lvl) != length(brk) + 1) {
          return("Levels count must equal breaks count + 1")
        }
        NULL
      }

      v$add_rule("sim_phi", .rule_number(function(x) is.numeric(x) && x >= 0 && x < 1,
                                           "Must be in [0, 1)"))
      v$add_rule("sim_n", .rule_number(function(x) as.integer(x) >= 10 && as.integer(x) <= 2000,
                                         "Must be an integer from 10 to 2000"))
      v$add_rule("sim_nsims", .rule_number(function(x) as.integer(x) >= 1,
                                             "Must be an integer >= 1"))
      v$add_rule("sim_nb", .rule_number(function(x) as.integer(x) >= 1,
                                          "Must be an integer >= 1"))
      v$add_rule("sim_maxp", .rule_number(function(x) as.integer(x) >= 1 && as.integer(x) <= 20,
                                            "Must be an integer from 1 to 20"))
      v$add_rule("sim_seed", .rule_optional_integer(min_val = 1L))

      v$add_rule("sim_norm_sd", .rule_when(
        function() identical(input$sim_innov_dist, "Normal"),
        .rule_number(function(x) x > 0, "Must be > 0")
      ))
      v$add_rule("sim_t_df", .rule_when(
        function() identical(input$sim_innov_dist, "Student's t"),
        .rule_number(function(x) x >= 1, "Must be >= 1")
      ))
      v$add_rule("sim_skt_df", .rule_when(
        function() identical(input$sim_innov_dist, "Skew-t"),
        .rule_number(function(x) x >= 3, "Must be >= 3")
      ))
      v$add_rule("sim_skt_alpha", .rule_when(
        function() identical(input$sim_innov_dist, "Skew-t"),
        .rule_number(function(x) is.finite(as.numeric(x)), "Must be numeric")
      ))
      v$add_rule("sim_ged_nu", .rule_when(
        function() identical(input$sim_innov_dist, "GED"),
        .rule_number(function(x) x > 0, "Must be > 0")
      ))
      v$add_rule("sim_ged_sd", .rule_when(
        function() identical(input$sim_innov_dist, "GED"),
        .rule_number(function(x) x > 0, "Must be > 0")
      ))
      v$add_rule("sim_lap_scale", .rule_when(
        function() identical(input$sim_innov_dist, "Laplace"),
        .rule_number(function(x) x > 0, "Must be > 0")
      ))
      v$add_rule("sim_unif_hw", .rule_when(
        function() identical(input$sim_innov_dist, "Uniform"),
        .rule_number(function(x) x > 0, "Must be > 0")
      ))
      v$add_rule("sim_mix_sd1", .rule_when(
        function() identical(input$sim_innov_dist, "Mixture Normal"),
        .rule_number(function(x) x > 0, "Must be > 0")
      ))
      v$add_rule("sim_mix_sd2", .rule_when(
        function() identical(input$sim_innov_dist, "Mixture Normal"),
        .rule_number(function(x) x > 0, "Must be > 0")
      ))
      v$add_rule("sim_mix_prob1", .rule_when(
        function() identical(input$sim_innov_dist, "Mixture Normal"),
        .rule_number(function(x) x > 0 && x < 1, "Must be between 0 and 1")
      ))

      v$add_rule("sim_garch_omega", .rule_when(
        function() identical(input$sim_innov_dist, "GARCH"),
        .rule_number(function(x) x > 0, "Must be > 0")
      ))
      v$add_rule("sim_garch_alpha", .rule_when(
        function() identical(input$sim_innov_dist, "GARCH"),
        .rule_csv_nonneg(required = TRUE)
      ))
      v$add_rule("sim_garch_beta", .rule_when(
        function() identical(input$sim_innov_dist, "GARCH"),
        .rule_csv_nonneg(required = FALSE)
      ))

      v$add_rule("sim_hetero_sd", .rule_when(
        function() identical(input$sim_innov_dist, "Heteroscedastic"),
        .rule_number(function(x) x > 0, "Must be > 0")
      ))
      v$add_rule("sim_hetero_from", .rule_when(
        function() {
          identical(input$sim_innov_dist, "Heteroscedastic") &&
            input$sim_hetero_shape %in% c("linear", "sqrt", "log", "power", "exp")
        },
        .rule_number(function(x) x > 0, "Must be > 0")
      ))
      v$add_rule("sim_hetero_to", .rule_when(
        function() {
          identical(input$sim_innov_dist, "Heteroscedastic") &&
            input$sim_hetero_shape %in% c("linear", "sqrt", "log", "power", "exp")
        },
        function(value) {
          to_val <- suppressWarnings(as.numeric(value))
          from_val <- suppressWarnings(as.numeric(input$sim_hetero_from))
          if (is.na(to_val) || to_val <= 0) return("Must be > 0")
          if (!is.na(from_val) && to_val < from_val) return("Must be >= From")
          NULL
        }
      ))
      v$add_rule("sim_hetero_power", .rule_when(
        function() identical(input$sim_innov_dist, "Heteroscedastic") &&
          identical(input$sim_hetero_shape, "power"),
        .rule_number(function(x) x > 0, "Must be > 0")
      ))
      v$add_rule("sim_hetero_breaks", .rule_when(
        function() identical(input$sim_innov_dist, "Heteroscedastic") &&
          identical(input$sim_hetero_shape, "step"),
        .rule_hetero_breaks
      ))
      v$add_rule("sim_hetero_levels", .rule_when(
        function() identical(input$sim_innov_dist, "Heteroscedastic") &&
          identical(input$sim_hetero_shape, "step"),
        .rule_hetero_levels
      ))
      v$add_rule("sim_hetero_base_w", .rule_when(
        function() identical(input$sim_innov_dist, "Heteroscedastic") &&
          identical(input$sim_hetero_shape, "periodic"),
        .rule_number(function(x) x > 0, "Must be > 0")
      ))
      v$add_rule("sim_hetero_period", .rule_when(
        function() identical(input$sim_innov_dist, "Heteroscedastic") &&
          identical(input$sim_hetero_shape, "periodic"),
        .rule_number(function(x) as.integer(x) >= 1, "Must be an integer >= 1")
      ))

      v$enable()
      adhoc_validator <- v
    }

    .format_make_gen_call <- function(dist, params, use_fast) {
      suffix <- if (isTRUE(use_fast)) "_fast" else ""
      switch(dist,
        "Normal" = {
          sd_val <- params$norm_sd
          if (is.null(sd_val) || is.na(sd_val) || sd_val <= 0) sd_val <- 1
          sprintf("make_gen_norm%s(sd = %s)", suffix, .fmt_log_value(sd_val))
        },
        "Student's t" = {
          df_val <- params$t_df
          if (is.null(df_val) || is.na(df_val) || df_val < 1) df_val <- 3
          scale_val <- isTRUE(params$t_scale)
          sprintf("make_gen_t%s(df = %s, scale = %s)",
                  suffix, .fmt_log_value(df_val), .fmt_log_value(scale_val))
        },
        "Skew-t" = {
          df_val <- params$skt_df
          if (is.null(df_val) || is.na(df_val) || df_val < 3) df_val <- 5
          alpha_val <- params$skt_alpha
          if (is.null(alpha_val) || is.na(alpha_val)) alpha_val <- 0
          scale_val <- isTRUE(params$skt_scale)
          sprintf("make_gen_skt%s(df = %s, alpha = %s, scale = %s)",
                  suffix, .fmt_log_value(df_val), .fmt_log_value(alpha_val), .fmt_log_value(scale_val))
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
        },
        sprintf("make_gen_<unknown>(label = %s)", .fmt_log_value(dist))
      )
    }

    .log_call_block <- function(dist, params, use_fast, n, phi, nb, maxp, criterion, bootadj, min_p) {
      .log_line("Calls for this run:")
      .log_line(sprintf("  %s", .format_make_gen_call(dist, params, use_fast = use_fast)))
      .log_line(sprintf(
        "  gen_aruma_flex(n = %s, phi = %s, innov_gen = <make_gen_*>, seed = NULL, plot = FALSE)$y",
        .fmt_log_value(as.integer(n)), .fmt_log_value(phi)
      ))
      .log_line(sprintf(
        "  wbg_boot_fast(x, nb = %s, maxp = %s, criterion = %s, bootadj = %s, min_p = %s)",
        .fmt_log_value(as.integer(nb)), .fmt_log_value(as.integer(maxp)),
        .fmt_log_value(criterion), .fmt_log_value(isTRUE(bootadj)), .fmt_log_value(as.integer(min_p))
      ))
    }

    # Helper: compute theoretical density for a distribution
    # Returns list(dfun, label, col) or NULL if not applicable
    .innov_theoretical_density <- function(dist, params) {
      switch(dist,
        "Normal" = {
          sd_val <- params$sd %||% 1
          list(dfun = function(x) dnorm(x, 0, sd_val),
               label = sprintf("N(0, %.2f)", sd_val), col = "red")
        },
        "Student's t" = {
          df_val <- params$df %||% 3
          scaled <- isTRUE(params$scale)
          if (scaled && df_val > 2) {
            c_sc <- sqrt(df_val / (df_val - 2))
            list(dfun = function(x) dt(x * c_sc, df = df_val) * c_sc,
                 label = sprintf("t(%d) scaled", df_val), col = "red")
          } else {
            list(dfun = function(x) dt(x, df = df_val),
                 label = sprintf("t(%d)", df_val), col = "red")
          }
        },
        "Skew-t" = {
          if (!requireNamespace("sn", quietly = TRUE)) return(NULL)
          df_val <- params$df %||% 5
          alpha_val <- params$alpha %||% 0
          scaled <- isTRUE(params$scale)
          theo_mean <- 0; theo_sd <- 1
          if (df_val > 1 && alpha_val != 0) {
            delta <- alpha_val / sqrt(1 + alpha_val^2)
            b_nu <- sqrt(df_val / pi) * exp(lgamma((df_val - 1) / 2) - lgamma(df_val / 2))
            theo_mean <- b_nu * delta
          }
          if (scaled && df_val > 2) {
            cumulants <- tryCatch(
              sn::st.cumulants(xi = 0, omega = 1, alpha = alpha_val, nu = df_val, n = 2),
              error = function(e) c(NA_real_, NA_real_))
            theo_var <- as.numeric(cumulants[2])
            if (is.finite(theo_var) && theo_var > 0) {
              theo_sd <- sqrt(theo_var)
              theo_mean <- as.numeric(cumulants[1])
            }
          }
          tm <- theo_mean; ts <- theo_sd; a <- alpha_val; d <- df_val
          list(dfun = function(x) sn::dst(x * ts + tm, xi = 0, omega = 1, alpha = a, nu = d) * ts,
               label = sprintf("Skew-t(df=%d, alpha=%.1f)", df_val, alpha_val), col = "red")
        },
        "GED" = {
          if (!requireNamespace("fGarch", quietly = TRUE)) return(NULL)
          nu_val <- params$nu %||% 2
          sd_val <- params$sd %||% 1
          list(dfun = function(x) fGarch::dged(x, mean = 0, sd = sd_val, nu = nu_val),
               label = sprintf("GED(nu=%.1f, sd=%.2f)", nu_val, sd_val), col = "red")
        },
        "Laplace" = {
          sc <- params$scale %||% (1 / sqrt(2))
          list(dfun = function(x) (1 / (2 * sc)) * exp(-abs(x) / sc),
               label = sprintf("Laplace(scale=%.3f)", sc), col = "red")
        },
        "Uniform" = {
          hw <- params$half_width %||% sqrt(3)
          list(dfun = function(x) dunif(x, min = -hw, max = hw),
               label = sprintf("Unif(%.2f, %.2f)", -hw, hw), col = "red")
        },
        "Mixture Normal" = {
          s1 <- params$sd1 %||% 1; s2 <- params$sd2 %||% 3; p1 <- params$prob1 %||% 0.9
          list(dfun = function(x) p1 * dnorm(x, 0, s1) + (1 - p1) * dnorm(x, 0, s2),
               label = sprintf("MixN(%.0f%%*N(0,%.1f) + %.0f%%*N(0,%.1f))",
                                p1 * 100, s1, (1 - p1) * 100, s2), col = "red")
        },
        NULL
      )
    }

    # Run simulation
    observeEvent(input$sim_run, {
      run_log_rv("")
      if (!is.null(adhoc_validator) && !isTRUE(adhoc_validator$is_valid())) {
        .log_line("Aborted: invalid input configuration.")
        showNotification("Please fix highlighted input errors before running.",
                         type = "warning", duration = 5)
        return()
      }

      dist <- input$sim_innov_dist
      use_fast <- isTRUE(input$sim_use_fast)

      # Check optional packages
      if (!use_fast && dist == "GARCH" && !requireNamespace("rugarch", quietly = TRUE)) {
        .log_line("Aborted: missing optional package 'rugarch' for GARCH innovations.")
        showNotification("Install the 'rugarch' package to use GARCH innovations.",
                         type = "error", duration = 8)
        return()
      }
      if (!use_fast && dist == "Skew-t" && !requireNamespace("sn", quietly = TRUE)) {
        .log_line("Aborted: missing optional package 'sn' for Skew-t innovations.")
        showNotification("Install the 'sn' package to use Skew-t innovations.",
                         type = "error", duration = 8)
        return()
      }
      if (!use_fast && dist == "GED" && !requireNamespace("fGarch", quietly = TRUE)) {
        .log_line("Aborted: missing optional package 'fGarch' for GED innovations.")
        showNotification("Install the 'fGarch' package to use GED innovations.",
                         type = "error", duration = 8)
        return()
      }

      # Read and clamp inputs
      phi <- input$sim_phi
      if (is.null(phi) || is.na(phi)) phi <- 0.95
      phi <- max(0, min(0.999, phi))

      n <- as.integer(input$sim_n)
      if (is.na(n) || n < 10) n <- 10L
      if (n > 2000) n <- 2000L

      nsims <- as.integer(input$sim_nsims)
      if (is.na(nsims) || nsims < 1) nsims <- 1L

      nb <- as.integer(input$sim_nb)
      if (is.na(nb) || nb < 1) nb <- 399L

      maxp <- as.integer(input$sim_maxp)
      if (is.na(maxp) || maxp < 1) maxp <- 5L
      maxp <- min(maxp, 20L)

      criterion <- match.arg(input$sim_criterion, c("aic", "aicc", "bic"))
      bootadj <- isTRUE(input$sim_bootadj)
      min_p <- if (isTRUE(input$sim_minp1)) 1L else 0L

      seed_val <- input$sim_seed
      has_seed <- !is.null(seed_val) && !is.na(seed_val)

      params_generic <- .adhoc_params_from_input(input)
      .log_line(sprintf(
        "Starting ad-hoc simulation: dist=%s, phi=%s, n=%s, nsims=%s, seed=%s",
        dist,
        .fmt_log_value(phi),
        .fmt_log_value(n),
        .fmt_log_value(nsims),
        if (has_seed) .fmt_log_value(as.integer(seed_val)) else "<random>"
      ))
      .log_call_block(dist, params_generic, use_fast, n, phi, nb, maxp, criterion, bootadj, min_p)

      # Build generator
      innov_gen <- tryCatch(.build_innov_gen(dist, input, use_fast = use_fast), error = function(e) {
        .log_line(sprintf("Generator error: %s", e$message))
        showNotification(paste("Error building generator:", e$message),
                         type = "error", duration = 8)
        return(NULL)
      })
      if (is.null(innov_gen)) return()

      innov_dist_str <- .build_innov_dist_str(dist, input)

      # Snapshot distribution parameters for theoretical density overlay
      innov_params <- switch(dist,
        "Normal" = {
          sd_v <- input$sim_norm_sd
          if (is.null(sd_v) || is.na(sd_v) || sd_v <= 0) sd_v <- 1
          list(sd = sd_v)
        },
        "Student's t" = {
          df_v <- input$sim_t_df
          if (is.null(df_v) || is.na(df_v) || df_v < 1) df_v <- 3
          sc_v <- isTRUE(input$sim_t_scale)
          if (df_v <= 2) sc_v <- FALSE
          list(df = df_v, scale = sc_v)
        },
        "Skew-t" = {
          df_v <- input$sim_skt_df
          if (is.null(df_v) || is.na(df_v) || df_v < 3) df_v <- 5
          al_v <- input$sim_skt_alpha
          if (is.null(al_v) || is.na(al_v)) al_v <- 0
          sc_v <- isTRUE(input$sim_skt_scale)
          if (df_v <= 2) sc_v <- FALSE
          list(df = df_v, alpha = al_v, scale = sc_v)
        },
        "GED" = {
          nu_v <- input$sim_ged_nu
          if (is.null(nu_v) || is.na(nu_v) || nu_v <= 0) nu_v <- 2
          sd_v <- input$sim_ged_sd
          if (is.null(sd_v) || is.na(sd_v) || sd_v <= 0) sd_v <- 1
          list(nu = nu_v, sd = sd_v)
        },
        "Laplace" = {
          sc_v <- input$sim_lap_scale
          if (is.null(sc_v) || is.na(sc_v) || sc_v <= 0) sc_v <- 1 / sqrt(2)
          list(scale = sc_v)
        },
        "Uniform" = {
          hw_v <- input$sim_unif_hw
          if (is.null(hw_v) || is.na(hw_v) || hw_v <= 0) hw_v <- sqrt(3)
          list(half_width = hw_v)
        },
        "Mixture Normal" = {
          s1 <- input$sim_mix_sd1; s2 <- input$sim_mix_sd2; p1 <- input$sim_mix_prob1
          if (is.null(s1) || is.na(s1) || s1 <= 0) s1 <- 1
          if (is.null(s2) || is.na(s2) || s2 <= 0) s2 <- 3
          if (is.null(p1) || is.na(p1) || p1 <= 0 || p1 >= 1) p1 <- 0.9
          list(sd1 = s1, sd2 = s2, prob1 = p1)
        },
        list()
      )

      # Run simulation loop
      results <- vector("list", nsims)
      failed_count <- 0L
      update_every <- max(1L, nsims %/% 20L)
      prog_id <- showNotification(
        sprintf("Running simulations... 0 / %d", nsims),
        duration = NULL, closeButton = FALSE, type = "message")
      on.exit(removeNotification(prog_id), add = TRUE)
      if (has_seed) set.seed(seed_val)
      t_start <- proc.time()[["elapsed"]]
      for (i in seq_len(nsims)) {
        results[[i]] <- tryCatch({
          x <- gen_aruma_flex(n, phi = phi, innov_gen = innov_gen,
                              seed = NULL, plot = FALSE)$y
          wbg_boot_fast(x, nb = nb, maxp = maxp,
                        criterion = criterion,
                        bootadj = bootadj, min_p = min_p)
        }, error = function(e) {
          failed_count <<- failed_count + 1L
          NULL
        })
        if (i %% update_every == 0L || i == nsims) {
          showNotification(
            sprintf("Running simulations... %d / %d", i, nsims),
            id = prog_id, duration = NULL, closeButton = FALSE, type = "message")
        }
      }
      elapsed_secs <- proc.time()[["elapsed"]] - t_start
      results <- Filter(Negate(is.null), results)
      if (length(results) == 0) {
        .log_line("All iterations failed. No results to display.")
        showNotification("All iterations failed. Try adjusting parameters.",
                         type = "error", duration = 8)
        return()
      }

      # Generate diagnostic samples (after loop so main results are reproducible)
      innov_sample <- innov_gen(n)
      sample_realization <- gen_aruma_flex(n, phi = phi, innov_gen = innov_gen,
                                           seed = NULL, plot = FALSE)$y
      is_time_dependent <- dist %in% c("GARCH", "Heteroscedastic")

      sim_results_rv(list(
        results = results,
        phi = phi,
        n = n,
        nsims = length(results),
        nb = nb,
        maxp = maxp,
        criterion = criterion,
        bootadj = bootadj,
        min_p = min_p,
        innov_dist = dist,
        innov_dist_str = innov_dist_str,
        seed = if (has_seed) seed_val else NULL,
        innov_sample = innov_sample,
        sample_realization = sample_realization,
        is_time_dependent = is_time_dependent,
        innov_params = innov_params,
        elapsed_secs = elapsed_secs
      ))

      showNotification(sprintf("Completed %d/%d simulations in %.1f s.",
                               length(results), nsims, elapsed_secs),
                       type = if (failed_count > 0) "warning" else "message",
                       duration = 4)
      .log_line(sprintf(
        "Completed ad-hoc simulation: %d/%d succeeded in %.1f s%s",
        length(results), nsims, elapsed_secs,
        if (failed_count > 0) sprintf(" (%d failed)", failed_count) else ""
      ))
    })

    # Summary tab
    output$sim_summary_ui <- renderUI({
      res <- sim_results_rv()
      if (is.null(res)) {
        return(div(class = "text-body-secondary", style = "padding: 20px;",
                   "Run a simulation to see results here."))
      }

      results <- res$results
      nsims <- res$nsims
      pvals <- vapply(results, function(r) r$pvalue, numeric(1))
      pvals_asymp <- vapply(results, function(r) r$pvalue_asymp, numeric(1))
      reject_boot <- mean(pvals < 0.05, na.rm = TRUE)
      reject_asymp <- mean(pvals_asymp < 0.05, na.rm = TRUE)

      # Build rejection rate rows with MC standard errors
      mc_se <- function(p, n) sqrt(p * (1 - p) / n)
      se_boot <- mc_se(reject_boot, nsims)
      se_asymp <- mc_se(reject_asymp, nsims)

      reject_rows <- tagList(
        tags$tr(tags$td(tags$strong("Bootstrap")),
                tags$td(sprintf("%.4f (SE: %.4f)", reject_boot, se_boot))),
        tags$tr(tags$td(tags$strong("Asymptotic")),
                tags$td(sprintf("%.4f (SE: %.4f)", reject_asymp, se_asymp)))
      )
      if (res$bootadj) {
        pvals_adj <- vapply(results, function(r) {
          if (is.null(r$pvalue_adj)) NA_real_ else r$pvalue_adj
        }, numeric(1))
        reject_adj <- mean(pvals_adj < 0.05, na.rm = TRUE)
        se_adj <- mc_se(reject_adj, nsims)
        reject_rows <- tagList(reject_rows,
          tags$tr(tags$td(tags$strong("COBA-Adjusted")),
                  tags$td(sprintf("%.4f (SE: %.4f)", reject_adj, se_adj)))
        )
      }

      innov_label <- sprintf("%s (%s)", res$innov_dist, res$innov_dist_str)

      tagList(
        h4("Simulation Summary"),
        hr(),
        fluidRow(
          column(6,
            wellPanel(
              h5("Configuration"),
              tags$table(class = "table table-sm",
                tags$tr(tags$td(tags$strong("phi")), tags$td(res$phi)),
                tags$tr(tags$td(tags$strong("n")), tags$td(res$n)),
                tags$tr(tags$td(tags$strong("Innovation")), tags$td(innov_label)),
                tags$tr(tags$td(tags$strong("n_sims")), tags$td(res$nsims)),
                tags$tr(tags$td(tags$strong("nb")), tags$td(res$nb)),
                tags$tr(tags$td(tags$strong("maxp")), tags$td(res$maxp)),
                tags$tr(tags$td(tags$strong("Criterion")), tags$td(res$criterion)),
                tags$tr(tags$td(tags$strong("COBA")),
                        tags$td(if (res$bootadj) "Yes" else "No")),
                tags$tr(tags$td(tags$strong("Min AR Order")), tags$td(res$min_p)),
                tags$tr(tags$td(tags$strong("Seed")),
                        tags$td(if (is.null(res$seed)) "random" else res$seed)),
                tags$tr(tags$td(tags$strong("Elapsed")),
                        tags$td(sprintf("%.1f s", res$elapsed_secs)))
              )
            )
          ),
          column(6,
            wellPanel(
              h5("Rejection Rates (alpha = 0.05)"),
              tags$table(class = "table table-sm",
                tags$thead(tags$tr(tags$th("Method"), tags$th("Rate (SE)"))),
                reject_rows
              )
            ),
            wellPanel(
              h5("P-Value Summary (Bootstrap)"),
              tags$table(class = "table table-sm",
                tags$tr(tags$td(tags$strong("Mean")),
                        tags$td(sprintf("%.4f", mean(pvals, na.rm = TRUE)))),
                tags$tr(tags$td(tags$strong("Median")),
                        tags$td(sprintf("%.4f", median(pvals, na.rm = TRUE))))
              )
            )
          )
        )
      )
    })

    # Individual Results table
    output$sim_results_table <- DT::renderDataTable({
      res <- sim_results_rv()
      if (is.null(res)) {
        return(DT::datatable(
          data.frame(Message = "Run a simulation to see individual results."),
          rownames = FALSE, selection = "none",
          options = list(dom = "t", ordering = FALSE)
        ))
      }

      results <- res$results
      df <- data.frame(
        Sim = seq_along(results),
        obs_stat = vapply(results, function(r) r$tco_obs, numeric(1)),
        pvalue = vapply(results, function(r) r$pvalue, numeric(1)),
        pvalue_asymp = vapply(results, function(r) r$pvalue_asymp, numeric(1)),
        AR_order = vapply(results, function(r) r$p, integer(1)),
        vara = vapply(results, function(r) r$vara, numeric(1))
      )

      round_cols <- c("obs_stat", "pvalue", "pvalue_asymp", "vara")
      display_names <- c("Sim", "Obs. Stat", "P-Value (Boot)",
                         "P-Value (Asymp)", "AR Order", "Innov. Var.")

      if (res$bootadj) {
        df$pvalue_adj <- vapply(results, function(r) {
          if (is.null(r$pvalue_adj)) NA_real_ else r$pvalue_adj
        }, numeric(1))
        round_cols <- c(round_cols, "pvalue_adj")
        display_names <- c(display_names, "P-Value (COBA)")
      }

      # formatRound uses original df column names, so pass colnames
      # as a simple character vector (positional) to avoid remapping
      DT::datatable(df, selection = "single", rownames = FALSE,
                    colnames = display_names,
                    options = list(pageLength = 50,
                                   lengthMenu = c(20, 50, 100, 200),
                                   scrollX = TRUE)) |>
        DT::formatRound(columns = round_cols, digits = 4)
    })

    # P-Value Distribution histogram
    output$sim_pval_hist <- renderPlot(bg = "transparent", {
      op <- viewer_base_par()
      on.exit(par(op), add = TRUE)
      res <- sim_results_rv()
      if (is.null(res)) {
        plot.new()
        text(0.5, 0.5, "Run a simulation to see the p-value distribution.",
             cex = 1.2, col = "grey50")
        return()
      }

      pvals <- vapply(res$results, function(r) r$pvalue, numeric(1))
      hist(pvals, breaks = seq(0, 1, by = 0.05), col = "steelblue",
           border = "white", main = "Bootstrap P-Value Distribution",
           xlab = "P-Value", ylab = "Frequency", xlim = c(0, 1))
      abline(v = 0.05, lty = 2, col = "red", lwd = 2)
      legend("topright", legend = sprintf("Reject at 0.05: %.1f%%",
             mean(pvals < 0.05, na.rm = TRUE) * 100),
             bty = "n", text.col = "red")
    })

    # Bootstrap Detail plot
    output$sim_boot_detail <- renderPlot(bg = "transparent", {
      op <- viewer_base_par()
      on.exit(par(op), add = TRUE)
      res <- sim_results_rv()
      sel <- input$sim_results_table_rows_selected

      if (is.null(res) || is.null(sel) || length(sel) == 0) {
        plot.new()
        text(0.5, 0.5, "Select a simulation from the Individual Results table.",
             cex = 1.2, col = "grey50")
        return()
      }

      r <- res$results[[sel]]
      boot_t <- r$boot_tstats
      obs_t <- r$tco_obs
      pval <- r$pvalue

      # Pre-compute ylim to fit both histogram and theoretical curve
      p_order <- r$p
      df_t <- res$n - 2 - p_order
      h <- hist(boot_t, breaks = 30, plot = FALSE)
      x_seq <- seq(min(boot_t), max(boot_t), length.out = 300)
      curve_max <- if (df_t > 0) max(dt(x_seq, df = df_t)) else 0
      y_max <- max(h$density, curve_max) * 1.05

      hist(boot_t, breaks = h$breaks, col = "lightblue", border = "white",
           main = sprintf("Bootstrap Distribution (Sim #%d)", sel),
           xlab = "t-statistic", ylab = "Density", probability = TRUE,
           ylim = c(0, y_max))

      # Theoretical t(n - 2 - p) overlay
      if (df_t > 0) {
        lines(x_seq, dt(x_seq, df = df_t), col = "darkgreen", lwd = 2, lty = 2)
      }

      abline(v = obs_t, col = "red", lwd = 2)

      leg_labels <- c(sprintf("Observed t = %.3f", obs_t),
                      sprintf("p-value = %.4f", pval))
      leg_cols <- c("red", NA)
      leg_lwd <- c(2, NA)
      leg_lty <- c(1, NA)
      if (df_t > 0) {
        leg_labels <- c(leg_labels, sprintf("t(%d) asymptotic", df_t))
        leg_cols <- c(leg_cols, "darkgreen")
        leg_lwd <- c(leg_lwd, 2)
        leg_lty <- c(leg_lty, 2)
      }
      legend("topright", legend = leg_labels,
             col = leg_cols, lwd = leg_lwd, lty = leg_lty, bty = "n")
    })

    # Sample realization plot on Summary tab
    output$sim_realization_plot <- renderPlot(bg = "transparent", {
      res <- sim_results_rv()
      if (is.null(res)) return(NULL)
      op <- viewer_base_par()
      on.exit(par(op), add = TRUE)
      y <- res$sample_realization
      plot(y, type = "l", col = "steelblue", lwd = 1.2,
           xlab = "Time", ylab = "Value",
           main = sprintf("Example AR(1) Realization (phi=%.3f, n=%d, %s)",
                          res$phi, res$n, res$innov_dist_str))
      abline(h = 0, lty = 3, col = "grey60")
    })

    # Innovation Diagnostics tab
    output$sim_innov_diag <- renderPlot(bg = "transparent", {
      res <- sim_results_rv()
      if (is.null(res)) {
        plot.new()
        text(0.5, 0.5, "Run a simulation to see innovation diagnostics.",
             cex = 1.2, col = "grey50")
        return()
      }

      innov <- res$innov_sample
      is_td <- res$is_time_dependent
      dist_label <- res$innov_dist_str

      if (is_td) {
        # Time-dependent: 3-panel layout (time series, histogram, squared innovations)
        op <- par(mfrow = c(3, 1), mar = c(5, 5.5, 3, 1), cex.axis = 1.3, cex.lab = 1.4)
        on.exit(par(op))

        # 1. Innovation time series
        plot(innov, type = "l", col = "darkorange", lwd = 1,
             xlab = "Time", ylab = "Innovation",
             main = sprintf("Innovation Time Series (%s)", dist_label))
        abline(h = 0, lty = 3, col = "grey60")

        # 2. Histogram with normal reference (GARCH/Hetero have no single marginal density)
        h <- hist(innov, breaks = 40, plot = FALSE)
        x_seq <- seq(min(innov), max(innov), length.out = 300)
        curve_vals <- dnorm(x_seq, mean(innov), sd(innov))
        y_max <- max(h$density, curve_vals) * 1.05
        hist(innov, breaks = h$breaks, col = "darkorange", border = "white",
             main = sprintf("Innovation Histogram (%s)", dist_label),
             xlab = "Value", ylab = "Density", probability = TRUE,
             ylim = c(0, y_max))
        lines(x_seq, curve_vals, col = "steelblue", lwd = 2, lty = 2)
        legend("topright", legend = "Normal reference",
               col = "steelblue", lwd = 2, lty = 2, bty = "n", cex = 0.9)

        # 3. Squared innovations (shows volatility structure)
        plot(innov^2, type = "l", col = "firebrick", lwd = 1,
             xlab = "Time", ylab = expression(epsilon^2),
             main = sprintf("Squared Innovations (%s) - Volatility Structure", dist_label))
        abline(h = mean(innov^2), lty = 2, col = "grey40")
        legend("topright", legend = sprintf("Mean = %.3f", mean(innov^2)),
               bty = "n", cex = 0.9)
      } else {
        # IID: 2-panel layout (histogram + QQ plot)
        op <- par(mfrow = c(2, 1), mar = c(5, 5.5, 3, 1), cex.axis = 1.3, cex.lab = 1.4)
        on.exit(par(op))

        # 1. Histogram with theoretical density overlay
        h <- hist(innov, breaks = 40, plot = FALSE)
        x_seq <- seq(min(innov), max(innov), length.out = 300)
        theo <- .innov_theoretical_density(res$innov_dist, res$innov_params)
        if (!is.null(theo)) {
          curve_vals <- theo$dfun(x_seq)
        } else {
          curve_vals <- dnorm(x_seq, mean(innov), sd(innov))
        }
        y_max <- max(h$density, curve_vals) * 1.05
        hist(innov, breaks = h$breaks, col = "steelblue", border = "white",
             main = sprintf("Innovation Histogram (%s, n=%d)", dist_label, length(innov)),
             xlab = "Value", ylab = "Density", probability = TRUE,
             ylim = c(0, y_max))
        if (!is.null(theo)) {
          lines(x_seq, curve_vals, col = "red", lwd = 2)
          legend("topright",
                 legend = c(sprintf("Mean=%.3f, SD=%.3f", mean(innov), sd(innov)),
                            theo$label),
                 col = c(NA, "red"), lwd = c(NA, 2), lty = c(NA, 1),
                 bty = "n", cex = 0.9)
        } else {
          lines(x_seq, curve_vals, col = "red", lwd = 2, lty = 2)
          legend("topright",
                 legend = c(sprintf("Mean=%.3f, SD=%.3f", mean(innov), sd(innov)),
                            "Normal reference"),
                 col = c(NA, "red"), lwd = c(NA, 2), lty = c(NA, 2),
                 bty = "n", cex = 0.9)
        }

        # 2. QQ plot against normal
        qqnorm(innov, main = sprintf("Normal QQ Plot (%s)", dist_label),
               col = "steelblue", pch = 16, cex = 0.6)
        qqline(innov, col = "red", lwd = 2)
      }
    })

    # Null Model Diagnostics (AR order, variance, rejection by order)
    output$sim_null_diag <- renderPlot(bg = "transparent", {
      op <- viewer_base_par()
      on.exit(par(op), add = TRUE)
      res <- sim_results_rv()
      if (is.null(res)) {
        plot.new()
        text(0.5, 0.5, "Run a simulation to see null model diagnostics.",
             cex = 1.2, col = "grey50")
        return()
      }

      results <- res$results
      nsims <- res$nsims
      fg <- viewer_plot_fg()
      reject_mode <- input$sim_null_reject_mode %||% "all"
      plot_null_model_diagnostics(results, nsims = nsims, maxp = res$maxp,
                                  min_p = res$min_p, fg = fg,
                                  reject_mode = reject_mode)
    })

    sim_null_df <- reactive({
      res <- sim_results_rv()
      if (is.null(res) || is.null(res$results) || length(res$results) == 0) {
        return(data.frame())
      }

      rows <- lapply(res$results, function(sim) {
        phi_hat <- if (is.null(sim$phi) || length(sim$phi) == 0) NA_real_ else sim$phi[1]
        data.frame(
          null_ar_order = sim$p,
          null_phi1 = phi_hat,
          obs_stat = sim$tco_obs %||% sim$obs_stat %||% NA_real_,
          pvalue = sim$pvalue,
          pvalue_asymp = sim$pvalue_asymp,
          pvalue_adj = if (is.null(sim$pvalue_adj)) NA_real_ else sim$pvalue_adj,
          n = as.integer(res$n %||% NA_integer_),
          phi = as.numeric(res$phi %||% NA_real_),
          innov_dist = as.character(res$innov_dist %||% "Ad-Hoc"),
          stringsAsFactors = FALSE
        )
      })

      do.call(rbind, rows)
    })

    output$sim_null_phi_dist <- renderPlot(bg = "transparent", {
      df <- sim_null_df()
      if (nrow(df) == 0) {
        plot.new(); text(0.5, 0.5, "No data", cex = 1.2, col = viewer_plot_fg()); return()
      }
      p <- plot_ar_coefficient_distribution(df)
      if (!is.null(p)) {
        print(p)
      } else {
        plot.new(); text(0.5, 0.5, "No AR(1)-selected fits for this scenario", cex = 1.1, col = viewer_plot_fg())
      }
    })

    output$sim_null_tstat <- renderPlot(bg = "transparent", {
      df <- sim_null_df()
      if (nrow(df) == 0 || all(is.na(df$obs_stat))) {
        plot.new(); text(0.5, 0.5, "No data", cex = 1.2, col = viewer_plot_fg()); return()
      }
      p <- plot_test_statistic_distribution(df[!is.na(df$obs_stat), ])
      if (!is.null(p)) print(p)
    })

    output$sim_null_convergence <- renderPlot(bg = "transparent", {
      df <- sim_null_df()
      pval_col <- input$sim_null_conv_method %||% "pvalue"
      if (nrow(df) == 0 || !pval_col %in% names(df)) {
        plot.new(); text(0.5, 0.5, "No data", cex = 1.2, col = viewer_plot_fg()); return()
      }
      method_labels <- c(pvalue = "Bootstrap (COB)",
                         pvalue_asymp = "Asymptotic (CO)",
                         pvalue_adj = "COBA")
      p <- plot_mc_convergence(df, pval_col = pval_col,
                               method_label = method_labels[pval_col])
      if (!is.null(p)) print(p)
    })

    # --- Performance Profile sub-tab ---
    adhoc_params <- reactive({
      dist <- input$sim_innov_dist
      use_fast <- isTRUE(input$sim_use_fast)
      innov_gen <- tryCatch(
        .build_innov_gen(dist, input, use_fast = use_fast),
        error = function(e) NULL)

      phi <- input$sim_phi
      if (is.null(phi) || is.na(phi)) phi <- 0.95
      phi <- max(0, min(0.999, phi))

      list(
        n         = as.integer(input$sim_n %||% 200),
        phi       = phi,
        innov_gen = innov_gen,
        use_fast  = use_fast,
        nb        = as.integer(input$sim_nb %||% 399),
        maxp      = min(as.integer(input$sim_maxp %||% 5), 20L),
        criterion = match.arg(input$sim_criterion %||% "aic",
                              c("aic", "aicc", "bic")),
        bootadj   = isTRUE(input$sim_bootadj),
        min_p     = if (isTRUE(input$sim_minp1)) 1L else 0L
      )
    })

    mod_adhoc_profile_server(input, output, session, adhoc_params = adhoc_params)

    output$sim_run_log <- renderText({
      log_text <- run_log_rv()
      if (!nzchar(log_text)) "No ad-hoc runs yet." else log_text
    })

  })
}
