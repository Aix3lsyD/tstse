# Module: Innovation Comparison (Tab 9)
# Compare how different innovation distributions affect AR(1) realizations.
# Features: 6 configurable slots, 3-state plot toggle (AR series / innovations / overlay).

# --- Module-level constants ---------------------------------------------------

.INNOV_CHOICES <- c("Normal", "Student's t", "Skew-t", "GED", "Laplace",
                     "Uniform", "Mixture Normal", "GARCH", "Heteroscedastic")

.INNOV_TYPE_CLASS <- c(
  "Normal" = "iid", "Student's t" = "iid", "Skew-t" = "iid",
  "GED" = "iid", "Laplace" = "iid", "Uniform" = "iid",
  "Mixture Normal" = "iid", "GARCH" = "time_dep", "Heteroscedastic" = "time_dep"
)

.INNOV_COLORS <- c("#e41a1c", "#ff7f00", "#984ea3", "#4daf4a", "#377eb8", "#a65628")

.SLOT_DEFAULTS <- c("Student's t", "GARCH", "Heteroscedastic",
                     "Laplace", "Mixture Normal", "Uniform")

.innov_picker_input <- function(inputId, label, choices, selected = NULL, options = list()) {
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

.innov_toggle_input <- function(inputId, label, value = FALSE) {
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

.innov_action_button <- function(inputId, label, icon_name = NULL,
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

# --- UI helper: per-slot parameter panels ------------------------------------

.innov_slot_param_ui <- function(ns, i) {
  sid <- as.character(i)  # slot id suffix
  type_id <- paste0("slot_type_", sid)

  tagList(
    # Normal
    conditionalPanel(
      condition = sprintf("input.%s == 'Normal'", type_id), ns = ns,
      numericInput(ns(paste0("slot_norm_sd_", sid)), "sd:", value = 1,
                   min = 0.001, step = 0.1)
    ),
    # Student's t
    conditionalPanel(
      condition = sprintf("input.%s == \"Student's t\"", type_id), ns = ns,
      numericInput(ns(paste0("slot_t_df_", sid)), "df:", value = 3,
                   min = 1, step = 1),
      checkboxInput(ns(paste0("slot_t_scale_", sid)),
                    "Scale to unit variance", value = FALSE)
    ),
    # Skew-t
    conditionalPanel(
      condition = sprintf("input.%s == 'Skew-t'", type_id), ns = ns,
      numericInput(ns(paste0("slot_skt_df_", sid)), "df:", value = 5,
                   min = 3, step = 1),
      numericInput(ns(paste0("slot_skt_alpha_", sid)), "Skewness (alpha):",
                   value = 0, step = 0.1),
      checkboxInput(ns(paste0("slot_skt_scale_", sid)),
                    "Scale to unit variance", value = FALSE)
    ),
    # GED
    conditionalPanel(
      condition = sprintf("input.%s == 'GED'", type_id), ns = ns,
      numericInput(ns(paste0("slot_ged_nu_", sid)), "Shape (nu):",
                   value = 2, min = 0.1, step = 0.1),
      numericInput(ns(paste0("slot_ged_sd_", sid)), "sd:",
                   value = 1, min = 0.001, step = 0.1)
    ),
    # Laplace
    conditionalPanel(
      condition = sprintf("input.%s == 'Laplace'", type_id), ns = ns,
      numericInput(ns(paste0("slot_lap_scale_", sid)), "Scale:",
                   value = 0.707, min = 0.001, step = 0.01)
    ),
    # Uniform
    conditionalPanel(
      condition = sprintf("input.%s == 'Uniform'", type_id), ns = ns,
      numericInput(ns(paste0("slot_unif_hw_", sid)), "Half-width:",
                   value = 1.732, min = 0.001, step = 0.1)
    ),
    # Mixture Normal
    conditionalPanel(
      condition = sprintf("input.%s == 'Mixture Normal'", type_id), ns = ns,
      numericInput(ns(paste0("slot_mix_sd1_", sid)), "sd1:",
                   value = 1, min = 0.001, step = 0.1),
      numericInput(ns(paste0("slot_mix_sd2_", sid)), "sd2:",
                   value = 3, min = 0.001, step = 0.1),
      numericInput(ns(paste0("slot_mix_prob1_", sid)),
                   "prob1 (weight of component 1):",
                   value = 0.9, min = 0.01, max = 0.99, step = 0.05)
    ),
    # GARCH
    conditionalPanel(
      condition = sprintf("input.%s == 'GARCH'", type_id), ns = ns,
      numericInput(ns(paste0("slot_garch_omega_", sid)), "omega:",
                   value = 0.1, min = 0.001, step = 0.01),
      textAreaInput(ns(paste0("slot_garch_alpha_", sid)),
                    "alpha (comma-separated):",
                    value = "0.2, 0.175, 0.15, 0.125, 0.1, 0.075, 0.05, 0.025",
                    rows = 2),
      textInput(ns(paste0("slot_garch_beta_", sid)),
                "beta (comma-separated, optional):", value = "")
    ),
    # Heteroscedastic
    conditionalPanel(
      condition = sprintf("input.%s == 'Heteroscedastic'", type_id), ns = ns,
      .innov_picker_input(ns(paste0("slot_hetero_shape_", sid)), "Weight Shape:",
                          choices = c("linear", "sqrt", "log", "power", "exp",
                                      "step", "periodic"),
                          selected = "linear"),
      # from/to for linear, sqrt, log, power, exp
      conditionalPanel(
        condition = sprintf(
          "input.%s == 'linear' || input.%s == 'sqrt' || input.%s == 'log' || input.%s == 'power' || input.%s == 'exp'",
          paste0("slot_hetero_shape_", sid), paste0("slot_hetero_shape_", sid),
          paste0("slot_hetero_shape_", sid), paste0("slot_hetero_shape_", sid),
          paste0("slot_hetero_shape_", sid)),
        ns = ns,
        numericInput(ns(paste0("slot_hetero_from_", sid)), "From:",
                     value = 1, min = 0.01, step = 0.1),
        numericInput(ns(paste0("slot_hetero_to_", sid)), "To:",
                     value = 10, min = 0.01, step = 0.1)
      ),
      # power exponent
      conditionalPanel(
        condition = sprintf("input.%s == 'power'",
                            paste0("slot_hetero_shape_", sid)),
        ns = ns,
        numericInput(ns(paste0("slot_hetero_power_", sid)), "Power:",
                     value = 2, min = 0.1, step = 0.1)
      ),
      # step breaks/levels
      conditionalPanel(
        condition = sprintf("input.%s == 'step'",
                            paste0("slot_hetero_shape_", sid)),
        ns = ns,
        textInput(ns(paste0("slot_hetero_breaks_", sid)),
                  "Breaks (comma-separated, 0-1):", value = "0.5"),
        textInput(ns(paste0("slot_hetero_levels_", sid)),
                  "Levels (comma-separated, SD weights):", value = "1, 5")
      ),
      # periodic params
      conditionalPanel(
        condition = sprintf("input.%s == 'periodic'",
                            paste0("slot_hetero_shape_", sid)),
        ns = ns,
        numericInput(ns(paste0("slot_hetero_base_w_", sid)), "Base Weight:",
                     value = 1, min = 0.01, step = 0.1),
        numericInput(ns(paste0("slot_hetero_amplitude_", sid)), "Amplitude:",
                     value = 0.5, step = 0.1),
        numericInput(ns(paste0("slot_hetero_period_", sid)),
                     "Period (observations):", value = 12, min = 1, step = 1)
      ),
      numericInput(ns(paste0("slot_hetero_sd_", sid)), "sd (base normal):",
                   value = 1, min = 0.001, step = 0.1)
    )
  )
}

# --- UI function --------------------------------------------------------------

mod_innovation_comp_ui <- function(id) {
  ns <- NS(id)

  tabPanel(
    "Innovation Comparison",
    br(),
    h4("Innovation Distribution Comparison"),
    p(class = "text-body-secondary", style = "margin-bottom: 15px;",
      "Compare how different innovation distributions affect AR(1) realizations. ",
      "All plots use the same seed, phi, and n so differences are purely due ",
      "to the innovation distribution. Click any plot to cycle views."),

    fluidRow(
      # --- Controls sidebar ---
      column(3, style = "max-height: 85vh; overflow-y: auto;",
        wellPanel(
          h5("Parameters"),
          numericInput(ns("innov_phi"), "AR(1) Coefficient (phi):",
                       value = 0.95, min = 0.01, max = 0.999, step = 0.01),
          numericInput(ns("innov_n"), "Sample Size (n):",
                       value = 200, min = 10, max = 2000, step = 10),
          numericInput(ns("innov_seed"), "Seed:", value = 42, min = 1, step = 1),
          .innov_toggle_input(ns("innov_compact_ref"),
                              "Compact reference plot", value = FALSE),
          .innov_toggle_input(ns("innov_global_innov_view"),
                              "Show all innovation views", value = FALSE),
          .innov_action_button(ns("innov_generate"), "Generate",
                               icon_name = "play", class = "btn-primary",
                               full_width = TRUE)
        ),

        # Slot configuration accordion
        bslib::accordion(
          id = ns("slot_accordion"),
          open = FALSE,
          multiple = TRUE,

          bslib::accordion_panel(
            "Slot 1",
            .innov_picker_input(ns("slot_type_1"), "Distribution:",
                                .INNOV_CHOICES, selected = .SLOT_DEFAULTS[1]),
            .innov_slot_param_ui(ns, 1)
          ),
          bslib::accordion_panel(
            "Slot 2",
            .innov_picker_input(ns("slot_type_2"), "Distribution:",
                                .INNOV_CHOICES, selected = .SLOT_DEFAULTS[2]),
            .innov_slot_param_ui(ns, 2)
          ),
          bslib::accordion_panel(
            "Slot 3",
            .innov_picker_input(ns("slot_type_3"), "Distribution:",
                                .INNOV_CHOICES, selected = .SLOT_DEFAULTS[3]),
            .innov_slot_param_ui(ns, 3)
          ),
          bslib::accordion_panel(
            "Slot 4",
            .innov_picker_input(ns("slot_type_4"), "Distribution:",
                                .INNOV_CHOICES, selected = .SLOT_DEFAULTS[4]),
            .innov_slot_param_ui(ns, 4)
          ),
          bslib::accordion_panel(
            "Slot 5",
            .innov_picker_input(ns("slot_type_5"), "Distribution:",
                                .INNOV_CHOICES, selected = .SLOT_DEFAULTS[5]),
            .innov_slot_param_ui(ns, 5)
          ),
          bslib::accordion_panel(
            "Slot 6",
            .innov_picker_input(ns("slot_type_6"), "Distribution:",
                                .INNOV_CHOICES, selected = .SLOT_DEFAULTS[6]),
            .innov_slot_param_ui(ns, 6)
          )
        )
      ),

      # --- Plots area ---
      column(9,
        # Reference normal plot: full-width or compact
        conditionalPanel(
          condition = "!input.innov_compact_ref", ns = ns,
          plotOutput(ns("innov_ref_plot_full"), height = "350px",
                     click = ns("innov_click_ref"))
        ),
        conditionalPanel(
          condition = "input.innov_compact_ref", ns = ns,
          fluidRow(
            column(6, plotOutput(ns("innov_ref_plot_compact"), height = "270px",
                               click = ns("innov_click_ref")))
          )
        ),

        # 3x2 comparison grid
        fluidRow(
          column(6, plotOutput(ns("innov_comp_plot_1"), height = "270px",
                               click = ns("innov_click_1"))),
          column(6, plotOutput(ns("innov_comp_plot_2"), height = "270px",
                               click = ns("innov_click_2")))
        ),
        fluidRow(
          column(6, plotOutput(ns("innov_comp_plot_3"), height = "270px",
                               click = ns("innov_click_3"))),
          column(6, plotOutput(ns("innov_comp_plot_4"), height = "270px",
                               click = ns("innov_click_4")))
        ),
        fluidRow(
          column(6, plotOutput(ns("innov_comp_plot_5"), height = "270px",
                               click = ns("innov_click_5"))),
          column(6, plotOutput(ns("innov_comp_plot_6"), height = "270px",
                               click = ns("innov_click_6")))
        )
      )
    )
  )
}

# --- Server helpers -----------------------------------------------------------

# Build an innovation generator from slot inputs
.build_slot_gen <- function(dist_type, input, slot_id) {
  sid <- as.character(slot_id)
  .get <- function(suffix, default = NULL) {
    val <- input[[paste0("slot_", suffix, "_", sid)]]
    if (is.null(val) || length(val) == 0 || is.na(val)) default else val
  }

  switch(dist_type,
    "Normal" = {
      sd_val <- .get("norm_sd", 1)
      if (sd_val <= 0) sd_val <- 1
      make_gen_norm(sd = sd_val)
    },
    "Student's t" = {
      df_val <- .get("t_df", 3)
      if (df_val < 1) df_val <- 3
      scale_val <- isTRUE(.get("t_scale", FALSE))
      make_gen_t(df = df_val, scale = scale_val)
    },
    "Skew-t" = {
      df_val <- .get("skt_df", 5)
      if (df_val < 3) df_val <- 5
      alpha_val <- .get("skt_alpha", 0)
      scale_val <- isTRUE(.get("skt_scale", FALSE))
      make_gen_skt(df = df_val, alpha = alpha_val, scale = scale_val)
    },
    "GED" = {
      nu_val <- .get("ged_nu", 2)
      if (nu_val <= 0) nu_val <- 2
      sd_val <- .get("ged_sd", 1)
      if (sd_val <= 0) sd_val <- 1
      make_gen_ged(nu = nu_val, sd = sd_val)
    },
    "Laplace" = {
      sc <- .get("lap_scale", 1 / sqrt(2))
      if (sc <= 0) sc <- 1 / sqrt(2)
      make_gen_laplace(scale = sc)
    },
    "Uniform" = {
      hw <- .get("unif_hw", sqrt(3))
      if (hw <= 0) hw <- sqrt(3)
      make_gen_unif(half_width = hw)
    },
    "Mixture Normal" = {
      sd1 <- .get("mix_sd1", 1)
      sd2 <- .get("mix_sd2", 3)
      p1 <- .get("mix_prob1", 0.9)
      if (sd1 <= 0) sd1 <- 1
      if (sd2 <= 0) sd2 <- 3
      if (p1 <= 0 || p1 >= 1) p1 <- 0.9
      make_gen_mixnorm(sd1 = sd1, sd2 = sd2, prob1 = p1)
    },
    "GARCH" = {
      omega <- .get("garch_omega", 0.1)
      if (omega <= 0) omega <- 0.1
      alpha_str <- input[[paste0("slot_garch_alpha_", sid)]]
      alpha <- tryCatch(
        as.numeric(trimws(strsplit(alpha_str, ",")[[1]])),
        warning = function(w) NULL, error = function(e) NULL)
      if (is.null(alpha) || any(is.na(alpha)) || length(alpha) == 0)
        alpha <- c(0.2, 0.175, 0.15, 0.125, 0.1, 0.075, 0.05, 0.025)
      beta_str <- trimws(input[[paste0("slot_garch_beta_", sid)]] %||% "")
      beta <- NULL
      if (nzchar(beta_str)) {
        beta <- tryCatch(
          as.numeric(trimws(strsplit(beta_str, ",")[[1]])),
          warning = function(w) NULL, error = function(e) NULL)
        if (!is.null(beta) && (any(is.na(beta)) || length(beta) == 0))
          beta <- NULL
      }
      make_gen_garch(omega = omega, alpha = alpha, beta = beta)
    },
    "Heteroscedastic" = {
      sd_val <- .get("hetero_sd", 1)
      if (sd_val <= 0) sd_val <- 1
      shape <- input[[paste0("slot_hetero_shape_", sid)]] %||% "linear"
      switch(shape,
        "linear" =, "sqrt" =, "log" =, "exp" = {
          from_val <- .get("hetero_from", 1)
          to_val <- .get("hetero_to", 10)
          make_gen_hetero(shape = shape, from = from_val, to = to_val, sd = sd_val)
        },
        "power" = {
          from_val <- .get("hetero_from", 1)
          to_val <- .get("hetero_to", 10)
          p_val <- .get("hetero_power", 2)
          make_gen_hetero(shape = "power", from = from_val, to = to_val,
                          power = p_val, sd = sd_val)
        },
        "step" = {
          breaks_str <- input[[paste0("slot_hetero_breaks_", sid)]] %||% "0.5"
          levels_str <- input[[paste0("slot_hetero_levels_", sid)]] %||% "1, 5"
          brk <- as.numeric(trimws(strsplit(breaks_str, ",")[[1]]))
          lvl <- as.numeric(trimws(strsplit(levels_str, ",")[[1]]))
          if (any(is.na(brk))) brk <- 0.5
          if (any(is.na(lvl))) lvl <- c(1, 5)
          make_gen_hetero(shape = "step", breaks = brk, levels = lvl, sd = sd_val)
        },
        "periodic" = {
          bw <- .get("hetero_base_w", 1)
          amp <- .get("hetero_amplitude", 0.5)
          per <- .get("hetero_period", 12)
          make_gen_hetero(shape = "periodic", base_w = bw, amplitude = amp,
                          period = per, sd = sd_val)
        }
      )
    },
    stop("Unknown distribution type: ", dist_type)
  )
}

# Build a human-readable label from slot inputs
.build_slot_label <- function(dist_type, input, slot_id) {
  sid <- as.character(slot_id)
  .get <- function(suffix, default = NULL) {
    val <- input[[paste0("slot_", suffix, "_", sid)]]
    if (is.null(val) || length(val) == 0 || is.na(val)) default else val
  }

  switch(dist_type,
    "Normal" = sprintf("Normal(sd=%.2f)", .get("norm_sd", 1)),
    "Student's t" = {
      df_val <- .get("t_df", 3)
      if (isTRUE(.get("t_scale", FALSE)))
        sprintf("Student's t(df=%d, scaled)", df_val)
      else
        sprintf("Student's t(df=%d)", df_val)
    },
    "Skew-t" = sprintf("Skew-t(df=%d, alpha=%.1f)", .get("skt_df", 5),
                        .get("skt_alpha", 0)),
    "GED" = sprintf("GED(nu=%.1f, sd=%.2f)", .get("ged_nu", 2),
                     .get("ged_sd", 1)),
    "Laplace" = sprintf("Laplace(scale=%.3f)", .get("lap_scale", 0.707)),
    "Uniform" = sprintf("Uniform(hw=%.2f)", .get("unif_hw", 1.732)),
    "Mixture Normal" = sprintf("MixNorm(sd1=%.1f,sd2=%.1f,p=%.2f)",
                                .get("mix_sd1", 1), .get("mix_sd2", 3),
                                .get("mix_prob1", 0.9)),
    "GARCH" = {
      alpha_str <- input[[paste0("slot_garch_alpha_", sid)]]
      alpha <- tryCatch(
        as.numeric(trimws(strsplit(alpha_str, ",")[[1]])),
        warning = function(w) NULL, error = function(e) NULL)
      if (is.null(alpha)) alpha <- rep(0.2, 8)
      beta_str <- trimws(input[[paste0("slot_garch_beta_", sid)]] %||% "")
      if (nzchar(beta_str)) sprintf("GARCH(%d,%d)", length(alpha), 1)
      else sprintf("ARCH(%d)", length(alpha))
    },
    "Heteroscedastic" = {
      shape <- input[[paste0("slot_hetero_shape_", sid)]] %||% "linear"
      sd_val <- .get("hetero_sd", 1)
      switch(shape,
        "linear" =, "sqrt" =, "log" =, "exp" =
          sprintf("Hetero(%s,%s-%s,sd=%.2f)", shape,
                  .get("hetero_from", 1), .get("hetero_to", 10), sd_val),
        "power" =
          sprintf("Hetero(power,%s-%s,p=%s,sd=%.2f)",
                  .get("hetero_from", 1), .get("hetero_to", 10),
                  .get("hetero_power", 2), sd_val),
        "step" = {
          breaks_str <- input[[paste0("slot_hetero_breaks_", sid)]] %||% "0.5"
          levels_str <- input[[paste0("slot_hetero_levels_", sid)]] %||% "1,5"
          sprintf("Hetero(step,%s|%s,sd=%.2f)", breaks_str, levels_str, sd_val)
        },
        "periodic" =
          sprintf("Hetero(periodic,bw=%s,amp=%s,per=%s,sd=%.2f)",
                  .get("hetero_base_w", 1), .get("hetero_amplitude", 0.5),
                  .get("hetero_period", 12), sd_val)
      )
    },
    dist_type
  )
}

# --- Server function ----------------------------------------------------------

mod_innovation_comp_server <- function(id) {
  moduleServer(id, function(input, output, session) {

    # View state per slot: 0 = AR series, 1 = innovations, 2 = overlay
    view_states <- lapply(1:6, function(i) reactiveVal(0L))
    # Reference plot: 0 = AR series, 1 = histogram (no overlay needed)
    ref_view_state <- reactiveVal(0L)

    .parse_num_list <- function(x) {
      if (is.null(x)) return(numeric(0))
      txt <- trimws(as.character(x))
      if (!nzchar(txt)) return(numeric(0))
      vals <- suppressWarnings(as.numeric(trimws(strsplit(txt, ",")[[1]])))
      if (length(vals) == 0 || any(is.na(vals))) return(NULL)
      vals
    }

    has_shinyvalidate <- requireNamespace("shinyvalidate", quietly = TRUE)
    innov_validator <- NULL

    if (has_shinyvalidate) {
      v <- shinyvalidate::InputValidator$new()

      .rule_number <- function(expr, msg) {
        force(expr)
        function(value) {
          if (is.null(value) || is.na(value) || !isTRUE(expr(as.numeric(value)))) msg else NULL
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

      v$add_rule("innov_phi", .rule_number(function(x) x >= 0.01 && x <= 0.999,
                                             "Must be between 0.01 and 0.999"))
      v$add_rule("innov_n", .rule_number(function(x) as.integer(x) >= 10 && as.integer(x) <= 2000,
                                           "Must be an integer from 10 to 2000"))
      v$add_rule("innov_seed", .rule_number(function(x) as.integer(x) >= 1,
                                              "Must be an integer >= 1"))

      for (i in 1:6) {
        local({
          idx <- i
          slot_type_id <- paste0("slot_type_", idx)
          .is_slot <- function(label) identical(input[[slot_type_id]], label)
          .is_hetero_shape <- function(shape) {
            identical(input[[slot_type_id]], "Heteroscedastic") &&
              identical(input[[paste0("slot_hetero_shape_", idx)]], shape)
          }

          v$add_rule(paste0("slot_norm_sd_", idx), .rule_when(
            function() .is_slot("Normal"),
            .rule_number(function(x) x > 0, "Must be > 0")
          ))
          v$add_rule(paste0("slot_t_df_", idx), .rule_when(
            function() .is_slot("Student's t"),
            .rule_number(function(x) x >= 1, "Must be >= 1")
          ))
          v$add_rule(paste0("slot_skt_df_", idx), .rule_when(
            function() .is_slot("Skew-t"),
            .rule_number(function(x) x >= 3, "Must be >= 3")
          ))
          v$add_rule(paste0("slot_skt_alpha_", idx), .rule_when(
            function() .is_slot("Skew-t"),
            .rule_number(function(x) is.finite(x), "Must be numeric")
          ))
          v$add_rule(paste0("slot_ged_nu_", idx), .rule_when(
            function() .is_slot("GED"),
            .rule_number(function(x) x > 0, "Must be > 0")
          ))
          v$add_rule(paste0("slot_ged_sd_", idx), .rule_when(
            function() .is_slot("GED"),
            .rule_number(function(x) x > 0, "Must be > 0")
          ))
          v$add_rule(paste0("slot_lap_scale_", idx), .rule_when(
            function() .is_slot("Laplace"),
            .rule_number(function(x) x > 0, "Must be > 0")
          ))
          v$add_rule(paste0("slot_unif_hw_", idx), .rule_when(
            function() .is_slot("Uniform"),
            .rule_number(function(x) x > 0, "Must be > 0")
          ))
          v$add_rule(paste0("slot_mix_sd1_", idx), .rule_when(
            function() .is_slot("Mixture Normal"),
            .rule_number(function(x) x > 0, "Must be > 0")
          ))
          v$add_rule(paste0("slot_mix_sd2_", idx), .rule_when(
            function() .is_slot("Mixture Normal"),
            .rule_number(function(x) x > 0, "Must be > 0")
          ))
          v$add_rule(paste0("slot_mix_prob1_", idx), .rule_when(
            function() .is_slot("Mixture Normal"),
            .rule_number(function(x) x > 0 && x < 1, "Must be between 0 and 1")
          ))

          v$add_rule(paste0("slot_garch_omega_", idx), .rule_when(
            function() .is_slot("GARCH"),
            .rule_number(function(x) x > 0, "Must be > 0")
          ))
          v$add_rule(paste0("slot_garch_alpha_", idx), .rule_when(
            function() .is_slot("GARCH"),
            .rule_csv_nonneg(required = TRUE)
          ))
          v$add_rule(paste0("slot_garch_beta_", idx), .rule_when(
            function() .is_slot("GARCH"),
            .rule_csv_nonneg(required = FALSE)
          ))

          v$add_rule(paste0("slot_hetero_sd_", idx), .rule_when(
            function() identical(input[[slot_type_id]], "Heteroscedastic"),
            .rule_number(function(x) x > 0, "Must be > 0")
          ))
          v$add_rule(paste0("slot_hetero_from_", idx), .rule_when(
            function() {
              identical(input[[slot_type_id]], "Heteroscedastic") &&
                input[[paste0("slot_hetero_shape_", idx)]] %in% c("linear", "sqrt", "log", "power", "exp")
            },
            .rule_number(function(x) x > 0, "Must be > 0")
          ))
          v$add_rule(paste0("slot_hetero_to_", idx), .rule_when(
            function() {
              identical(input[[slot_type_id]], "Heteroscedastic") &&
                input[[paste0("slot_hetero_shape_", idx)]] %in% c("linear", "sqrt", "log", "power", "exp")
            },
            function(value) {
              to_val <- suppressWarnings(as.numeric(value))
              from_val <- suppressWarnings(as.numeric(input[[paste0("slot_hetero_from_", idx)]]))
              if (is.na(to_val) || to_val <= 0) return("Must be > 0")
              if (!is.na(from_val) && to_val < from_val) return("Must be >= From")
              NULL
            }
          ))
          v$add_rule(paste0("slot_hetero_power_", idx), .rule_when(
            function() .is_hetero_shape("power"),
            .rule_number(function(x) x > 0, "Must be > 0")
          ))
          v$add_rule(paste0("slot_hetero_breaks_", idx), .rule_when(
            function() .is_hetero_shape("step"),
            function(value) {
              brk <- .parse_num_list(value)
              if (is.null(brk) || length(brk) == 0) return("Provide comma-separated numeric breaks")
              if (any(brk <= 0 | brk >= 1)) return("Breaks must be between 0 and 1")
              if (is.unsorted(brk, strictly = TRUE)) return("Breaks must be strictly increasing")
              NULL
            }
          ))
          v$add_rule(paste0("slot_hetero_levels_", idx), .rule_when(
            function() .is_hetero_shape("step"),
            function(value) {
              lvl <- .parse_num_list(value)
              if (is.null(lvl) || length(lvl) == 0) return("Provide comma-separated numeric levels")
              if (any(lvl <= 0)) return("Levels must be > 0")
              brk <- .parse_num_list(input[[paste0("slot_hetero_breaks_", idx)]])
              if (!is.null(brk) && length(brk) > 0 && length(lvl) != length(brk) + 1) {
                return("Levels count must equal breaks count + 1")
              }
              NULL
            }
          ))
          v$add_rule(paste0("slot_hetero_base_w_", idx), .rule_when(
            function() .is_hetero_shape("periodic"),
            .rule_number(function(x) x > 0, "Must be > 0")
          ))
          v$add_rule(paste0("slot_hetero_period_", idx), .rule_when(
            function() .is_hetero_shape("periodic"),
            .rule_number(function(x) as.integer(x) >= 1, "Must be an integer >= 1")
          ))
        })
      }

      v$enable()
      innov_validator <- v
    }

    # Click handler for reference plot (toggle between 0 and 1)
    observeEvent(input$innov_click_ref, {
      ref_view_state((ref_view_state() + 1L) %% 2L)
    })

    # Click handlers to cycle view state for comparison slots
    local({
      for (idx in 1:6) {
        local({
          i <- idx
          observeEvent(input[[paste0("innov_click_", i)]], {
            current <- view_states[[i]]()
            view_states[[i]]((current + 1L) %% 3L)
          })
        })
      }
    })

    # Global toggle: switch all to innovation view (1) or back to AR (0)
    observeEvent(input$innov_global_innov_view, {
      target <- if (isTRUE(input$innov_global_innov_view)) 1L else 0L
      ref_view_state(target)
      for (i in 1:6) view_states[[i]](target)
    }, ignoreInit = TRUE)

    # Main reactive: generate all realizations on button click
    innov_realizations <- eventReactive(input$innov_generate, {
      if (!is.null(innov_validator) && !isTRUE(innov_validator$is_valid())) {
        showNotification("Please fix highlighted input errors before generating.",
                         type = "warning", duration = 5)
        return(NULL)
      }

      phi <- input$innov_phi
      n <- as.integer(input$innov_n)
      if (is.na(n) || n < 10) n <- 10L
      if (n > 2000) n <- 2000L
      seed <- input$innov_seed
      if (is.null(phi) || is.na(phi)) phi <- 0.95
      phi <- max(0.01, min(0.999, phi))

      withProgress(message = "Generating realizations...", value = 0, {
        # Reference: normal
        incProgress(1 / 7, detail = "Normal (reference)")
        ref_gen <- make_gen_norm(sd = 1)
        ref_y <- gen_aruma_flex(n, phi = phi, innov_gen = ref_gen,
                                seed = seed, plot = FALSE)$y
        set.seed(seed)
        ref_innov <- ref_gen(n)

        # 6 configurable slots
        slots <- vector("list", 6)
        for (i in 1:6) {
          incProgress(1 / 7, detail = sprintf("Slot %d", i))
          dist_type <- input[[paste0("slot_type_", i)]]
          if (is.null(dist_type)) dist_type <- .SLOT_DEFAULTS[i]

          slots[[i]] <- tryCatch({
            gen <- .build_slot_gen(dist_type, input, i)
            y <- gen_aruma_flex(n, phi = phi, innov_gen = gen,
                                seed = seed, plot = FALSE)$y
            # Raw innovations (separate draw with same seed)
            set.seed(seed)
            raw_innov <- gen(n)

            list(y = y, innov = raw_innov,
                 label = .build_slot_label(dist_type, input, i),
                 color = .INNOV_COLORS[i], error = NULL,
                 gen_type = unname(.INNOV_TYPE_CLASS[dist_type]))
          }, error = function(e) {
            list(y = NULL, innov = NULL,
                 label = dist_type, color = .INNOV_COLORS[i],
                 error = e$message, gen_type = "iid")
          })
        }
      })

      list(ref_y = ref_y, ref_innov = ref_innov,
           slots = slots, phi = phi, n = n, seed = seed)
    })

    # Helper: render the reference normal plot (2-state: AR series or histogram)
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

      state <- ref_view_state()
      state_labels <- c("AR series", "Innovations")

      if (state == 0L) {
        y <- data$ref_y
        plot(y, type = "l", col = "steelblue", lwd = 1.2,
             xlab = "Time", ylab = "Value", ylim = range(y),
             main = sprintf("Normal (Reference)  |  phi=%.2f, n=%d, seed=%d",
                            data$phi, data$n, data$seed))
        abline(h = 0, lty = 3, col = "grey60")
      } else {
        innov <- data$ref_innov
        hist(innov, breaks = 30,
             col = adjustcolor("steelblue", alpha.f = 0.6),
             border = "white", probability = TRUE,
             main = sprintf("Normal (Reference) [innovations]  |  n=%d, seed=%d",
                            data$n, data$seed),
             xlab = "Value", ylab = "Density")
        x_seq <- seq(min(innov), max(innov), length.out = 200)
        lines(x_seq, dnorm(x_seq, 0, 1),
              col = "grey40", lwd = 1.5, lty = 2)
      }

      mtext(paste0("[", state_labels[state + 1L], "] click to cycle"),
            side = 1, adj = 1, line = -1, cex = 0.65, col = "grey50")
    }

    # Reference plots (full and compact)
    output$innov_ref_plot_full <- renderPlot(bg = "transparent", {
      .render_innov_ref()
    })
    output$innov_ref_plot_compact <- renderPlot(bg = "transparent", {
      .render_innov_ref()
    })

    # Comparison plots (slots 1-6) with 3-state view toggle
    local({
      for (idx in 1:6) {
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
              text(0.5, 0.5, slot$error, cex = 0.9, col = "grey50")
              return()
            }

            state <- view_states[[i]]()
            state_labels <- c("AR series", "Innovations", "Overlay")

            if (state == 0L) {
              # State 0: AR-filtered time series (default)
              y <- slot$y
              plot(y, type = "l", col = slot$color, lwd = 1.2,
                   xlab = "Time", ylab = "Value", ylim = range(y),
                   main = slot$label)
              abline(h = 0, lty = 3, col = "grey60")

            } else if (state == 1L) {
              # State 1: Innovation view
              if (slot$gen_type == "iid") {
                # Histogram for IID distributions
                hist(slot$innov, breaks = 30,
                     col = adjustcolor(slot$color, alpha.f = 0.6),
                     border = "white", probability = TRUE,
                     main = paste0(slot$label, "  [innovations]"),
                     xlab = "Value", ylab = "Density")
                # Normal reference density curve
                x_seq <- seq(min(slot$innov), max(slot$innov), length.out = 200)
                lines(x_seq, dnorm(x_seq, 0, 1),
                      col = "grey40", lwd = 1.5, lty = 2)
              } else {
                # Time plot for time-dependent distributions
                innov <- slot$innov
                plot(innov, type = "l", col = slot$color, lwd = 1.0,
                     xlab = "Time", ylab = "Innovation",
                     ylim = range(innov),
                     main = paste0(slot$label, "  [innovations]"))
                abline(h = 0, lty = 3, col = "grey60")
              }

            } else if (state == 2L) {
              # State 2: Overlay with normal reference
              y_range <- range(c(data$ref_y, slot$y))
              plot(data$ref_y, type = "l",
                   col = adjustcolor("grey50", alpha.f = 0.4), lwd = 1.5,
                   xlab = "Time", ylab = "Value", ylim = y_range,
                   main = paste0(slot$label, "  [overlay]"))
              lines(slot$y, col = slot$color, lwd = 1.2)
              abline(h = 0, lty = 3, col = "grey60")
              legend("topright",
                     legend = c("Normal ref", slot$label),
                     col = c(adjustcolor("grey50", alpha.f = 0.4), slot$color),
                     lwd = c(1.5, 1.2), bty = "n", cex = 0.8)
            }

            # View state indicator
            mtext(paste0("[", state_labels[state + 1L], "] click to cycle"),
                  side = 1, adj = 1, line = -1, cex = 0.65, col = "grey50")
          })
        })
      }
    })

  })
}
