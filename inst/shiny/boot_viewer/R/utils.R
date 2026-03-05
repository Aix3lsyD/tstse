# Shared helper functions for boot_viewer modules

# Read root dark-mode toggle from any module/server reactive context.
viewer_dark_mode <- function(session = shiny::getDefaultReactiveDomain()) {
  if (is.null(session)) return(NA)
  root_session <- tryCatch(session$rootScope(), error = function(e) session)
  dark_val <- tryCatch(root_session$input$dark_mode, error = function(e) NULL)
  if (is.null(dark_val) || length(dark_val) == 0 || is.na(dark_val)) return(NA)
  if (is.character(dark_val)) {
    mode <- tolower(trimws(dark_val[[1]]))
    if (mode %in% c("dark", "light")) return(identical(mode, "dark"))
  }
  if (is.logical(dark_val)) return(isTRUE(dark_val[[1]]))
  if (is.numeric(dark_val)) return(isTRUE(as.logical(dark_val[[1]])))
  NA
}

# Resolve a thematic color option, with fallback for non-thematic contexts.
resolve_thematic_color <- function(color_opt, fallback) {
  if (is.null(color_opt) || length(color_opt) == 0) return(fallback)
  if (inherits(color_opt, "thematic_auto")) return(fallback)
  color_val <- as.character(color_opt[[1]])
  if (is.na(color_val) || color_val == "" || identical(color_val, "auto")) {
    return(fallback)
  }
  color_val
}

# Resolve the active foreground color for plot text in light/dark mode.
viewer_plot_fg <- function(session = shiny::getDefaultReactiveDomain()) {
  dark_mode <- viewer_dark_mode(session)
  if (isTRUE(dark_mode)) return("#e9ecef")
  if (identical(dark_mode, FALSE)) return("#212529")
  resolve_thematic_color(thematic::thematic_get_option("fg"), "#212529")
}

# Shared base-graphics defaults so non-ggplot charts use larger, consistent text.
viewer_base_par <- function(scale = 1.3) {
  old <- graphics::par(no.readonly = TRUE)
  graphics::par(
    cex = scale,
    cex.axis = scale,
    cex.lab = scale * 1.05,
    cex.main = scale * 1.12,
    cex.sub = scale
  )
  old
}

# Shared ggplot theme that respects thematic_shiny() dark/light palettes.
viewer_plot_theme <- function(base_size = 14, session = shiny::getDefaultReactiveDomain()) {
  effective_base_size <- base_size + 4
  fg <- viewer_plot_fg(session = session)
  grid_major <- grDevices::adjustcolor(fg, alpha.f = 0.16)
  grid_minor <- grDevices::adjustcolor(fg, alpha.f = 0.08)
  axis_tick <- grDevices::adjustcolor(fg, alpha.f = 0.45)

  ggplot2::theme_minimal(base_size = effective_base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = fg),
      axis.text = ggplot2::element_text(colour = fg, size = effective_base_size * 0.9),
      axis.title = ggplot2::element_text(colour = fg, size = effective_base_size * 1.02),
      plot.title = ggplot2::element_text(colour = fg, size = effective_base_size * 1.16),
      plot.subtitle = ggplot2::element_text(colour = fg),
      plot.caption = ggplot2::element_text(colour = fg),
      strip.text = ggplot2::element_text(colour = fg, size = effective_base_size, face = "bold"),
      legend.title = ggplot2::element_text(colour = fg, size = effective_base_size * 0.95),
      legend.text = ggplot2::element_text(colour = fg, size = effective_base_size * 0.9),
      panel.background = ggplot2::element_rect(fill = "transparent", colour = NA),
      plot.background = ggplot2::element_rect(fill = "transparent", colour = NA),
      legend.background = ggplot2::element_rect(fill = "transparent", colour = NA),
      legend.key = ggplot2::element_rect(fill = "transparent", colour = NA),
      panel.grid.major = ggplot2::element_line(colour = grid_major),
      panel.grid.minor = ggplot2::element_line(colour = grid_minor),
      axis.ticks = ggplot2::element_line(colour = axis_tick)
    )
}

# Format a grid data frame as an interactive DT with color-coded rate columns
format_grid_dt <- function(df, row_label_col, row_label_name) {
  if (nrow(df) == 0) {
    return(datatable(data.frame(Message = "No data for this combination"),
                     rownames = FALSE, options = list(dom = "t")))
  }

  display_df <- data.frame(
    Label    = df[[row_label_col]],
    n_sims   = df$n_sims,
    CO       = df$reject_asymp_05,
    CO_SE    = df$reject_asymp_05_se,
    COB      = df$reject_05,
    COB_SE   = df$reject_05_se,
    COBA     = df$reject_adj_05,
    COBA_SE  = df$reject_adj_05_se
  )

  col_names <- c(row_label_name, "Sims", "CO Rate", "CO SE",
                 "COB Rate", "COB SE", "COBA Rate", "COBA SE")

  dt <- datatable(
    display_df,
    colnames = col_names,
    rownames = FALSE,
    options = list(
      dom = "t",
      paging = FALSE,
      searching = FALSE,
      ordering = FALSE,
      columnDefs = list(
        list(targets = c(3, 5, 7), className = "dt-body-right")
      )
    )
  )

  rate_cols <- c("CO", "COB", "COBA")
  se_cols <- c("CO_SE", "COB_SE", "COBA_SE")

  dt <- formatRound(dt, columns = rate_cols, digits = 4)
  dt <- formatRound(dt, columns = se_cols, digits = 4)

  for (col in rate_cols) {
    dt <- formatStyle(dt, col,
      backgroundColor = styleInterval(
        c(0.03, 0.07),
        c("#fff3cd", "#d4edda", "#f8d7da")
      )
    )
  }

  dt
}

# =============================================================================
# Shared Innovation Generator Helpers
# =============================================================================

# Build innov_dist string from distribution label + params list.
# params uses generic keys (no "sim_" prefix).
# Returns a self-describing string for DuckDB storage.
build_innov_dist_str <- function(label, params = list()) {
  switch(label,
    "Normal" = "norm",
    "Student's t" = sprintf("t(%s)", params$t_df %||% 3),
    "Skew-t" = sprintf("skt(%s,%s)", params$skt_df %||% 5, params$skt_alpha %||% 0),
    "GED" = sprintf("ged(%s)", params$ged_nu %||% 2),
    "Laplace" = "laplace",
    "Uniform" = "unif",
    "Mixture Normal" = sprintf("mixnorm(%s,%s,%s)",
                                params$mix_sd1 %||% 1, params$mix_sd2 %||% 3,
                                params$mix_prob1 %||% 0.9),
    "GARCH" = {
      alpha <- params$garch_alpha
      beta <- params$garch_beta
      if (is.character(alpha)) {
        alpha <- tryCatch(
          as.numeric(trimws(strsplit(alpha, ",")[[1]])),
          warning = function(w) NULL, error = function(e) NULL)
      }
      has_beta <- !is.null(beta)
      if (is.character(beta)) {
        beta_str <- trimws(beta)
        has_beta <- nzchar(beta_str) && beta_str != ""
        if (has_beta) {
          beta <- tryCatch(as.numeric(trimws(strsplit(beta_str, ",")[[1]])),
                           warning = function(w) NULL, error = function(e) NULL)
          has_beta <- !is.null(beta)
        }
      }
      if (has_beta) {
        sprintf("garch(%s,%s)", params$garch_omega %||% 0.1,
                paste(c(alpha, beta), collapse = ","))
      } else {
        sprintf("arch(%s)", paste(alpha, collapse = ","))
      }
    },
    "Heteroscedastic" = {
      shape <- params$hetero_shape %||% "linear"
      sd_val <- params$hetero_sd %||% 1
      switch(shape,
        "linear" =, "sqrt" =, "log" =, "exp" =
          sprintf("hetero(%s,%s-%s,sd=%s)", shape,
                  params$hetero_from %||% 1, params$hetero_to %||% 10, sd_val),
        "power" =
          sprintf("hetero(power,%s-%s,p=%s,sd=%s)",
                  params$hetero_from %||% 1, params$hetero_to %||% 10,
                  params$hetero_power %||% 2, sd_val),
        "step" =
          sprintf("hetero(step,%s|%s,sd=%s)",
                  params$hetero_breaks %||% "0.5",
                  params$hetero_levels %||% "1,5", sd_val),
        "periodic" =
          sprintf("hetero(periodic,bw=%s,amp=%s,per=%s,sd=%s)",
                  params$hetero_base_w %||% 1,
                  params$hetero_amplitude %||% 0.5,
                  params$hetero_period %||% 12, sd_val)
      )
    },
    "unknown"
  )
}

# Build innov_gen function from distribution label + params list.
# params uses generic keys (no "sim_" prefix).
# Returns a function(n) -> numeric(n) innovation generator.
build_innov_gen <- function(label, params = list(), use_fast = FALSE) {
  switch(label,
    "Normal" = {
      sd_val <- params$norm_sd
      if (is.null(sd_val) || is.na(sd_val) || sd_val <= 0) sd_val <- 1
      if (isTRUE(use_fast)) make_gen_norm_fast(sd = sd_val) else make_gen_norm(sd = sd_val)
    },
    "Student's t" = {
      df_val <- params$t_df
      if (is.null(df_val) || is.na(df_val) || df_val < 1) df_val <- 3
      scale_val <- isTRUE(params$t_scale)
      if (isTRUE(use_fast)) make_gen_t_fast(df = df_val, scale = scale_val) else make_gen_t(df = df_val, scale = scale_val)
    },
    "Skew-t" = {
      df_val <- params$skt_df
      if (is.null(df_val) || is.na(df_val) || df_val < 3) df_val <- 5
      alpha_val <- params$skt_alpha
      if (is.null(alpha_val) || is.na(alpha_val)) alpha_val <- 0
      scale_val <- isTRUE(params$skt_scale)
      if (isTRUE(use_fast)) make_gen_skt_fast(df = df_val, alpha = alpha_val, scale = scale_val) else make_gen_skt(df = df_val, alpha = alpha_val, scale = scale_val)
    },
    "GED" = {
      nu_val <- params$ged_nu
      if (is.null(nu_val) || is.na(nu_val) || nu_val <= 0) nu_val <- 2
      sd_val <- params$ged_sd
      if (is.null(sd_val) || is.na(sd_val) || sd_val <= 0) sd_val <- 1
      if (isTRUE(use_fast)) make_gen_ged_fast(nu = nu_val, sd = sd_val) else make_gen_ged(nu = nu_val, sd = sd_val)
    },
    "Laplace" = {
      sc <- params$lap_scale
      if (is.null(sc) || is.na(sc) || sc <= 0) sc <- 1 / sqrt(2)
      if (isTRUE(use_fast)) make_gen_laplace_fast(scale = sc) else make_gen_laplace(scale = sc)
    },
    "Uniform" = {
      hw <- params$unif_hw
      if (is.null(hw) || is.na(hw) || hw <= 0) hw <- sqrt(3)
      if (isTRUE(use_fast)) make_gen_unif_fast(half_width = hw) else make_gen_unif(half_width = hw)
    },
    "Mixture Normal" = {
      sd1 <- params$mix_sd1
      sd2 <- params$mix_sd2
      p1 <- params$mix_prob1
      if (is.null(sd1) || is.na(sd1) || sd1 <= 0) sd1 <- 1
      if (is.null(sd2) || is.na(sd2) || sd2 <= 0) sd2 <- 3
      if (is.null(p1) || is.na(p1) || p1 <= 0 || p1 >= 1) p1 <- 0.9
      if (isTRUE(use_fast)) make_gen_mixnorm_fast(sd1 = sd1, sd2 = sd2, prob1 = p1) else make_gen_mixnorm(sd1 = sd1, sd2 = sd2, prob1 = p1)
    },
    "GARCH" = {
      omega <- params$garch_omega
      if (is.null(omega) || is.na(omega) || omega <= 0) omega <- 0.1
      alpha <- params$garch_alpha
      if (is.character(alpha)) {
        alpha <- tryCatch(
          as.numeric(trimws(strsplit(alpha, ",")[[1]])),
          warning = function(w) NULL, error = function(e) NULL)
      }
      if (is.null(alpha) || any(is.na(alpha)) || length(alpha) == 0)
        alpha <- c(0.2, 0.175, 0.15, 0.125, 0.1, 0.075, 0.05, 0.025)
      beta <- params$garch_beta
      if (is.character(beta)) {
        beta_str <- trimws(beta)
        if (!nzchar(beta_str) || beta_str == "") {
          beta <- NULL
        } else {
          beta <- tryCatch(as.numeric(trimws(strsplit(beta_str, ",")[[1]])),
                           warning = function(w) NULL, error = function(e) NULL)
          if (!is.null(beta) && (any(is.na(beta)) || length(beta) == 0)) beta <- NULL
        }
      }
      if (isTRUE(use_fast)) {
        if (is.null(beta)) beta <- numeric(0)
        make_gen_garch_fast(omega = omega, alpha = alpha, beta = beta)
      } else {
        make_gen_garch(omega = omega, alpha = alpha, beta = beta)
      }
    },
    "Heteroscedastic" = {
      # Support legacy w interface (vector or function)
      if (!is.null(params$hetero_w)) {
        if (isTRUE(use_fast)) {
          make_gen_hetero_fast(w = params$hetero_w, sd = params$hetero_sd %||% 1)
        } else {
          make_gen_hetero(w = params$hetero_w, sd = params$hetero_sd %||% 1)
        }
      } else {
        # Shape interface
        sd_val <- params$hetero_sd
        if (is.null(sd_val) || is.na(sd_val) || sd_val <= 0) sd_val <- 1
        shape <- params$hetero_shape %||% "linear"
        switch(shape,
          "linear" =, "sqrt" =, "log" =, "exp" = {
            from_val <- params$hetero_from %||% 1
            to_val <- params$hetero_to %||% 10
            if (isTRUE(use_fast)) {
              make_gen_hetero_fast(shape = shape, from = from_val, to = to_val, sd = sd_val)
            } else {
              make_gen_hetero(shape = shape, from = from_val, to = to_val, sd = sd_val)
            }
          },
          "power" = {
            from_val <- params$hetero_from %||% 1
            to_val <- params$hetero_to %||% 10
            p_val <- params$hetero_power %||% 2
            if (isTRUE(use_fast)) {
              make_gen_hetero_fast(shape = "power", from = from_val, to = to_val,
                                   power = p_val, sd = sd_val)
            } else {
              make_gen_hetero(shape = "power", from = from_val, to = to_val,
                              power = p_val, sd = sd_val)
            }
          },
          "step" = {
            breaks_str <- params$hetero_breaks %||% "0.5"
            levels_str <- params$hetero_levels %||% "1, 5"
            if (is.character(breaks_str)) {
              brk <- as.numeric(trimws(strsplit(breaks_str, ",")[[1]]))
            } else brk <- breaks_str
            if (is.character(levels_str)) {
              lvl <- as.numeric(trimws(strsplit(levels_str, ",")[[1]]))
            } else lvl <- levels_str
            if (any(is.na(brk))) brk <- 0.5
            if (any(is.na(lvl))) lvl <- c(1, 5)
            if (isTRUE(use_fast)) {
              make_gen_hetero_fast(shape = "step", breaks = brk, levels = lvl, sd = sd_val)
            } else {
              make_gen_hetero(shape = "step", breaks = brk, levels = lvl, sd = sd_val)
            }
          },
          "periodic" = {
            bw <- params$hetero_base_w %||% 1
            amp <- params$hetero_amplitude %||% 0.5
            per <- params$hetero_period %||% 12
            if (isTRUE(use_fast)) {
              make_gen_hetero_fast(shape = "periodic", base_w = bw, amplitude = amp,
                                   period = per, sd = sd_val)
            } else {
              make_gen_hetero(shape = "periodic", base_w = bw, amplitude = amp,
                              period = per, sd = sd_val)
            }
          }
        )
      }
    },
    stop("Unrecognized innovation distribution label: ", label)
  )
}

# Extract params list from adhoc module input (converts sim_* prefixed keys).
.adhoc_params_from_input <- function(input) {
  list(
    norm_sd       = input$sim_norm_sd,
    t_df          = input$sim_t_df,
    t_scale       = input$sim_t_scale,
    skt_df        = input$sim_skt_df,
    skt_alpha     = input$sim_skt_alpha,
    skt_scale     = input$sim_skt_scale,
    ged_nu        = input$sim_ged_nu,
    ged_sd        = input$sim_ged_sd,
    lap_scale     = input$sim_lap_scale,
    unif_hw       = input$sim_unif_hw,
    mix_sd1       = input$sim_mix_sd1,
    mix_sd2       = input$sim_mix_sd2,
    mix_prob1     = input$sim_mix_prob1,
    garch_omega   = input$sim_garch_omega,
    garch_alpha   = input$sim_garch_alpha,
    garch_beta    = input$sim_garch_beta,
    hetero_shape  = input$sim_hetero_shape,
    hetero_from   = input$sim_hetero_from,
    hetero_to     = input$sim_hetero_to,
    hetero_sd     = input$sim_hetero_sd,
    hetero_power  = input$sim_hetero_power,
    hetero_breaks = input$sim_hetero_breaks,
    hetero_levels = input$sim_hetero_levels,
    hetero_base_w = input$sim_hetero_base_w,
    hetero_amplitude = input$sim_hetero_amplitude,
    hetero_period = input$sim_hetero_period
  )
}

# Build grouped horizontal bar chart comparing CO / COB / COBA rejection rates.
# df must contain reject_asymp_05, reject_05, reject_adj_05 and label_col.
# Returns a ggplot object (or NULL if df is empty).
build_rejection_barchart <- function(df, label_col = "innov_dist", title = NULL) {
  if (nrow(df) == 0) return(NULL)

  bar_rows <- list()
  for (i in seq_len(nrow(df))) {
    lbl <- df[[label_col]][i]
    if ("reject_asymp_05" %in% names(df)) {
      bar_rows <- c(bar_rows, list(data.frame(
        label = lbl, method = "Asymptotic (CO)",
        rate = df$reject_asymp_05[i], stringsAsFactors = FALSE)))
    }
    if ("reject_05" %in% names(df)) {
      bar_rows <- c(bar_rows, list(data.frame(
        label = lbl, method = "Bootstrap (COB)",
        rate = df$reject_05[i], stringsAsFactors = FALSE)))
    }
    if ("reject_adj_05" %in% names(df)) {
      bar_rows <- c(bar_rows, list(data.frame(
        label = lbl, method = "COBA",
        rate = df$reject_adj_05[i], stringsAsFactors = FALSE)))
    }
  }
  bar_df <- do.call(rbind, bar_rows)
  bar_df$method <- factor(bar_df$method,
                           levels = c("Asymptotic (CO)", "Bootstrap (COB)", "COBA"))

  ggplot(bar_df, aes(x = rate, y = label, fill = method)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
    geom_vline(xintercept = 0.05, linetype = "dashed", linewidth = 0.8) +
    scale_fill_manual(values = c("Asymptotic (CO)" = "#e41a1c",
                                  "Bootstrap (COB)" = "#377eb8",
                                  "COBA" = "#4daf4a")) +
    labs(x = "Rejection Rate", y = NULL, fill = "Method",
         title = title %||% "Rejection Rates") +
    viewer_plot_theme(base_size = 13) +
    theme(legend.position = "bottom")
}

# =============================================================================
# Extracted Plot-Building Functions (from standalone modules)
# =============================================================================
# These are pure rendering functions: take data.frames/vectors, return ggplot
# objects (or base plots). No reactive dependencies.

# Standard rate column label map used by several plot functions
.rate_labels <- c(
  reject_05       = "Bootstrap (COB)",
  reject_asymp_05 = "Asymptotic (CO)",
  reject_adj_05   = "COBA"
)

# --- From mod_plots.R ---------------------------------------------------------

# Power curve: rejection rate vs phi, colored by innov_dist
plot_power_curve <- function(df, rate_col = "reject_05", facet_by = "none",
                             base_size = 14) {
  if (nrow(df) == 0) return(NULL)
  df$phi <- as.numeric(df$phi)

  p <- ggplot2::ggplot(df, ggplot2::aes(
    x = phi, y = .data[[rate_col]],
    color = factor(innov_dist), group = factor(innov_dist)
  )) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_hline(yintercept = 0.05, linetype = "dashed") +
    ggplot2::labs(
      x = expression(phi), y = "Rejection Rate",
      color = "Innovation", title = .rate_labels[rate_col]
    ) +
    viewer_plot_theme(base_size = base_size) +
    ggplot2::theme(legend.position = "bottom")

  se_col <- paste0(rate_col, "_se")
  if (se_col %in% names(df)) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data[[rate_col]] - .data[[se_col]],
                   ymax = .data[[rate_col]] + .data[[se_col]]),
      width = 0.02, alpha = 0.5)
  }

  if (facet_by == "n") {
    p <- p + ggplot2::facet_wrap(~ n, scales = "free_x",
                                  labeller = ggplot2::labeller(
                                    n = function(x) paste0("n = ", x)))
  } else if (facet_by == "innov_dist") {
    p <- p + ggplot2::facet_wrap(~ innov_dist)
  }
  p
}

# Heatmap: n vs phi grid, fill = rejection rate, faceted by innov_dist
plot_heatmap_rejection <- function(df, rate_col = "reject_05", base_size = 14,
                                   session = shiny::getDefaultReactiveDomain()) {
  if (nrow(df) == 0) return(NULL)
  df$n <- factor(df$n)
  df$phi <- factor(df$phi)
  df$label <- sprintf("%.3f", df[[rate_col]])
  df$label_col <- ifelse(
    is.na(df[[rate_col]]), "#212529",
    ifelse(df[[rate_col]] <= 0.02 | df[[rate_col]] >= 0.08, "#f8f9fa", "#212529"))

  fg <- viewer_plot_fg(session)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = phi, y = n, fill = .data[[rate_col]])) +
    ggplot2::geom_tile(
      color = grDevices::adjustcolor(fg, alpha.f = 0.20), linewidth = 0.25) +
    ggplot2::scale_fill_gradient2(
      low = "#2b8cbe", mid = "#fee8c8", high = "#d7301f",
      midpoint = 0.05, name = "Rejection\nRate") +
    ggplot2::geom_text(ggplot2::aes(label = label, colour = label_col),
                       size = 4.2, fontface = "bold") +
    ggplot2::scale_colour_identity() +
    ggplot2::labs(x = expression(phi), y = "Sample Size (n)",
                  title = .rate_labels[rate_col]) +
    viewer_plot_theme(base_size = base_size) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   legend.position = "right") +
    ggplot2::facet_wrap(~ innov_dist)
  p
}

# Deviation from nominal: dot plot of (rate - 0.05) for CO/COB/COBA
plot_deviation_from_nominal <- function(df, base_size = 13) {
  if (nrow(df) == 0) return(NULL)
  dev_rows <- list()
  for (i in seq_len(nrow(df))) {
    config_label <- sprintf("%s / n=%s / phi=%s", df$innov_dist[i], df$n[i], df$phi[i])
    if ("reject_asymp_05" %in% names(df)) {
      dev_rows <- c(dev_rows, list(data.frame(
        config = config_label, method = "CO",
        deviation = df$reject_asymp_05[i] - 0.05,
        se = if ("reject_asymp_05_se" %in% names(df)) df$reject_asymp_05_se[i] else NA_real_,
        stringsAsFactors = FALSE)))
    }
    if ("reject_05" %in% names(df)) {
      dev_rows <- c(dev_rows, list(data.frame(
        config = config_label, method = "COB",
        deviation = df$reject_05[i] - 0.05,
        se = if ("reject_05_se" %in% names(df)) df$reject_05_se[i] else NA_real_,
        stringsAsFactors = FALSE)))
    }
    if ("reject_adj_05" %in% names(df)) {
      dev_rows <- c(dev_rows, list(data.frame(
        config = config_label, method = "COBA",
        deviation = df$reject_adj_05[i] - 0.05,
        se = if ("reject_adj_05_se" %in% names(df)) df$reject_adj_05_se[i] else NA_real_,
        stringsAsFactors = FALSE)))
    }
  }
  dev_df <- do.call(rbind, dev_rows)
  dev_df$method <- factor(dev_df$method, levels = c("CO", "COB", "COBA"))

  p <- ggplot2::ggplot(dev_df, ggplot2::aes(x = deviation, y = config, color = method)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = c(-0.02, 0.02), linetype = "dotted") +
    ggplot2::geom_point(size = 2.5, position = ggplot2::position_dodge(width = 0.6))

  if (!all(is.na(dev_df$se))) {
    p <- p + ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = deviation - 1.96 * se, xmax = deviation + 1.96 * se),
      height = 0.2, position = ggplot2::position_dodge(width = 0.6), alpha = 0.5)
  }

  p + ggplot2::scale_color_manual(
      values = c(CO = "#e41a1c", COB = "#377eb8", COBA = "#4daf4a")) +
    ggplot2::labs(x = "Deviation from Nominal (Rate - 0.05)", y = NULL,
                  color = "Method", title = "Deviation from Nominal 5% Rate") +
    viewer_plot_theme(base_size = base_size) +
    ggplot2::theme(legend.position = "bottom")
}

# Rate vs sample size: line plot of rate vs n
plot_rate_vs_n <- function(df, rate_col = "reject_05", facet_by = "none",
                           base_size = 14) {
  if (nrow(df) == 0) return(NULL)
  df$n <- as.numeric(df$n)

  p <- ggplot2::ggplot(df, ggplot2::aes(
    x = n, y = .data[[rate_col]],
    color = factor(innov_dist), group = factor(innov_dist)
  )) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_hline(yintercept = 0.05, linetype = "dashed") +
    ggplot2::labs(
      x = "Sample Size (n)", y = "Rejection Rate",
      color = "Innovation",
      title = paste(.rate_labels[rate_col], "vs Sample Size")) +
    viewer_plot_theme(base_size = base_size) +
    ggplot2::theme(legend.position = "bottom")

  se_col <- paste0(rate_col, "_se")
  if (se_col %in% names(df)) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data[[rate_col]] - .data[[se_col]],
                   ymax = .data[[rate_col]] + .data[[se_col]]),
      width = 10, alpha = 0.5)
  }

  if (facet_by == "innov_dist") {
    p <- p + ggplot2::facet_wrap(~ innov_dist)
  } else if (facet_by == "phi") {
    p <- p + ggplot2::facet_wrap(~ phi,
                                  labeller = ggplot2::labeller(
                                    phi = function(x) paste0("phi = ", x)))
  }
  p
}

# Method forest plot: CO/COB/COBA rates with 95% CI per config
plot_forest_methods <- function(df, base_size = 13) {
  if (nrow(df) == 0) return(NULL)
  forest_rows <- list()
  for (i in seq_len(nrow(df))) {
    config_label <- sprintf("%s / n=%s / phi=%s", df$innov_dist[i], df$n[i], df$phi[i])
    if ("reject_asymp_05" %in% names(df)) {
      forest_rows <- c(forest_rows, list(data.frame(
        config = config_label, method = "CO", rate = df$reject_asymp_05[i],
        se = if ("reject_asymp_05_se" %in% names(df)) df$reject_asymp_05_se[i] else NA_real_,
        stringsAsFactors = FALSE)))
    }
    if ("reject_05" %in% names(df)) {
      forest_rows <- c(forest_rows, list(data.frame(
        config = config_label, method = "COB", rate = df$reject_05[i],
        se = if ("reject_05_se" %in% names(df)) df$reject_05_se[i] else NA_real_,
        stringsAsFactors = FALSE)))
    }
    if ("reject_adj_05" %in% names(df)) {
      forest_rows <- c(forest_rows, list(data.frame(
        config = config_label, method = "COBA", rate = df$reject_adj_05[i],
        se = if ("reject_adj_05_se" %in% names(df)) df$reject_adj_05_se[i] else NA_real_,
        stringsAsFactors = FALSE)))
    }
  }
  forest_df <- do.call(rbind, forest_rows)
  forest_df$method <- factor(forest_df$method, levels = c("CO", "COB", "COBA"))

  p <- ggplot2::ggplot(forest_df, ggplot2::aes(x = rate, y = config, color = method)) +
    ggplot2::geom_vline(xintercept = 0.05, linetype = "dashed") +
    ggplot2::geom_point(size = 2.5, position = ggplot2::position_dodge(width = 0.6))

  if (!all(is.na(forest_df$se))) {
    p <- p + ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = rate - 1.96 * se, xmax = rate + 1.96 * se),
      height = 0.2, position = ggplot2::position_dodge(width = 0.6), alpha = 0.5)
  }

  p + ggplot2::scale_color_manual(
      values = c(CO = "#e41a1c", COB = "#377eb8", COBA = "#4daf4a")) +
    ggplot2::labs(x = "Rejection Rate", y = NULL, color = "Method",
                  title = "Method Comparison: Rejection Rate with 95% CI") +
    viewer_plot_theme(base_size = base_size) +
    ggplot2::theme(legend.position = "bottom")
}

# --- From mod_pvalue.R --------------------------------------------------------

# P-value histogram: side-by-side histograms for each available method
plot_pvalue_histogram <- function(df, base_size = 13,
                                  session = shiny::getDefaultReactiveDomain()) {
  methods <- c(pvalue = "COB (Bootstrap)", pvalue_asymp = "CO (Asymptotic)",
               pvalue_adj = "COBA (Adjusted)")
  long_rows <- list()
  for (col in names(methods)) {
    if (col %in% names(df) && !all(is.na(df[[col]]))) {
      vals <- df[[col]][!is.na(df[[col]])]
      long_rows <- c(long_rows, list(data.frame(
        pval = vals, method = methods[col], stringsAsFactors = FALSE)))
    }
  }
  if (length(long_rows) == 0) return(NULL)
  pval_long <- do.call(rbind, long_rows)

  fg <- viewer_plot_fg(session)
  p <- ggplot2::ggplot(pval_long, ggplot2::aes(x = pval)) +
    ggplot2::geom_histogram(bins = 20, fill = "#4292c6", color = NA, alpha = 0.8,
                            boundary = 0) +
    ggplot2::facet_wrap(~ method, scales = "free_y") +
    ggplot2::labs(x = "P-value", y = "Count",
                  title = "P-value Distributions (uniform = well-calibrated)") +
    viewer_plot_theme(base_size = base_size) +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold", colour = fg))

  n_per_method <- table(pval_long$method)
  ref_df <- data.frame(method = names(n_per_method),
                        yint = as.numeric(n_per_method) / 20,
                        stringsAsFactors = FALSE)
  p + ggplot2::geom_hline(data = ref_df, ggplot2::aes(yintercept = yint),
                           color = "red", linetype = "dashed", linewidth = 0.6)
}

# QQ plot: p-values vs Uniform(0,1) quantiles
plot_pvalue_qq <- function(df, base_size = 13,
                           session = shiny::getDefaultReactiveDomain()) {
  methods <- c(pvalue = "COB (Bootstrap)", pvalue_asymp = "CO (Asymptotic)",
               pvalue_adj = "COBA (Adjusted)")
  qq_rows <- list()
  for (col in names(methods)) {
    if (col %in% names(df)) {
      vals <- df[[col]][!is.na(df[[col]])]
      if (length(vals) > 0) {
        n_vals <- length(vals)
        theoretical <- (seq_len(n_vals) - 0.5) / n_vals
        qq_rows <- c(qq_rows, list(data.frame(
          theoretical = theoretical, observed = sort(vals),
          method = methods[col], stringsAsFactors = FALSE)))
      }
    }
  }
  if (length(qq_rows) == 0) return(NULL)
  qq_df <- do.call(rbind, qq_rows)
  fg <- viewer_plot_fg(session)

  ggplot2::ggplot(qq_df, ggplot2::aes(x = theoretical, y = observed)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    ggplot2::geom_point(size = 0.8, alpha = 0.5, color = "#2171b5") +
    ggplot2::facet_wrap(~ method) +
    ggplot2::labs(x = "Theoretical Quantiles (Uniform)",
                  y = "Observed P-value Quantiles",
                  title = "QQ Plot: P-values vs Uniform(0, 1)") +
    ggplot2::coord_equal() +
    viewer_plot_theme(base_size = base_size) +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold", colour = fg))
}

# P-value scatter: COB vs CO, COB vs COBA
plot_pvalue_scatter <- function(df, base_size = 13,
                                session = shiny::getDefaultReactiveDomain()) {
  scatter_rows <- list()
  if ("pvalue_asymp" %in% names(df)) {
    idx <- !is.na(df$pvalue) & !is.na(df$pvalue_asymp)
    if (any(idx)) {
      scatter_rows <- c(scatter_rows, list(data.frame(
        boot_pval = df$pvalue[idx], other_pval = df$pvalue_asymp[idx],
        comparison = "COB vs CO (Asymptotic)", stringsAsFactors = FALSE)))
    }
  }
  if ("pvalue_adj" %in% names(df)) {
    idx <- !is.na(df$pvalue) & !is.na(df$pvalue_adj)
    if (any(idx)) {
      scatter_rows <- c(scatter_rows, list(data.frame(
        boot_pval = df$pvalue[idx], other_pval = df$pvalue_adj[idx],
        comparison = "COB vs COBA (Adjusted)", stringsAsFactors = FALSE)))
    }
  }
  if (length(scatter_rows) == 0) return(NULL)
  scatter_df <- do.call(rbind, scatter_rows)
  fg <- viewer_plot_fg(session)

  ggplot2::ggplot(scatter_df, ggplot2::aes(x = boot_pval, y = other_pval)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    ggplot2::geom_point(size = 0.6, alpha = 0.2, color = "#2171b5") +
    ggplot2::facet_wrap(~ comparison) +
    ggplot2::labs(x = "Bootstrap P-value (COB)", y = "Comparison Method P-value",
                  title = "P-value Agreement Between Methods") +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    viewer_plot_theme(base_size = base_size) +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold", colour = fg))
}

# --- From mod_bootstrap_dist.R ------------------------------------------------

# Bootstrap distribution histogram with observed t-stat overlay
plot_bootstrap_distribution <- function(boot_dist, obs_stat, pval,
                                        title = "Bootstrap Distribution",
                                        base_size = 14) {
  hist_data <- data.frame(t_stat = boot_dist)
  ggplot2::ggplot(hist_data, ggplot2::aes(x = t_stat)) +
    ggplot2::geom_histogram(bins = 40, fill = "#4292c6", color = NA, alpha = 0.8) +
    ggplot2::geom_vline(xintercept = obs_stat, color = "red",
                        linewidth = 1.2, linetype = "solid") +
    ggplot2::annotate("text", x = obs_stat, y = Inf, vjust = 2, hjust = -0.1,
                      label = sprintf("T_obs = %.3f\np = %.4f", obs_stat, pval),
                      color = "red", size = 4, fontface = "bold") +
    ggplot2::labs(x = "Bootstrap t-statistic", y = "Count", title = title) +
    viewer_plot_theme(base_size = base_size)
}

# --- From mod_diagnostics.R ---------------------------------------------------

# AR order distribution bar chart, optionally faceted by config
plot_ar_order_distribution <- function(df, base_size = 14) {
  if (nrow(df) == 0) return(NULL)
  df$null_ar_order <- factor(df$null_ar_order)
  df$config <- sprintf("n=%s, phi=%s, %s", df$n, df$phi, df$innov_dist)
  n_configs <- length(unique(df$config))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = null_ar_order)) +
    ggplot2::geom_bar(fill = "#4292c6", color = NA, alpha = 0.8) +
    ggplot2::labs(x = "Fitted AR Order", y = "Count",
                  title = "Distribution of Fitted AR Orders Under the Null") +
    viewer_plot_theme(base_size = base_size)

  if (n_configs > 1 && n_configs <= 12) {
    p <- p + ggplot2::facet_wrap(~ config, scales = "free_y")
  } else if (n_configs > 12) {
    p <- p + ggplot2::facet_wrap(~ config, scales = "free_y", ncol = 4)
  }
  p
}

# Estimated null AR(1) coefficient distribution
plot_ar_coefficient_distribution <- function(df, base_size = 14) {
  if (nrow(df) == 0 || all(is.na(df$null_phi1))) return(NULL)
  if (!"null_ar_order" %in% names(df)) return(NULL)
  df <- df[!is.na(df$null_phi1) & df$null_ar_order == 1, ]
  if (nrow(df) == 0) return(NULL)
  df$true_phi_label <- paste0("True phi = ", df$phi)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = null_phi1)) +
    ggplot2::geom_histogram(bins = 50, fill = "#4292c6", color = NA, alpha = 0.8) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = phi), color = "red",
                        linewidth = 1, linetype = "dashed") +
    ggplot2::labs(x = "Estimated First AR Coefficient", y = "Count",
                  title = "Distribution of Estimated Null AR(1) Coefficient") +
    viewer_plot_theme(base_size = base_size)

  n_configs <- length(unique(df$true_phi_label))
  if (n_configs > 1) p <- p + ggplot2::facet_wrap(~ true_phi_label, scales = "free")
  p
}

# MC convergence: cumulative rejection rate vs simulation number
plot_mc_convergence <- function(df, pval_col = "pvalue", method_label = NULL,
                                base_size = 13,
                                session = shiny::getDefaultReactiveDomain()) {
  if (nrow(df) == 0) return(NULL)
  if (is.null(method_label)) method_label <- .rate_labels[pval_col] %||% pval_col
  df$config <- sprintf("%s / n=%s / phi=%s", df$innov_dist, df$n, df$phi)
  configs <- unique(df$config)

  conv_rows <- list()
  for (cfg in configs) {
    sub <- df[df$config == cfg, ]
    pvals <- sub[[pval_col]]
    rejected <- pvals < 0.05
    cum_rate <- cumsum(rejected) / seq_along(rejected)
    conv_rows <- c(conv_rows, list(data.frame(
      sim_num = seq_along(cum_rate), cum_reject_rate = cum_rate,
      config = cfg, stringsAsFactors = FALSE)))
  }
  conv_df <- do.call(rbind, conv_rows)
  fg <- viewer_plot_fg(session)

  ggplot2::ggplot(conv_df, ggplot2::aes(
    x = sim_num, y = cum_reject_rate, color = config)) +
    ggplot2::geom_line(alpha = 0.7, linewidth = 0.6) +
    ggplot2::geom_hline(yintercept = 0.05, linetype = "dashed") +
    ggplot2::labs(x = "Simulation Number", y = "Cumulative Rejection Rate",
                  color = "Configuration",
                  title = paste0("MC Convergence: ", method_label)) +
    viewer_plot_theme(base_size = base_size) +
    ggplot2::theme(legend.position = "bottom",
                   legend.text = ggplot2::element_text(size = 8, colour = fg))
}

# Batch consistency: boxplot + jitter of per-batch rejection rates
plot_batch_consistency <- function(df, rate_col = "reject_05", base_size = 13,
                                   session = shiny::getDefaultReactiveDomain()) {
  if (nrow(df) == 0 || !rate_col %in% names(df)) return(NULL)
  rate_label <- .rate_labels[rate_col] %||% rate_col
  df$config <- sprintf("%s / n=%s / phi=%s", df$innov_dist, df$n, df$phi)
  fg <- viewer_plot_fg(session)

  ggplot2::ggplot(df, ggplot2::aes(x = config, y = .data[[rate_col]])) +
    ggplot2::geom_hline(yintercept = 0.05, linetype = "dashed") +
    ggplot2::geom_boxplot(fill = "#5fad8c", alpha = 0.6, outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "#2171b5") +
    ggplot2::labs(x = NULL, y = "Rejection Rate",
                  title = paste0("Batch-to-Batch Consistency: ", rate_label)) +
    viewer_plot_theme(base_size = base_size) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(
      angle = 45, hjust = 1, size = 9, colour = fg))
}

# Test statistic distribution histogram, optionally faceted by config
plot_test_statistic_distribution <- function(df, base_size = 14) {
  if (nrow(df) == 0) return(NULL)
  df$config <- sprintf("%s / n=%s / phi=%s", df$innov_dist, df$n, df$phi)
  n_configs <- length(unique(df$config))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = obs_stat)) +
    ggplot2::geom_histogram(bins = 50, fill = "#4292c6", color = NA, alpha = 0.8) +
    ggplot2::labs(x = "Observed Test Statistic", y = "Count",
                  title = "Distribution of Observed Test Statistics") +
    viewer_plot_theme(base_size = base_size)

  if (n_configs > 1 && n_configs <= 12) {
    p <- p + ggplot2::facet_wrap(~ config, scales = "free")
  } else if (n_configs > 12) {
    p <- p + ggplot2::facet_wrap(~ config, scales = "free", ncol = 4)
  }
  p
}

# --- From mod_parallel_coords.R -----------------------------------------------

# Generic parallel coordinates plot builder
build_parcoord <- function(df, cols, color_col, color_label, title = "",
                           rate_cols = c("reject_asymp_05", "reject_05",
                                         "reject_adj_05"),
                           base_size = 13) {
  if (nrow(df) == 0) return(NULL)

  col_names <- names(cols)
  axis_labels <- unname(cols)
  n_axes <- length(col_names)
  df$config_id <- seq_len(nrow(df))

  # Determine shared rate scale
  rate_in_cols <- intersect(col_names, rate_cols)
  if (length(rate_in_cols) > 0) {
    all_rates <- unlist(lapply(rate_in_cols, function(c) as.numeric(df[[c]])))
    rate_max <- max(max(all_rates, na.rm = TRUE) * 1.1, 0.5)
  } else {
    rate_max <- 1
  }

  # Scale columns to [0,1]
  ranges <- list()
  scaled <- data.frame(config_id = df$config_id)
  for (col in col_names) {
    vals <- as.numeric(df[[col]])
    if (col %in% rate_cols) {
      ranges[[col]] <- c(0, rate_max)
      scaled[[col]] <- vals / rate_max
    } else {
      rng <- range(vals, na.rm = TRUE)
      ranges[[col]] <- rng
      scaled[[col]] <- if (rng[1] == rng[2]) 0.5 else (vals - rng[1]) / (rng[2] - rng[1])
    }
  }

  # Build long-format
  long_rows <- vector("list", nrow(df) * n_axes)
  idx <- 1
  for (i in seq_len(nrow(df))) {
    for (j in seq_len(n_axes)) {
      long_rows[[idx]] <- data.frame(
        config_id = i, axis_num = j, axis_label = axis_labels[j],
        value = scaled[[col_names[j]]][i],
        color_var = as.character(df[[color_col]][i]),
        stringsAsFactors = FALSE)
      idx <- idx + 1
    }
  }
  long_df <- do.call(rbind, long_rows)
  long_df$axis_num <- as.integer(long_df$axis_num)

  # Tick labels (min/max at each axis)
  tick_labels <- unlist(lapply(col_names, function(col) {
    rng <- ranges[[col]]
    c(format(rng[1], digits = 3), format(rng[2], digits = 3))
  }))
  tick_df <- data.frame(
    axis_num = rep(seq_len(n_axes), each = 2),
    y = rep(c(0, 1), times = n_axes),
    label = tick_labels, stringsAsFactors = FALSE)

  # 0.05 reference for rate axes
  ref_lines <- data.frame()
  for (j in seq_along(col_names)) {
    if (col_names[j] %in% rate_cols) {
      ref_lines <- rbind(ref_lines, data.frame(
        axis_num = j, y = 0.05 / rate_max, stringsAsFactors = FALSE))
    }
  }

  p <- ggplot2::ggplot(long_df, ggplot2::aes(
    x = axis_num, y = value, group = config_id, color = color_var)) +
    ggplot2::geom_line(alpha = 0.6, linewidth = 0.7) +
    ggplot2::geom_point(size = 1.5, alpha = 0.7) +
    ggplot2::scale_x_continuous(breaks = seq_len(n_axes), labels = axis_labels,
                                limits = c(0.5, n_axes + 0.5)) +
    ggplot2::scale_y_continuous(limits = c(-0.05, 1.05),
                                breaks = seq(0, 1, by = 0.25)) +
    ggplot2::geom_text(data = tick_df,
                       ggplot2::aes(x = axis_num, y = y, label = label),
                       inherit.aes = FALSE, size = 2.8,
                       hjust = 1.2, fontface = "italic") +
    ggplot2::geom_vline(xintercept = seq_len(n_axes), linewidth = 0.3) +
    ggplot2::labs(x = NULL, y = NULL, color = color_label, title = title) +
    viewer_plot_theme(base_size = base_size) +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "bottom")

  if (nrow(ref_lines) > 0) {
    p <- p + ggplot2::geom_point(data = ref_lines,
                                  ggplot2::aes(x = axis_num, y = y),
                                  inherit.aes = FALSE, shape = 4, size = 3,
                                  stroke = 1.2)
  }
  p
}

# DGP parameter scatter: rate vs phi, colored by n, shaped by innov_dist
plot_pc_dgp_scatter <- function(df, base_size = 13) {
  if (nrow(df) == 0) return(NULL)
  rate_map <- c(reject_asymp_05 = "CO", reject_05 = "COB", reject_adj_05 = "COBA")
  long_rows <- list()
  for (i in seq_len(nrow(df))) {
    for (rate_col in names(rate_map)) {
      if (rate_col %in% names(df)) {
        long_rows <- c(long_rows, list(data.frame(
          n = df$n[i], phi = df$phi[i], innov_dist = df$innov_dist[i],
          method = rate_map[rate_col], rate = df[[rate_col]][i],
          stringsAsFactors = FALSE)))
      }
    }
  }
  long_df <- do.call(rbind, long_rows)
  long_df$method <- factor(long_df$method, levels = c("CO", "COB", "COBA"))

  ggplot2::ggplot(long_df, ggplot2::aes(
    x = phi, y = rate, color = factor(n), shape = innov_dist)) +
    ggplot2::geom_point(size = 2.5, alpha = 0.7) +
    ggplot2::geom_hline(yintercept = 0.05, linetype = "dashed") +
    ggplot2::facet_wrap(~ method) +
    ggplot2::labs(x = expression(phi), y = "Rejection Rate",
                  color = "n", shape = "Innovation",
                  title = "DGP Parameter Relationships: Rate vs phi by n") +
    viewer_plot_theme(base_size = base_size) +
    ggplot2::theme(legend.position = "bottom")
}

# --- From mod_adhoc_sim.R -----------------------------------------------------

# Theoretical density overlay for innovation distributions
innov_theoretical_density <- function(dist, params) {
  switch(dist,
    "Normal" = {
      sd_val <- params$sd %||% 1
      list(dfun = function(x) stats::dnorm(x, 0, sd_val),
           label = sprintf("N(0, %.2f)", sd_val), col = "red")
    },
    "Student's t" = {
      df_val <- params$df %||% 3
      scaled <- isTRUE(params$scale)
      if (scaled && df_val > 2) {
        c_sc <- sqrt(df_val / (df_val - 2))
        list(dfun = function(x) stats::dt(x * c_sc, df = df_val) * c_sc,
             label = sprintf("t(%d) scaled", df_val), col = "red")
      } else {
        list(dfun = function(x) stats::dt(x, df = df_val),
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
        b_nu <- sqrt(df_val / pi) * exp(lgamma((df_val - 1) / 2) -
                                          lgamma(df_val / 2))
        theo_mean <- b_nu * delta
      }
      if (scaled && df_val > 2) {
        cumulants <- tryCatch(
          sn::st.cumulants(xi = 0, omega = 1, alpha = alpha_val,
                           nu = df_val, n = 2),
          error = function(e) c(NA_real_, NA_real_))
        theo_var <- as.numeric(cumulants[2])
        if (is.finite(theo_var) && theo_var > 0) {
          theo_sd <- sqrt(theo_var)
          theo_mean <- as.numeric(cumulants[1])
        }
      }
      tm <- theo_mean; ts <- theo_sd; a <- alpha_val; d <- df_val
      list(dfun = function(x) sn::dst(x * ts + tm, xi = 0, omega = 1,
                                       alpha = a, nu = d) * ts,
           label = sprintf("Skew-t(df=%d, alpha=%.1f)", df_val, alpha_val),
           col = "red")
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
      list(dfun = function(x) stats::dunif(x, min = -hw, max = hw),
           label = sprintf("Unif(%.2f, %.2f)", -hw, hw), col = "red")
    },
    "Mixture Normal" = {
      s1 <- params$sd1 %||% 1; s2 <- params$sd2 %||% 3
      p1 <- params$prob1 %||% 0.9
      list(dfun = function(x) p1 * stats::dnorm(x, 0, s1) +
                               (1 - p1) * stats::dnorm(x, 0, s2),
           label = sprintf("MixN(%.0f%%*N(0,%.1f) + %.0f%%*N(0,%.1f))",
                           p1 * 100, s1, (1 - p1) * 100, s2), col = "red")
    },
    NULL
  )
}

# Innovation diagnostics: base graphics multi-panel
# For IID: 2 panels (histogram + QQ); for time-dependent: 3 panels
plot_innovation_diagnostics <- function(innov, is_time_dependent, dist_label,
                                        innov_dist = NULL, innov_params = list(),
                                        fg = "#212529") {
  if (is_time_dependent) {
    op <- graphics::par(mfrow = c(3, 1), mar = c(5, 5.5, 3, 1),
                        cex.axis = 1.3, cex.lab = 1.4,
                        col.axis = fg, col.lab = fg, col.main = fg, fg = fg)
    on.exit(graphics::par(op))

    # Time series
    plot(innov, type = "l", col = "darkorange", lwd = 1,
         xlab = "Time", ylab = "Innovation",
         main = sprintf("Innovation Time Series (%s)", dist_label))
    graphics::abline(h = 0, lty = 3, col = "grey60")

    # Histogram
    h <- graphics::hist(innov, breaks = 40, plot = FALSE)
    x_seq <- seq(min(innov), max(innov), length.out = 300)
    curve_vals <- stats::dnorm(x_seq, mean(innov), stats::sd(innov))
    y_max <- max(h$density, curve_vals) * 1.05
    graphics::hist(innov, breaks = h$breaks, col = "darkorange", border = "white",
                   main = sprintf("Innovation Histogram (%s)", dist_label),
                   xlab = "Value", ylab = "Density", probability = TRUE,
                   ylim = c(0, y_max))
    graphics::lines(x_seq, curve_vals, col = "steelblue", lwd = 2, lty = 2)
    graphics::legend("topright", legend = "Normal reference",
                     col = "steelblue", lwd = 2, lty = 2, bty = "n", cex = 0.9,
                     text.col = fg)

    # Squared innovations
    plot(innov^2, type = "l", col = "firebrick", lwd = 1,
         xlab = "Time", ylab = expression(epsilon^2),
         main = sprintf("Squared Innovations (%s) - Volatility Structure",
                        dist_label))
    graphics::abline(h = mean(innov^2), lty = 2, col = "grey40")
    graphics::legend("topright", legend = sprintf("Mean = %.3f", mean(innov^2)),
                     bty = "n", cex = 0.9, text.col = fg)
  } else {
    op <- graphics::par(mfrow = c(2, 1), mar = c(5, 5.5, 3, 1),
                        cex.axis = 1.3, cex.lab = 1.4,
                        col.axis = fg, col.lab = fg, col.main = fg, fg = fg)
    on.exit(graphics::par(op))

    h <- graphics::hist(innov, breaks = 40, plot = FALSE)
    x_seq <- seq(min(innov), max(innov), length.out = 300)
    theo <- innov_theoretical_density(innov_dist, innov_params)
    if (!is.null(theo)) {
      curve_vals <- theo$dfun(x_seq)
    } else {
      curve_vals <- stats::dnorm(x_seq, mean(innov), stats::sd(innov))
    }
    y_max <- max(h$density, curve_vals) * 1.05
    graphics::hist(innov, breaks = h$breaks, col = "steelblue", border = "white",
                   main = sprintf("Innovation Histogram (%s, n=%d)",
                                  dist_label, length(innov)),
                   xlab = "Value", ylab = "Density", probability = TRUE,
                   ylim = c(0, y_max))
    if (!is.null(theo)) {
      graphics::lines(x_seq, curve_vals, col = "red", lwd = 2)
      graphics::legend("topright",
                       legend = c(sprintf("Mean=%.3f, SD=%.3f",
                                          mean(innov), stats::sd(innov)),
                                  theo$label),
                       col = c(NA, "red"), lwd = c(NA, 2), lty = c(NA, 1),
                       bty = "n", cex = 0.9, text.col = fg)
    } else {
      graphics::lines(x_seq, curve_vals, col = "red", lwd = 2, lty = 2)
      graphics::legend("topright",
                       legend = c(sprintf("Mean=%.3f, SD=%.3f",
                                          mean(innov), stats::sd(innov)),
                                  "Normal reference"),
                       col = c(NA, "red"), lwd = c(NA, 2), lty = c(NA, 2),
                       bty = "n", cex = 0.9, text.col = fg)
    }

    # QQ plot
    stats::qqnorm(innov, main = sprintf("Normal QQ Plot (%s)", dist_label),
                  col = "steelblue", pch = 16, cex = 0.6)
    stats::qqline(innov, col = "red", lwd = 2)
  }
}

# Null model diagnostics: 3-panel base graphics (AR order, variance, rejection by order)
plot_null_model_diagnostics <- function(results, nsims, maxp, min_p = 1L,
                                        fg = "#212529",
                                        reject_mode = c("all", "cob", "co", "coba")) {
  reject_mode <- match.arg(reject_mode)
  ar_orders <- vapply(results, function(r) r$p, integer(1))
  varas <- vapply(results, function(r) r$vara, numeric(1))
  pvals_boot <- vapply(results, function(r) r$pvalue, numeric(1))
  pvals_asymp <- vapply(results, function(r) r$pvalue_asymp, numeric(1))
  pvals_adj <- vapply(results, function(r) {
    if (is.null(r$pvalue_adj) || length(r$pvalue_adj) == 0 || is.na(r$pvalue_adj)) {
      NA_real_
    } else {
      as.numeric(r$pvalue_adj)
    }
  }, numeric(1))

  graphics::par(
    mfrow = c(3, 1),
    mar = c(5.5, 5.3, 3.8, 1.8),
    mgp = c(2.6, 0.85, 0),
    col.axis = fg,
    col.lab = fg,
    col.main = fg,
    fg = fg,
    cex = 1.02,
    cex.axis = 1.0,
    cex.lab = 1.04,
    cex.main = 1.08
  )

  # Panel 1: AR Order Distribution
  order_tbl <- table(factor(ar_orders, levels = 0:maxp))
  base_cols <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f")
  order_cols <- grDevices::colorRampPalette(base_cols)(maxp + 1)
  graphics::barplot(order_tbl, col = order_cols, border = NA,
                    main = sprintf("Fitted Null AR Order  (n=%d sims, min_p=%d, maxp=%d)",
                                   nsims, min_p, maxp),
                    xlab = "AR Order", ylab = "Count")
  modal_order <- as.integer(names(which.max(order_tbl)))
  modal_pct <- max(order_tbl) / nsims * 100
  graphics::mtext(sprintf("Mode: AR(%d) = %.0f%%", modal_order, modal_pct),
                  side = 3, adj = 1, line = 0, cex = 1.0, col = fg)

  # Panel 2: Innovation Variance Distribution
  graphics::hist(varas, breaks = 30,
                 col = grDevices::adjustcolor("steelblue", alpha.f = 0.6),
                 border = "white",
                 main = "Estimated Innovation Variance (null_vara)",
                 xlab = expression(hat(sigma)[a]^2), ylab = "Frequency")
  graphics::abline(v = mean(varas), col = "#e41a1c", lwd = 2, lty = 1)
  graphics::abline(v = stats::median(varas), col = "#ff7f00", lwd = 2, lty = 2)
  graphics::legend("topright",
                   legend = c(sprintf("Mean = %.4f", mean(varas)),
                              sprintf("Median = %.4f", stats::median(varas)),
                              sprintf("SD = %.4f", stats::sd(varas))),
                   col = c("#e41a1c", "#ff7f00", NA),
                   lwd = c(2, 2, NA), lty = c(1, 2, NA),
                   bty = "n", cex = 1.1, text.col = fg)

  # Panel 3: Rejection Rate by AR Order
  unique_orders <- sort(unique(ar_orders))
  if (length(unique_orders) > 0) {
    boot_rates <- numeric(length(unique_orders))
    asymp_rates <- numeric(length(unique_orders))
    adj_rates <- numeric(length(unique_orders))
    counts <- integer(length(unique_orders))
    for (j in seq_along(unique_orders)) {
      idx <- ar_orders == unique_orders[j]
      counts[j] <- sum(idx)
      boot_rates[j] <- mean(pvals_boot[idx] < 0.05, na.rm = TRUE)
      asymp_rates[j] <- mean(pvals_asymp[idx] < 0.05, na.rm = TRUE)
      adj_rates[j] <- mean(pvals_adj[idx] < 0.05, na.rm = TRUE)
    }

    if (identical(reject_mode, "cob")) {
      rate_mat <- rbind(boot_rates)
      row_labs <- c("COB")
      fill_cols <- c("#377eb8")
    } else if (identical(reject_mode, "co")) {
      rate_mat <- rbind(asymp_rates)
      row_labs <- c("CO")
      fill_cols <- c("#e41a1c")
    } else if (identical(reject_mode, "coba")) {
      rate_mat <- rbind(adj_rates)
      row_labs <- c("COBA")
      fill_cols <- c("#4daf4a")
    } else {
      if (all(is.na(adj_rates))) {
        rate_mat <- rbind(boot_rates, asymp_rates)
        row_labs <- c("COB", "CO")
        fill_cols <- c("#377eb8", "#e41a1c")
      } else {
        rate_mat <- rbind(boot_rates, asymp_rates, adj_rates)
        row_labs <- c("COB", "CO", "COBA")
        fill_cols <- c("#377eb8", "#e41a1c", "#4daf4a")
      }
    }

    rownames(rate_mat) <- row_labs
    colnames(rate_mat) <- paste0("AR(", unique_orders, ")\nn=", counts)

    if (all(is.na(rate_mat))) {
      graphics::plot.new()
      graphics::text(0.5, 0.5,
                     "No rejection-rate values available for selected method(s).",
                     cex = 1.05, col = fg)
      return(invisible(NULL))
    }

    max_rate <- suppressWarnings(max(rate_mat, na.rm = TRUE))
    if (!is.finite(max_rate)) max_rate <- 0.15
    graphics::barplot(rate_mat, beside = TRUE,
                      col = fill_cols, border = NA,
                      width = 0.55,
                      main = "Rejection Rate by Fitted AR Order",
                      xlab = "AR Order (with sim count)", ylab = "Rejection Rate",
                      ylim = c(0, max(0.15, max_rate * 1.2)))
    graphics::abline(h = 0.05, lty = 2, col = "grey50", lwd = 1.5)

    legend_cols <- c(fill_cols, "grey50")
    legend_pch <- c(rep(15, length(fill_cols)), NA)
    legend_lwd <- c(rep(NA, length(fill_cols)), 1.5)
    legend_lty <- c(rep(NA, length(fill_cols)), 2)
    legend_labs <- c(row_labs, "Nominal 0.05")
    graphics::legend("topright",
                     legend = legend_labs,
                     col = legend_cols,
                     lwd = legend_lwd,
                     lty = legend_lty,
                     pch = legend_pch,
                     pt.cex = 2,
                     bty = "n", cex = 1.1, text.col = fg)
  }
}
