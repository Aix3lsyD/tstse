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

# Build parameterized SQL query with optional filters
build_query <- function(view_name, n_vals = NULL, phi_vals = NULL, innov_vals = NULL) {
  sql <- paste0("SELECT * FROM ", view_name, " WHERE 1=1")
  params <- list()

  if (length(n_vals) > 0) {
    placeholders <- paste(rep("?", length(n_vals)), collapse = ", ")
    sql <- paste0(sql, " AND n IN (", placeholders, ")")
    params <- c(params, as.list(as.integer(n_vals)))
  }
  if (length(phi_vals) > 0) {
    placeholders <- paste(rep("?", length(phi_vals)), collapse = ", ")
    sql <- paste0(sql, " AND phi IN (", placeholders, ")")
    params <- c(params, as.list(as.numeric(phi_vals)))
  }
  if (length(innov_vals) > 0) {
    placeholders <- paste(rep("?", length(innov_vals)), collapse = ", ")
    sql <- paste0(sql, " AND innov_dist IN (", placeholders, ")")
    params <- c(params, as.list(innov_vals))
  }

  sql <- paste0(sql, " ORDER BY innov_dist, n, phi")
  list(sql = sql, params = params)
}

# Query grid data with optional single-value filters
grid_query <- function(con, innov_dist = NULL, phi = NULL, n = NULL,
                       order_col = "n", batch_id = NULL) {
  if (!is.null(batch_id) && batch_id != "all") {
    sql <- "SELECT * FROM v_rejection_rates_by_batch WHERE 1=1"
  } else {
    sql <- "SELECT * FROM v_rejection_rates WHERE 1=1"
  }
  params <- list()

  if (!is.null(batch_id) && batch_id != "all") {
    sql <- paste0(sql, " AND batch_id = ?")
    params <- c(params, list(as.integer(batch_id)))
  }
  if (!is.null(innov_dist)) {
    sql <- paste0(sql, " AND innov_dist = ?")
    params <- c(params, list(innov_dist))
  }
  if (!is.null(phi)) {
    sql <- paste0(sql, " AND phi = ?")
    params <- c(params, list(as.numeric(phi)))
  }
  if (!is.null(n)) {
    sql <- paste0(sql, " AND n = ?")
    params <- c(params, list(as.integer(n)))
  }

  sql <- paste0(sql, " ORDER BY ", order_col)
  tryCatch(
    dbGetQuery(con, sql, params = params),
    error = function(e) data.frame()
  )
}

# Get batch choices matching a config filter
grid_batch_choices <- function(con, innov_dist = NULL, phi = NULL, n = NULL) {
  sql <- paste0(
    "SELECT DISTINCT b.batch_id, b.label ",
    "FROM simulations s JOIN batches b ON s.batch_id = b.batch_id WHERE 1=1"
  )
  params <- list()
  if (!is.null(innov_dist)) {
    sql <- paste0(sql, " AND s.innov_dist = ?")
    params <- c(params, list(innov_dist))
  }
  if (!is.null(phi)) {
    sql <- paste0(sql, " AND s.phi = ?")
    params <- c(params, list(as.numeric(phi)))
  }
  if (!is.null(n)) {
    sql <- paste0(sql, " AND s.n = ?")
    params <- c(params, list(as.integer(n)))
  }
  sql <- paste0(sql, " ORDER BY b.batch_id")
  batches <- tryCatch(dbGetQuery(con, sql, params = params),
                      error = function(e) data.frame())
  if (nrow(batches) == 0) return(c("All (Pooled)" = "all"))
  choices <- setNames(
    as.character(batches$batch_id),
    ifelse(is.na(batches$label) | batches$label == "",
           paste("Batch", batches$batch_id),
           paste0("Batch ", batches$batch_id, ": ", batches$label))
  )
  c("All (Pooled)" = "all", choices)
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
