# Silence R CMD check notes for data frame column references in gt/tidyverse code
utils::globalVariables(c(
  "p_val", "Estimate", "Std_Error", "t_value", "Parameter", "p_display", "Sig",
  "ARCH", "GARCH", "WLB1", "WLB2", "WLB3", "Nyblom", "Nyblom_crit", "Nyblom_pass",
  "SignBias", "n_sig", "n_coef", "CoefSig_pass", "AIC", "AICC", "BIC",
  "WLB1_fmt", "WLB2_fmt", "WLB3_fmt", "Nyblom_fmt", "SignBias_fmt", "CoefSig", "Model"
))

#' GARCH Model Display Functions
#'
#' Functions for displaying GARCH model comparison results and coefficient
#' tables in various formats.
#'
#' @name display_garch
#' @keywords internal
NULL


# ==============================================================================
# Internal Utility Functions
# ==============================================================================

#' Format p-value for display
#' @noRd
.fmt_pval <- function(p, digits = 3) {
  if (is.na(p)) {
    return("NA")
  }
  if (p < 0.001) {
    "<0.001"
  } else {
    sprintf(paste0("%.", digits, "f"), p)
  }
}


#' Find indices of best values in comparison data frame
#' @noRd
.find_best_indices <- function(df) {
  # Safe wrappers that handle all-NA columns
  safe_which_min <- function(x) {
    if (all(is.na(x))) return(NA_integer_)
    which.min(x)
  }
  safe_which_max <- function(x) {
    if (all(is.na(x))) return(NA_integer_)
    which.max(x)
  }

  list(
    aic = safe_which_min(df$AIC),
    bic = safe_which_min(df$BIC),
    aicc = safe_which_min(df$AICC),
    wlb1 = safe_which_max(df$WLB1),
    wlb2 = safe_which_max(df$WLB2),
    wlb3 = safe_which_max(df$WLB3),
    nyblom = safe_which_min(df$Nyblom),
    signbias = safe_which_max(df$SignBias)
  )
}


# ==============================================================================
# Base R Fallback Functions
# ==============================================================================

#' Base R fallback for GARCH comparison table
#' @noRd
.table_garch_base <- function(results, title = "GARCH Model Comparison") {
  cat("\n", title, "\n")
  cat(paste(rep("-", nchar(title)), collapse = ""), "\n")
  cat("Distribution:", results$distribution, "\n")
  cat("Sample size:", results$n, "\n\n")

  # Select display columns (exclude internal columns)
  df <- results$comparison
  display_cols <- c("Model", "AIC", "AICC", "BIC")

  # Add WLB columns ONLY if not all NA (hide when WeightedPortTest missing)
  if (!all(is.na(df$WLB1))) {
    display_cols <- c(display_cols, "WLB1", "WLB2", "WLB3")
  }

  # Add other diagnostics
  display_cols <- c(display_cols, "Nyblom", "SignBias")

  # Print with formatting
  print(df[, display_cols, drop = FALSE], row.names = FALSE, digits = 3)

  # Best models
  cat("\nBest by AIC: ", df$Model[which.min(df$AIC)], "\n")
  cat("Best by BIC: ", df$Model[which.min(df$BIC)], "\n")

  invisible(df)
}


#' Base R fallback for coefficient table
#' @noRd
.table_coef_base <- function(fit, title = NULL) {
  coef_mat <- fit@fit$matcoef
  df <- as.data.frame(coef_mat)
  colnames(df) <- c("Estimate", "Std.Error", "t value", "Pr(>|t|)")

  arch_order <- fit@model$modelinc["alpha"]
  garch_order <- fit@model$modelinc["beta"]

  if (is.null(title)) {
    title <- paste0(.garch_label(arch_order, garch_order), " Coefficients")
  }

  cat("\n", title, "\n")
  cat(paste(rep("-", nchar(title)), collapse = ""), "\n")
  cat("Distribution:", fit@model$modeldesc$distribution, "\n\n")
  print(df, digits = 4)

  invisible(df)
}


# ==============================================================================
# gt Table Display
# ==============================================================================

#' Display GARCH Comparison Table with gt
#'
#' Creates a publication-ready comparison table using the \pkg{gt} package
#' with color-coded highlighting for best values and diagnostic failures.
#'
#' @param results A \code{garch_comparison} object from
#'   \code{\link{compare_garch}}.
#' @param title Character string for table title. Default is
#'   "GARCH Model Comparison".
#' @param color_best Color for best/passing values. Default is "seagreen".
#' @param color_fail Color for failing values. Default is "tomato".
#'
#' @return A \code{gt} table object.
#'
#' @details
#' The table displays:
#' \itemize{
#'   \item \strong{Information Criteria}: AIC, AICc, BIC (lower is better,
#'     minimum highlighted in green)
#'   \item \strong{Weighted Ljung-Box}: P-values at lags 1-3 (higher is better,
#'     red if < 0.05)
#'   \item \strong{Nyblom}: Stability test statistic (lower is better,
#'     red if >= critical value)
#'   \item \strong{Sign Bias}: Joint test p-value (higher is better,
#'     red if < 0.05)
#'   \item \strong{Coef}: Significant/total coefficients (green if all
#'     significant)
#' }
#'
#' @seealso \code{\link{compare_garch}}, \code{\link{table_garch_cli}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' garch_gen <- make_gen_garch(omega = 0.1, alpha = 0.15, beta = 0.8)
#' y <- garch_gen(1000)
#' results <- compare_garch(y)
#'
#' # Display table
#' table_garch_gt(results)
#'
#' # Custom colors
#' table_garch_gt(results, color_best = "darkgreen", color_fail = "red")
#' }
table_garch_gt <- function(results,
                           title = "GARCH Model Comparison",
                           color_best = "seagreen",
                           color_fail = "tomato") {

  if (!requireNamespace("gt", quietly = TRUE)) {
    warning("Package 'gt' not installed. Using basic text output.",
            "\n  Install gt for publication-ready tables: install.packages('gt')",
            call. = FALSE)
    return(.table_garch_base(results, title = title))
  }

  if (!inherits(results, "garch_comparison")) {
    stop("results must be a 'garch_comparison' object from compare_garch()")
  }

  df <- results$comparison

  if (nrow(df) == 0) {
    stop("No models to display (comparison table is empty)")
  }

  # Check if WLB data is available
  has_wlb <- !all(is.na(df$WLB1))

  # Build display dataframe with formatted columns
  display_df <- df
  display_df$WLB1_fmt <- sapply(df$WLB1, .fmt_pval)
  display_df$WLB2_fmt <- sapply(df$WLB2, .fmt_pval)
  display_df$WLB3_fmt <- sapply(df$WLB3, .fmt_pval)
  display_df$Nyblom_fmt <- ifelse(is.na(df$Nyblom), "NA", sprintf("%.2f", df$Nyblom))
  display_df$Nyblom_pass <- !is.na(df$Nyblom) & !is.na(df$Nyblom_crit) & df$Nyblom < df$Nyblom_crit
  display_df$SignBias_fmt <- sapply(df$SignBias, .fmt_pval)
  display_df$CoefSig <- paste0(df$n_sig, "/", df$n_coef)
  display_df$CoefSig_pass <- df$n_sig == df$n_coef

  tbl <- gt::gt(display_df)
  tbl <- gt::tab_header(tbl,
    title = title,
    subtitle = paste0("Distribution: ", results$distribution)
  )
  tbl <- gt::cols_hide(tbl, columns = c(
    ARCH, GARCH, WLB1, WLB2, WLB3, Nyblom, Nyblom_crit,
    Nyblom_pass, SignBias, n_sig, n_coef, CoefSig_pass
  ))
  tbl <- gt::cols_label(tbl,
    Model = "Model",
    AIC = "AIC",
    AICC = "AICc",
    BIC = "BIC",
    WLB1_fmt = "Lag 1",
    WLB2_fmt = "Lag 2",
    WLB3_fmt = "Lag 3",
    Nyblom_fmt = "Nyblom",
    SignBias_fmt = "Sign Bias",
    CoefSig = "Coef"
  )
  tbl <- gt::tab_spanner(tbl, label = "Info Criteria", columns = c(AIC, AICC, BIC))
  tbl <- gt::tab_spanner(tbl, label = "Weighted Ljung-Box", columns = c(WLB1_fmt, WLB2_fmt, WLB3_fmt))
  tbl <- gt::tab_spanner(tbl, label = "Diagnostics", columns = c(Nyblom_fmt, SignBias_fmt))
  tbl <- gt::fmt_number(tbl, columns = c(AIC, AICC, BIC), decimals = 2)
  tbl <- gt::cols_align(tbl, align = "right", columns = c(AIC, AICC, BIC, WLB1_fmt, WLB2_fmt, WLB3_fmt, Nyblom_fmt, SignBias_fmt))
  tbl <- gt::cols_align(tbl, align = "center", columns = CoefSig)
  tbl <- gt::cols_align(tbl, align = "left", columns = Model)
  tbl <- gt::tab_style(tbl, style = gt::cell_text(weight = "bold"), locations = gt::cells_body(columns = Model))

  # Hide WLB columns if all NA (WeightedPortTest not installed)
  if (!has_wlb) {
    tbl <- gt::cols_hide(tbl, columns = c(WLB1_fmt, WLB2_fmt, WLB3_fmt))
  }

  # Highlight best IC (green bold) - with NA handling
  tbl <- gt::tab_style(tbl,
    style = gt::cell_text(color = color_best, weight = "bold"),
    locations = gt::cells_body(columns = AIC, rows = !is.na(AIC) & AIC == min(AIC, na.rm = TRUE))
  )
  tbl <- gt::tab_style(tbl,
    style = gt::cell_text(color = color_best, weight = "bold"),
    locations = gt::cells_body(columns = AICC, rows = !is.na(AICC) & AICC == min(AICC, na.rm = TRUE))
  )
  tbl <- gt::tab_style(tbl,
    style = gt::cell_text(color = color_best, weight = "bold"),
    locations = gt::cells_body(columns = BIC, rows = !is.na(BIC) & BIC == min(BIC, na.rm = TRUE))
  )

  # Highlight best WLB p-values (green bold) and failures (red) - only if WLB available
  if (has_wlb) {
    tbl <- gt::tab_style(tbl,
      style = gt::cell_text(color = color_best, weight = "bold"),
      locations = gt::cells_body(columns = WLB1_fmt, rows = !is.na(WLB1) & WLB1 == max(WLB1, na.rm = TRUE))
    )
    tbl <- gt::tab_style(tbl,
      style = gt::cell_text(color = color_best, weight = "bold"),
      locations = gt::cells_body(columns = WLB2_fmt, rows = !is.na(WLB2) & WLB2 == max(WLB2, na.rm = TRUE))
    )
    tbl <- gt::tab_style(tbl,
      style = gt::cell_text(color = color_best, weight = "bold"),
      locations = gt::cells_body(columns = WLB3_fmt, rows = !is.na(WLB3) & WLB3 == max(WLB3, na.rm = TRUE))
    )
    tbl <- gt::tab_style(tbl,
      style = gt::cell_text(color = color_fail),
      locations = gt::cells_body(columns = WLB1_fmt, rows = !is.na(WLB1) & WLB1 < 0.05)
    )
    tbl <- gt::tab_style(tbl,
      style = gt::cell_text(color = color_fail),
      locations = gt::cells_body(columns = WLB2_fmt, rows = !is.na(WLB2) & WLB2 < 0.05)
    )
    tbl <- gt::tab_style(tbl,
      style = gt::cell_text(color = color_fail),
      locations = gt::cells_body(columns = WLB3_fmt, rows = !is.na(WLB3) & WLB3 < 0.05)
    )
  }

  # Nyblom: lower is better, red if unstable - with NA handling
  tbl <- gt::tab_style(tbl,
    style = gt::cell_text(color = color_best, weight = "bold"),
    locations = gt::cells_body(columns = Nyblom_fmt, rows = !is.na(Nyblom) & Nyblom == min(Nyblom, na.rm = TRUE))
  )
  tbl <- gt::tab_style(tbl,
    style = gt::cell_text(color = color_fail),
    locations = gt::cells_body(columns = Nyblom_fmt, rows = !is.na(Nyblom) & !is.na(Nyblom_crit) & Nyblom >= Nyblom_crit)
  )

  # Sign Bias: higher p-value is better - with NA handling
  tbl <- gt::tab_style(tbl,
    style = gt::cell_text(color = color_best, weight = "bold"),
    locations = gt::cells_body(columns = SignBias_fmt, rows = !is.na(SignBias) & SignBias == max(SignBias, na.rm = TRUE))
  )
  tbl <- gt::tab_style(tbl,
    style = gt::cell_text(color = color_fail),
    locations = gt::cells_body(columns = SignBias_fmt, rows = !is.na(SignBias) & SignBias < 0.05)
  )

  # Coefficient significance
  tbl <- gt::tab_style(tbl,
    style = gt::cell_text(color = color_best, weight = "bold"),
    locations = gt::cells_body(columns = CoefSig, rows = CoefSig_pass == TRUE)
  )
  tbl <- gt::tab_style(tbl,
    style = gt::cell_text(color = color_fail),
    locations = gt::cells_body(columns = CoefSig, rows = CoefSig_pass == FALSE)
  )

  # Footnotes
  tbl <- gt::tab_source_note(tbl, "Green = best/pass. Red = fail. Weighted LB lags are model-dependent.")
  tbl <- gt::tab_source_note(tbl, gt::html("Nyblom: H<sub>0</sub> = parameters stable. Green = stable (&lt; 5% critical value)."))
  tbl <- gt::tab_source_note(tbl, gt::html("Sign Bias: H<sub>0</sub> = symmetric effects. If rejected, consider EGARCH/GJR-GARCH."))
  tbl <- gt::tab_options(tbl,
    table.font.size = gt::px(12),
    heading.title.font.size = gt::px(16),
    heading.subtitle.font.size = gt::px(13),
    column_labels.font.size = gt::px(11),
    column_labels.padding = gt::px(8),
    data_row.padding = gt::px(7),
    source_notes.font.size = gt::px(10),
    table.width = gt::pct(100)
  )

  tbl
}


# ==============================================================================
# CLI Table Display
# ==============================================================================

#' Display GARCH Comparison Table in Console
#'
#' Creates a formatted console table using the \pkg{cli} package with
#' color-coded highlighting.
#'
#' @param results A \code{garch_comparison} object from
#'   \code{\link{compare_garch}}.
#' @param show_signbias Logical. Include Sign Bias column? Default is
#'   \code{TRUE}.
#'
#' @return Invisibly returns the comparison data frame.
#'
#' @details
#' Colors in terminal output:
#' \itemize{
#'   \item \strong{Green}: Best value for that criterion, or passing diagnostic
#'   \item \strong{Red}: Failing diagnostic (p < 0.05 for tests, or
#'     Nyblom >= critical value)
#' }
#'
#' @seealso \code{\link{compare_garch}}, \code{\link{table_garch_gt}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' garch_gen <- make_gen_garch(omega = 0.1, alpha = 0.15, beta = 0.8)
#' y <- garch_gen(1000)
#' results <- compare_garch(y)
#'
#' # Display in console
#' table_garch_cli(results)
#' }
table_garch_cli <- function(results, show_signbias = TRUE) {

  if (!requireNamespace("cli", quietly = TRUE)) {
    warning("Package 'cli' not installed. Using basic text output.",
            "\n  Install cli for colored terminal output: install.packages('cli')",
            call. = FALSE)
    return(.table_garch_base(results, title = "GARCH Model Comparison"))
  }

  if (!inherits(results, "garch_comparison")) {
    stop("results must be a 'garch_comparison' object from compare_garch()")
  }

  df <- results$comparison

  if (nrow(df) == 0) {
    cli::cli_alert_warning("No models to display (comparison table is empty)")
    return(invisible(df))
  }

  # Check if WLB data is available
  has_wlb <- !all(is.na(df$WLB1))

  best <- .find_best_indices(df)

  # Column widths
  w_model <- 12
  w_ic <- 10
  w_wlb <- 7
  w_nyb <- 7
  w_sb <- 9
  w_coef <- 5
  sep <- " "

  # Helper to check if index matches best (handles NA)
  is_best <- function(i, best_idx) {
    !is.na(best_idx) && i == best_idx
  }

  # Formatting helpers (apply styling after sprintf)
  fmt_ic <- function(val, is_best_val) {
    if (is.na(val)) {
      s <- sprintf("%*s", w_ic, "NA")
      return(s)
    }
    s <- sprintf("%*.2f", w_ic, val)
    if (is_best_val) cli::col_green(cli::style_bold(s)) else s
  }

  fmt_pval_cli <- function(val, is_best_val, width = w_wlb) {
    if (is.na(val)) {
      s <- sprintf("%*s", width, "NA")
      return(s)
    }
    if (val < 0.001) {
      s <- sprintf("%*s", width, "<.001")
      cli::col_red(s)
    } else if (val < 0.05) {
      s <- sprintf("%*.3f", width, val)
      cli::col_red(s)
    } else if (is_best_val) {
      s <- sprintf("%*.3f", width, val)
      cli::col_green(cli::style_bold(s))
    } else {
      sprintf("%*.3f", width, val)
    }
  }

  fmt_nyblom <- function(stat, crit, is_best_val) {
    if (is.na(stat)) {
      s <- sprintf("%*s", w_nyb, "NA")
      return(s)
    }
    s <- sprintf("%*.2f", w_nyb, stat)
    if (is_best_val) {
      cli::col_green(cli::style_bold(s))
    } else if (!is.na(crit) && stat >= crit) {
      cli::col_red(s)
    } else {
      s
    }
  }

  fmt_coef <- function(n_sig, n_coef) {
    s <- sprintf("%*s", w_coef, paste0(n_sig, "/", n_coef))
    if (n_sig == n_coef) cli::col_green(s) else cli::col_red(s)
  }

  fmt_model <- function(name) {
    sprintf("%-*s", w_model, name)
  }

  # Header
  cli::cli_h2("GARCH Model Comparison")
  cli::cli_text("Distribution: {.val {results$distribution}}")
  cli::cli_text("")

  # Build header string (consistent order: AIC, AICc, BIC)
  header <- paste0(
    sprintf("%-*s", w_model, "Model"), sep,
    sprintf("%*s", w_ic, "AIC"), sep,
    sprintf("%*s", w_ic, "AICc"), sep,
    sprintf("%*s", w_ic, "BIC")
  )

  # Add WLB columns only if available
  if (has_wlb) {
    header <- paste0(
      header, sep,
      sprintf("%*s", w_wlb, "WLB1"), sep,
      sprintf("%*s", w_wlb, "WLB2"), sep,
      sprintf("%*s", w_wlb, "WLB3")
    )
  }

  header <- paste0(header, sep, sprintf("%*s", w_nyb, "Nyblom"))

  if (show_signbias) {
    header <- paste0(header, sep, sprintf("%*s", w_sb, "SignBias"))
  }
  header <- paste0(header, sep, sprintf("%*s", w_coef, "Coef"))

  cat(cli::style_bold(header), "\n")
  cli::cli_rule()

  # Data rows
  for (i in seq_len(nrow(df))) {
    r <- df[i, ]

    row <- paste0(
      fmt_model(r$Model), sep,
      fmt_ic(r$AIC, is_best(i, best$aic)), sep,
      fmt_ic(r$AICC, is_best(i, best$aicc)), sep,
      fmt_ic(r$BIC, is_best(i, best$bic))
    )

    # Add WLB columns only if available
    if (has_wlb) {
      row <- paste0(
        row, sep,
        fmt_pval_cli(r$WLB1, is_best(i, best$wlb1)), sep,
        fmt_pval_cli(r$WLB2, is_best(i, best$wlb2)), sep,
        fmt_pval_cli(r$WLB3, is_best(i, best$wlb3))
      )
    }

    row <- paste0(row, sep, fmt_nyblom(r$Nyblom, r$Nyblom_crit, is_best(i, best$nyblom)))

    if (show_signbias) {
      row <- paste0(row, sep, fmt_pval_cli(r$SignBias, is_best(i, best$signbias), width = w_sb))
    }
    row <- paste0(row, sep, fmt_coef(r$n_sig, r$n_coef))

    cat(row, "\n")
  }

  cli::cli_rule()
  cli::cli_text(cli::col_silver("Green = best/pass | Red = fail | Nyblom: green = stable"))

  invisible(df)
}


# ==============================================================================
# Coefficient Table
# ==============================================================================

#' Display GARCH Coefficient Table
#'
#' Creates a formatted table of coefficient estimates, standard errors,
#' t-values, and p-values for a fitted GARCH model.
#'
#' @param fit A \code{uGARCHfit} object from \pkg{rugarch}, typically
#'   extracted from the \code{fits} element of a \code{garch_comparison}
#'   object.
#' @param title Optional custom title. If \code{NULL}, auto-generated from
#'   model specification.
#' @param color_sig Color for significant p-values (< 0.05). Default is
#'   "seagreen".
#' @param color_nonsig Color for non-significant p-values. Default is
#'   "tomato".
#'
#' @return A \code{gt} table object.
#'
#' @details
#' Significance codes in the table:
#' \itemize{
#'   \item \code{***}: p < 0.001
#'   \item \code{**}: p < 0.01
#'   \item \code{*}: p < 0.05
#'   \item \code{.}: p < 0.1
#' }
#'
#' @seealso \code{\link{compare_garch}}, \code{\link{table_garch_gt}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' garch_gen <- make_gen_garch(omega = 0.1, alpha = 0.15, beta = 0.8)
#' y <- garch_gen(1000)
#' results <- compare_garch(y)
#'
#' # Get coefficient table for best model
#' best_fit <- results$fits[["GARCH(1,1)"]]
#' table_coef_gt(best_fit)
#' }
table_coef_gt <- function(fit,
                          title = NULL,
                          color_sig = "seagreen",
                          color_nonsig = "tomato") {

  if (!requireNamespace("gt", quietly = TRUE)) {
    warning("Package 'gt' not installed. Using basic text output.",
            "\n  Install gt for publication-ready tables: install.packages('gt')",
            call. = FALSE)
    return(.table_coef_base(fit, title = title))
  }

  if (!inherits(fit, "uGARCHfit")) {
    stop("fit must be a 'uGARCHfit' object from rugarch::ugarchfit()")
  }

  # Extract model orders for title
  arch_order <- fit@model$modelinc["alpha"]
  garch_order <- fit@model$modelinc["beta"]

  if (is.null(title)) {
    title <- paste0(.garch_label(arch_order, garch_order), " Coefficient Estimates")
  }

  # Extract coefficient matrix
  coef_mat <- fit@fit$matcoef

  df <- as.data.frame(coef_mat)
  colnames(df) <- c("Estimate", "Std_Error", "t_value", "p_val")
  df$Parameter <- rownames(coef_mat)

  df <- df[, c("Parameter", "Estimate", "Std_Error", "t_value", "p_val")]

  df$p_display <- sapply(df$p_val, function(p) {
    if (is.na(p)) return("NA")
    if (p < 0.001) "<0.001" else sprintf("%.4f", p)
  })

  df$Sig <- sapply(df$p_val, function(p) {
    if (is.na(p)) return("")
    if (p < 0.001) "***"
    else if (p < 0.01) "**"
    else if (p < 0.05) "*"
    else if (p < 0.1) "."
    else ""
  })

  tbl <- gt::gt(df)
  tbl <- gt::tab_header(tbl,
    title = title,
    subtitle = paste0("Distribution: ", fit@model$modeldesc$distribution)
  )
  tbl <- gt::cols_hide(tbl, columns = c(p_val))
  tbl <- gt::cols_label(tbl,
    Parameter = "Parameter",
    Estimate = "Estimate",
    Std_Error = "Std. Error",
    t_value = "t value",
    p_display = "p-value",
    Sig = ""
  )
  tbl <- gt::fmt_number(tbl, columns = c(Estimate, Std_Error, t_value), decimals = 4)
  tbl <- gt::cols_align(tbl, align = "left", columns = Parameter)
  tbl <- gt::cols_align(tbl, align = "right", columns = c(Estimate, Std_Error, t_value, p_display))
  tbl <- gt::cols_align(tbl, align = "center", columns = Sig)
  tbl <- gt::tab_style(tbl,
    style = gt::cell_text(weight = "bold"),
    locations = gt::cells_body(columns = Parameter)
  )
  tbl <- gt::tab_style(tbl,
    style = gt::cell_text(color = color_sig, weight = "bold"),
    locations = gt::cells_body(columns = p_display, rows = !is.na(p_val) & p_val < 0.05)
  )
  tbl <- gt::tab_style(tbl,
    style = gt::cell_text(color = color_nonsig),
    locations = gt::cells_body(columns = p_display, rows = !is.na(p_val) & p_val >= 0.05)
  )
  tbl <- gt::tab_source_note(tbl, "Significance: *** p < 0.001, ** p < 0.01, * p < 0.05, . p < 0.1")
  tbl <- gt::tab_options(tbl,
    table.font.size = gt::px(12),
    heading.title.font.size = gt::px(15),
    heading.subtitle.font.size = gt::px(12),
    column_labels.padding = gt::px(8),
    data_row.padding = gt::px(6),
    source_notes.font.size = gt::px(10)
  )

  tbl
}
