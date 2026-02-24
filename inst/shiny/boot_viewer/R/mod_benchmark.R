# Module: Benchmark tab
# Runtime benchmark across tswge and tstse WBG implementations.

mod_benchmark_ui <- function(id) {
  ns <- NS(id)

  tabPanel(
    "Benchmark",
    br(),
    h4("WBG Performance Benchmark"),
    p(
      class = "text-body-secondary",
      "Compare runtime on matched AR(1) simulations across tswge and tstse implementations."
    ),

    fluidRow(
      column(
        4,
        wellPanel(
          h5("Scenario"),
          numericInput(ns("bench_n"), "Series length (n)", value = 500, min = 20, step = 10),
          numericInput(ns("bench_nb"), "Bootstrap replicates (nb)", value = 199, min = 19, step = 10),
          numericInput(ns("bench_phi"), "AR(1) phi", value = 0.9, min = -0.99, max = 0.99, step = 0.01),
          numericInput(ns("bench_nsims"), "Num simulations (Y)", value = 10, min = 1, max = 500, step = 1),
          numericInput(ns("bench_seed"), "Seed", value = 42, min = 1, step = 1),
          hr(),
          h5("Methods"),
          checkboxGroupInput(
            ns("bench_methods"),
            label = NULL,
            choices = c(
              "tswge::wbg.boot.wge()" = "tswge",
              "tstse::wbg_boot()" = "wbg_boot",
              "tstse::wbg_boot(cores = 0)" = "wbg_boot_cores0",
              "tstse::wbg_boot_fast()" = "wbg_boot_fast"
            ),
            selected = c("tswge", "wbg_boot", "wbg_boot_cores0", "wbg_boot_fast")
          ),
          p(class = "text-body-secondary", style = "margin-top: -4px; margin-bottom: 8px;",
            "Select at least two methods."),
          p(class = "text-body-secondary", style = "margin-top: 10px; margin-bottom: 8px;",
            "Innovation distribution: Normal(0, 1); maxp = 5; criterion = AIC"),
          actionButton(
            ns("bench_run"),
            "Run Benchmark",
            icon = icon("play"),
            class = "btn-primary btn-block"
          ),
          hr(),
          textOutput(ns("bench_status"))
        )
      ),
      column(
        8,
        plotOutput(ns("bench_plot"), height = "380px"),
        DTOutput(ns("bench_table"))
      )
    )
  )
}

mod_benchmark_server <- function(id, con = NULL) {
  moduleServer(id, function(input, output, session) {
    scenario_maxp <- 5L

    .safe_int <- function(x, default, min_val = NULL, max_val = NULL) {
      val <- suppressWarnings(as.integer(x))
      if (is.na(val)) val <- as.integer(default)
      if (!is.null(min_val)) val <- max(val, as.integer(min_val))
      if (!is.null(max_val)) val <- min(val, as.integer(max_val))
      val
    }

    .safe_num <- function(x, default, min_val = NULL, max_val = NULL) {
      val <- suppressWarnings(as.numeric(x))
      if (is.na(val)) val <- as.numeric(default)
      if (!is.null(min_val)) val <- max(val, as.numeric(min_val))
      if (!is.null(max_val)) val <- min(val, as.numeric(max_val))
      val
    }

    .generate_series <- function(n, phi, series_seed) {
      tryCatch(
        gen_aruma_flex(
          n = n,
          phi = phi,
          innov_gen = make_gen_norm(sd = 1),
          seed = series_seed,
          plot = FALSE
        )$y,
        error = function(e) {
          set.seed(series_seed)
          as.numeric(stats::arima.sim(list(ar = phi), n = n))
        }
      )
    }

    bench_results <- eventReactive(input$bench_run, {
      n <- .safe_int(input$bench_n, default = 500L, min_val = 20L)
      nb <- .safe_int(input$bench_nb, default = 199L, min_val = 19L)
      phi <- .safe_num(input$bench_phi, default = 0.9, min_val = -0.99, max_val = 0.99)
      nsims <- .safe_int(input$bench_nsims, default = 10L, min_val = 1L, max_val = 500L)
      seed <- .safe_int(input$bench_seed, default = 42L, min_val = 1L)

      # Required by wbg_boot_fast() guard.
      if (n <= scenario_maxp + 10L) n <- scenario_maxp + 11L

      set.seed(seed)
      series_seeds <- sample.int(.Machine$integer.max, nsims)
      boot_seeds <- sample.int(.Machine$integer.max, nsims)

      all_methods <- list(
        tswge = list(
          label = "tswge::wbg.boot.wge()",
          run = function(x, run_seed) {
            if (!suppressPackageStartupMessages(requireNamespace("tswge", quietly = TRUE))) {
              stop("Package 'tswge' is not installed", call. = FALSE)
            }
            tswge::wbg.boot.wge(x, nb = nb, pvalue = TRUE, sn = run_seed)
          }
        ),
        wbg_boot = list(
          label = "tstse::wbg_boot()",
          run = function(x, run_seed) {
            wbg_boot(
              x,
              nb = nb,
              maxp = scenario_maxp,
              method = "burg",
              type = "aic",
              use_fast = TRUE,
              cores = 1L,
              seed = run_seed
            )
          }
        ),
        wbg_boot_cores0 = list(
          label = "tstse::wbg_boot(cores = 0)",
          run = function(x, run_seed) {
            wbg_boot(
              x,
              nb = nb,
              maxp = scenario_maxp,
              method = "burg",
              type = "aic",
              use_fast = TRUE,
              cores = 0L,
              seed = run_seed
            )
          }
        ),
        wbg_boot_fast = list(
          label = "tstse::wbg_boot_fast()",
          run = function(x, run_seed) {
            wbg_boot_fast(
              x,
              nb = nb,
              maxp = scenario_maxp,
              criterion = "aic",
              seed = run_seed
            )
          }
        )
      )

      selected_ids <- input$bench_methods
      if (is.null(selected_ids)) selected_ids <- character(0)
      selected_ids <- intersect(names(all_methods), selected_ids)

      if (length(selected_ids) < 2L) {
        showNotification("Select at least two methods to run the benchmark.", type = "error")
        return(
          list(
            meta = list(n = n, nb = nb, phi = phi, nsims = nsims, seed = seed),
            results = data.frame(),
            error = "Select at least two methods."
          )
        )
      }

      methods <- all_methods[selected_ids]

      out_rows <- vector("list", length(methods))
      total_steps <- length(methods) * nsims

      withProgress(message = "Running benchmark...", value = 0, {
        for (i in seq_along(methods)) {
          m <- methods[[i]]

          elapsed_total <- 0
          user_total <- 0
          system_total <- 0
          n_ok <- 0L
          first_error <- ""

          for (j in seq_len(nsims)) {
            incProgress(
              1 / total_steps,
              detail = sprintf("%s [%d/%d]", m$label, j, nsims)
            )

            x <- .generate_series(n = n, phi = phi, series_seed = series_seeds[j])
            run_seed <- boot_seeds[j]

            gc(FALSE)
            run_result <- tryCatch({
              t <- system.time({
                m$run(x, run_seed)
              })
              list(
                status = "ok",
                elapsed = as.numeric(t["elapsed"]),
                user = as.numeric(t["user.self"]),
                system = as.numeric(t["sys.self"]),
                note = ""
              )
            }, error = function(e) {
              list(
                status = "error",
                elapsed = 0,
                user = 0,
                system = 0,
                note = conditionMessage(e)
              )
            })

            if (identical(run_result$status, "ok")) {
              n_ok <- n_ok + 1L
              elapsed_total <- elapsed_total + run_result$elapsed
              user_total <- user_total + run_result$user
              system_total <- system_total + run_result$system
            } else if (!nzchar(first_error)) {
              first_error <- sprintf("Sim %d: %s", j, run_result$note)
            }
          }

          status <- if (n_ok == nsims) {
            "ok"
          } else if (n_ok > 0L) {
            "partial"
          } else {
            "error"
          }

          out_rows[[i]] <- data.frame(
            method = m$label,
            total_elapsed_sec = elapsed_total,
            avg_elapsed_sec = elapsed_total / nsims,
            user_sec = user_total,
            system_sec = system_total,
            n_simulations = nsims,
            n_success = n_ok,
            status = status,
            note = first_error,
            stringsAsFactors = FALSE
          )
        }
      })

      list(
        meta = list(n = n, nb = nb, phi = phi, nsims = nsims, seed = seed),
        results = do.call(rbind, out_rows)
      )
    })

    output$bench_status <- renderText({
      obj <- bench_results()
      if (is.null(obj) || is.null(obj$results) || nrow(obj$results) == 0) return("Ready")
      if (!is.null(obj$error) && nzchar(obj$error)) return(obj$error)
      df <- obj$results
      n_ok <- sum(df$status == "ok", na.rm = TRUE)
      n_partial <- sum(df$status == "partial", na.rm = TRUE)
      n_err <- sum(df$status == "error", na.rm = TRUE)
      sprintf("Methods: %d ok, %d partial, %d failed", n_ok, n_partial, n_err)
    })

    output$bench_plot <- renderPlot(bg = "transparent", {
      obj <- bench_results()
      if (is.null(obj) || is.null(obj$results) || nrow(obj$results) == 0) {
        plot.new()
        msg <- if (!is.null(obj$error) && nzchar(obj$error)) obj$error else "Click Run Benchmark to execute."
        text(0.5, 0.5, msg, cex = 1.2, col = "grey50")
        return()
      }

      df <- obj$results
      ok_df <- df[df$n_success > 0 & is.finite(df$total_elapsed_sec), , drop = FALSE]
      if (nrow(ok_df) == 0) {
        plot.new()
        text(0.5, 0.5, "No successful benchmark runs.", cex = 1.2, col = "grey50")
        return()
      }

      ok_df$method <- factor(ok_df$method, levels = ok_df$method[order(ok_df$total_elapsed_sec)])

      p <- ggplot(ok_df, aes(x = method, y = total_elapsed_sec, fill = method)) +
        geom_col(width = 0.72, show.legend = FALSE, alpha = 0.9) +
        geom_text(
          aes(label = sprintf("total %.2fs\navg %.2fs", total_elapsed_sec, avg_elapsed_sec)),
          hjust = -0.05, size = 4.5, lineheight = 0.95
        ) +
        coord_flip(clip = "off") +
        labs(
          x = NULL,
          y = "Total Time (seconds)",
          title = "Runtime by Method",
          subtitle = sprintf(
            "n=%d, nb=%d, phi=%.2f, Y=%d, normal innovations",
            obj$meta$n, obj$meta$nb, obj$meta$phi, obj$meta$nsims
          )
        ) +
        viewer_plot_theme(base_size = 14) +
        theme(plot.margin = margin(12, 100, 12, 12))

      print(p)
    })

    output$bench_table <- renderDT({
      obj <- bench_results()
      if (is.null(obj) || is.null(obj$results) || nrow(obj$results) == 0) {
        msg <- if (!is.null(obj$error) && nzchar(obj$error)) obj$error else "Click Run Benchmark to execute."
        return(
          datatable(
            data.frame(Message = msg),
            rownames = FALSE,
            options = list(dom = "t")
          )
        )
      }

      df <- obj$results
      dt <- datatable(
        df,
        rownames = FALSE,
        options = list(
          pageLength = 10,
          autoWidth = TRUE,
          scrollX = TRUE
        )
      )

      dt <- formatRound(
        dt,
        columns = c("total_elapsed_sec", "avg_elapsed_sec", "user_sec", "system_sec"),
        digits = 3
      )
      dt
    })
  })
}
