# Console tab module

mod_console_ui <- function(id) {
  ns <- NS(id)

  tabPanel("Console",
    br(),
    fluidRow(
      column(4,
        wellPanel(
          h4("R Console"),
          p(class = "text-body-secondary",
            "Run multiline R code in a persistent session environment. Shift+Enter runs all; Ctrl/Cmd+Enter runs current line."),
          textAreaInput(
            ns("code"), "Code",
            value = "# Example\nsummary(cars)\nplot(cars)",
            width = "100%", height = "320px",
            placeholder = "Type R code here"
          ),
          tags$script(HTML(sprintf(
            "(function() {
               var codeId = '%s';
               var runId = '%s';
               var runLineId = '%s';
               var el = document.getElementById(codeId);
               if (!el || el.dataset.ctrlEnterBound === '1') return;
               el.dataset.ctrlEnterBound = '1';
               el.addEventListener('keydown', function(e) {
                 if (e.shiftKey && e.key === 'Enter') {
                   e.preventDefault();
                   var runBtn = document.getElementById(runId);
                   if (runBtn) runBtn.click();
                   return;
                 }
                 if ((e.ctrlKey || e.metaKey) && e.key === 'Enter') {
                   e.preventDefault();
                   var start = el.selectionStart || 0;
                   var val = el.value || '';
                   var lineStart = val.lastIndexOf('\\n', Math.max(0, start - 1));
                   lineStart = (lineStart === -1) ? 0 : lineStart + 1;
                   var lineEnd = val.indexOf('\\n', start);
                   lineEnd = (lineEnd === -1) ? val.length : lineEnd;
                   var line = val.slice(lineStart, lineEnd);
                   if (window.Shiny && typeof Shiny.setInputValue === 'function') {
                     Shiny.setInputValue(runLineId, { code: line, nonce: Date.now() }, { priority: 'event' });
                   }
                   var nextStart = lineEnd < val.length ? lineEnd + 1 : val.length;
                   el.selectionStart = nextStart;
                   el.selectionEnd = nextStart;
                   el.focus();
                 }
               });
             })();",
            ns("code"), ns("run"), ns("run_line")
          ))),
          fluidRow(
            column(4,
              actionButton(ns("run"), "Run", icon = icon("play"),
                           class = "btn-primary", width = "100%")
            ),
            column(4,
              actionButton(ns("clear"), "Clear Output", icon = icon("eraser"),
                           class = "btn-outline-secondary", width = "100%")
            ),
            column(4,
              actionButton(ns("reset"), "Reset Env", icon = icon("rotate-left"),
                           class = "btn-outline-danger", width = "100%")
            )
          ),
          br(),
          p(class = "text-body-secondary", style = "margin-bottom: 0;",
            "Full session mode: code runs with normal R capabilities.")
        )
      ),
      column(8,
        wellPanel(
          h5("Console Output"),
          verbatimTextOutput(ns("stdout"), placeholder = TRUE)
        ),
        wellPanel(
          h5("Latest Plot"),
          plotOutput(ns("plot"), height = "420px")
        )
      )
    )
  )
}

mod_console_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    console_env <- reactiveVal(new.env(parent = .GlobalEnv))
    console_text <- reactiveVal("Console ready. Enter R code and click Run.")
    console_plot <- reactiveVal(NULL)

    .evaluate_code <- function(code, append = FALSE) {
      if (is.null(code) || !nzchar(trimws(code))) {
        if (!isTRUE(append)) showNotification("Enter R code to run.", type = "warning")
        return()
      }

      exprs <- tryCatch(parse(text = code), error = function(e) e)
      if (inherits(exprs, "error")) {
        msg <- paste0("Parse error: ", conditionMessage(exprs))
        if (isTRUE(append) && nzchar(console_text())) {
          console_text(paste(console_text(), msg, sep = "\n"))
        } else {
          console_text(msg)
          console_plot(NULL)
        }
        return()
      }

      env <- console_env()
      out <- character(0)

      grDevices::pdf(NULL)
      grDevices::dev.control(displaylist = "enable")
      on.exit({
        while (grDevices::dev.cur() > 1L) grDevices::dev.off()
      }, add = TRUE)

      for (expr in exprs) {
        src <- paste(deparse(expr, width.cutoff = 500L), collapse = "\n")
        src_lines <- strsplit(src, "\n", fixed = TRUE)[[1]]
        out <- c(out, paste0("> ", src_lines[1]))
        if (length(src_lines) > 1) {
          out <- c(out, paste0("+ ", src_lines[-1]))
        }

        chunk <- withCallingHandlers(
          tryCatch(
            capture.output({
              vis <- withVisible(eval(expr, envir = env))
              if (isTRUE(vis$visible)) print(vis$value)
            }),
            error = function(e) paste0("Error: ", conditionMessage(e))
          ),
          warning = function(w) {
            out <<- c(out, paste0("Warning: ", conditionMessage(w)))
            invokeRestart("muffleWarning")
          },
          message = function(m) {
            out <<- c(out, paste0("Message: ", conditionMessage(m)))
            invokeRestart("muffleMessage")
          }
        )

        if (length(chunk) > 0) out <- c(out, chunk)
        out <- c(out, "")
      }

      recorded_plot <- tryCatch(grDevices::recordPlot(), error = function(e) NULL)
      has_plot <- !is.null(recorded_plot) && !is.null(recorded_plot[[1]]) &&
        length(recorded_plot[[1]]) > 0

      new_text <- paste(out, collapse = "\n")
      if (isTRUE(append) && nzchar(console_text())) {
        console_text(paste(console_text(), new_text, sep = "\n"))
      } else {
        console_text(new_text)
      }

      if (isTRUE(has_plot)) {
        console_plot(recorded_plot)
      } else if (!isTRUE(append)) {
        console_plot(NULL)
      }
    }

    observeEvent(input$run, {
      .evaluate_code(input$code, append = FALSE)
    })

    observeEvent(input$run_line, {
      payload <- input$run_line
      code <- if (is.list(payload) && !is.null(payload$code)) payload$code else payload
      .evaluate_code(code, append = TRUE)
    })

    observeEvent(input$clear, {
      console_text("Output cleared.")
      console_plot(NULL)
    })

    observeEvent(input$reset, {
      console_env(new.env(parent = .GlobalEnv))
      console_text("Environment reset. All console-created objects were removed.")
      console_plot(NULL)
      showNotification("Console environment reset.", type = "message")
    })

    output$stdout <- renderText({
      console_text()
    })

    output$plot <- renderPlot(bg = "transparent", {
      rp <- console_plot()
      if (is.null(rp)) {
        plot.new()
        text(0.5, 0.5, "No plot yet", cex = 1.2, col = "grey50")
      } else {
        replayPlot(rp)
      }
    })
  })
}
