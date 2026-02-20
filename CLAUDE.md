# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

<!-- RSHINY-AGENTS-MD-START -->[R Shiny Docs Index]|root: ./.rshiny-docs|IMPORTANT: Prefer retrieval-led reasoning over pre-training-led reasoning for any R Shiny tasks.|root-files:{ExtendedTask.Rd,MockShinySession.Rd,NS.Rd,Progress.Rd,absolutePanel.Rd,actionButton.Rd,app-state.R,app_template.R,applyInputHandlers.Rd,bind-cache.R,bind-event.R,bindCache.Rd,bindEvent.Rd,bookmark-state-local.R,bookmark-state.R,bookmarkButton.Rd,bootstrap-deprecated.R,bootstrap-layout.R,bootstrap.R,bootstrapLib.Rd,bootstrapPage.Rd,brushOpts.Rd,brushedPoints.Rd,busy-indicators-spinners.R,busy-indicators.R,busyIndicatorOptions.Rd,cache-utils.R,callModule.Rd,checkboxGroupInput.Rd,checkboxInput.Rd,clickOpts.Rd,column.Rd,conditionalPanel.Rd,conditions.R,createRenderFunction.Rd,createWebDependency.Rd,dateInput.Rd,dateRangeInput.Rd,debounce.Rd,deprecated.R,devmode.R,devmode.Rd,diagnose.R,diskCache.Rd,domains.Rd,downloadButton.Rd,downloadHandler.Rd,enableBookmarking.Rd,exportTestValues.Rd,exprToFunction.Rd,extended-task.R,fileInput.Rd,fileupload.R,fillPage.Rd,fillRow.Rd,fixedPage.Rd,flowLayout.Rd,fluidPage.Rd,freezeReactiveValue.Rd,getCurrentOutputInfo.Rd,getCurrentTheme.Rd,getQueryString.Rd,globals.R,graph.R,headerPanel.Rd,helpText.Rd,history.R,hooks.R,html-deps.R,htmlOutput.Rd,httpResponse.Rd,icon.Rd,image-interact-opts.R,image-interact.R,imageutils.R,input-action.R,input-checkbox.R,input-checkboxgroup.R,input-date.R,input-daterange.R,input-file.R,input-numeric.R,input-password.R,input-radiobuttons.R,input-select.R,input-slider.R,input-submit.R,input-text.R,input-textarea.R,input-utils.R,inputPanel.Rd,insert-tab.R,insert-ui.R,insertTab.Rd,insertUI.Rd,invalidateLater.Rd,is.reactivevalues.Rd,isRunning.Rd,isTruthy.Rd,isolate.Rd,jqueryui.R,knitr.R,knitr_methods.Rd,loadSupport.Rd,makeReactiveBinding.Rd,map.R,markOutputAttrs.Rd,markRenderFunction.Rd,markdown.Rd,maskReactiveContext.Rd,memoryCache.Rd,middleware-shiny.R,middleware.R,mock-session.R,modal.R,modalDialog.Rd,moduleServer.Rd,modules.R,navbarPage.Rd,navlistPanel.Rd,notifications.R,numericInput.Rd,observe.Rd,observeEvent.Rd,onBookmark.Rd,onFlush.Rd,onStop.Rd,otel-attr-srcref.R,otel-collect.R,otel-enable.R,otel-error.R,otel-label.R,otel-reactive-update.R,otel-session.R,otel-shiny.R,otel-with.R,outputOptions.Rd,pageWithSidebar.Rd,parseQueryString.Rd,passwordInput.Rd,plotOutput.Rd,plotPNG.Rd,priorityqueue.R,progress.R,radioButtons.Rd,react.R,reactive-domains.R,reactive.Rd,reactiveConsole.Rd,reactiveFileReader.Rd,reactivePoll.Rd,reactiveTimer.Rd,reactiveVal.Rd,reactiveValues.Rd,reactiveValuesToList.Rd,reactives.R,reactlog.Rd,reexports.R,reexports.Rd,registerInputHandler.Rd,registerThemeDependency.Rd,removeInputHandler.Rd,render-cached-plot.R,render-plot.R,render-table.R,renderCachedPlot.Rd,renderDataTable.Rd,renderImage.Rd,renderPlot.Rd,renderPrint.Rd,renderTable.Rd,renderUI.Rd,repeatable.Rd,req.Rd,resourcePaths.Rd,restoreInput.Rd,run-url.R,runApp.Rd,runExample.Rd,runGadget.Rd,runTests.Rd,runUrl.Rd,runapp.R,safeError.Rd,selectInput.Rd,serializers.R,server-input-handlers.R,server-resource-paths.R,server.R,serverInfo.Rd,session.Rd,setBookmarkExclude.Rd,setSerializer.Rd,shiny-options.R,shiny-package.R,shiny-package.Rd,shiny.R,shiny.appobj.Rd,shinyApp.Rd,shinyAppTemplate.Rd,shinyDeprecated.Rd,shinyOptions.Rd,shinyServer.Rd,shinyUI.Rd,shinyapp.R,shinyui.R,shinywrappers.R,showBookmarkUrlModal.Rd,showModal.Rd,showNotification.Rd,showTab.Rd,showcase.R,sidebarLayout.Rd,sizeGrowthRatio.Rd,sliderInput.Rd,snapshot.R,snapshotExclude.Rd,snapshotPreprocessInput.Rd,snapshotPreprocessOutput.Rd,splitLayout.Rd,stacktrace.Rd,staticimports.R,stopApp.Rd,submitButton.Rd,tabPanel.Rd,tabsetPanel.Rd,tar.R,test-export.R,test-server.R,test.R,testServer.Rd,textAreaInput.Rd,textInput.Rd,textOutput.Rd,timer.R,titlePanel.Rd,update-input.R,updateActionButton.Rd,updateCheckboxGroupInput.Rd,updateCheckboxInput.Rd,updateDateInput.Rd,updateDateRangeInput.Rd,updateNumericInput.Rd,updateQueryString.Rd,updateRadioButtons.Rd,updateSelectInput.Rd,updateSliderInput.Rd,updateTabsetPanel.Rd,updateTextAreaInput.Rd,updateTextInput.Rd,urlModal.Rd,useBusyIndicators.Rd,utils-lang.R,utils-tags.R,utils.R,validate.Rd,varSelectInput.Rd,version_bs_date_picker.R,version_ion_range_slider.R,version_jquery.R,version_jqueryui.R,version_selectize.R,version_strftime.R,verticalLayout.Rd,viewer.R,viewer.Rd,wellPanel.Rd,withMathJax.Rd,withOtelCollect.Rd,withProgress.Rd}|01_hello:{app.R}|02_text:{app.R}|03_reactivity:{app.R}|04_mpg:{app.R}|05_sliders:{app.R}|06_tabsets:{app.R}|07_widgets:{app.R}|09_upload:{app.R}|10_download:{app.R}|11_timer:{app.R}<!-- RSHINY-AGENTS-MD-END -->

## Package Overview

**tstse** is an R package for time series statistical estimation focused on ARMA/ARIMA/GARCH modeling. It provides simulation, estimation, model selection, and diagnostic tools with C++ acceleration via Rcpp.

**Current Status**: Phase 1 (tswge rewrite) is complete. Phase 2 focuses on extended features: flexible generators, generalized bootstrap, GARCH tools, and performance optimization.

## Build & Development Commands

```bash
# Load package for development (without installing)
R -e "devtools::load_all()"

# Regenerate documentation from Roxygen2 comments
R -e "devtools::document()"

# Run all tests
R -e "devtools::test()"

# Run a single test file
R -e "testthat::test_file('tests/testthat/test-est_ar.R')"

# Full package check
R CMD check --no-manual .

# Install from source
R CMD INSTALL .

# Regenerate Rcpp exports after modifying C++ code
R -e "Rcpp::compileAttributes()"
```

## Architecture

### R/C++ Integration Pattern

The package uses Rcpp for performance-critical operations. After modifying C++ code in `src/`, run `Rcpp::compileAttributes()` to regenerate exports.

**Current C++ functions** (in `src/`):
- `backcast_cpp()` - Backcasting residuals for ARMA models
- `gen_arch_cpp()` - ARCH process generation
- `gen_garch_cpp()` - GARCH process generation
- `sigplusnoise_signal_cpp()` - Signal component generation
- `gegenb_cpp()` - Gegenbauer coefficient computation
- `convolve_truncated_cpp()` - Truncated polynomial convolution

**Auto-generated files** (do not edit manually):
- `src/RcppExports.cpp` - C++ side
- `R/RcppExports.R` - R side

### Parallel Processing

Platform-aware parallelization pattern used throughout:
- Unix/Mac: `mclapply()` (fork-based, preferred)
- Windows: `parLapply()` (socket-based fallback)
- Control via `options(tstse.cores = N)` or function arguments


### S3 Classes

Output classes with custom print/plot methods:
- `est_ar`, `est_arma` - Estimation results
- `aic_ar`, `aic_arma` - Model selection results
- `ljung_box_test` - Diagnostic test results
- `fore_aruma`, `fore_arima`, `fore_arma` - Forecast results
- `aruma` - Generator output (Phase 2)
- `garch.comparison.tse` - GARCH comparison results (Phase 2)

### MA Coefficient Convention

MA coefficients use **negated sign** from R's `arima()` to match textbook (ATSA) conventions. This is intentional and documented in `est_arma()`.

## DuckDB Monte Carlo Storage

Bootstrap simulation results are stored in DuckDB for efficient columnar storage and querying. The schema is defined in `R/duckdb_utils.R` and initialized idempotently via `mc_db_init()`.

### Schema

**2 tables** with a flat design optimized for a single capstone study. A `batch` groups simulation iterations from one execution run; the same config can be run multiple times across batches and results pooled or compared.

```
batches
├── batch_id    INTEGER PK (auto-increment)
├── label       VARCHAR
└── created_at  TIMESTAMP DEFAULT now

simulations
├── sim_id          INTEGER PK (auto-increment)
├── batch_id        INTEGER FK → batches
├── iteration       INTEGER NOT NULL
├── n               INTEGER NOT NULL        -- observation length
├── phi             DOUBLE NOT NULL         -- AR(1) coefficient
├── innov_dist      VARCHAR NOT NULL        -- 'norm','t(3)','arch(0.2)'
├── obs_stat        DOUBLE NOT NULL         -- CO t-statistic
├── boot_dist       DOUBLE[]                -- full bootstrap distribution
├── pvalue          DOUBLE                  -- two-sided (abs-based)
├── pvalue_upper    DOUBLE                  -- upper-tail
├── pvalue_lower    DOUBLE                  -- lower-tail
├── pvalue_asymp    DOUBLE                  -- asymptotic
├── pvalue_adj      DOUBLE                  -- COBA adjusted
├── null_ar_order   INTEGER                 -- selected AR order
├── null_ar_phi     DOUBLE[]                -- fitted AR coefficients
├── null_vara       DOUBLE                  -- innovation variance
└── created_at      TIMESTAMP DEFAULT now
```

**Indexes**: `idx_sim_config` on (n, phi, innov_dist), `idx_sim_batch` on (batch_id).

### Views

- **`v_rejection_rates`** - Pooled rejection rates across all batches. Groups by (n, phi, innov_dist). Computes `reject_05` (bootstrap), `reject_asymp_05` (asymptotic), `reject_adj_05` (COBA) with binomial standard errors.
- **`v_rejection_rates_by_batch`** - Per-batch rejection rates for comparing individual runs. Same columns plus batch label and timestamp.

### Key Functions (`R/duckdb_utils.R`)

| Function | Purpose |
|----------|---------|
| `mc_db_connect(path, read_only)` | Open DuckDB connection |
| `mc_db_init(con)` | Create schema (idempotent) |
| `mc_db_write_batch(con, results, batch_id, n, phi, innov_dist)` | Batch write `wbg_boot_fast` results in a transaction |
| `mc_db_query(con, n, phi, innov_dist, batch_id, limit)` | Query simulations with optional filters |
| `mc_db_rejection_rates(con, by_batch, n, phi, innov_dist)` | Get pooled or per-batch rejection rates |
| `mc_db_recalc_pvalue(con, sim_id)` | Recompute p-values from stored boot_dist |
| `mc_study(path, label)` | High-level closure: `$save_batch()`, `$query()`, `$rejection_rates()`, `$end()` |

### Shiny Viewer

`boot_db_viewer(db_path)` launches an interactive Shiny app (`inst/shiny/boot_viewer/app.R`) for exploring results. Features: cascading filters (n, phi, innov_dist), pooled vs per-batch toggle, color-coded rejection rate tables, power curve plots, heatmaps, bootstrap distribution histograms, CSV/PNG export. DB path can also be set via `options(tstse.viewer_db = "path")`.

### Design Patterns

- **Auto-increment integers**: batch_id and sim_id use DuckDB sequences (no UUIDs)
- **Transactional batch writes**: All inserts wrapped in DBI transactions with rollback on failure
- **innov_dist as self-describing string**: `"norm"`, `"t(3)"`, `"arch(0.2)"` -- no separate params column
- **Direct field mapping**: `wbg_boot_fast` fields mapped 1:1 (tco_obs → obs_stat, boot_tstats → boot_dist, phi → null_ar_phi)

## Phase 2: Extended Features


## Code Style & Conventions

### Naming
- **Functions**: `snake_case` (e.g., `est_arma`, `gen_garch`, `true_acf`)
- **Factory functions**: `make_*` prefix (e.g., `make_gen_t`, `make_stat_co`)
- **Arguments**: `snake_case` (e.g., `lag_max`, `n_back`, not `lag.max`)
- **Internal helpers**: prefix with `.` (e.g., `.arma_variance`, `.parallel_map`)

### Function Arguments
- Use `seed = NULL` for RNG control (NULL = random, integer = reproducible)
- Use proper logicals `plot = TRUE` (not string `"TRUE"`)
- Use `match.arg()` for string arguments with fixed choices
- Factory functions should use `force()` on captured variables

### Return Values
- Use `invisible(result)` when function primarily plots
- Return S3 objects with print/plot methods for complex results
- Include input parameters in output for reproducibility

### Dependencies
- Core dependencies go in `Imports` (Rcpp)
- Optional packages go in `Suggests` (rugarch, sn, gt, cli, Kendall, sandwich, lmtest, nlme)
- Use `requireNamespace("pkg", quietly = TRUE)` for conditional use
- Fail fast with helpful install instructions

### Plotting
- Always save and restore `par()` via `on.exit(par(op))`
- Keep plot helpers as internal functions (`.plot_*`)
- Consider ggplot2 for new visualization functions

## Testing

Tests use testthat (Edition 3). One test file per major function in `tests/testthat/`.

**Testing priorities for Phase 2**:
1. Factory functions return callable functions with correct signature
2. Parallel and sequential results match (given same seeds)
3. Edge cases: empty phi, single observation, near-unit-root
4. Stationarity checks catch non-stationary models
5. GARCH comparison handles convergence failures gracefully
6. For existing functions that change, we must be able to still pass equilivance tests against tswge where applicable. 

## Performance Guidelines

### When to use Rcpp
- Inner loops with >10,000 iterations
- Recursive computations (AR/MA filtering)
- Operations that can't be vectorized in R

### When to use mclapply
- Bootstrap replications (100+ iterations)
- Grid search over model orders
- Rolling window computations
- Always provide sequential fallback for Windows

### Profiling
```r
# Profile bottlenecks
Rprof("profile.out")
# ... code to profile ...
Rprof(NULL)
summaryRprof("profile.out")

# Or use profvis
profvis::profvis({ ... })
```

## File Organization

```
R/
├── RcppExports.R          # Auto-generated
├── est_*.R                # Estimation functions
├── fore_*.R               # Forecasting functions
├── gen_*.R                # Generation functions
├── aic*.R                 # Model selection
├── duckdb_utils.R         # DuckDB schema, connection, read/write utilities
├── boot_db_viewer.R       # Shiny viewer launcher
├── generators.R           # Innovation factory functions (Phase 2)
├── stat_generators.R      # Statistic factory functions (Phase 2)
├── wbg_boot_test.R        # Generalized bootstrap (Phase 2)
├── compare_garch.R        # GARCH comparison (Phase 2)
└── display_garch.R        # GARCH display utilities (Phase 2)

src/
├── RcppExports.cpp        # Auto-generated
├── backcast.cpp
├── gen_arch.cpp
├── gen_garch.cpp
├── gen_ar_fast.cpp        # (Phase 2 candidate)
└── ...

inst/shiny/boot_viewer/
└── app.R                  # DuckDB results viewer Shiny app

R_extended/                # Prototypes (not installed)
├── generators_tse.R
├── stat_generators_tse.R
├── wbg_boot_test_tse.R
├── gen_aruma_tse.R
├── gen_ar_fast_tse.R
├── ar_mle_tse.R
├── compare_garch_tse.R
└── display_garch_tse.R

tests/testthat/
├── test-*.R
└── test-duckdb_utils.R    # DuckDB schema and utility tests
```

