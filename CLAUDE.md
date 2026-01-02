# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

**tstse** is an R package for time series statistical estimation focused on ARMA/ARIMA/GARCH modeling. It provides simulation, estimation, model selection, and diagnostic tools with C++ acceleration via Rcpp.

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
```

## Architecture

### R/C++ Integration Pattern
The package uses Rcpp for performance-critical operations. Six C++ functions are exposed:
- `backcast_cpp()` - Backcasting residuals for ARMA models (src/backcast.cpp)
- `gen_arch_cpp()` - ARCH process generation (src/gen_arch.cpp)
- `gen_garch_cpp()` - GARCH process generation (src/gen_garch.cpp)
- `sigplusnoise_signal_cpp()` - Signal component generation (src/gen_sigplusnoise.cpp)
- `gegenb_cpp()` - Gegenbauer coefficient computation (src/gegenb.cpp)
- `convolve_truncated_cpp()` - Truncated polynomial convolution (src/convolve_truncated.cpp)

R wrappers handle input validation, then call C++ via `.Call()`. Auto-generated files:
- `src/RcppExports.cpp` - C++ side (do not edit manually)
- `R/RcppExports.R` - R side (do not edit manually)

After modifying C++ code, run `Rcpp::compileAttributes()` to regenerate exports.

### S3 Classes
Output classes with custom print methods:
- `est_ar`, `est_arma` - Estimation results
- `aic_ar`, `aic_arma` - Model selection results
- `ljung_box_test` - Diagnostic test results
- `fore_aruma`, `fore_arima`, `fore_arma` - Forecast results

### Parallel Processing
- Platform-aware: uses `mclapply` (Unix) or `parLapply` (Windows)
- Controlled via `options(tstse.cores = N)` or `getOption("tstse.cores", 1L)`
- Used in: `aic()`, `aic_ar()` (grid search), `wbg_boot()` (bootstrap), `roll_win_rmse()` (rolling windows)

### MA Coefficient Convention
MA coefficients use **negated sign** from R's `arima()` to match textbook conventions. This is intentional and documented in `est_arma()`.

## Key Function Groups

- **Generation**: `gen_arch()`, `gen_arima()`, `gen_arma()`, `gen_aruma()`, `gen_garch()`, `gen_geg()`, `gen_sigplusnoise()`
- **Estimation**: `est_ar()` (MLE/Burg/Yule-Walker), `est_arma()` (MLE), `slr()` (signal+noise)
- **Forecasting**: `fore_aruma()`, `fore_arima()`, `fore_arma()`, `fore_farma()`, `fore_sigplusnoise()` (see hierarchy below)
- **Model Selection**: `aic()`, `aic_ar()`, `aic5()`, `aic_burg()`
- **Diagnostics**: `ljung_box()`, `backcast()`
- **Trend Testing**: `co()` (Cochrane-Orcutt), `wbg_boot()` (WBG bootstrap)
- **Evaluation**: `roll_win_rmse()`, `roll_win_rmse_nn()`
- **Theoretical Properties**: `true_acf()`, `true_garma_acf()`, `true_spec()`, `psi_weights()`, `pi_weights()`
- **Spectral**: `parzen()`, `period()`, `sample_spec()`
- **Visualization**: `plotts()`, `plotts_sample()`, `plotts_true()`, `plotts_parzen()`, `plotts_dwt()`, `plotts_mra()`, `unit_circle()`
- **Filtering**: `butterworth()`, `expsmooth()`, `kalman()`, `kalman_miss()`, `ma_smooth()`, `ma_pred()`
- **G-Lambda Transform**: `is_glambda()`, `is_sample()`, `trans_to_dual()`, `trans_to_original()`
- **Utilities**: `factor()`, `factor_comp()`, `artrans()`, `mult()`, `gegenb()`, `macoef_geg()`, `hilbert()`, `wv()`

### Forecasting Function Hierarchy

```
fore_arma  ──▶  fore_arima  ──▶  fore_aruma
(d=0,s=0)       (lambda=0)       (full impl)
```

- `fore_aruma(x, phi, theta, d, s, lambda, ...)` — Full ARUMA implementation
- `fore_arima(x, phi, theta, d, s, ...)` — Wrapper with lambda=0
- `fore_arma(x, phi, theta, ...)` — Wrapper with d=0, s=0, adds RMSE/MAD for holdout mode

## Testing

Tests use testthat (Edition 3). One test file per major function in `tests/testthat/`.

## Code Style & Conventions

### Naming
- **Functions**: `snake_case` (e.g., `est_arma`, `gen_garch`, `true_acf`)
- **Arguments**: `snake_case` (e.g., `lag_max`, `n_back`, not `lag.max`)
- **Internal helpers**: prefix with `.` (e.g., `.arma_variance`, `.plot_acf_panel`)

### Function Arguments
- Use `seed = NULL` for RNG control (NULL = random, integer = reproducible)
- Use proper logicals `plot = TRUE` (not string `"TRUE"`)
- Use `match.arg()` for string arguments with fixed choices

### Return Values
- Use `invisible(result)` when function primarily plots
- Return S3 objects with print methods for estimation/test results

### Dependencies
- Core dependencies go in `Imports` (Rcpp)
- Optional packages go in `Suggests` (astsa, PolynomF, ggplot2, zoo)
- Use `requireNamespace("pkg", quietly = TRUE)` for conditional use

### Plotting Functions
- Always save and restore `par()` via `on.exit(par(op))`
- Keep plot helpers as internal functions (`.plot_*`)

## Migration from tswge

This package rewrites functions from the `tswge` package. When adding new functions:

1. Review original for bugs/issues before rewriting
2. Use our existing functions as dependencies (e.g., `backcast()`, `est_ar()`)
3. Write tests comparing against original with appropriate tolerance
4. Document any intentional deviations from original behavior

### Completed Functions (60)
See `R/` directory and PROGRESS.md for full list. Key mappings:
- `artrans.wge` → `artrans()`
- `backcast.wge` → `backcast()`
- `co.wge` → `co()`
- `est.arma.wge` → `est_arma()`
- `est.ar.wge` → `est_ar()`
- `fore.aruma.wge` → `fore_aruma()`
- `fore.arima.wge` → `fore_arima()`
- `fore.arma.wge` → `fore_arma()`
- `fore.sigplusnoise.wge` → `fore_sigplusnoise()`
- `gen.arch.wge` → `gen_arch()`
- `gen.arima.wge` → `gen_arima()`
- `gen.aruma.wge` → `gen_aruma()`
- `gen.geg.wge` → `gen_geg()`
- `fore.farma.wge` → `fore_farma()`
- `est.farma.wge` → `est_farma()`
- `is.glambda.wge` → `is_glambda()`
- `is.sample.wge` → `is_sample()`
- `unit.circle.wge` → `unit_circle()`
- `aic.wge` → `aic()`
- `roll.win.rmse.wge` → `roll_win_rmse()`
- `trans.to.dual.wge` → `trans_to_dual()`
- `trans.to.original.wge` → `trans_to_original()`
- `true.arma.aut.wge` → `true_acf()`
- `true.arma.spec.wge` → `true_spec()`
- `wbg.boot.wge` → `wbg_boot()`

### Known Original Bugs Fixed
- `gen.garch.wge`: Crashes when p0 > q0 (array bounds)
- `true.arma.spec.wge`: Returns scalar `f` instead of vector when plot=FALSE
- `true.arma.aut.wge`: `plot=="TRUE"` string comparison instead of logical
- `fore.arima.wge`: Dead lambda/Box-Cox code (hardcoded to 0)
- `fore.arma.wge`, `fore.arima.wge`: `plot=='TRUE'` string comparison bugs
