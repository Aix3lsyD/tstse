# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

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

## Key Function Groups

### Core (Phase 1 - Complete)
- **Generation**: `gen_arch()`, `gen_arima()`, `gen_arma()`, `gen_aruma()`, `gen_garch()`, `gen_geg()`, `gen_sigplusnoise()`
- **Estimation**: `est_ar()` (MLE/Burg/Yule-Walker), `est_arma()` (MLE), `slr()` (signal+noise)
- **Forecasting**: `fore_aruma()`, `fore_arima()`, `fore_arma()`, `fore_farma()`, `fore_sigplusnoise()`
- **Model Selection**: `aic()`, `aic_ar()`, `aic5()`, `aic_burg()`
- **Diagnostics**: `ljung_box()`, `backcast()`
- **Trend Testing**: `co()` (Cochrane-Orcutt), `wbg_boot()` (WBG bootstrap)

### Extended (Phase 2 - In Progress)
See `R_extended/` for prototypes.

## Phase 2: Extended Features

### Priority 1: Flexible ARUMA Generator

**Goal**: ARUMA generation with non-normal/non-iid innovations (t-distributed, skew-t, GARCH, mixture).

**Prototype**: `gen_aruma_tse.R`, `generators_tse.R`

**Design Pattern** - Factory functions for innovation generators:
```r
# Factory creates a generator function
make_gen_t <- function(df, scale = FALSE) {
  force(df); force(scale)
  function(n) {
    x <- rt(n, df = df)
    if (scale && df > 2) x <- x * sqrt((df - 2) / df)
    x
  }
}

# Usage
t_gen <- make_gen_t(df = 5, scale = TRUE)
result <- gen_aruma(n = 500, phi = 0.7, innov_gen = t_gen)
```

**Innovation generators to implement**:
| Generator | File | Status | Notes |
|-----------|------|--------|-------|
| `make_gen_norm()` | generators.R | Prototype | Standard normal, configurable sd |
| `make_gen_t()` | generators.R | Prototype | Student's t, optional scaling |
| `make_gen_skt()` | generators.R | Prototype | Skew-t via `sn` package |
| `make_gen_unif()` | generators.R | Prototype | Uniform, unit variance by default |
| `make_gen_mixnorm()` | generators.R | Prototype | Mixture of normals (outliers) |
| `make_gen_garch()` | generators.R | Prototype | GARCH process via `rugarch` |

**Key implementation notes**:
- Precompute theoretical moments at factory time (not per-call)
- Use `force()` to capture closure variables
- Single call to `innov_gen(total_n)` preserves dependence (critical for GARCH)
- Handle burn-in properly: `n_start` for ARMA stabilization, `spin` for nonstationary operators

### Priority 2: Generalized WBG Bootstrap

**Goal**: Bootstrap framework that accepts any test statistic function, not just Cochrane-Orcutt.

**Prototype**: `wbg_boot_test_tse.R`, `stat_generators_tse.R`

**Design Pattern** - Statistic factory functions:
```r
# Factory creates a statistic function
make_stat_co <- function(maxp = 5, ar_method = "mle") {
  function(x) {
    co(x, maxp = maxp, ar_method = ar_method)$tco
  }
}

# Usage
stat_fn <- make_stat_co(maxp = 5)
result <- wbg_boot_test(x, stat_fn = stat_fn, nb = 399, parallel = TRUE)
```

**Statistic generators to implement**:
| Generator | Purpose | Dependencies |
|-----------|---------|--------------|
| `make_stat_co()` | Cochrane-Orcutt t-stat | Internal |
| `make_stat_ols_t()` | OLS trend t-stat | None |
| `make_stat_ols_slope()` | OLS slope estimate | None |
| `make_stat_mk()` | Mann-Kendall S | `Kendall` |
| `make_stat_spearman()` | Spearman correlation | None |
| `make_stat_sen()` | Sen's slope | None |
| `make_stat_hac()` | Newey-West t-stat | `sandwich`, `lmtest` |
| `make_stat_bn()` | Bloomfield-Nychka t-stat | None |
| `make_stat_lr()` | Likelihood ratio | None |
| `make_stat_gls()` | GLS t-stat | `nlme` |

**Performance requirements**:
- Use `gen_ar_fast()` (not `gen_arma()`) for AR bootstrap - 30-50x speedup
- Support parallel execution via `mclapply`/`parLapply`
- Reproducible via explicit `boot_seeds` vector

### Priority 3: Fast AR Generator

**Goal**: Optimized AR generation for bootstrap (adaptive burn-in, C-implemented filter).

**Prototype**: `gen_ar_fast_tse.R`

**Key optimizations**:
```r
gen_ar_fast <- function(n, phi, vara = 1, seed = NULL) {
  # Adaptive burn-in based on persistence
  persistence <- sum(abs(phi))
  burn <- if (persistence >= 0.999) {
    max(50, round(n * 0.5), 500)
  } else if (persistence < 0.01) {
    50
  } else {
    max(50, round(3 * abs(log(0.001) / log(persistence))))
  }
  burn <- min(burn, 2000)
  
  # Use stats::filter (C implementation)
  a <- rnorm(n + burn, sd = sqrt(vara))
  x <- stats::filter(a, filter = phi, method = "recursive")
  as.numeric(x[(burn + 1):(n + burn)])
}
```

**Rcpp candidate**: If R version is bottleneck, implement in C++.

### Priority 4: GARCH Comparison Tools

**Goal**: Systematic GARCH model comparison with diagnostics and display.

**Prototype**: `compare_garch_tse.R`, `display_garch_tse.R`

**Components**:
- `compare_garch()` - Fit grid of ARCH/GARCH models, compute IC and diagnostics
- `table_garch_gt()` - Publication-ready gt table with color coding
- `table_garch_cli()` - Terminal display via cli package
- `table_coef_gt()` - Coefficient table for single model

**Diagnostics computed**:
- Information criteria: AIC, AICc, BIC
- Weighted Ljung-Box on squared standardized residuals
- Nyblom stability test
- Sign bias test
- Coefficient significance counts

### Priority 5: MLE AR with Stationarity Check

**Goal**: MLE estimation with fallback through top-N candidates if non-stationary.

**Prototype**: `ar_mle_tse.R`

**Key features**:
- `check_stationary()` - Verify roots outside unit circle
- `aic_ar_mle()` - MLE with IC selection
- Fallback: check top `n_best` models for stationarity
- Option to disable stationarity enforcement for simulation studies

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
└── test-*.R
```

## Migration Notes

### From tswge
Phase 1 rewrote all tswge functions. Key mappings documented in git history.

### From Prototypes to Production
When moving from `R_extended/` to `R/`:
1. Remove `.tse` suffix from function names
2. Update internal function calls to use package functions
3. Add roxygen2 documentation
4. Add testthat tests
5. Update NAMESPACE via `devtools::document()`
6. Check `R/` for existing functions we can improve on.  We might not need to make a whole new one. 
